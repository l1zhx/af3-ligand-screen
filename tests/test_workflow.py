"""Input contracts, Slurm submission recovery, and output validation."""
import argparse
from contextlib import nullcontext
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import unittest
from unittest.mock import patch
import uuid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import demo
import workflow as w


class WorkflowTest(unittest.TestCase):
    def setUp(self):
        self.root = ROOT / 'tests' / ('scratch_' + uuid.uuid4().hex[:6])
        self.root.mkdir()
        self.addCleanup(shutil.rmtree, self.root)

    def bundle(self):
        source = self.root / 'msa.json'
        source.write_text(json.dumps(dict(sequences=[dict(protein=dict(id='A', sequence='AAA',
                          unpairedMsa='>arbitrary_query_header\nAA\nA\n', pairedMsa='old', templates=[]))])))
        manifest = self.root / 'ligands.csv'
        manifest.write_text('id,inchikey,smiles\n1,DEMO_KEY_1,CC\n2,DEMO_KEY_2,CO\n')
        bundle = self.root / 'bundle'
        subprocess.run([sys.executable, str(ROOT/'scripts/build_jobs.py'), '--msa-json', str(source),
                        '--ligands-csv', str(manifest), '--out-dir', str(bundle)], check=True, capture_output=True)
        return bundle

    def test_shared_msa_input_names_and_plan(self):
        bundle = self.bundle()
        jobs = sorted((bundle/'jobs').glob('*/*.json'))
        self.assertEqual(len(jobs), 2)
        job = json.loads(jobs[0].read_text())
        protein = job['sequences'][0]['protein']
        self.assertEqual(protein['pairedMsa'], '')
        self.assertEqual(job['name'], 'bcat_lig_000001_seed_1')
        self.assertTrue((jobs[0].parent/protein['unpairedMsaPath']).is_file())
        w.plan(argparse.Namespace(work=bundle, round='seed_1', start_id=1, end_id=None,
                                 out=self.root/'plan', chunk_size=1, batch_tasks=1, samples_per_seed=5))
        data = json.loads((self.root/'plan/plan.json').read_text())
        self.assertEqual(len(data['batches']), 2)
        self.assertEqual(len(data['shared_sha256']), 1)
        cmd = w.sbatch_command(self.root/'plan/plan.json', 1, data, 12, '1234')
        self.assertIn('afterany:1234', cmd)
        self.assertIn('0-0%12', cmd)
        self.assertNotIn('--resolve_msa_overlaps', ' '.join(cmd))

    def test_complete_check_all_seeds_and_root_copy(self):
        demo.make_fixture(self.root/'fixture')
        folder = self.root/'fixture/results/000/lig_000001/seed_1/attempt_100_0'
        self.assertTrue(w.valid_attempt(folder, [1], 5))
        self.assertFalse(w.valid_attempt(folder, [1, 2], 5))
        next(folder.rglob('model.cif')).write_text('')
        self.assertFalse(w.valid_attempt(folder, [1], 5))

    @patch('workflow.lock', side_effect=lambda path: nullcontext())
    def test_submit_ledger_resumes_and_unknown_blocks(self, mocked_lock):
        bundle = self.bundle()
        w.plan(argparse.Namespace(work=bundle, round='seed_1', start_id=1, end_id=None,
                                 out=self.root/'plan', chunk_size=1, batch_tasks=1, samples_per_seed=5))
        args = argparse.Namespace(plan=self.root/'plan/plan.json', max_active=4, max_batches=1, dry_run=False)
        with patch('workflow.subprocess.run', return_value=subprocess.CompletedProcess([], 0, '1234;cluster\n', '')) as run:
            w.submit(args)
            self.assertEqual(run.call_count, 1)
        with patch('workflow.subprocess.run', return_value=subprocess.CompletedProcess([], 0, '1235\n', '')) as run:
            w.submit(args)
            self.assertIn('afterany:1234', run.call_args.args[0])
        with patch('workflow.subprocess.run') as run:
            w.submit(args)
            run.assert_not_called()
        ledger_path = self.root/'plan/submission.json'
        ledger = json.loads(ledger_path.read_text())
        ledger['inflight'] = dict(batch=2)
        w.save(ledger_path, ledger)
        with self.assertRaisesRegex(SystemExit, 'unknown status'):
            w.submit(args)

    @patch('workflow.lock', side_effect=lambda path: nullcontext())
    def test_worker_marker_skip_changed_input(self, mocked_lock):
        bundle = self.bundle()
        w.plan(argparse.Namespace(work=bundle, round='seed_1', start_id=1, end_id=1,
                                 out=self.root/'plan', chunk_size=1, batch_tasks=1, samples_per_seed=5))
        af3 = self.root/'af3'
        af3.mkdir()
        (af3/'run_alphafold.py').write_text('# Fake AF3 runner; subprocess mocked by test\n')
        model = self.root/'model'
        model.mkdir()
        env = dict(WORK_DIR=str(bundle), AF3_DIR=str(af3), AF3_PYTHON=sys.executable, MODEL_DIR=str(model), SLURM_JOB_ID='900')
        demo.make_fixture(self.root/'fixture')
        source = self.root/'fixture/results/000/lig_000001/seed_1/attempt_100_0'
        def fake_run(cmd, **kwargs):
            output = Path(next(x.split('=', 1)[1] for x in cmd if x.startswith('--output_dir=')))
            shutil.copytree(source, output, dirs_exist_ok=True)
            return subprocess.CompletedProcess(cmd, 0)
        args = argparse.Namespace(plan=self.root/'plan/plan.json', batch=0, task=0)
        with patch.dict(os.environ, env), patch('workflow.subprocess.run', side_effect=fake_run) as call:
            w.run(args)
            w.run(args)
            self.assertEqual(call.call_count, 1)
        marker = bundle/'results/000/lig_000001/seed_1/complete.json'
        self.assertTrue(marker.is_file())
        job = next((bundle/'jobs').glob('*/*.json'))
        job.write_text(job.read_text() + ' ')
        with patch.dict(os.environ, env), self.assertRaisesRegex(SystemExit, 'Failed inputs'):
            w.run(args)

    @unittest.skipUnless(os.name == 'posix', 'fcntl locking targets Linux')
    def test_exclusive_lock(self):
        with w.lock(self.root/'lock'):
            with self.assertRaisesRegex(RuntimeError, 'Already locked'):
                with w.lock(self.root/'lock'):
                    pass

    @unittest.skipUnless(importlib.util.find_spec('rdkit'), 'RDKit optional for analysis-only installation')
    def test_library_inchikey_dedup_preserves_record_audit(self):
        from rdkit import Chem
        import prepare
        source = self.root/'sdf'
        source.mkdir()
        writer = Chem.SDWriter(str(source/'tiny.sdf'))
        for smiles in ['CCO', 'OCC', 'CC(=O)O']:
            writer.write(Chem.MolFromSmiles(smiles))
        writer.close()
        prepare.library(argparse.Namespace(sdf_dir=source, out=self.root/'prepared'))
        audit = json.loads((self.root/'prepared/preparation_audit.json').read_text())
        self.assertEqual((audit['records'], audit['unique'], audit['duplicates']), (3, 2, 1))


if __name__ == '__main__':
    unittest.main()
