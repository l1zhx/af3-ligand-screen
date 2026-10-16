#!/usr/bin/env python3
"""Plan bounded Slurm arrays, submit with a ledger, and validate each AF3 run."""
import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys

REPO = Path(__file__).resolve().parents[1]
JOB = re.compile(r'bcat_(lig_(\d{6}))_([A-Za-z0-9_]+)\.json$')
SAMPLE = re.compile(r'seed-(\d+)_sample-(\d+)$')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def save(path, data):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')
    tmp.replace(path)


@contextmanager
def lock(path):
    # flock locks are released when a killed worker exits. Verify support on your shared FS.
    import fcntl
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a+') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError(f'Already locked by another process: {path}')
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def input_record(path):
    match = JOB.fullmatch(path.name)
    if not match:
        raise ValueError(f'Unexpected job name: {path}')
    data = json.loads(path.read_text(encoding='utf-8'))
    if data['name'] != path.stem:
        raise ValueError(f'JSON name and filename differ: {path}')
    proteins = [s['protein'] for s in data['sequences'] if 'protein' in s]
    ligands = [s['ligand'] for s in data['sequences'] if 'ligand' in s]
    if len(proteins) != 1 or len(ligands) != 1 or proteins[0]['id'] != 'A' or ligands[0]['id'] != 'B':
        raise ValueError('This workflow requires one protein A and one ligand B')
    protein = proteins[0]
    if protein.get('pairedMsa') != '' or 'pairedMsaPath' in protein:
        raise ValueError('pairedMsa must be explicitly empty')
    refs = []
    if 'unpairedMsaPath' in protein:
        refs.append((path.parent / protein['unpairedMsaPath']).resolve())
    elif not protein.get('unpairedMsa'):
        raise ValueError('Missing unpaired MSA')
    for template in protein.get('templates', []):
        if 'mmcifPath' in template:
            refs.append((path.parent / template['mmcifPath']).resolve())
    seeds = data['modelSeeds']
    if not seeds or any(type(s) is not int or s < 0 for s in seeds) or len(set(seeds)) != len(seeds):
        raise ValueError('Invalid model seeds')
    return dict(path=str(path.resolve()), sha256=digest(path), name=data['name'],
                ligand=match[1], ligand_number=int(match[2]), round=match[3],
                shard=f'{(int(match[2])-1)//1000:03d}', seeds=seeds, refs=[str(p) for p in refs])


def plan(args):
    work = args.work.resolve()
    entries, shared, seen = [], {}, set()
    for path in sorted((work / 'jobs').glob('*/*.json')):
        match = JOB.fullmatch(path.name)
        if not match or match[3] != args.round or int(match[2]) < args.start_id:
            continue
        if args.end_id is not None and int(match[2]) > args.end_id:
            continue
        rec = input_record(path)
        if rec['ligand'] in seen:
            raise ValueError(f'Duplicate input for {rec["ligand"]}')
        seen.add(rec['ligand'])
        for ref in rec['refs']:
            if ref not in shared:
                shared[ref] = digest(ref)
        entries.append(rec)
    if not entries:
        raise SystemExit('No matching jobs')
    args.out.mkdir(parents=True, exist_ok=False)
    width = args.chunk_size * args.batch_tasks
    batches = [entries[i:i+width] for i in range(0, len(entries), width)]
    save(args.out / 'plan.json', dict(schema=1, work=str(work), round=args.round,
         chunk_size=args.chunk_size, samples_per_seed=args.samples_per_seed, shared_sha256=shared,
         created=datetime.now(timezone.utc).isoformat(), batches=batches))
    print(f'Planned {len(entries)} inputs in {len(batches)} arrays; plan={args.out.resolve() / "plan.json"}')


def batch_tasks(batch, chunk):
    return math.ceil(len(batch) / chunk)


def sbatch_command(plan_path, batch_index, data, throttle, previous):
    logs = Path(data['work']) / 'logs'
    cmd = ['sbatch', '--parsable', '--chdir', str(REPO), '--job-name=bcat-af3',
           '--partition', os.environ.get('GPU_PARTITION', 'gpu-c'),
           '--cpus-per-task', os.environ.get('CPUS_PER_TASK', '4'),
           '--mem', os.environ.get('MEMORY', '64G'), '--time', os.environ.get('TIME_LIMIT', '2-00:00:00'),
           '--array', f'0-{batch_tasks(data["batches"][batch_index], data["chunk_size"])-1}%{throttle}',
           '--output', str(logs / 'af3-%A_%a.out')]
    if previous:
        cmd += ['--dependency', f'afterany:{previous}']
    # PLAN and BATCH_INDEX are passed in the environment; avoids comma escaping in --export.
    return cmd + ['--export=ALL', str(REPO / 'slurm/infer.slurm')]


def submit(args):
    plan_path = args.plan.resolve()
    data = json.loads(plan_path.read_text())
    ledger_path = plan_path.parent / 'submission.json'
    (Path(data['work']) / 'logs').mkdir(exist_ok=True)
    with lock(plan_path.parent / 'submit.lock'):
        ledger = json.loads(ledger_path.read_text()) if ledger_path.exists() else dict(plan_sha256=digest(plan_path), accepted=[], inflight=None)
        if ledger['plan_sha256'] != digest(plan_path):
            raise SystemExit('Plan changed after submission')
        if ledger.get('inflight') is not None:
            raise SystemExit('An earlier submission has unknown status. Reconcile submission.json with squeue/sacct before continuing; see docs/HPC.md.')
        start = len(ledger['accepted'])
        previous = ledger['accepted'][-1]['job_id'] if start else None
        stop = len(data['batches']) if args.max_batches == 0 else min(len(data['batches']), start + args.max_batches)
        for i in range(start, stop):
            cmd = sbatch_command(plan_path, i, data, args.max_active, previous)
            if args.dry_run:
                print(json.dumps(dict(batch=i, argv=cmd, env=dict(PLAN=str(plan_path), BATCH_INDEX=str(i)))))
                previous = f'<accepted-array-{i}>'
                continue
            ledger['inflight'] = dict(batch=i, submitted_at=datetime.now(timezone.utc).isoformat())
            save(ledger_path, ledger)
            env = dict(os.environ, PLAN=str(plan_path), PLAN_SHA256=digest(plan_path), BATCH_INDEX=str(i))
            proc = subprocess.run(cmd, env=env, text=True, capture_output=True)
            response = proc.stdout.strip().split(';')[0]
            if proc.returncode or not response.isdigit():
                print(proc.stdout + proc.stderr, file=sys.stderr)
                raise SystemExit('Submission stopped; inflight retained to prevent an uncertain duplicate submission. See docs/HPC.md.')
            previous = response
            ledger['accepted'].append(dict(batch=i, job_id=response))
            ledger['inflight'] = None
            save(ledger_path, ledger)
            print(f'Submitted batch {i}: {response}', flush=True)
        if start == len(data['batches']):
            print('All batches in this plan have already been submitted. Check failures separately.')


def valid_attempt(root, seeds, count):
    """Require exactly one complete sample group for every requested seed."""
    from af3_postprocess import companion_files, summary_prefix, number
    samples = {}
    for path in Path(root).rglob('*.json'):
        match = SAMPLE.fullmatch(path.parent.name)
        if not match or summary_prefix(path.name) is None:
            continue
        key = (int(match[1]), int(match[2]))
        if key in samples:
            return False
        try:
            data = json.loads(path.read_text())
            if number(data.get('ranking_score')) is None:
                return False
            if any(not p.is_file() or p.stat().st_size == 0 for p in companion_files(path)):
                return False
        except (OSError, ValueError, TypeError):
            return False
        samples[key] = path
    return set(samples) == {(s, i) for s in seeds for i in range(count)}


def run(args):
    if os.environ.get('PLAN_SHA256') and digest(args.plan) != os.environ['PLAN_SHA256']:
        raise SystemExit('Plan changed after submission')
    data = json.loads(args.plan.read_text())
    batch = data['batches'][args.batch]
    first = args.task * data['chunk_size']
    selected = batch[first:first + data['chunk_size']]
    if not selected:
        raise SystemExit('Array task outside plan')
    work = Path(data['work'])
    if Path(os.environ['WORK_DIR']).resolve() != work:
        raise SystemExit('Configured WORK_DIR differs from submitted plan')
    af3 = Path(os.environ['AF3_DIR'])
    python = Path(os.environ['AF3_PYTHON'])
    model = Path(os.environ['MODEL_DIR'])
    if not (af3 / 'run_alphafold.py').is_file() or not python.is_file() or not model.is_dir():
        raise SystemExit('Check AF3_DIR/run_alphafold.py, AF3_PYTHON and MODEL_DIR')
    verified_refs, failed = set(), []
    for rec in selected:
        result = work / 'results' / rec['shard'] / rec['ligand'] / rec['round']
        try:
            result.mkdir(parents=True, exist_ok=True)
            with lock(result / 'worker.lock'):
                if digest(rec['path']) != rec['sha256']:
                    raise ValueError('Input JSON changed after plan creation')
                for ref in rec['refs']:
                    if ref not in verified_refs:
                        if digest(ref) != data['shared_sha256'][ref]:
                            raise ValueError(f'Shared input changed: {ref}')
                        verified_refs.add(ref)
                signature = hashlib.sha256(json.dumps(dict(input=rec, samples=data['samples_per_seed'],
                    shared={r: data['shared_sha256'][r] for r in rec['refs']}), sort_keys=True).encode()).hexdigest()
                marker = result / 'complete.json'
                if marker.exists():
                    old = json.loads(marker.read_text())
                    if old['signature'] != signature:
                        raise ValueError('Existing completion belongs to different inputs; choose a new round')
                    if valid_attempt(old['attempt'], rec['seeds'], data['samples_per_seed']):
                        print(f'SKIP complete {rec["ligand"]}', flush=True)
                        continue
                job_id = os.environ.get('SLURM_JOB_ID', 'local')
                attempt = result / f'attempt_{job_id}_{args.task}'
                # Never overwrite a failed attempt, including a Slurm requeue of the same task.
                attempt.mkdir(exist_ok=False)
                cmd = [str(python), str(af3 / 'run_alphafold.py'), f'--json_path={rec["path"]}',
                       f'--output_dir={attempt}', f'--model_dir={model}', '--run_data_pipeline=false',
                       '--flash_attention_implementation=' + os.environ.get('FLASH_ATTENTION_IMPLEMENTATION', 'xla')]
                print(f'{datetime.now(timezone.utc).isoformat()} START {rec["ligand"]}', flush=True)
                save(attempt / 'run_record.json', dict(command=cmd, input_signature=signature,
                     slurm_job_id=job_id, hostname=os.environ.get('SLURMD_NODENAME'),
                     cuda_visible_devices=os.environ.get('CUDA_VISIBLE_DEVICES'),
                     plan_sha256=digest(args.plan)))
                subprocess.run(cmd, cwd=af3, check=True)
                if not valid_attempt(attempt, rec['seeds'], data['samples_per_seed']):
                    raise ValueError('Incomplete sample group; no completion marker written')
                save(marker, dict(signature=signature, attempt=str(attempt), completed=datetime.now(timezone.utc).isoformat()))
                print(f'COMPLETE {rec["ligand"]}', flush=True)
        except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
            failed.append(rec['ligand'])
            print(f'FAILED {rec["ligand"]}: {exc}', file=sys.stderr, flush=True)
    if failed:
        raise SystemExit(f'Failed inputs: {failed}')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    a = sub.add_parser('plan')
    a.add_argument('--work', type=Path, required=True)
    a.add_argument('--out', type=Path, required=True)
    a.add_argument('--round', default='seed_1')
    a.add_argument('--start-id', type=int, default=1)
    a.add_argument('--end-id', type=int)
    a.add_argument('--chunk-size', type=int, default=4)
    a.add_argument('--batch-tasks', type=int, default=500)
    a.add_argument('--samples-per-seed', type=int, default=5)
    a.set_defaults(func=plan)
    a = sub.add_parser('submit')
    a.add_argument('--plan', type=Path, required=True)
    a.add_argument('--max-active', type=int, default=4)
    a.add_argument('--max-batches', type=int, default=5, help='0 submits all remaining arrays')
    a.add_argument('--dry-run', action='store_true')
    a.set_defaults(func=submit)
    a = sub.add_parser('run')
    a.add_argument('--plan', type=Path, required=True)
    a.add_argument('--batch', type=int, required=True)
    a.add_argument('--task', type=int, required=True)
    a.set_defaults(func=run)
    args = p.parse_args()
    for key in ('start_id', 'end_id', 'chunk_size', 'batch_tasks', 'samples_per_seed', 'max_active'):
        if getattr(args, key, None) is not None and getattr(args, key) < 1:
            p.error(f'{key} must be positive')
    for key in ('max_batches', 'batch', 'task'):
        if getattr(args, key, 0) < 0:
            p.error(f'{key} must be nonnegative')
    args.func(args)


if __name__ == '__main__':
    main()
