# Validation record

## Local validation scope

The portable repository was assembled from the working campaign scripts and tested on Windows with Python 3.12 and NumPy. Synthetic outputs exercise the analysis path without invoking AF3. The actual validation counts and environment are recorded below after the local test run.

## Automated checks

- Bare plural, bare singular and job-prefixed sample summary filenames.
- Exclusion of root best copies and accidental higher-scoring repeated attempts.
- Complete, partial, missing and unrecognized result groups.
- Explicit chain order and directional PAE extraction.
- Protein rigid transform invariance and a displaced ligand's RMSD.
- Shared MSA paths, wrapped query sequence, unique job names and empty paired MSA.
- Requirement for every planned seed in a successful inference attempt.
- Slurm command construction, continuation of a submission ledger and blocking uncertain submissions.
- Simulated AF3 worker completion, skipping verified outputs and rejecting changed inputs.
- Linux advisory lock exclusion (requires Linux).
- InChIKey deduplication and source-record provenance (requires RDKit).

Submission and worker tests use mocked subprocesses; they never submit jobs or execute AF3. This verifies orchestration behavior without representing a real cluster run. The Linux CI workflow additionally exercises advisory locks and RDKit when dependencies are available. Real AF3 compatibility, distributed filesystem locking and resource requests must still be checked by the four-ligand cluster smoke test.

The full-confidence geometry parser is tested with minimal synthetic atom tables. Before scaling scientific interpretation, visually inspect representative real CIFs and compare numbering/chain labels to the original receptor and chemical records.

## Recorded local run

- Date: 2026-10-15.
- Environment: Windows, Python 3.12.14, NumPy 2.3.5.
- Unit tests: 8 discovered; 6 passed; 2 skipped (native Linux advisory locks and RDKit availability).
- CPU demo: passed; 2 complete ligand groups, 10 retained samples, 1 excluded duplicate, 1 partial, 1 missing, 1 exported structure.
- Python source compilation: passed.
- Submission and worker subprocesses were mocked; no cluster jobs were launched during this verification.
- Bash syntax, native Linux locks, RDKit preparation/Glide matching and GPU execution were not executed in this local environment. Linux CI is configured for shell syntax, lock and RDKit preparation checks after upload; Glide matching still needs a real-data integration check.
