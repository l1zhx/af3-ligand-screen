# Methods and interpretation

## Biological task and input contract

The case study asks how to prioritize small-molecule poses against a fixed β-catenin construct using a pretrained complex-structure model and conventional docking outputs. The supplied protein sequence is 556 amino acids. Protein chain A and ligand chain B are fixed throughout this workflow. Generalizing to multiple proteins, covalent ligands, metals or alternate chain layouts needs deliberate changes to input validation and structural analysis.

The example FASTA is the original project construct. Its working mapping is `reference residue = AF3 label_seq_id + 130`, corresponding to reference positions 131–686. This offset was checked against the project's reference alignment, but a prepared docking receptor may use different numbering. Validate the sequence alignment before using the offset with any new receptor. `inspect` defaults to offset 0 and accepts explicit pocket residue lists through `--pockets`.

## Molecular identity

Preparation sanitizes SDF records with RDKit, removes explicit hydrogens through the supplier, obtains canonical isomeric SMILES, and deduplicates on the full **standard InChIKey**. The first occurrence supplies the retained state and name. All source records, including duplicates and rejects, remain in `record_audit.csv` and SQLite. Standard InChI has its own chemical normalization and tautomer conventions, so this is not an assertion of exact atom-graph equality. Salts and charge are not explicitly neutralized or stripped during this step.

CSV IDs are sequential, and job names include ligand ID and round. The job bundle hashes its source MSA JSON and manifest. The execution plan records hashes of each JSON and referenced MSA/template file. These detect changed inputs before a worker runs. Save dependency versions and AF3/model provenance alongside the campaign; matching input hashes alone does not certify identical model environments.

## MSA and sampling

A protein-only AF3 data pipeline computes unpaired MSA and templates once. Complex inputs reuse those data and explicitly set `pairedMsa` to an empty string. No second protein partner exists for a paired protein MSA in this workflow. Seeds are specified per input; the round name separates experiments.

Reusing MSA avoids redundant database searches. It does not remove feature preparation, JAX compilation or diffusion inference. Five diffusion samples from one seed are not five independent external validations. Different seeds can probe sampling variability, while the computational cost increases. Current JSON construction uses AF3 input version 2 and relative MSA/template paths; check compatibility with the installed AF3 release.

## Completeness and repeated attempts

The scanner reads sample directories `seed-N_sample-M`. Accepted summary filenames include bare or job-prefixed `summary_confidences.json`, with singular `summary_confidence.json` compatibility. Companion model and full confidence files must exist and be nonempty. Root best structures and their summary copies are excluded from the sample count.

Each retained seed group must contain exactly the configured sample indices, normally 0–4. Repeated complete attempts for the same ligand/round/seed are resolved by the lowest numeric attempt job ID, then path order. This approximates original submission order and avoids selecting the highest score across accidental reruns. It is not a timestamp-based rule across different clusters or recycled Slurm IDs. Originals are retained and exclusions recorded.

`coverage.csv` marks a ligand complete when **at least one seed group** is complete. It does not certify that every planned seed/round completed. The new execution worker validates all requested seeds before writing its completion record; historical multi-seed studies additionally need expected-seed reconciliation. A zero-recognized-file scan is labeled unresolved and exits with diagnostics rather than treating the entire library as failed predictions.

## Confidence ranking

The representative sample is chosen in this order:

1. Explicitly `has_clash=False`, followed by unknown clash status, then `True`.
2. Available A–B `chain_pair_iptm`, descending.
3. `ranking_score`, descending; stable path order breaks remaining ties.

The matrix chain order comes from `chain_ids` when supplied, or explicitly from `--chain-order A,B` for this two-chain input. Verify that fallback against the AF3 input and output; an incorrect chain order changes the meaning of pair metrics. Global `iptm`, `ptm`, per-seed counts, interface mean/SD and the maximum ranking-score model are also retained. The selected interface representative may differ from AF3's root best model.

These are prioritization heuristics. AF3 confidence describes model certainty, not binding free energy, ligand specificity or inhibition. A protein–protein ipTM threshold should not be treated as a universally calibrated small-molecule cutoff. See the [AF3 confidence output documentation](https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md).

## Structural inspection

Only the chosen candidate subset is inspected with full CIF/confidence data. Default settings are deliberately simple review aids:

| Quantity | Calculation | Interpretation limit |
| --- | --- | --- |
| Ligand pLDDT | Mean atom pLDDT for chain B | Local confidence, not affinity |
| Directional cross PAE | Mean/minimum over protein-to-ligand and ligand-to-protein blocks | Whole-protein block includes distant residues |
| Contacts | Protein–ligand heavy-atom distances ≤4 Å | No chemical typing of interactions |
| Very short contacts | Heavy-atom pairs <1.5 Å | Coarse flag, not a van der Waals clash model |
| Pose agreement | Ligand atom RMSD after protein pocket CA alignment, support ≤2 Å | No symmetry correction; atom names must agree |
| Residue frequencies | Fraction of inspected samples contacting a residue | Conditional on selected samples and candidates |

The reference pocket uses protein residues within 10 Å of the representative ligand. Kabsch fitting requires at least three noncollinear matching CA atoms. The ligand is never fitted to itself. Missing comparable atom sets or insufficient alignment points produce unavailable RMSD, not zero. The reference sample counts toward its own support fraction.

The CIF reader targets AF3 atom-site loops and is not a general-purpose crystallographic parser. `geometry_triage_pass` requires known absence of AF3 global clash, at least one 4 Å contact and no <1.5 Å contacts. It is not a comprehensive chemical-validity test. Exported Top 100 structures are review candidates, with original filenames and confidence files preserved.

## Docking comparison and diversity

AF3 and Glide rows are linked after `RemoveHs → FragmentParent → Uncharger → Canonicalize tautomer → canonical isomeric SMILES`. This follows the original docking clustering convention while preserving stereochemical information supported by the standardization. It can merge protonation/tautomer states; exact-state coordinate comparison needs explicit atom mapping and is not implemented here.

For each standardized parent, retain the lowest Glide score and the preferred AF3 nonclashing/interface representative. Keep all matched AF3 states in a separate table. Compare average-tie ranks with Spearman correlation within the matched, nonclashing parent cohort. An exploratory consensus selects both rankings in their top 10% of that cohort. Review representatives require geometry triage and are limited to one per existing docking-cluster/pocket combination.

If the docking input was already restricted to Top 10,000, correlation and consensus are conditional on that subset. Unmatched molecules are not negatives. Imported docking clusters do not cover new AF3-only chemotypes, and this repository does not recompute their chemical clusters. Report those candidates separately for further diversity analysis.

## Next scientific validation

Inspect chemical geometry, chirality, strain, plausible interactions and receptor-pocket placement; consider a tool such as [PoseBusters](https://github.com/maabuu/posebusters). For SMILES-defined ligands, recover bond orders from the original chemical record with explicit atom correspondence before exporting formats or performing energy calculations. Cross-method pose comparison should align the protein and handle ligand symmetry. Follow-up docking/rescoring, additional seeds and experimental binding/function measurements address questions that confidence ranking alone cannot resolve.
