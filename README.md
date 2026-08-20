# SQRCD Evidence Qualification Framework

This repository contains the public computational materials for the Evidence Qualification Framework (EQF) study. EQF evaluates whether evidence assembled for a traditional medicine formula is sufficiently specified for formula-level mechanistic interpretation. SQRCD is the case study, and fourteen independently specified formulae form a cross-formula portability benchmark. The manuscript is currently unpublished; this repository is a journal-independent code, data-product and reproducibility archive.

## Linked research records

- Time-stamped registration: [OSF n6gjy](https://osf.io/n6gjy/)
- Versioned research archive: [Zenodo DOI 10.5281/zenodo.21873862](https://doi.org/10.5281/zenodo.21873862)
- Public code repository: [Jeerryyoung/SQRCD_EQF](https://github.com/Jeerryyoung/SQRCD_EQF)

The three services have different roles. OSF preserves the registered analysis intent and amendments; GitHub provides version-controlled code and lightweight derived outputs; Zenodo archives immutable research snapshots cited by DOI. Future manuscript versions should cite a Git tag and the matching Zenodo version DOI, not a moving branch alone.

## Public contents

- `config/`: locked SQRCD identity, preparation metadata and analysis protocol.
- `external_benchmark/`: harmonized formula definitions, portability scripts, tests and summary outputs.
- `scripts/`: case-study, sensitivity and evidence-state reproducibility code.
- `results/`: de-identified aggregate outputs used in the manuscript and figures.
- `figures/`: final figure images and their directly corresponding summary tables.
- `journal_neutral_release/`: compact journal-independent release package, including the readable supplementary PDF, core machine-readable outputs ZIP, de-identified paired assessor ratings CSV, and data/code manifest.
- `docs/`: data-access, archive-linkage and release instructions.
- `provenance/`: release checksums and provenance metadata.

## Deliberately excluded

The public repository does not contain manuscript drafts, authors' correspondence, raw second-rater workbooks, evaluator initials or dates, credentials, local absolute paths, or third-party database exports whose redistribution terms are unclear. GEO data are referenced by accession and should be obtained from the original repository. Licensed or externally hosted pharmacology data must be retrieved from the named source under its current terms.

De-identified agreement summaries and the disagreement taxonomy are provided because they reproduce the reported framework-reproducibility results without releasing the assessors' source workbooks.

## Reproduction scope

The cross-formula portability module can be rerun with the scripts under `external_benchmark/scripts/`, subject to source availability and current database terms. The repository also provides the frozen derived summaries used for manuscript figures. Some case-study reconstruction steps depend on source exports that cannot be redistributed; those steps are represented by source identifiers, version records, hashes and derived audit tables rather than by copied database dumps.

## Environment

Python 3.11 or later is recommended.

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

## Citation

Use the metadata in `CITATION.cff`. For any future manuscript version, replace provisional software/archive identifiers with the final Git tag and matching Zenodo version DOI.

## License

Repository-authored code is released under GPL-3.0-only. Original derived tables and figures are released under CC BY 4.0 as described in `LICENSE-DATA`. These licenses do not override the terms of third-party databases or source publications.
