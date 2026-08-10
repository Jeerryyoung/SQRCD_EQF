# Reproducibility notes

This file summarizes the minimal reproducibility information for the PLOS ONE upload package.

## Software

- Python: 3.11 or later recommended
- Core packages: pandas, numpy, scipy, statsmodels, matplotlib, networkx
- Repository requirements: `requirements.txt`

## Main data products

- Complete SQRCD non-empty subformula enumeration: `subformula/subformula_retention_complete.csv`
- FDI and retention summaries: `subformula/formula_fdi_summary.csv`, `subformula/retention_summary_by_herb_count.csv`
- Exact Shapley attribution: `subformula/shapley_contributions.csv`
- Annotation perturbation and null summaries: `subformula/edge_thinning_summary.csv`, `subformula/degree_preserving_null_summary.csv`
- Cross-formula benchmark: `cross_formula/cross_formula_benchmark_summary.csv`
- De-identified paired assessment matrix: `S3 Data. De-identified paired assessor ratings.csv`

## Expected checks

- SQRCD complete subformula rows: 65,535
- Paired assessor rating rows: 750
- Formulae in portability benchmark: 15 total, including SQRCD and 14 external formulae
- Supporting Information PDF size: under 20 MB
- Each PLOS upload file size: under 20 MB

## Versioned links

- Zenodo: https://doi.org/10.5281/zenodo.21694013
- GitHub: https://github.com/Jeerryyoung/SQRCD_EQF
- OSF: https://osf.io/n6gjy/
