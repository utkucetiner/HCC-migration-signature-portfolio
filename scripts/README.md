# Analysis scripts

Browse the [stage-by-stage source index](../docs/CODE_INDEX.md). The 49 original scripts retain their filenames; two shared numerical modules and one R runtime helper support the refactor. Data and intermediate results are kept under an ignored `work/` directory. They are not mixed with source files.

Use the root launcher for one documented stage:

```sh
python run_analysis.py patient_uncertainty.py --workspace /path/to/local-workspace
python run_analysis.py run_infercnv_no_signature.R --workspace /path/to/local-workspace --patient HCC03
```

The second command is computationally expensive and requires prepared CNV inputs, inferCNV and JAGS. Neither command downloads missing data. See [reproduction scope](../docs/REPRODUCING.md) before executing it.

Historical exploratory scripts remain visible for provenance. `stability.py` includes cell-level bootstrap intervals; the final uncertainty analysis uses `patient_uncertainty.py`. `cnv_preview.py` is an HCC06 preview. `report_cnv_no_signature.py` is interim reporting; `finalize_cnv_no_signature.py` checks the completed cohort. CopyKAT is not implemented in this release.
