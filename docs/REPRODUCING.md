# Reproduction scope and execution

## What can be checked from this repository alone

Use Python 3.11. The system Python 3.9 installation was unsuitable for the pinned scientific packages. The independent environment used Python 3.11.0 on macOS arm64.

```sh
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r environment/requirements-lock.txt
python -m pip check
python check_environment.py
python -m unittest discover -s tests -v
```

Six tests exercise exact patient resampling, equal-patient influence, invalid inputs, known weighting contrasts and outcome-shift invariance. They use synthetic data. Passing them does not rerun the biological study.

`requirements.txt` pins direct Python dependencies. `requirements-lock.txt` records every installed Python package from the successful clean-environment installation. It is a version snapshot for the tested platform, without wheel hashes or an operating-system lock. Other platforms remain untested.

## Source, data and output paths

Source files are under `scripts/`. Runtime inputs and caches use `work/`; the original main count file uses `data/external/`. Generated reports use `outputs/`. These directories are ignored by Git.

The launcher chooses the current Python interpreter, finds Rscript through PATH (or `HCC_RSCRIPT`), sets the source root and selects a separate working-data directory. Example:

```sh
python run_analysis.py patient_uncertainty.py --workspace /path/to/local-workspace
```

That workspace must contain `work/expanded_control_results.json`, `work/balance_results.json`, `work/cnv_signature_results.json` and `work/GSE189175_score_results.json`. The stage recomputes summaries from fixed recorded patient contrasts. It does not refit cell annotations, CNV labels or patient weights.

For covariate weighting, `balance_analysis.py` requires four caches: `hepatocyte_scoring_inputs.npz`, `expanded_control_scores.npz`, `GSE189175_cluster_inputs.npz` and `GSE189175_scores.npz`, all under `work/`. The independent snRNA inputs must be prepared before this joint reporting script runs; calculations remain cohort-specific.

[INPUT_PATHS.json](INPUT_PATHS.json) inventories literal data, cache and output paths in source. It includes both reads and writes, plus path prefixes; it is an audit aid, not a complete machine-readable dependency graph.

## Biological input preparation

- Main data: [GSE149614](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614), with the count file at `data/external/GSE149614_HCC.scRNAseq.S71915.count.txt.gz` and metadata at `work/data/GSE149614_HCC.metadata.updated.txt.gz`.
- Independent nuclear cohort: [GSE189175](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189175), with matrix, feature, barcode and metadata filenames expected by `gse189175_cluster.py`.
- Genomic ordering: [GENCODE v38](https://www.gencodegenes.org/human/release_38.html), expected at `work/data/gencode.v38.annotation.gtf.gz`.
- Expanded controls require `work/data/s_genes.txt` and `g2m_genes.txt`. Scoring and annotation stages also depend on preparation plans such as `expanded_control_plan.json` and `GSE189175_analysis_plan.json`.

These sources, preparation plans and caches are not all distributed here. Automated retrieval and complete public reconstruction of preparation inputs remain unfinished. Consult the [code index](CODE_INDEX.md) and original publications. Do not invent missing metadata or substitute a different gene list silently.

## R and CNV

Recorded biological-analysis runtime: R 4.4.1, inferCNV 1.22.0, JAGS 4.3.2; Matrix 1.7.0 was available in the original runtime. The R environment has not been reconstructed from an empty package library. Follow [official inferCNV documentation](https://github.com/broadinstitute/infercnv/wiki/Running-InferCNV) for its dependencies.

`HCC_R_LIBRARY` optionally selects a local R package library; `HCC_JAGS_MODULE_DIR` optionally selects installed JAGS modules. No machine-specific JAGS path is embedded in the updated source. Install these dependencies separately before running CNV stages.

Original CNV preparation precedes both reference workflows, candidate classification and candidate contrasts. Signature exclusion uses `prepare_cnv_no_signature.py`, then both exclusion reference workflows and classification. The preparation script no longer generates or overwrites source files; its three R entry points are tracked in Git. The driver skips completed stages by their markers. A completed marker must not be reused after changing inputs without explicitly preparing a fresh workspace.

## What the clean-environment check actually did

See [VALIDATION.md](VALIDATION.md). Python dependencies were freshly installed, six synthetic tests passed, 13 patient-summary rows matched exactly, and 11 real patient weighting results were recomputed from local private caches. Full raw-count processing, clustering, all 16 inferCNV objects, candidate relabelling and CopyKAT were not rerun. This is partial computational reproduction, not end-to-end reproduction from public inputs.
