# Computational validation — 10 October 2026

A new Python 3.11.0 virtual environment was created with system-site packages disabled. All eight direct dependencies were installed afresh from the package registry; transitive versions are recorded in `environment/requirements-lock.txt`. `pip check` reported no broken requirements. This was a clean Python environment on the same macOS arm64 host, not a container or a separate machine.

| Check | Observed outcome | Scope |
|---|---|---|
| Independent synthetic numerical tests | Six passed | Patient resampling, influence and weighting functions |
| Original patient-uncertainty stage | All 13 summary rows; complete result JSON matched exactly | Recomputed from four local aggregate-result inputs |
| Separate Monte Carlo and LOO verification | Passed | 200,000 resampling draws and independent three-patient enumeration |
| Real cohort covariate weighting | Eight scRNA and three snRNA patient fits completed; maximum raw/overlap/entropy contrast error 0.0 | Recomputed from four existing private NPZ caches |
| Python syntax | All public Python files parsed | Syntax only for unexecuted stages |
| R syntax | 13 R source files parsed | Does not load inferCNV or JAGS |
| Independent base-R HCC03 overlap check | glm converged; contrast 0.0100092941486234 | Fresh vanilla R session using locally prepared model input |

The private NPZ and model inputs were copied to a separate validation workspace only. No cells, barcodes, raw counts or large CNV objects were uploaded by this update. Aggregate comparisons check refactoring consistency; they do not provide independent biological validation.

The original eight-patient CNV runs and the previously recorded R/Python candidate checks are historical analysis results. They were not rerun in this new environment. A full rebuild from downloaded counts, clustering and CNV inference remains unverified. The three-patient nuclear cohort retains its original statistical limitations.

Machine-readable outcome: [VALIDATION.json](VALIDATION.json). Source provenance and original/refactored checksums: [SOURCE_MANIFEST.json](SOURCE_MANIFEST.json).
