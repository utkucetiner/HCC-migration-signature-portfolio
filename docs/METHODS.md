# Analysis scope and implementation

## Study units and signature

The main analysis used eight paired GSE149614 patients, HCC03–HCC10: 1,788 adjacent hepatocytes and 11,974 tumour hepatocytes. GSE189175 contributed three separately analysed paired patients and 25,977 hepatocyte-like nuclei. Cohorts were not pooled.

The source human migration list contained 50 genes. Forty-seven were mapped, including verified CYTOR/LINC00152 and PRAG1/SGK223 aliases. Three unmapped features were omitted rather than treated as zero expression. The score is the arithmetic mean across 47 genes of ln(1 + 10,000 × gene count / total feature counts).

## Adjustment and weighting

Five covariates were used: log1p total counts, continuous S and G2/M expression scores, hepatocyte identity and portal-minus-central zonation. These are expression scores, not Seurat cell-phase assignments. Within-patient regression gave equal total weight to the two tissues.

Unpenalized logistic propensity models excluded the outcome signature. Tumour weights were 1−e(X); adjacent weights were e(X), normalized separately by tissue. No one-to-one propensity matching was performed. Entropy balancing used a symmetric midpoint target with convex-hull feasibility checks; infeasible targets were not reported as successful fits. Quadratic propensity sensitivity added squared covariates without interactions. Effective sample sizes and weighted distribution diagnostics accompanied mean balance.

Conditioning on identity, cycle or zonation may remove biological disease variation as well as technical variation. The contrasts are descriptive conditional associations, not total or direct causal effects.

## RNA-CNV candidate selection

inferCNV 1.22.0 used GENCODE v38 genomic coordinates. Primary references were adjacent T/NK, B and endothelial cells, up to 100 each per patient, selected deterministically. Main settings included cutoff 0.1, minimum three expressing cells per gene, window 101, denoising, i6 HMM, Leiden subclusters and two threads. Additional Bayesian filtering was disabled with BayesMaxPNormal=0.

Candidate criteria required non-neutral i6 segments of at least 10 Mb on at least two autosomes, supported by a second adjacent-hepatocyte-reference workflow. The second workflow did not use HMM. It required same-direction segment signal exceeding the larger of 0.05 and the adjacent-reference absolute segment 99th percentile, at least ten genes per segment and at least 20 supported hepatocytes in the HMM group. These are research criteria, not a clinically validated malignancy test. Indeterminate cells were not called benign. Adjacent tissue was not DNA-verified benign.

For the exclusion sensitivity, all 47 signature genes were removed from raw CNV matrices and gene-order inputs. All eight patients were rerun using both reference workflows and the existing criteria. The original signature score and its normalization were retained. Jaccard agreement compared candidate sets before and after exclusion.

## Patient uncertainty

Each recorded within-patient contrast received equal patient weight. The patient-resampling distribution was enumerated exactly using multinomial count vectors; eight patients produce 6,435 vectors. Each vector has probability n!/(product of count factorials × n^n). Intervals use the weighted 2.5th and 97.5th percentiles.

Within-patient fits, labels and weights remained fixed during this resampling. Student-t intervals were additional sensitivity estimates under an unverified approximate-normality assumption. Bootstrap calculations were cross-checked by separate 200,000-replicate Monte Carlo calculations. Small patient samples limit coverage, particularly n=3.

Leave-one-patient-out analysis recomputed the mean of the remaining recorded patient contrasts. It is influence sensitivity, not prediction cross-validation. Entropy summaries of feasible subsets represent different patient sets from full-cohort summaries.

## Recorded runtime and source order

R 4.4.1; inferCNV 1.22; Seurat 5.4; JAGS 4.3.2. Python 3.11 with NumPy, SciPy, h5py, scikit-learn, igraph, leidenalg and Matplotlib. This records the working environment; it is not a tested installation recipe.

Suggested dependency order is input and metadata audit, gene mapping and scoring, expanded controls, balance diagnostics, CNV preparation, both reference runs, candidate classification, candidate contrasts, R verification and final reporting. The separate snRNA clustering and scoring precede its models. Signature-excluded CNV reruns precede the new candidate uncertainty summaries.

The refactored source contains 49 original scripts plus shared numerical and runtime helpers. See [CODE_INDEX.md](CODE_INDEX.md) for stages, [REPRODUCING.md](REPRODUCING.md) for environment and input requirements, and [VALIDATION.md](VALIDATION.md) for the clean Python checks. A complete raw-data-to-CNV rerun remains pending.
