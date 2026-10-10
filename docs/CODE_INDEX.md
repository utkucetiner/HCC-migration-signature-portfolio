# Analysis source index

Run scripts from the repository root, or use `python run_analysis.py SCRIPT --workspace WORKSPACE`. Runtime data and caches remain under `work/`; source is under `scripts/`. This separation preserves the original data paths.

Stages below describe dependencies; they are not an automatic complete pipeline. Script-specific literal paths are listed in `INPUT_PATHS.json`. Some preparation plans and external data are not distributed.

## 01 Input audit and scoring

- [audit_matrix.py](../scripts/audit_matrix.py)
- [expanded_controls.py](../scripts/expanded_controls.py)
- [qc_audit.py](../scripts/qc_audit.py)
- [random_controls.py](../scripts/random_controls.py)
- [score_descriptive.py](../scripts/score_descriptive.py)

## 02 Independent snRNA cohort

- [assay_limit_report.py](../scripts/assay_limit_report.py)
- [gse189175_cluster.py](../scripts/gse189175_cluster.py)
- [gse189175_deliver.py](../scripts/gse189175_deliver.py)
- [gse189175_score.py](../scripts/gse189175_score.py)
- [sharma_input_audit.py](../scripts/sharma_input_audit.py)

## 03 Covariate adjustment

- [balance_R_verify.py](../scripts/balance_R_verify.py)
- [balance_analysis.py](../scripts/balance_analysis.py)
- [balance_diagnostics.py](../scripts/balance_diagnostics.py)
- [balance_figures.py](../scripts/balance_figures.py)
- [build_balance_report.py](../scripts/build_balance_report.py)
- [common_support.py](../scripts/common_support.py)
- [parallel_models_report.py](../scripts/parallel_models_report.py)
- [quadratic_overlap.py](../scripts/quadratic_overlap.py)
- [stability.py](../scripts/stability.py) — historical exploratory cell bootstrap; not the final patient uncertainty method
- [verify_overlap.R](../scripts/verify_overlap.R) — R verification

## 04 Original RNA-CNV

- [build_gene_order.py](../scripts/build_gene_order.py)
- [classify_cnv.R](../scripts/classify_cnv.R)
- [cnv_all_figures.py](../scripts/cnv_all_figures.py)
- [cnv_consensus_fastpath.R](../scripts/cnv_consensus_fastpath.R)
- [cnv_preview.py](../scripts/cnv_preview.py) — HCC06 display preview; not full-cohort evidence
- [cnv_report.py](../scripts/cnv_report.py)
- [cnv_signature_sensitivity.py](../scripts/cnv_signature_sensitivity.py)
- [cnv_verify_final.py](../scripts/cnv_verify_final.py)
- [prepare_infercnv.py](../scripts/prepare_infercnv.py)
- [run_infercnv.R](../scripts/run_infercnv.R)
- [run_infercnv_hepref.R](../scripts/run_infercnv_hepref.R)
- [verify_candidate_overlap.R](../scripts/verify_candidate_overlap.R) — R verification
- [verify_candidate_overlap.py](../scripts/verify_candidate_overlap.py)
- [verify_cnv_fastpath.R](../scripts/verify_cnv_fastpath.R) — R verification

## 05 Signature-excluded RNA-CNV

- [classify_cnv_no_signature.R](../scripts/classify_cnv_no_signature.R)
- [cnv_signature_sensitivity_no_signature.py](../scripts/cnv_signature_sensitivity_no_signature.py)
- [compare_cnv_no_signature.py](../scripts/compare_cnv_no_signature.py)
- [drive_cnv_no_signature.py](../scripts/drive_cnv_no_signature.py)
- [finalize_cnv_no_signature.py](../scripts/finalize_cnv_no_signature.py)
- [plot_cnv_no_signature.py](../scripts/plot_cnv_no_signature.py)
- [prepare_cnv_no_signature.py](../scripts/prepare_cnv_no_signature.py)
- [report_cnv_no_signature.py](../scripts/report_cnv_no_signature.py) — interim reporting; use finalization for completed-cohort reports
- [run_infercnv_hepref_no_signature.R](../scripts/run_infercnv_hepref_no_signature.R)
- [run_infercnv_no_signature.R](../scripts/run_infercnv_no_signature.R)
- [verify_candidate_overlap_no_signature.R](../scripts/verify_candidate_overlap_no_signature.R) — R verification
- [verify_candidate_overlap_no_signature.py](../scripts/verify_candidate_overlap_no_signature.py)
- [verify_cnv_no_signature_final.R](../scripts/verify_cnv_no_signature_final.R) — R verification

## 06 Patient uncertainty

- [patient_uncertainty.py](../scripts/patient_uncertainty.py)
- [verify_patient_uncertainty.py](../scripts/verify_patient_uncertainty.py)
