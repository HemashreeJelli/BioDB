# Pipeline Status & Runbook

**Project:** Pharmacogenomic Variant Frequency Mapping Across Global Populations  
**Last updated:** 2026-10-25  
**Status:** Steps 1–6 complete. Step 7 (validation & deliverables) to do.

---

## 1. Pipeline Execution Order

| Step | Script | Input | Output | Status | Notes |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1. PharmGKB Core** | `01_build_core_list.py` | `summaryAnnotations/summary_annotations.tsv` | `core_variants.csv`, `rsid_list.txt` | DONE | 2,604 rsIDs, 4,016 rows; duplicates keep the strongest evidence level |
| **2a. gnomAD Fetch** | `02_fetch_gnomad.py` | `rsid_list.txt` | `rsid_to_variant_ids.csv`, `gnomad_candidates.csv` | DONE | Resumable; re-run to fill gaps |
| **2b. Allele Alignment** | `03_match_alleles.py` | `gnomad_candidates.csv`, `rsid_to_variant_ids.csv`, `variantAnnotations/` | `site_allele_frequencies.csv`, `annotation_allele_matches.csv` | DONE | |
| **2c. Resolution v1** | `04_resolve_matches.py` | `annotation_allele_matches.csv` | `annotation_resolution.csv` | DONE | |
| **2d. Resolution v2** | `05_refine_strand_notes.py` | `annotation_resolution.csv` | `annotation_resolution_v2.csv` | DONE | |
| **2e. Core Usability** | `06_core_usability.py` | `core_variants.csv`, `annotation_resolution_v2.csv` | `core_rows_usable_v3.csv`, `annotation_resolution_v3.csv`, `manual_review_list_v3.csv` | DONE | 3,834 usable rows |
| **3. ClinVar** | `07_process_clinvar.py` | `variant_summary.txt.gz`, `rsid_list.txt` | `clinvar_matched.csv`, `clinvar_all_records.csv` | DONE | |
| **4–5. Merge & Skew** | `08_merge_and_skew.py` | step 2 & 3 outputs | `merged_table.csv`, `sanity_checks.txt` | DONE | Skew on minor allele; AF=0 with max AF ≥1% counts as high skew |
| **6. Classification** | `09_classify_and_report.py` | `merged_table.csv` | `final_table.csv`, `results/*.csv` | DONE | |
| **7. Validation & Deliverables** | — | `final_table.csv` | validation log, workbook, figures, report | TODO | |

All data files above live in `data_clean/`. Superseded outputs are kept in `data_clean/archive/`.

## 2. Change log

- **2026-10-25:** Fixed zero-AF variants being classed as low skew; skew now measured on the minor allele; ClinVar `has_pathogenic` no longer matches "Conflicting classifications of pathogenicity"; core dedup keeps the strongest evidence level; restored the two resolution steps as active steps. Previous outputs are in `data_clean/archive/pre_skewfix_2026-10-25/` and `data_clean/archive/pre_dedup_clinvar_fix_2026-10-25/`.
- **2026-10-25 (cleanup):** Scripts renamed into a 01–09 sequence with uniform headers; dead code and superseded side outputs removed; `run_all.py` added. Fixed a random tie-break in effect-allele selection (step 08) that made results change between runs. Full re-run passed all sanity checks. Pre-cleanup scripts and outputs are in `data_clean/archive/pre_cleanup_2026-10-25/`.
