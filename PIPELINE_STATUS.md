# Pipeline Status & Runbook

**Project:** Pharmacogenomic Variant Frequency Mapping Across Global Populations  
**Current Date:** 2026-10-12  
**Status:** Step 1 & Step 2 Complete (3,834 usable core relationships resolved with gnomAD). Ready for ClinVar (Step 3).

---

## 1. Directory Structure

```
pharmacogenomic-project/
├── data_raw/
│   ├── gnomad/
│   └── pharmgkb/
├── data_clean/
│   ├── archive/                       # Superseded/intermediate files from earlier passes
│   ├── core_variants.csv              # Step 1: Filtered PharmGKB core variant annotations
│   ├── rsid_list.txt                  # Step 1: 2,604 unique rsIDs in scope
│   ├── gnomad_candidates.csv          # Step 2: gnomAD GraphQL raw candidate variants
│   ├── site_allele_frequencies.csv    # Step 2: Allele-level counts (AC, AN, AF) for gnomAD populations
│   ├── core_rows_usable_v3.csv        # Canonical Step 2 output: 3,834 usable core relationships
│   ├── annotation_resolution_v3.csv   # Per-annotation resolution classification
│   └── manual_review_list_v3.csv      # 87 edge cases flagged for manual review
├── scripts/
│   ├── archive/                       # Archived intermediate and pilot scripts
│   ├── 01_build_core_list.py          # Step 1: Build core PharmGKB variant list
│   ├── 02d_fetch_all_candidates.py    # Step 2: Batch fetch gnomAD candidate variants & frequencies
│   ├── 03c_match_per_annotation.py    # Step 2: Per-annotation allele matching & site frequency table
│   └── 03f_fix_core_matching.py       # Step 2: Refined core relationship matching (v3)
└── PIPELINE_STATUS.md                 # This runbook
```

---

## 2. Pipeline Execution Order

| Step | Script | Input | Output | Status | Notes |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1. PharmGKB Core** | `01_build_core_list.py` | `data_raw/pharmgkb/` | `data_clean/core_variants.csv`, `rsid_list.txt` | **DONE** | 2,604 unique rsIDs, 4,016 core relationships |
| **2a. gnomAD Fetch** | `02d_fetch_all_candidates.py` | `rsid_list.txt` | `gnomad_candidates.csv` | **DONE** | 2,560 rsIDs matched (98.3%), 3,419 candidate variants fetched |
| **2b. Allele Alignment** | `03c_match_per_annotation.py` | `gnomad_candidates.csv`, `variantAnnotations/` | `site_allele_frequencies.csv`, `annotation_allele_matches.csv` | **DONE** | Resolved multiallelic and strand orientations |
| **2c. Core Usability** | `03f_fix_core_matching.py` | `core_variants.csv`, `annotation_resolution_v2.csv` | `core_rows_usable_v3.csv`, `manual_review_list_v3.csv` | **DONE** | **3,834 usable core rows (95.5%)**; 87 flagged for review |
| **3. ClinVar** | *Next script* | NCBI ClinVar `variant_summary.txt.gz` | `data_clean/clinvar_matched.csv` | **TODO** | Person B responsibility |
| **4. Merge** | *Upcoming* | `core_rows_usable_v3.csv`, `site_allele_frequencies.csv`, `clinvar_matched.csv` | `data_clean/merged_table.csv` | **TODO** | Person A responsibility |
| **5. Skew Calculation** | *Upcoming* | `merged_table.csv` | Analysis columns (max_AF, min_AF, ratio, 3x/5x/10x) | **TODO** | Person A responsibility |
| **6. 4-Group Classification** | *Upcoming* | `merged_table.csv` | `results/final_table.csv`, `results/group_counts.csv` | **TODO** | Person A responsibility |
| **7. Validation & Deliverables** | *Upcoming* | `final_table.csv` | `validation_log.csv`, Excel workbook, figures, report | **TODO** | Person B & C responsibilities |
