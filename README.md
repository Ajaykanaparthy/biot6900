# BIOT 6900 coursework
Name: Ajay Kanaparthy
This repository contains my BIOT 6900 Module 1 computational lab work.
The biot6900 conda environment is working with Python, Git, Jupyter, and BioPython. Parts A-C are completed. Part D will be completed after the assigned dataset/accession is provided in class
## Module 2 Assignment 2 — Independent Multi-Omics Target Discovery

**Disease:** Lung Adenocarcinoma (LUAD)

### Data Sources
- Transcriptomics: CPTAC pan-cancer LUAD RNA-seq — LinkedOmics — Open access
- Proteomics: CPTAC pan-cancer LUAD proteomics — LinkedOmics — Open access
- Genomics: CPTAC pan-cancer LUAD gene-level somatic mutation data — LinkedOmics — Open access

Source: https://www.linkedomics.org/data_download/CPTAC-pancan-LUAD/

### Results
- 6,758 genes survived the three-way join
- 4,918 genes were RNA-protein concordant
- Top 15 targets exported to `targets_luad.csv`

### Files
- `BIOT6900_Module2_Assignment2_LUAD.ipynb`
- `targets_luad.csv`
- `LUAD_Multi_Omics_Report.md`

## Module 3 — Computational Lab 3

Name: Ajay Kanaparthy
Disease working on: Lung Adenocarcinoma (LUAD)

### Results
- 15 LUAD targets successfully mapped
- 2 targets removed to prevent data leakage
- Family-split AUROC: 0.726
- Family-split AP: 0.537
- Re-ranked 15 targets for Week 4

### Files
- BIOT6900_Module3_Assignment.ipynb
- targets_luad_w3.csv
- Assignment 3- Target Triage with Protein Language Models.pdf

### Issues
Initial scikit-learn installation and file-path errors were resolved.


## Module 4 – Structure-Based Druggability

Name: Ajay Kanaparthy
Disease working on: Lung Adenocarcinoma (LUAD)

### Targets
- SCN7A (Q01118)
- MASP1 (P48740)

### Results
- SCN7A RMSD: 5.35 Å
- MASP1 RMSD: 1.68 Å
- Best pocket scores: SCN7A 0.982, MASP1 0.687
- Structure assessment completed

### Files
- BIOT6900_Module4_Assignment.ipynb
- structure_assessment_luad.csv
- BIOT6900_Module4_LUAD_Report.pdf

### Issues
fpocket was unavailable on Windows. I ran it on Northeastern Explorer and transferred the results back to VS Code.


## Module 5 – Variant Interpretation

**Name:** Ajay Kanaparthy  
**Gene:** CYP2C19

- Analyzed 80 CYP2C19 variants.
- Evaluated CADD, REVEL, and AlphaMissense.
- Applied allele frequency filtering.
- Calibration was unreliable due to insufficient labelled variants.
- All 33 VUS remained uncertain.
- AlphaGenome Atlas was skipped because no API key was configured.
- Completed discordance analysis and submitted the report.

**Files:**
- BIOT6900_Module5_Part3_Assignment.ipynb
- pgx_classified_CYP2C19.csv
- Module5_CYP2C19_Report.pdf
