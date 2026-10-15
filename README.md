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
Disease: Lung Adenocarcinoma (LUAD)

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