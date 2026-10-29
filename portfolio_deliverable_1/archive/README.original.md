# BIOT 6900 coursework

Name: Ajay Kanaparthy
This repository contains my BIOT 6900 coursework, organized by module.

## Module guide

| Folder | Topic | Contents |
| --- | --- | --- |
| [module_1](module_1/) | Environment setup | Setup and scratch notebooks |
| [module_2](module_2/) | Multi-omics target discovery | CPTAC/LUAD notebooks, raw downloads, `data/`, targets, and report |
| [module_3](module_3/) | Protein language models | Toy and completed LUAD assignment notebooks, assignment PDF, `data/`, and re-ranked targets |
| [module_4](module_4/) | Structure-based druggability | Lab and assignment notebooks, assignment PDF, `structures/`, fpocket results, and visualization scripts |
| [module_5](module_5/) | Variant interpretation | Toy, real-data, and assignment notebooks, assignment PDF, CYP2C19 data and classified variants |
| [module_6](module_6/README.md) | Polygenic Risk Score Portability (LUAD) | PGS003393, real 1000 Genomes genotypes, harmonized PRS analysis, figures, and disclosed simulated outcomes |
| [Module_7](Module_7/) | Knowledge graphs | Neo4j real-model notebook |

Each module keeps its own notebooks and supporting files. `module_1/`, `module_2/`, and `module_3/` also contain `archive/root_copies/`, preserving the extra files previously stored at the repository root. These include an unfinished Module 3 template; the completed LUAD notebook remains directly in `module_3/`. No original copies were discarded.

## Running the notebooks

Activate the `biot6900` environment and run each notebook with its module folder as the working directory. For example, from the repository root:

```powershell
conda activate biot6900
cd module_2
jupyter notebook
```

When using an IDE, set the notebook's working directory to its module folder so relative `data/` and `structures/` paths resolve correctly. Module 3 reads its input directly from `../module_2/targets_luad.csv` and writes the re-ranked CSV into `module_3/`.

For the Module 4 PyMOL script, use `module_4/` as the working directory. The ChimeraX script already points to the structures in that folder. Module 6's graded LUAD notebook reads the real genotype and published score inputs from `module_6/data/` with `USE_DEMO = False`; ages and LUAD outcomes are simulated and disclosed in the Module 6 README.

## Module 1 - Environment setup

Notebook: [module1_setup.ipynb](module_1/module1_setup.ipynb)

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

- [LUAD assignment notebook](module_2/BIOT6900_Module2_Assignment2_LUAD.ipynb)
- [Starter notebook](module_2/BIOT6900_Module2_Starter.ipynb)
- [Real CPTAC data loader](module_2/BIOT6900_Module2_Part2_LoadRealCPTAC.ipynb)
- [LUAD data preparation](module_2/LUAD_Prepare_Data.ipynb)
- [targets_luad.csv](module_2/targets_luad.csv)
- [LUAD multi-omics report](module_2/LUAD_Multi_Omics_Report.md)
- [Module 2 data](module_2/data/)

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

- [Completed LUAD assignment notebook](module_3/BIOT6900_Module3_Assignment.ipynb)
- [Toy protein language model](module_3/BIOT6900_Module3_Part1_ToyModel.ipynb)
- [targets_luad_w3.csv](module_3/targets_luad_w3.csv)
- [Assignment instructions](module_3/Assignment%203-%20Target%20Triage%20with%20Protein%20Language%20Models.pdf)
- [Module 3 data pack](module_3/data/)

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

- [Assignment notebook](module_4/BIOT6900_Module4_Assignment.ipynb)
- [structure_assessment_luad.csv](module_4/structure_assessment_luad.csv)
- [Assignment instructions](module_4/Assignment%204-b%20Structure-Based%20Druggability%20Assessment.pdf)
- [Structure and fpocket files](module_4/structures/)
- [PyMOL script](module_4/demo_pymol.pml) and [ChimeraX script](module_4/demo_chimerax.cxc)
- `BIOT6900_Module4_LUAD_Report.pdf` is referenced in the coursework notes but is not present in this repository.

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

- [Assignment notebook](module_5/BIOT6900_Module5_Part3_Assignment.ipynb)
- [Real-data notebook](module_5/BIOT6900_Module5_Part2_RealData.ipynb)
- [pgx_classified_CYP2C19.csv](module_5/pgx_classified_CYP2C19.csv)
- [Assignment instructions](module_5/Module%205%20Assignment-%20Variant%20Interpretation%20of%20CYP2C19.pdf)
- [Module 5 data](module_5/data/)
- `Module5_CYP2C19_Report.pdf` is referenced in the coursework notes but is not present in this repository.

## Module 6 — Polygenic Risk Score Portability (LUAD)

- Published **PGS003393 (best_LUAD)** lung adenocarcinoma score.
- Real **1000 Genomes Phase 3 genotypes from 2,504 participants**, using GRCh37 coordinates.
- **64 harmonized variants**; ten strand-ambiguous variants intentionally excluded.
- Ancestry-specific PRS distributions, thresholds, AUC, and calibration.
- **Simulated ages and LUAD outcomes, clearly disclosed**; these teaching outcomes cannot validate clinical LUAD prediction or real-world ancestry portability.
- [Module 6 analysis, results, and run instructions](module_6/README.md).

### Files

- [Original graded notebook](module_6/BIOT6900_Module6_Part3_Graded_Starter.ipynb)
- [Toy PRS notebook](module_6/BIOT6900_Module6_Toy_PRS.ipynb)
- [Real PRS notebook](module_6/BIOT6900_Module6_Part2_RealPRS.ipynb)

## Module 7 - Knowledge graphs

- [Neo4j real-model notebook](Module_7/BIOT6900_Module7_RealModel_Neo4j.ipynb)

## Additional Data Sources and Limitations

### Module 3 – Protein Language Models

**Data sources:** My Module 2 LUAD target list (`targets_luad.csv`) and the course-provided protein training panel, human proteome embeddings, clinical target labels, and anchor variant dataset.

**Work completed:** Mapped all 15 targets, trained small-molecule and antibody tractability models, validated with a family-aware split, and re-ranked the LUAD targets.

**Limitations:** The model estimates similarity to previously drugged proteins. It does not prove drug binding. The family-split AUROC was 0.726 and average precision was 0.537.

**Issues:** Initial installation and file-path problems were resolved.

### Module 4 – Structure-Based Druggability

**Data sources:** AlphaFold Protein Structure Database (SCN7A Q01118 and MASP1 P48740), Protein Data Bank (7TJ8 and 4IGD), and the course p53 structural reference.

**Work completed:** Analyzed structure confidence, compared predictions with experimental structures, identified binding pockets using fpocket, and assessed potential druggability.

**Limitations:** Predicted structures do not fully represent protein movement, membrane environments, or interactions with binding partners.

**Issues:** fpocket was unavailable on Windows. I ran it using Northeastern's Explorer cluster and analyzed the results locally.

### Module 5 – Variant Interpretation

**Data sources:** MyVariant.info, ClinVar annotations, CADD, REVEL, AlphaMissense, and gnomAD allele-frequency data.

**Work completed:** Analyzed 80 CYP2C19 variants, evaluated predictor coverage, applied frequency filtering, assessed calibration, identified discordant variants, and exported provisional classifications.

**Limitations:** Only five benign/likely benign variants and no pathogenic/likely pathogenic variants were available for calibration. The evidence was insufficient to reclassify VUS.

**Issues:** AlphaGenome Atlas was skipped because an API key was not configured. No Atlas results were claimed.
