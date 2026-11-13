# BIOT 6900 — AI-Powered Biotech Workflows

**Student:** Ajay Kanaparthy  
**Program:** M.S. Biotechnology, Northeastern University  
**Primary disease focus:** Lung adenocarcinoma (LUAD); Module 5 separately covers CYP2C19 pharmacogenomics.

This repository contains my BIOT 6900 coursework, organized by module, including notebooks, supporting data, reports, and reproducibility records.

[Portfolio Deliverable 1: audited LUAD evidence from Modules 2-4](portfolio_deliverable_1/README.md) brings together verified figures, target prioritization, reproducibility records, and scientific limitations for the separate portfolio report.

## Module guide

| Assignment | Topic | Main contents |
| --- | --- | --- |
| [Module 1 notebook](module1_setup.ipynb) | Environment setup | Setup notebook at repository root |
| [Module 2](module_2/) | LUAD multi-omics discovery | CPTAC data, target shortlist, figures, report and README |
| [Module 3](module_3/) | Protein language models | LUAD target triage, model validation, rankings and README |
| [Module 4](module_4/) | Structure-based druggability | SCN7A/MASP1 analysis, pocket figures and README |
| [Module 5](module_5/) | CYP2C19 variant interpretation | Completed assignment, classified variants and README |
| [Module 6](module_6/) | LUAD polygenic risk scores | Real genotypes, simulated outcomes, figures and report |
| [Module 7](Module_7/) | LUAD knowledge graphs | Graded Neo4j notebook, exported graph, biomarker rankings and report |
| [Module 9](Module_9/) | LUAD clinical NLP | Completed graded notebook, D1–D4 answers and leakage analysis |
| [Portfolio Deliverable 1](portfolio_deliverable_1/) | Modules 2–4 synthesis | Verified evidence and prioritization audit; final PDF not in repository |

The original Module 1 notebook is stored at the repository root, not in a `module_1/` directory. Historical working copies for Modules 2–4 are preserved in their respective archive folders.

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

Notebook: [module1_setup.ipynb](module1_setup.ipynb)

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
- 4,918 genes were RNA-protein concordant; 1,840 were discordant
- Top 15 targets exported to `targets_luad.csv`, with the original ranking preserved
- Thresholded RNA/protein/genomic evidence overlap: 26 genes meet all three declared thresholds

### Files

- [LUAD assignment notebook](module_2/BIOT6900_Module2_Assignment2_LUAD.ipynb)
- [Starter notebook](module_2/BIOT6900_Module2_Starter.ipynb)
- [Real CPTAC data loader](module_2/BIOT6900_Module2_Part2_LoadRealCPTAC.ipynb)
- [LUAD data preparation](module_2/LUAD_Prepare_Data.ipynb)
- [targets_luad.csv](module_2/targets_luad.csv)
- [LUAD multi-omics report](module_2/LUAD_Multi_Omics_Report.md)
- [Module 2 data](module_2/data/)
- [Verified methods, figure captions, and run instructions](module_2/README.md)
- [Multi-omics integration figure](module_2/Figures/luad_multiomics_integration.png) and [three-set Venn](module_2/Figures/luad_evidence_venn.png)
- [Integrated 6,758-gene table](module_2/data/luad_integrated_gene_level.csv)

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
- [Verified methods and limitations](module_3/README.md)
- [Independent mapping, leakage, family-split, anchor, and ranking audit](module_3/validation/validation_report.json)

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
- Best kept-pocket fpocket Druggability Scores: SCN7A 0.982, MASP1 0.687
- Structure assessment completed using verified retained pocket outputs; these scores are computational evidence

### Files

- [Assignment notebook](module_4/BIOT6900_Module4_Assignment.ipynb)
- [structure_assessment_luad.csv](module_4/structure_assessment_luad.csv)
- [Verified structural methods, confidence rules, and limitations](module_4/README.md)
- [SCN7A structural figure](module_4/Figures/scn7a_structure_confidence_pockets.png) and [MASP1 structural figure](module_4/Figures/masp1_structure_confidence_pockets.png)
- [All 215 pocket assessments](module_4/results/pocket_quality_luad.csv)
- [Portable retained structural inputs](module_4/reproducibility/luad_structure_inputs.zip)
- [Assignment instructions](module_4/Assignment%204-b%20Structure-Based%20Druggability%20Assessment.pdf)
- [Preserved structures and fpocket inputs](module_4/reproducibility/luad_structure_inputs.zip)
- [PyMOL script](module_4/demo_pymol.pml) and [ChimeraX script](module_4/demo_chimerax.cxc)
- `BIOT6900_Module4_LUAD_Report.pdf` is referenced in the coursework notes but is not present in this repository.

### Issues

fpocket was unavailable on Windows. I ran it on Northeastern Explorer and transferred the results back to VS Code. The portfolio audit verifies the retained pockets against the exact source-model atoms and recomputes pocket confidence; it does not claim a fresh fpocket run. The historical job log, command, and version are unavailable.


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
- [pgx_classified_CYP2C19.csv](module_5/pgx_classified_CYP2C19.csv)
- [Assignment instructions](module_5/Module%205%20Assignment-%20Variant%20Interpretation%20of%20CYP2C19.pdf)
- [Module 5 README](module_5/README.md)
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
- [Final LUAD PRS report](module_6/Module%206-%20Polygenic%20Risk%20Score%20Portability%20in%20Lung%20Adenocarcinoma%20%28LUAD%29.pdf)
- [PRS distribution figure](module_6/Figures/prs_distribution.png) and [calibration figure](module_6/Figures/calibration.png)

## Module 7 — LUAD knowledge graphs

I built and queried a Neo4j graph joining LUAD targets to pathways, phenotypes, variants, related diseases and drugs. The graded notebook contains two-hop drug and shared-phenotype queries, biomarker prioritization, and community detection. The analysis summary lists **ACTN2, COL6A6, SCN7A, TNR, and ZNF536** as its top five combined candidates.

- [Completed graded notebook](Module_7/BIOT6900_Module7_Starter.ipynb)
- [Graph visualization](Module_7/my_knowledge_graph.png) and [focused paths](Module_7/luad_focused_paths.png)
- [Nodes](Module_7/nodes.csv), [relationships](Module_7/edges.csv), and [ranked biomarkers](Module_7/biomarkers.csv)
- [Four-page interpretation report](Module_7/Module7_LUAD_Knowledge_Graph_Report.pdf) and [README](Module_7/README.md)

**Consistency issue to resolve:** The module README and `analysis_summary.json` say **27 nodes / 30 relationships**, but the committed CSV exports contain **29 nodes / 29 relationships**. These are conflicting saved versions and must be reconciled before treating the counts as final.

## Module 9 — Clinical NLP on LUAD

The completed graded notebook generates **400 synthetic patient records**, extracts clinical concepts, and identifies outcome leakage. For post-obstructive pneumonia, naive keyword matching counted **155** patients versus **71** with assertion-aware detection.

The cross-validated AUROC was **0.939** with all features, **0.839** after removing declared leakage, and **0.769** after removing all identified leakage. The hidden leak was `renal_replacement`, while `lactate_not_clearing` was retained as an early signal. Answers to **D1–D4** are included in the graded notebook.

- [Completed graded notebook](Module_9/BIOT6900_Module9_Notebook3_Graded.ipynb)
- [Module 9 ReadME.md](Module_9/ReadME.md)

**Limitation and outstanding documentation:** These are simulated records and outcomes, not clinical validation. As of this update, `Module_9/ReadME.md` is on GitHub but **empty (0 bytes)**; it must be populated with the disease, data sources, work and limitations.

## Portfolio Deliverable 1 — Modules 2–4 synthesis

The [portfolio evidence folder](portfolio_deliverable_1/) combines the audits and supporting analyses from Modules 2–4, including a [prioritization summary](portfolio_deliverable_1/target_prioritization_summary.csv) and [methods documentation](portfolio_deliverable_1/README.md). The separate required **5–6 page final portfolio PDF is not present** in the committed folder; these support files do not replace it.

## Additional Data Sources and Limitations

### Module 3 – Protein Language Models

**Data sources:** My Module 2 LUAD target list (`targets_luad.csv`) and the course-provided protein training panel, human proteome embeddings, clinical target labels, and anchor variant dataset.

**Work completed:** Mapped all 15 targets, trained small-molecule and antibody tractability models, validated with a family-aware split, and re-ranked the LUAD targets.

**Limitations:** The model estimates similarity to previously drugged proteins. It does not prove drug binding. The family-split AUROC was 0.726 and average precision was 0.537.

**Issues:** Initial installation and file-path problems were resolved.

### Module 4 – Structure-Based Druggability

**Data sources:** AlphaFold Protein Structure Database (SCN7A Q01118 and MASP1 P48740), Protein Data Bank (7TJ8 and 4IGD), and the course p53 structural reference.

**Work completed:** Analyzed structure confidence, compared predictions with experimental structures, identified binding pockets using fpocket, and assessed potential druggability.

**Limitations:** Predicted structures do not fully represent protein movement, membrane environments, or interactions with binding partners. SCN7A inter-component PAE limits relative placement, and MASP1 structure 4IGD covers an engineered catalytic-region fragment rather than the full wild-type protein. Kept pockets are not experimentally validated binding sites.

**Issues:** fpocket was unavailable on Windows. I ran it using Northeastern's Explorer cluster and analyzed the results locally.

### Module 5 – Variant Interpretation

**Data sources:** MyVariant.info, ClinVar annotations, CADD, REVEL, AlphaMissense, and gnomAD allele-frequency data.

**Work completed:** Analyzed 80 CYP2C19 variants, evaluated predictor coverage, applied frequency filtering, assessed calibration, identified discordant variants, and exported provisional classifications.

**Limitations:** Only five benign/likely benign variants and no pathogenic/likely pathogenic variants were available for calibration. The evidence was insufficient to reclassify VUS.

**Issues:** AlphaGenome Atlas was skipped because an API key was not configured. No Atlas results were claimed.
