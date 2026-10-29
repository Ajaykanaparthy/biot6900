# Portfolio Deliverable 1: Verified LUAD Evidence from Modules 2-4

This directory brings together the existing LUAD analyses for the separate 5-6 page portfolio PDF. The three assignment notebooks were repaired where needed and executed in fresh kernels. Original coursework is preserved in each module's `archive/portfolio_deliverable_1_originals/` directory. The final portfolio PDF has not been created.

## Readiness and evidence

| Required evidence | Verified artifact |
|---|---|
| Multi-omics integration figure | [Module 2 integration figure](../module_2/Figures/luad_multiomics_integration.png), with caption in the [module README](../module_2/README.md) |
| Three-layer evidence overlap | [Three-set Venn](../module_2/Figures/luad_evidence_venn.png), [declared thresholds](../module_2/data/luad_evidence_thresholds.json), and [gene membership](../module_2/data/luad_evidence_set_membership.csv) |
| Joined genes and concordance | [6,758-gene integrated table](../module_2/data/luad_integrated_gene_level.csv) and [verification](../module_2/data/luad_portfolio_verification.json) |
| Tractability and final target ranking | [Module 3 ranking](../module_3/targets_luad_w3.csv), [validation report](../module_3/validation/validation_report.json), and the [cross-module summary](target_prioritization_summary.md) |
| AlphaFold confidence and pockets | [SCN7A figure](../module_4/Figures/scn7a_structure_confidence_pockets.png), [MASP1 figure](../module_4/Figures/masp1_structure_confidence_pockets.png), [captions](../module_4/Figures/captions.md), and [all pocket assessments](../module_4/results/pocket_quality_luad.csv) |
| Known LUAD genes and a discordant hit | [Known/discordant gene table](../module_2/data/luad_known_and_discordant_genes.csv), including TP53, KRAS, EGFR, and ZNF536 |
| Limitations and experimental validation | Candidate-specific limits in the [summary CSV](target_prioritization_summary.csv), the module READMEs, and the discussion below |

The listed computational evidence is ready to synthesize into the portfolio report. It supports hypothesis generation, not claims of experimentally established LUAD therapeutic targets.

The complete [file-change list](file_changes.csv) identifies modified originals and new artifacts. The [audit manifest](audit_manifest.json) records measured hashes, archived originals, execution checks, and preservation of unrelated coursework and Git state before publication.

## Independently reproduced results

### Module 2: multi-omics integration

The original five CPTAC-LUAD matrices from [LinkedOmics](https://www.linkedomics.org/data_download/CPTAC-pancan-LUAD/) reproduce **6,758 joined genes**, **4,918 concordant**, and **1,840 discordant**. The canonical 15-target CSV is byte-for-byte unchanged. Paired tumor-minus-normal RNA/protein effects use shared sample IDs; mutation frequencies use 108 tumors. The original score remains an equal mean of rank percentiles of absolute RNA effect, absolute protein effect, and mutation frequency.

The Venn represents positive evidence, rather than matrix availability. RNA/protein require layer-wide BH q <= 0.05, absolute paired log2 effect >= 1, and at least 10 valid pairs; genomic evidence requires mutation frequency >= 5%. The recurrence cutoff is descriptive, not a driver-significance test. Set sizes are RNA **1,443**, protein **371**, and genomic **383**; **26** genes meet all three thresholds and **4,990** meet none. Thresholds do not alter the original ranking.

One current MyGene mapping was missing relative to the completed analysis. The existing course annotation and a retained official Ensembl response establish **ENSG00000215695 -> RSC1A1**, restoring the original universe without inventing a symbol. Cached mappings remove this source of annotation drift. The preceding BRCA worked example remains course scaffolding, separate from the LUAD assignment.

TP53, KRAS, and EGFR retain ranks **572**, **1,624**, and **1,807**, with observed mutation frequencies **59.26%**, **32.41%**, and **36.11%**. Their lower integrated ranks reflect the chosen expression-and-recurrence score, not low biological importance. ZNF536 remains the notable rank-2 discordant candidate: RNA **-2.716**, protein **+2.712** log2 units, with protein evidence from **15** valid pairs. Its discordance has no demonstrated mechanism in this analysis.

### Module 3: protein language model triage

All 15 targets map uniquely. Removing SCN7A and GDF10 from the 1,201-protein panel leaves **1,199** training proteins. Five family folds have no family overlap between training and test sets. Independent calculations reproduce family-split **AUROC 0.726455** and **average precision 0.537090**; random-split values are 0.843971 and 0.724053. All 45 nearest-neighbor entries and four required anchor scores/percentiles were verified.

The original CSV regenerates exactly. Its existing rule is:

```text
final = 0.6 * percentile_rank(Module_2_score)
      + 0.4 * percentile_rank(max(p_sm, p_ab))
```

No new ranking or optimized weight was introduced. Course-provided embeddings and substitution scores were reused and audited; language-model inference was not repeated. The cross-disease APOE/TREM2 anchors remain assignment controls, while the target analysis remains LUAD.

### Module 4: structural druggability

| Measurement | SCN7A | MASP1 |
|---|---:|---:|
| Mean model pLDDT | 73.74 | 90.46 |
| Residues below pLDDT 50 | 14.74% | 2.72% |
| Experimental comparison | 7TJ8 | 4IGD |
| Confident-region RMSD | 5.35 Angstrom | 1.68 Angstrom |
| Best kept pocket | 98 | 1 |
| fpocket Score | -0.081 | 0.341 |
| fpocket Druggability Score | 0.982 | 0.687 |
| Pocket mean pLDDT | 85.06 | 70.91 |
| Pocket mean within-PAE | 3.88 Angstrom | 5.89 Angstrom |
| Professor's confidence verdict | keep | keep |

All **215** pocket descriptors and **5,105** pocket-lining atom entries were checked against the retained models; pLDDT, PAE, and verdicts were independently recomputed. Reject pockets with mean lining pLDDT <70; otherwise mark suspect if mean lining PAE >10 Angstrom; otherwise keep. SCN7A pocket 149 is shown among the top five but rejected despite its 0.926 druggability score. The legacy structural CSV is unchanged.

These are verified retained fpocket-format outputs, not a fresh fpocket run. Exact source-model coordinates support reuse; the original Explorer job log, command, and software version remain unavailable. A portable [452-file input ZIP](../module_4/reproducibility/luad_structure_inputs.zip) preserves the existing models, PAE, experimental structures, and pocket files without relying on current model downloads.

## Final prioritization and scientific limits

The [15-row summary CSV](target_prioritization_summary.csv) contains both original ranks, RNA/protein effects and support, mutation frequency, concordance, modality scores, and available structure/pocket evidence. Its [readable table](target_prioritization_summary.md) follows the existing Module 3 order. Structural values for the other 13 candidates are blank and explicitly unassessed.

- **COL6A6: Module 2 rank 3 -> Module 3 rank 1.** Concordant decreases and combined within-list percentiles support its existing priority, but its small-molecule and antibody scores are only 0.280420 and 0.030515. Its supplied embedding is truncated and no structural pocket assessment was performed. High rank does not establish that inhibiting a downregulated gene would help.
- **SCN7A: rank 5 -> rank 3.** Small-molecule score 0.644801 and a locally confident pocket support computational follow-up. Inter-component PAE **20.72 Angstrom** limits confidence in relative placement; the confidence-component heuristic does not establish biological domain boundaries. Pocket geometry and membrane-associated lipid/glycan contacts do not establish inhibitor binding or selectivity.
- **MASP1: rank 11 -> rank 8.** Scores 0.604960 and 0.549767 support follow-up for both modalities. Its best pocket only narrowly passes the mean-pLDDT rule. [4IGD](https://www.rcsb.org/structure/4IGD) is the **R448Q zymogen catalytic-region fragment, residues 298-699**, so its comparison does not validate full-length wild-type MASP1.

Strong multi-omics evidence measures dysregulation and recurrence. It does not establish tumor dependency, target accessibility, ligand affinity, or safety. Protein coverage varies: ZNF536 has 15 valid pairs, MYH2 12, and ADAMTS16 46, versus 101 for SCN7A/MASP1. The retained >=50-pair sensitivity analysis gives a different list; it does not silently replace the canonical ranking. COL6A6 and TNXB embeddings are truncated. Clinical-precedence classifier scores are not calibrated LUAD drug-success probabilities. Mutation recurrence may include passenger mutations and does not control for gene-specific background mutation rates.

Experimental validation remains future work: replicate the expression and discordance findings with orthogonal measurements in matched samples; test direction-appropriate perturbation and rescue in LUAD models with normal-tissue controls; and test pocket binding, functional modulation, selectivity, and toxicity. None of these experiments is claimed to have been performed.

## Execution and reproducibility

| Existing assignment notebook | Fresh-kernel code cells | Saved error outputs |
|---|---:|---:|
| [Module 2](../module_2/BIOT6900_Module2_Assignment2_LUAD.ipynb) | 17 | 0 |
| [Module 3](../module_3/BIOT6900_Module3_Assignment.ipynb) | 11 | 0 |
| [Module 4](../module_4/BIOT6900_Module4_Assignment.ipynb) | 11 | 0 |

Execution records are retained for [Module 2](../module_2/data/luad_notebook_execution.json), [Module 3](../module_3/validation/validation_report.json), and [Module 4](../module_4/results/notebook_validation.json). Original notebooks and report/CSV backups were preserved before modification. Both target CSVs and the original structural summary remain byte-identical.

Use the existing Python 3.11 `biot6900` environment; exact installed package versions are recorded in [requirements.txt](requirements.txt). Run each notebook from its own module directory using **Kernel -> Restart & Run All**. The three processed BRCA teaching tables are needed for the complete Module 2 notebook, in addition to its genuine LUAD source matrices. The structural notebook verifies/restores missing ignored inputs from its ZIP.

From the repository root, these existing-data commands reproduce the analysis artifacts and audits:

```powershell
conda activate biot6900
python -B module_2/recover_luad_portfolio.py
python -B module_3/scripts/validate_portfolio_module3.py
python -B module_4/scripts/execute_assignment.py
python -B portfolio_deliverable_1/scripts/build_summary.py
```

The Module 2 recovery command regenerates its tables/figures; use the notebook itself for full Module 2 execution. Other environments can differ in floating-point tie behavior; preserved rankings and versioned evidence provide the comparison reference. The summary builder refuses to replace different existing outputs.

The scoped root `.gitattributes` preserves audited file bytes, including the original CRLF source matrices and CSVs, so Git's platform-specific newline conversion cannot invalidate their checksums. It applies to the root README, Modules 2-4, and this portfolio directory; Modules 5-7 retain their existing attributes.

## Primary references for interpretation

- [TCGA LUAD, Nature 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4231481/): established LUAD genomic context.
- [Gillette et al., Cell 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7373300/): the CPTAC-LUAD proteogenomic study.
- [Ma et al., Frontiers in Oncology 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC7968342/): bioinformatic COL6A6 expression/prognostic associations, not experimental proof for these portfolio candidates.
- [RCSB 7TJ8](https://www.rcsb.org/structure/7TJ8) and [4IGD](https://www.rcsb.org/structure/4IGD): experimental structure identities and coverage.
