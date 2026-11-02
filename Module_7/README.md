# BIOT 6900 Module 7 — LUAD

Ajay Kanaparthy. The professor's original `BIOT6900_Module7_Starter.ipynb` contains the completed analysis; the reference notebooks are not submission replacements.

## Outputs

- `BIOT6900_Module7_Starter.ipynb`
- `my_knowledge_graph.png`
- `luad_focused_paths.png`
- `nodes.csv`
- `edges.csv`
- `biomarkers.csv`
- `communities.csv`
- `query_evidence.csv`
- `analysis_summary.json`
- `sources.csv`
- `example_cleanup_audit.json`
- `Module7_LUAD_Knowledge_Graph_Report.pdf`
- `README.md`

The report has 4 pages including references. It uses the executed graph, query results, candidate rankings and community output.

## Graph and actual results

The assignment graph contains **27 nodes and 30 relationships**. Counts describe this analysis rather than unrelated data already in AuraDB.

| Node type | Count |
|---|---:|
| disease | 2 |
| drug | 2 |
| gene | 17 |
| pathway | 3 |
| phenotype | 2 |
| variant | 1 |

| Relationship type | Count |
|---|---:|
| ASSOCIATED_WITH | 17 |
| HAS_PHENOTYPE | 4 |
| OCCURS_IN | 1 |
| PARTICIPATES_IN | 6 |
| TARGETS | 2 |

The original 15 candidates are ranked separately from contextual EGFR and TOP1 genes. The top five combined priorities are **ACTN2, COL6A6, SCN7A, TNR, ZNF536**. Directed PageRank tied across all 15 cohort candidates: **True**. Priority is `0.4 * centrality + 0.4 * within-list omics percentile + 0.2 * within-list max(p_sm, p_ab) percentile`; a tied centrality component is set equally to 0.5. Both tractability priors remain separate in `biomarkers.csv`. Weights are descriptive teaching choices, not clinically optimized parameters.

The two-hop query returned: Osimertinib targets EGFR (FDA-approved in specified EGFR-mutant NSCLC settings; 2024 label).

One actual shared-phenotype path has five nodes and four relationships (the final edge is traversed in reverse):

- `Topotecan` —`TARGETS`→ `TOP1`
- `TOP1` —`ASSOCIATED_WITH`→ `Small-cell lung carcinoma`
- `Small-cell lung carcinoma` —`HAS_PHENOTYPE`→ `Cough`
- `Lung adenocarcinoma (LUAD)` —`HAS_PHENOTYPE`→ `Cough`

This topotecan path suggests a testable hypothesis; shared symptoms do not establish LUAD treatment efficacy. TOP1-to-SCLC represents therapeutic-target context. Greedy modularity found **5 communities**. The interpreted community contains: COL11A1; COL6A6; Extracellular matrix organization; TNR; TNXB. Full membership appears in `communities.csv`.

## Run the original notebook

1. Select the Python `biot6900` kernel with the required packages installed.
2. Store Aura credentials locally in `Module_7/.env`, using `NEO4J_URI`, `NEO4J_USERNAME` and `NEO4J_PASSWORD`; quoted values are supported. This file is git-ignored. Environment variables override local configuration. A hidden password prompt is used only if the URI and username are available and the password is missing. Never publish credentials or paste them into outputs.
3. Use the Aura home database by default. Set optional `NEO4J_DATABASE` only when its actual name is known; do not assume the database is named `neo4j`.
4. Run **Kernel → Restart & Run All**. The connection cell verifies connectivity before database work; missing or invalid settings stop execution.
5. Inspect the query tables, figure, CSVs and final four-page PDF. Output paths resolve to `Module_7` when launched from the repository or module directory.

Before `MERGE`, the notebook checks the exact original 13-node/13-relationship teaching fixture. It removes that fixture only when all labels, properties, directed relationships and isolation match. Any mismatch stops and rolls back. Verified relationships are deleted explicitly before plain node deletion; there is no database-wide reset or broad deletion. `example_cleanup_audit.json` records the removal and preservation checks. All non-example data are preserved, and later runs skip cleanup when the fixture is absent.

## Existing data and limitations

- Module 2: `module_2/targets_luad.csv` and `module_2/data/luad_integrated_gene_level.csv` provide measured molecular associations and harmonized identifiers. GWAS p-values are unavailable and remain missing.
- Module 3: `module_3/targets_luad_w3.csv` provides separate model-derived `p_sm` and `p_ab` tractability priors and the previous ranking. These are not calibrated probabilities of LUAD treatment success; COL6A6 and TNXB have truncated embeddings.
- Module 4: `module_4/results/verified_structure_summary.csv` provides SCN7A/MASP1 structural context. Predicted pockets do not prove binding or functional LUAD dependence.
- Module 5: CYP2C19 variants concern pharmacogenomics and are not labeled as LUAD driver mutations. The graph's EGFR L858R is independently sourced external context, not an asserted cohort variant.
- Module 6: Genotypes are real, while ages and disease labels are simulated. The optional patient layer is empty.
- Curated pathway membership and observational abundance associations do not establish causality. Shared cough/dyspnea are nonspecific; PageRank and communities depend on graph annotation and edge directions. Qualitative confidence documents support for individual edges rather than clinical certainty.
- AI assisted source review, identifier harmonization, coding and report drafting. Results were checked against original data, official sources, executed Neo4j queries and exported relationships.

## Sources

Every relationship records `source_url`, confidence and its rationale in `edges.csv`; the source catalog appears in `sources.csv`. Confidence criteria: high: Direct curated identifier/pathway assertion or regulatory drug mechanism/indication.; low: An indirect hypothesis without direct experimental support; not used as an asserted molecular edge.; moderate: Measured observational association in one cohort, or general clinical symptom description.. These categories describe edge support rather than causal certainty. The report provides APA-style citations and references.

- Non-Small Cell Lung Cancer Treatment (PDQ) - Patient Version: https://www.cancer.gov/types/lung/patient/non-small-cell-lung-treatment-pdq
- Small Cell Lung Cancer Treatment (PDQ) - Patient Version: https://www.cancer.gov/types/lung/patient/small-cell-lung-treatment-pdq
- Tagrisso (osimertinib) tablets: Prescribing information: https://www.accessdata.fda.gov/drugsatfda_docs/label/2024/208065s033lbl.pdf
- Topotecan injection: Prescribing information: https://www.accessdata.fda.gov/drugsatfda_docs/label/2019/200582s003lbl.pdf
- Extracellular matrix organization (R-HSA-1474244): https://reactome.org/content/detail/R-HSA-1474244
- Lectin pathway of complement activation (R-HSA-166662): https://reactome.org/content/detail/R-HSA-166662
- Signaling by EGFR (R-HSA-177929): https://reactome.org/content/detail/R-HSA-177929
- Cough (HP:0012735): https://hpo.jax.org/browse/term/HP:0012735
- Dyspnea (HP:0002094): https://hpo.jax.org/browse/term/HP:0002094
- NM_005228.5(EGFR):c.2573T>G (p.Leu858Arg), Variation ID 16609: https://www.ncbi.nlm.nih.gov/clinvar/variation/16609/
- HGNC gene symbol reports and REST API: https://www.genenames.org/
- Module 2: LUAD integrated multi-omics results: https://github.com/Ajaykanaparthy/biot6900/blob/main/module_2/data/luad_integrated_gene_level.csv
- Module 3: LUAD tractability predictions and target ranking: https://github.com/Ajaykanaparthy/biot6900/blob/main/module_3/targets_luad_w3.csv
- Module 4: Verified structural pocket summary: https://github.com/Ajaykanaparthy/biot6900/blob/main/module_4/results/verified_structure_summary.csv
