# Module 3: LUAD target triage with protein language models

The completed assignment is [BIOT6900_Module3_Assignment.ipynb](BIOT6900_Module3_Assignment.ipynb). It uses the existing 15 LUAD targets from [Module 2](../module_2/targets_luad.csv). The shorter notebook in `archive/root_copies/` is an unfinished historical starter, not the completed analysis.

## Reproduced results

The original assignment was executed in a fresh kernel from beginning to end. All 11 code cells passed, including the modality plot, mapping, leakage exclusion, validation, nearest-neighbor assessment, anchor variants, and export. The professor's source cells and instructions were preserved; existing TODO comments label implemented calculations.

- All 15 target genes map uniquely to reviewed-protein rows in the supplied proteome table; none are unmapped.
- GDF10 and SCN7A were removed from the 1,201-protein training panel before model fitting, leaving 1,199 proteins. No target remains in the training panel.
- Five-fold family validation has no family-group overlap between training and test sets.
- The regenerated `targets_luad_w3.csv` is byte-for-byte identical to the existing CSV. All original target ranks and scores are retained.

| Validation | AUROC | Average precision |
|---|---:|---:|
| Random split | 0.844 | 0.724 |
| Family split | 0.726 | 0.537 |

The positive-label prevalence, and hence chance average precision, is 0.334. Unrounded metrics and out-of-fold predictions are retained in [validation/](validation/).

| Candidate | Module 2 rank | Module 3 rank | Small-molecule prior | Antibody prior |
|---|---:|---:|---:|---:|
| COL6A6 | 3 | 1 | 0.280420 | 0.030515 |
| SCN7A | 5 | 3 | 0.644801 | 0.000398 |
| MASP1 | 11 | 8 | 0.604960 | 0.549767 |

## Data, model, and ranking

The course data pack supplies protein tables, ESM C 300M main embeddings, clinical-precedence labels, and ESM-2 anchor substitution scores. Its [manifest](data/manifest.json) records `dry_run=false`; this audit reuses those supplied inputs and does not rerun the protein language models. Table/embedding dimensions, unique target mapping, finite embedding values, and input hashes were checked.

The provided classifier is `StandardScaler` followed by `LogisticRegression(C=0.05, max_iter=5000)`, fitted separately to small-molecule and antibody clinical-precedence labels. These outputs are **tractability priors, not calibrated probabilities of clinical success in LUAD**. The three nearest training proteins with clinical precedence were independently verified by cosine similarity for every target.

The existing weighting method is unchanged:

```text
tractability_prior = max(p_sm, p_ab)
final = 0.6 * percentile_rank(Module_2_score)
      + 0.4 * percentile_rank(tractability_prior)
rank = descending_rank(final, method="first")
```

Percentiles are calculated within the 15 targets, with larger values preferred. The 0.6 weight retains a majority contribution from LUAD multi-omics evidence while using the prior to triage feasibility. It is a chosen teaching weight, not an optimized or clinically validated parameter. The professor's arithmetic and tie handling are retained; apparently equal rounded final scores can differ at floating-point precision.

COL6A6 rises from rank 3 to 1 through the combined within-list percentiles; its classifier priors remain modest. ZNF536 falls from rank 2 to 7 because its tractability priors are low despite strong Module 2 evidence. SCN7A has a much larger small-molecule than antibody prior. MASP1 has substantial priors for both modalities and remains the second existing structural follow-up candidate alongside SCN7A. Neither the ranking nor sequence similarity establishes binding, disease causality, or therapeutic efficacy.

## Required anchor variants

| Gene | Variant | Supplied score | Within-protein percentile |
|---|---|---:|---:|
| TP53 | R175H | -5.974362 | 19.3 |
| APOE | C130R | 7.207705 | 100.0 |
| APOE | R176C | -8.559011 | 11.8 |
| TREM2 | R47H | -5.155936 | 47.7 |

Each score matches its full supplied substitution-score table; wild-type residues match the precursor protein sequence. Empirical percentiles reproduce the supplied values at one decimal place. APOE uses precursor numbering, including the 18-residue signal peptide. These required cross-disease anchors are retained even though the project disease is LUAD. Their scores and percentile ranks alone do not establish clinical pathogenicity.

## Reproducibility and audit files

Use the `biot6900` environment and `module_3/` as the notebook working directory. `TARGETS_FILE = "../module_2/targets_luad.csv"` deliberately resolves the existing Module 2 file without duplicating it; `DISEASE="luad"` and `TOP_N=15` remain set.

```powershell
conda activate biot6900
cd module_3
jupyter notebook BIOT6900_Module3_Assignment.ipynb
```

Use **Kernel -> Restart & Run All**. To rerun the fresh-kernel execution and independent audit:

```powershell
python -B scripts/validate_portfolio_module3.py
```

- [Mapping, scores, and ranking](validation/target_mapping_and_scores.csv)
- [Unrounded nearest-neighbor evidence](validation/nearest_drug_targets.csv)
- [Anchor audit](validation/anchor_variant_audit.csv)
- [Validation metrics](validation/validation_metrics.csv), [family-fold audit](validation/family_fold_audit.csv), and [removed training overlaps](validation/training_overlap_removed.csv)
- [Out-of-fold predictions](validation/out_of_fold_predictions.csv)
- [Exact regenerated target CSV](validation/reproduced_targets_luad_w3.csv)
- [Machine-readable validation report](validation/validation_report.json), including versions and genuine input hashes

Original notebook and CSV copies are preserved in [archive/portfolio_deliverable_1_originals/](archive/portfolio_deliverable_1_originals/); earlier historical files are retained separately in `archive/root_copies/`.

## Limitations and next step

Family validation measures prediction of the supplied clinical-precedence labels, not independent LUAD therapeutic efficacy. The data pack flags COL6A6 and TNXB main embeddings as truncated, so their full protein sequences are not represented. Sequence-based priors do not establish tumor selectivity, target accessibility, pocket geometry, toxicity, or response to an intervention. The existing SCN7A and MASP1 follow-up requires the separate Module 4 structural and pocket evidence. Experimental work would be needed to test expression/localization, binding, target perturbation, and LUAD-specific biological effects. No candidate is presented as an experimentally proven LUAD therapeutic target.
