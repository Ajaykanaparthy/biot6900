# BIOT 6900 Module 6: LUAD Polygenic Risk Scores

## Objective

Evaluate a lung adenocarcinoma (LUAD) polygenic risk score (PRS) using real genotype calls, examining allele harmonization, score distributions, a pooled risk threshold, discrimination, calibration, and value beyond age across ancestry groups.

**Ages and LUAD phenotype labels are simulated. Participants are not documented LUAD patients; all performance results below use synthetic teaching outcomes.**

## Data sources

- **Published score:** [PGS003393 (best_LUAD)](https://www.pgscatalog.org/score/PGS003393/), Namba S et al., *Cancer Research* (2022), [doi:10.1158/0008-5472.CAN-22-1492](https://doi.org/10.1158/0008-5472.CAN-22-1492). The original GRCh37 score contains 74 variants; published alleles and beta weights are preserved in `data/score.csv`.
- **Genotypes:** [1000 Genomes Phase 3, release 20130502](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/), v5b VCFs on GRCh37/hs37d5. Targeted indexed retrieval supplied genuine individual genotype calls for 2,504 participants.
- **Ancestry:** Real super-population and population labels from the [official Phase 3 sample panel](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel), matched by sample ID.

The four inputs are `data/genotypes.csv`, `data/variants.csv`, `data/score.csv`, and `data/meta.csv`. See the [data README](data/README.md) for source provenance and preparation methods. The Alzheimer's datasets remain preserved in their backup folders.

## Variant harmonization

All 74 published positions were found in the source VCFs. Ten A/T or C/G variants were intentionally excluded because independent evidence could not resolve their strand orientation. The analysis therefore uses a **64-variant subset** of PGS003393.

Harmonization: **64 used, 0 flipped, 0 dropped** among the retained genotype columns. Every retained effect allele matches VCF ALT, and every other allele matches REF. Exported dosages count ALT copies (0/1/2) and were verified against source VCF calls. The ten excluded variants are absent from the genotype table and may appear as "not found" during notebook matching; they were not missing from the source data. Details are recorded in the [variant match audit](data/variant_match_audit.csv).

## Simulated ages and phenotypes

As documented in the [data README](data/README.md), simulation uses seed **6900** with NumPy `Generator(PCG64)`. Ages are drawn from Normal(60, 8), clipped to 18-90 years, and rounded to three decimal places. The partial PRS is standardized across all participants using population SD (`ddof=0`).

```text
logit(p) = log(0.10 / 0.90) + 0.5 * partial_score_z + 0.03 * (age - 60)
```

After the age draws, one uniform draw per participant assigns phenotype 1 when the draw is below p, otherwise 0. Both fields are explicitly flagged as simulated in `meta.csv`. Parameters were fixed without tuning to performance metrics.

## Results

AUC measures discrimination against the simulated phenotype.

| Group | Participants | PRS AUC |
|---|---:|---:|
| AFR | 661 | 0.629 |
| AMR | 347 | 0.585 |
| EAS | 504 | 0.686 |
| EUR | 503 | 0.605 |
| SAS | 489 | 0.662 |
| **Overall** | **2,504** | **0.639** |

| Model | Cross-validated AUC |
|---|---:|
| Age only | 0.609 |
| Age + PRS | 0.675 |

Figures:

- [PRS distribution by ancestry](Figures/prs_distribution.png)
- [Calibration by ancestry](Figures/calibration.png)

## Run the graded notebook

Activate the course environment and open the professor's original notebook from `module_6`:

```powershell
conda activate biot6900
cd C:\Users\kanap\OneDrive\Desktop\biot6900\module_6
jupyter notebook BIOT6900_Module6_Part3_Graded_Starter.ipynb
```

Keep **`USE_DEMO = False`**, ensure the working directory is `module_6`, and use **Kernel -> Restart & Run All**. The loader reads the four CSVs from `./data/`.

SHA-256 fingerprints (first eight hexadecimal characters):

| File | SHA-256 prefix |
|---|---|
| `genotypes.csv` | `62fbc28e` |
| `variants.csv` | `5349e113` |
| `score.csv` | `e8d74cb7` |
| `meta.csv` | `fcf43773` |

Full checksums are in [data/SHA256SUMS.txt](data/SHA256SUMS.txt).

## Limitations

Simulated outcomes cannot validate clinical LUAD prediction or real-world ancestry portability. The outcome-generation formula includes the evaluated PRS, so discrimination and its apparent value beyond age are partly built into the simulation. The same score coefficient applies across ancestries and does not model ancestry-specific LD effects. Results also concern a 64-variant subset, rather than the full published score. Independent cohorts with observed ages, verified LUAD outcomes, and relevant covariates such as smoking would be needed for clinical validation.
