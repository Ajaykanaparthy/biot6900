# LUAD Task 1 data: PGS003393 and 1000 Genomes

**The ages and phenotype labels are SIMULATED. These are not observed LUAD patients or diagnoses.**
The genotype calls, REF/ALT alleles, ancestry labels, and published score weights are real.
This document describes the LUAD files in this directory. The pre-existing
`../README.md` describes the separate Alzheimer's teaching cohort and does not describe these files.

## Sources

- Score: [PGS003393 / best_LUAD](https://www.pgscatalog.org/score/PGS003393/),
  Namba S et al., Cancer Research (2022), DOI [10.1158/0008-5472.CAN-22-1492](https://doi.org/10.1158/0008-5472.CAN-22-1492).
  [Original scoring file](https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS003393/ScoringFiles/PGS003393.txt.gz): 74 SNPs, GRCh37, beta weights.
  The publisher's compressed-file MD5 was downloaded and verified. All 74 published
  allele pairs and exact weight strings are retained in `score.csv`.
- Genotype calls: [1000 Genomes Phase 3 release 20130502](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/),
  per-autosome `ALL.chrN.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`
  and corresponding Tabix indices. Each VCF declares `hs37d5.fa.gz`, which is GRCh37
  plus decoy sequence. No GRCh38 coordinates or liftover are used.
- Ancestry: [official Phase 3 sample panel](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel),
  `super_pop` and `pop`, joined by exact sample ID. No Alzheimer's genotype or metadata row was reused.
- Official callset and known-issues READMEs are preserved in `../data_sources/`.

## Files and cohort

| File | Schema / contents |
|---|---|
| genotypes.csv | sample_id, then 64 usable variant columns; 2,504 rows; VCF ALT allele counts 0/1/2 |
| variants.csv | variant,ref,alt for those 64 genotype columns, in the same order |
| score.csv | variant,effect_allele,other_allele,weight; all 74 original score rows |
| meta.csv | sample_id,ancestry,age,phenotype plus population and explicit simulation flags/source columns |
| variant_match_audit.csv | All 74 targets, their source coordinates/alleles, VCF matches, inclusion status, and exclusion reasons |
| data_preparation_manifest.json | Source URLs, method, counts, simulation specification, range-log paths, and actual SHA-256 hashes |
| SHA256SUMS.txt | Full SHA-256 checksums of the four required CSVs |

All 2,504 participants are included in original VCF sample order, identical in
`genotypes.csv` and `meta.csv`; no ancestry balancing or participant sampling was performed.

| Ancestry code | Participants |
|---|---:|
| AFR | 661 |
| AMR | 347 |
| EAS | 504 |
| EUR | 503 |
| SAS | 489 |
| Total | 2504 |

AFR, AMR, EAS, EUR, and SAS are study super-population labels, not precise
individual ancestry fractions. Several source populations are admixed.

## Matching and exclusions

Variant IDs have the form `chrN:position`, using the original 1-based GRCh37
coordinates. The source score does not supply rsIDs; the v5b release intentionally
removed rsIDs. No identifiers were invented to look like rsIDs.

Match exact chromosome, position, and unordered two-allele pair. Only PASS,
nonpalindromic biallelic SNPs with complete diploid 0/1 calls are included.
Genotypes count allele index 1 (VCF ALT); score alleles and weights remain as published.
Effect-allele orientation is recorded in the audit for later notebook harmonization.
No missing calls are imputed during this preparation. Phase 3 source calls themselves
were statistically genotyped and phased as described in the official release README.

The pipeline excludes multiallelic/marked-MULTI_ALLELIC matches rather than collapsing
third alleles. A/T and C/G pairs are excluded conservatively because this source
does not provide independent forward-strand orientation or effect-allele frequencies.
Coordinate harmonization alone cannot resolve that ambiguity. No strand complements
are guessed. Every found, absent, incompatible, or excluded target is retained in the audit.

| Audit status | Variants |
|---|---:|
| excluded_strand_ambiguous | 10 |
| usable | 64 |

The resulting usable set is a **64-variant subset of PGS003393**, not the complete
published 74-variant score. Published validation metrics for the full score do not apply.

## Simulated metadata (explicitly user-authorized)

`age_is_simulated=True` and `phenotype_is_simulated=True` for every participant.
The source columns label ages as synthetic teaching values and phenotypes as synthetic
teaching outcomes, not observed LUAD. Phenotype 1 is a simulated case label; 0 is a
simulated control label. Real ancestry/population fields were not simulated.

Fixed seed 6900; NumPy 2.4.6, `Generator(PCG64)`. In original VCF sample order:

1. Draw age from Normal(mean=60, SD=8), clip to [18,90], and round to three decimal places.
2. Form a partial score from the 64 usable real variants and unchanged published
   weights, counting effect-allele copies. Standardize across all 2,504 samples with population SD (`ddof=0`).
3. Set `logit(p) = log(0.10/0.90) + 0.5 * partial_score_z + 0.03 * (age - 60)`.
4. After all age draws, draw one Uniform(0,1) per sample from the same RNG; label 1 when the draw is below p.

All coefficients and distributions above are arbitrary teaching choices, fixed before
any performance evaluation. No AUC, calibration, or ancestry-specific metric was used
to tune them. The intercept does not establish a measured LUAD prevalence.

**Scientific limits:** The outcome mechanism includes the score later evaluated,
so apparent discrimination or value beyond age is partly built into the simulation.
The same score coefficient is used in every ancestry; this simulation does not model
ancestry-specific LD tagging or effect differences. It cannot validate clinical LUAD
risk, establish real portability failure, or show threshold fairness for real patients.
Actual ages, verified LUAD status, relevant covariates (including smoking), and suitable
independent cohorts would be needed for disease-performance validation. The score's
Catalog development/evaluation samples were European; these reference samples are not
an independent LUAD case-control cohort.

## Provenance and verification

Only Tabix indices and bounded HTTP 206 byte ranges around the 74 positions were
retrieved; 69,206,016 compressed VCF bytes in total.
No whole chromosome/genome VCF was downloaded. Every successful range has its returned
Content-Range, byte count, and computed SHA-256 recorded in `../data_sources/extracted/`.
The raw target VCF records, including excluded variants, are preserved there as small
gzip subsets. A range hash or subset hash is not a whole-VCF checksum; none is claimed.

Preparation code: `../scripts/prepare_luad_data.py` and `../scripts/remote_tabix.py`.
Independent validation: `../scripts/validate_luad_data.py`. Preparation refuses to
overwrite existing required CSVs. Download caches and source subsets are retained.

All pre-existing files were checked against the initial SHA-256 baseline. The original
professor notebook remains byte-for-byte unchanged, `USE_DEMO=False`, and was not run.
The Alzheimer's data remain in their original misspelled backup and are also copied
unchanged to `../data_alzheimers_backup/`. No files were deleted and no Git staging,
commit, or push was performed.

Task 1 stops here. No notebook TODOs, figures, performance results, or final report
have been completed. The original notebook still contains Alzheimer's/APOE-specific
instructions that need a later, explicitly approved LUAD adaptation.
