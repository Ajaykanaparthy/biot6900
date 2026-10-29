# Verified cross-module LUAD target summary

Rows follow the existing Module 3 ranking. No new composite ranking was invented.
RNA and protein values are paired tumor-minus-normal log2 effects; mutation frequency is a cohort proportion.
Tractability outputs are course classifier scores, not validated probabilities of clinical success.
Displayed classifier scores are rounded to three decimals; exact nonzero values are retained in the CSV.
A dash means the structure/pocket was not assessed, rather than a score of zero.

| Gene | M2 rank | M3 rank | RNA log2FC | Protein log2FC | Mutation % | Concordant | SM score | Antibody score | Mean pLDDT | Pocket druggability / confidence |
|---|---:|---:|---:|---:|---:|---|---:|---:|---:|---|
| COL6A6 | 3 | 1 | -3.042 | -1.885 | 8.33 | yes | 0.280 | 0.031 | - | - |
| ACTN2 | 1 | 2 | -3.447 | -1.731 | 12.04 | yes | 0.220 | 0.000 | - | - |
| SCN7A | 5 | 3 | -2.936 | -1.407 | 9.26 | yes | 0.645 | 0.000 | 73.74 | 0.982 / keep |
| TNR | 4 | 4 | -3.763 | -1.266 | 11.11 | yes | 0.140 | 0.236 | - | - |
| TNXB | 8 | 5 | -2.952 | -1.456 | 5.56 | yes | 0.044 | 0.279 | - | - |
| ABCA8 | 6 | 6 | -3.257 | -1.581 | 5.56 | yes | 0.180 | 0.004 | - | - |
| ZNF536 | 2 | 7 | -2.716 | 2.712 | 12.04 | no | 0.003 | 0.000 | - | - |
| MASP1 | 11 | 8 | -3.403 | -1.135 | 6.48 | yes | 0.605 | 0.550 | 90.46 | 0.687 / keep |
| NOS1 | 9 | 9 | -3.691 | -1.069 | 7.41 | yes | 0.204 | 0.001 | - | - |
| ADAMTS16 | 7 | 10 | 2.558 | 1.234 | 10.19 | yes | 0.099 | 0.099 | - | - |
| SEMA5A | 14 | 11 | -2.431 | -1.096 | 12.04 | yes | 0.100 | 0.251 | - | - |
| MYH2 | 10 | 12 | -2.726 | -1.611 | 5.56 | yes | 0.047 | 0.003 | - | - |
| COL11A1 | 13 | 13 | 6.157 | 0.871 | 20.37 | yes | 0.102 | 0.000 | - | - |
| ERICH3 | 12 | 14 | -3.340 | -0.989 | 9.26 | yes | 0.000 | 0.005 | - | - |
| GDF10 | 15 | 15 | -3.710 | -1.354 | 4.63 | yes | 0.142 | 0.032 | - | - |

The [CSV](target_prioritization_summary.csv) includes all pair counts, p/q-values, thresholded evidence flags, UniProt mappings, embedding truncation flags, separate fpocket Score and Druggability Score, and candidate-specific rationales/limitations.

A high Module 2 rank measures relative dysregulation and mutation recurrence. It does not demonstrate ligand binding, tumor dependency, benefit from inhibiting a downregulated gene, or an acceptable safety margin.
