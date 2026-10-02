# Independent Multi-Omics Target Discovery in Lung Adenocarcinoma (LUAD)

## Disease and Data Choice

For this project, I chose **lung adenocarcinoma (LUAD)**. I selected LUAD because there are good-quality public datasets available for all three omics layers required for this assignment: transcriptomics, proteomics, and genomics. LUAD is also a good disease for this type of analysis because several important genes, including **TP53, KRAS, EGFR, STK11, KEAP1, MET, and BRAF**, are already well studied. This gave me a way to check whether my analysis was recovering known LUAD biology while also looking for less obvious genes that may be worth further investigation (The Cancer Genome Atlas Research Network, 2014).

I obtained the data from the **CPTAC pan-cancer LUAD dataset available through LinkedOmics**. The datasets were openly available and did not require an application or data-use agreement (LinkedOmics, n.d.). I used five raw files: RNA-seq tumor and normal data, proteomics tumor and normal data, and a gene-level somatic mutation dataset.

The RNA dataset contained **110 tumor samples and 101 normal samples**. The proteomics dataset also contained **110 tumor samples and 101 normal samples**, while the mutation dataset contained **108 samples**. These data came from the CPTAC LUAD proteogenomic cohort, which was designed to study lung adenocarcinoma across multiple molecular layers (Gillette et al., 2020).

The raw RNA and protein files were gene-by-sample matrices, so I first reduced them to the gene-level format required for the assignment. For each gene, I calculated the average tumor-versus-normal log2 fold change and a paired t-test p-value using the overlapping tumor and normal sample IDs. This produced the required structure:

`gene | log2fc | pval`

For the genomics layer, the mutation dataset was already a binary gene-by-sample matrix showing whether a gene was mutated in each tumor sample. I calculated the **mutation frequency** for each gene and used this as the genomic association-strength measure. This was appropriate for this dataset because the assignment allowed the genomics layer to use `neglog10p` or another association-strength measure.

The raw files used Ensembl gene IDs, so I removed the version numbers and converted the Ensembl IDs to standard gene symbols before integrating the three datasets. Even though CPTAC provides matched multi-omics data, my final analysis was done using **gene-level summary values**. Therefore, the final integration should be interpreted at the gene level rather than as a patient-by-patient analysis.

## Harmonization, Concordance, and Multi-Evidence Scoring

After preparing the three datasets, I joined the transcriptomic, proteomic, and genomic tables using the gene symbol. A total of **6,758 genes survived the three-way join**. This meant that all 6,758 genes had information available from RNA, protein, and mutation data.

The next step was RNA-protein concordance. A gene was considered concordant when the RNA and protein log2 fold changes moved in the same direction. For example, if both RNA and protein were increased in tumor compared with normal, the gene was considered concordant. The same was true if both were decreased.

Out of the 6,758 genes, **4,918 genes were concordant**, which is about **72.8%** of the genes in the final integrated dataset. The remaining genes showed disagreement between RNA and protein direction.

For the multi-evidence score, I used the absolute RNA log2 fold change as the transcriptomic signal, the absolute protein log2 fold change as the proteomic signal, and mutation frequency as the genomic signal. Each layer was rank-normalized before being combined.

I kept the default **equal weighting** for all three layers. I felt this was the most reasonable choice because each layer provides different information. RNA data show changes at the transcriptional level, protein data show what is actually observed at the protein level, and mutation data show how often a gene is genetically altered in the tumors.

I did not have a strong reason to make one layer more important than another. Using equal weights also made the final score easier to interpret because transcriptomics, proteomics, and genomics each contributed equally to the ranking. If the goal of the study were different, such as focusing specifically on genetic drivers, it may make sense to give more weight to the mutation layer, but I did not do that here.

## Top Targets and Biological Interpretation

After calculating the multi-evidence score, I ranked the genes from highest to lowest and exported the top 15 genes to `targets_luad.csv`.

The top 15 genes were:

**ACTN2, ZNF536, COL6A6, TNR, SCN7A, ABCA8, ADAMTS16, TNXB, NOS1, MYH2, MASP1, ERICH3, COL11A1, SEMA5A, and GDF10.**

The highest-ranked gene was **ACTN2**, followed by **ZNF536** and **COL6A6**. I would not automatically consider every top-ranked gene to be a proven LUAD driver. The score only shows that these genes had strong combined evidence across the RNA, protein, and mutation layers used in this analysis.

One result that stood out to me was **COL6A6**, which ranked third. In my analysis, COL6A6 had an RNA log2 fold change of approximately **-3.04** and a protein log2 fold change of approximately **-1.89**. Both RNA and protein therefore showed a strong decrease in tumor compared with normal tissue.

This result is interesting because previous research has also reported reduced COL6A6 expression in lung adenocarcinoma. Ma et al. (2021) reported that COL6A6 was downregulated in LUAD and also found associations between COL6A6 expression and clinical outcomes. This gives some biological support to the pattern found in my analysis (Ma et al., 2021).

I also checked whether known LUAD genes were recovered in the complete ranked list. Several important LUAD genes were present, including **TP53, MET, KRAS, ALK, EGFR, ERBB2, STK11, NF1, KEAP1, and BRAF**.

For example, **TP53 ranked 572** and had a mutation frequency of approximately **59.3%**. **KRAS ranked 1,624** with a mutation frequency of about **32.4%**, while **EGFR ranked 1,807** with a mutation frequency of about **36.1%**.

These genes are well established in LUAD. Large genomic studies have identified recurrent alterations involving TP53, KRAS, EGFR, BRAF, STK11, KEAP1, NF1, MET, and ERBB2 in lung adenocarcinoma (The Cancer Genome Atlas Research Network, 2014). The CPTAC LUAD study also showed important proteogenomic effects related to known driver alterations such as KRAS, EGFR, and ALK (Gillette et al., 2020).

It may seem surprising that genes such as TP53, KRAS, and EGFR were not in my top 15. However, my score does not depend only on mutation frequency. It combines evidence from RNA, protein, and mutation data equally. A gene can therefore have a high mutation frequency but still receive a lower final score if its RNA or protein changes are not as strong.

Because of this, the final ranking should not be interpreted as a list of the most important LUAD genes overall. It is better to think of it as a **multi-omics prioritization list based on the specific datasets and scoring method used in this project**.

## Discordant Gene

One of the most interesting genes in my results was **ZNF536**, which ranked second overall.

ZNF536 showed a strong disagreement between RNA and protein. Its RNA log2 fold change was approximately **-2.72**, while its protein log2 fold change was approximately **+2.71**. Its mutation frequency was around **0.12**, or 12%.

This means that ZNF536 RNA was strongly decreased in tumor samples, while the protein level showed the opposite pattern and was increased. Because the RNA and protein changes moved in different directions, ZNF536 was classified as a discordant gene.

This result shows why using both transcriptomics and proteomics can be useful. If I had only looked at RNA data, I would have concluded that ZNF536 was strongly decreased. However, the protein data showed a completely different pattern.

There are several possible biological explanations for this type of disagreement. RNA levels do not always directly determine protein levels. Translation can be regulated after transcription, and proteins can also have different rates of degradation and different levels of stability. Post-transcriptional regulation can therefore cause protein abundance to behave differently from RNA abundance. Proteogenomic analysis in LUAD has shown that protein-level information can reveal biological effects that cannot always be seen from genomic or transcriptomic data alone (Gillette et al., 2020).

However, my analysis does not prove why ZNF536 showed this disagreement. I can only say that the RNA and protein layers showed opposite effects in this dataset. More experimental work would be required to determine whether the difference is caused by translation, protein stability, degradation, cell-type composition, or another mechanism.

Because ZNF536 ranked second overall and showed strong RNA-protein discordance, I think it would be an interesting non-obvious gene to study further.

## Limitations and What the Analysis Can Tell Us

The main limitation of this project is that the final integration was performed using **gene-level summary values**.

Although the original CPTAC LUAD cohort contains matched multi-omics samples, I first summarized the RNA and protein datasets into tumor-versus-normal effects and summarized the mutation dataset into gene-level mutation frequencies. I then combined these summary values for each gene.

Because of this, the final score tells me which genes show relatively strong combined evidence across the three layers, but it does not show exactly what happened in an individual patient.

For example, I cannot say that a mutation in one patient's tumor directly caused an RNA or protein change in that same patient. A true patient-level matched analysis would need to compare the RNA, protein, and genomic measurements within each individual patient.

Another limitation is that a high multi-evidence score does not prove that a gene causes lung adenocarcinoma. It also does not prove that the gene will be a useful therapeutic target. Some highly ranked genes may reflect tumor-related changes, tissue composition, extracellular matrix changes, or other biological processes rather than being direct cancer drivers.

The results should therefore be considered **hypothesis-generating**. The ranked list helps reduce thousands of genes into a smaller group of candidates that may be worth investigating further, but experimental validation would still be required.

Overall, this project showed the advantage of combining multiple types of molecular data. Concordant genes provide stronger support when RNA and protein agree, while discordant genes can reveal additional biology that may not be visible from transcriptomics alone. The final target list provides a useful starting point for further biological and experimental validation.

## References

Gillette, M. A., Satpathy, S., Cao, S., Dhanasekaran, S. M., Vasaikar, S. V., Krug, K., Petralia, F., Li, Y., Liang, W. W., Reva, B., Krek, A., Ji, J., Song, X., Liu, W., Hong, R., Yao, L., Blumenberg, L., Savage, S. R., Wendl, M. C., ... Clinical Proteomic Tumor Analysis Consortium. (2020). Proteogenomic characterization reveals therapeutic vulnerabilities in lung adenocarcinoma. *Cell, 182*(1), 200–225.e35. https://doi.org/10.1016/j.cell.2020.06.013

LinkedOmics. (n.d.). *CPTAC pan-cancer LUAD data download*. https://www.linkedomics.org/data_download/CPTAC-pancan-LUAD/

Ma, Y., Qiu, M., Guo, H., Chen, H., Li, J., Li, X., & Yang, F. (2021). Comprehensive analysis of the immune and prognostic implication of COL6A6 in lung adenocarcinoma. *Frontiers in Oncology, 11*, 633420. https://doi.org/10.3389/fonc.2021.633420

The Cancer Genome Atlas Research Network. (2014). Comprehensive molecular profiling of lung adenocarcinoma. *Nature, 511*, 543–550. https://doi.org/10.1038/nature13385
