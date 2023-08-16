# procap-network
Active transcriptional regulatory elements are bidirectionally transcribed, the level of transcription correlating with regulatory activity. Using capped nascent RNA sequencing (PRO-cap) in 69 human lymphoblastoid cell lines, we previously identified over 80 thousand transcribed transcriptional regulatory elements (tTREs), both enhancers and promoters. Co-expression analyses, comparing tTRE expression across individuals, reveal interactions between tTREs. We find evidence of TFs influencing mid- to long-range interactions and of strand-dependent cooperativity at closely-spaced tTREs.

Procap-Networking

<img src ="https://github.com/sl2665/procap-network/assets/42938330/1b3e1813-5b49-445b-9182-68d53e0f7b1c" style="height: 500px;"/>

## Identification of tTREs and selection of variably expressed tTREs

Transcribed Transcriptional Regulatory Elements (tTREs) were identified from 76 partially replicated PRO-cap data in Lymphoblastoid Cell Lines (LCLs) from 67 individuals from the Yoruban population (YRI) as described by Kristjánsdóttir et al21 (Supplementary Table 1). Briefly, we merged all the PRO-cap reads (~ 1.4 billion unique molecular identifiers separated reads) from the dataset that were mapped to the hg19 reference genome. The reads were scanned along the genome to pick out the local maxima within the 300 base window in a strand specific manner. The local maxima peaks were then matched to another local maxima of the opposite strand between 50 - 250 bases upstream on the antisense direction, so that the paired PRO-cap peaks on both strands form a divergent bidirectional transcription pattern. Single strand peaks without a divergent bidirectional pair were discarded. If there were multiple single strand peaks within the 150 base pair window, we selected the peak with the highest amount of PRO-cap reads. As a result, the closest elements are at least 150 base pairs apart, which is well above our resolution for distinguishing two nascent transcription start sites. We excluded tTREs with ambiguous start sites within 150 bp. The 150 bp cutoff corresponds to one nucleosome distance, which has a structural rationale that enhancers with accessible chromatin have at least one nucleosome removed to create an open chromatin, and eRNA transcription occurs in a bidirectional manner around the boundaries of this region. 89.5% of the tTREs we used were at least 300 bp apart from each other. This process identified 76,730 tTREs that were bidirectionally transcribed.

To identify tTREs that are variably expressed, we used normalized reads-per-million (RPM) normalized PRO-cap read count data containing partial replicates as described previously21 . Briefly, we used the q-value method described by Storey et al22. We used partially replicated samples as the level of technical variation and used the variation in technical variation as the reference to calculate p-values of pairwise differences between non-replicated different individuals. For each tTRE, we calculated the deviation from the mean of the normalized read counts between replicates and between different samples. We then used a one-sided Wilcoxon’s rank sum test to test the alternative hypothesis that the differences between samples were greater than between the replicates for each tTRE, and calculated p-values. We estimated the number of variably expressed tTREs by analyzing the complete distribution of the p-values as described previously22. Under the null hypothesis, p-values should have a uniform distribution with a density of 1, but the observed p-values are only uniformly distributed only for large p-values. The density of the portion of the p-value distribution that is uniform is ~0.281, indicating that up to ~71.9% of tTREs can be considered variably expressed. Using FDR < 0.2, we identified 29,694 variably expressed tTREs (40% variably expressed). 
## Distance-dependent pair-wise co-expression analysis of tTREs. 
We used the variably expressed tTREs (n=29,694; promoter - 4,006; enhancer - 25,688 using the CAGE based criteria), and calculated correlation coefficients of the 75 individual normalized read counts (67 individuals + 8 replicates) between two tTREs within 5 Mb distance (2,249,839 pairs). For the distance analysis, we binned the correlation coefficients by the distance between 2 tTREs from all the tTREs. The bins are generated based on fixed distance intervals up to 1,024 kb (0 - 1 kb, 1 - 2 kb, 2 - 4, kb, etc.), or a fixed number of tTRE pairs (1,000 pairs per bin) with variable distance intervals. All the tTRE pairs were grouped into distance groups, and the box plots were generated to display the median, 25’th and 75’th percentiles of the correlation coefficients. All of the distance groups in the box plot analyses contained at least 200 tTRE pairs, and the comparison of their means could be expected to follow a Gaussian statistic in student t-tests.

With the variable interval bins of 1,000 pairs, we generated plots for the 5th percentile, median and the 95th percentiles of the correlation coefficient distributions within each bin of 1,000 pairs along the distance (Supplementary Fig. 1d-e). These percentiles were used to generate a step function of correlation coefficient percentiles as a function of distance. The top 5’th percentile is the 50’th highest correlation coefficient per bin (of 1,000 elements) in this step function analysis, and using the false discovery analysis of correlation p-values, the 5’th percentile corresponded to FDR < 0.006. The correlation percentile step function was superimposed on the color density scatterplot of the correlation coefficients for visualization (traced scatterplot) and Area Under the Curve (∆AUC) analysis.

## Estimation of the correlation bias of the spurious correlation and the genotype linkage
To obtain an estimate of spurious correlations, we used 2 million random inter-chromosomal correlations and correlations that are more distant than 1 megabase away as the background distribution of PRO-cap correlations. Both distributions are fitted to the Gaussian distribution with the standard deviation of the correlation distribution and mean of 0, and tested with QQ-plot to indicate their normality using the “qqnorm” function in the R statistical package.

We performed an independent genotype analysis to exclude the possibility that some of the tTRE variation is genetically driven by SNPs and to ensure that the observed correlations were not confounded by the genetic association of the SNPs in the YRI population. To exclude this possibility, we performed the same correlation analysis on tTREs that are not genetically associated (discrete Pearson’s correlation of the genotypes labeled as 0, 1, or 2 - reference, heterozygous, alternative alleles). If there were no variable SNPs within the enhancer region, these regions were excluded from the genotype association and considered genotype-independent. Minor allele frequency selection (greater than 0.05) was applied only to the enhancers with variable SNPs to exclude tTREs whose eRNA expression levels were suspected to be genetically associated with SNPs. From this analysis, approximately 5% of tTRE-tTRE pairs were significantly genetically associated by this criterion (FDR<0.1) and were removed from the analysis. At least 70% of the tTRE-tTRE pairs did not contain any SNPs associated in the population (discrete Pearson correlation less than 0.05, n=129,660), allowing for more rigorous cut-off of genotype-independence to exclude that genetic linkage in the population confounded the co-expression patterns.

## Co-expression analysis between tTRE nascent transcription and mRNA expression. 
We used the RNA-seq expression data from Pickrell et al.7 which included normalized RNA levels in 161 replicated LCL datasets from the same 67 YRI individuals that we used in PRO-cap. We selected 13,002 genes with the mean expression levels greater than 1 RPKM. 275,660 pairs of tTREs and annotated mRNA TSS within 1 Mb were tested, and the correlation coefficients of the 67 individual samples were calculated. To assign mRNA gene positions, we selected mRNA TSS positions according to the following criteria, in contrast to Kristjánsdóttir et al.21 which considered all annotated TSSs for each mRNA. First, we selected annotated mRNA TSSs that overlapped with a promoter tTRE within 250 base pairs of distance. Then, we further selected the mRNA TSS with the highest PRO-cap expression level within the same mRNA transcripts. Therefore, 1 representative mRNA TSS position was selected for each gene. Correlation coefficients were calculated as described above for tTRE-tTRE co-expression. 

## Analysis of the ENCODE factor dependencies in tTRE co-expression.
The top 5 percentile, as well as the median, of the correlation coefficients in the variable interval bins were used as the indicators of tTRE co-expression. This correlation coefficient percentile serves as a step function of the distance, and we refer to it as the correlation decay plot. The correlation decay plots were generated by subsetting the tTRE pairs by their distance with 1,000 tTRE pairs per bin as described in the ‘Distance-dependent pair-wise co-expression analysis of tTREs’ section. The top 5 percentile of the correlation coefficient was generated from the pairwise correlation coefficients between any pair of tTREs within 200 kb distance, resulting in a total of 192.4 thousand total pairs. The Area Under the Curve was calculated for the 5’th percentile step function we generated between 0 kb and 200 kb range. The same ∆AUC was generated using the median step function to validate the ∆AUC based on the 5'th percentile step function to evaluate the consistency of the metric and liability to noises.
∆AUC values were calculated to evaluate whether TFs located either between the tTREs or at the tTREs affected the correlation decay plots. We used the ENCODE FACTORBOOK binding sites in the representative LCL GM12878- to separate the 192.4 thousand tTRE pairs between 0 kb and 200 kb distances into 2 groups based on the respective TF binding or intersection statues, and calculated the difference between the area under the correlation decay curves in the 0 - 200 kb distance range (∆AUC). Both the 5’th percentile and the median correlation decay step functions were used, and ∆AUCs from the 5’th percentile and the median values correlated well across TFs (Supplementary Fig. 4a-b). For the TF intersection analysis, we compared 0 or 1 TF intersection against 2 or more TF intersections between the tTRE pairs, as there were not enough tTRE pairs with 0 TF intersection in longer distance ranges for statistical comparisons. 
The statistical significance of ∆AUC values on TF intersection and occupancy was estimated using a bootstrapping randomization strategy. We also considered whether the number of factor binding sites affects the dispersion of ∆AUC values, by using a different number of mock sites to simulate the effects of different numbers of TF binding sites on the ∆AUC calculation, especially for the TFs with fewer or higher binding sites. The expected mean and standard deviation of ∆AUC with the randomized sets allowed us to estimate the p-values and FDR of ∆AUC in the TF intersection and occupancy plots.
Specifically, the significance of the ∆AUC was estimated by comparison with the background distribution of ∆AUC which was generated by randomly shuffling the genomic locations of the factor binding sites. We calculated the background ∆AUC distributions by 1,000 permutations of 5,000 to 70,000 randomly shuffled ENCODE factor binding sites, maintaining the same tTRE expression vector levels across individuals. The background ∆AUC distributions followed a normal distribution and were dependent on the number of binding sites. We approximated the expected mean and standard deviation (sd) of the ∆AUC as a function of the actual number of binding sites for the specific factor, and used this to generate a z-score and p-value of the ∆AUC between correlation decay curves enriched or depleted with the factor binding sites (Supplementary Fig. 4c-d). For example, CTCF contains 41,465 ChIP-seq binding sites28, and the ∆AUC = +27.8 between more than 2 CTCF sites or less than 1 site intersecting a pair of tTREs. The background distribution of ∆AUC with 40,000 sites intersecting tTREs is +6.41 ± 2.56 (mean ± sd), and we obtained the z-score = +8.34 and p = 7.26 x 10-17 for CTCF intersections affecting correlation decay plots. The intersection and the occupancy scores of all other 76 factors available for GM12878 in FACTORBOOK were calculated in this way. Only factors with a sufficient number of binding sites in each category were reported. These p-values and z-scores were used in the clustering analysis to cluster the TFs. 
Hi-C analysis
As a representative LCL, GM12878 Hi-C data were obtained from the 4D Nucleome Consortium public database (accession # 4DNFIXP4QG5B). To extract the Hi-C contact frequency data, we first converted our tTRE coordinates from the hg19 to the hg38 reference genome used by the 4D Nucleome Consortium using NCBI liftover. Contact frequencies from the regions 500 kb upstream and downstream of the tTRE positions were extracted from the Hi-C mcool data using the cooler package (https://cooler.readthedocs.io/en/latest/index.html) at 5 kb resolution. The extracted contact frequency table was used to query the exact contact frequency between tTRE pairs within 500 kb distance. The distance-dependent decay, ENCODE TF dependency, and the ∆AUC analyses were performed in the same way as the co-expression correlation analyses, by using the contact frequency values and the log of contact frequencies instead of the co-expression coefficients, and using the same tTRE pair classifications. Background subtraction was performed on the linear values of Hi-C contact frequencies prior to log transformation.
Strand-specific co-expression analysis for tTRE nascent transcription at adjacent sites
For each tTRE, the nearest adjacent downstream tTRE was identified and a linear regression between its strands was computed across the samples according to the following categories: Sense (plus-plus, minus-minus), Antisense convergent (plus-minus), and Antisense divergent (minus-plus). Pairs of tTREs were binned based on the distance between them, and the distributions of Pearson correlation coefficients were compared between the categories. Pairs within 250 bp were not included to avoid the possibility of counting the same reading in both tTREs. A total of 21,486 pairs of adjacent tTREs (16,857 for enhancers only) within 10 kb distances were compared.  To compare distance dependence, LOESS fits were generated by taking 200 local data points on the distance axis for each data point to calculate the local polynomial regression curves (Supplementary Fig. 6a). 
For the adjacent vs interleaved tTRE analysis, two sets of tTRE pairs were selected, those that were immediately neighboring (adjacent) or those that contained one other tTRE in between (interleaved). The LOESS fits of their PRO-cap co-expression correlation coefficients or Hi-C contact frequencies in 1 kb resolution were generated.



