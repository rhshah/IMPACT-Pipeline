================
File Description
================


**The following files are output by the CMO bioinformatics pipeline and are organized in the data delivery drives as outlined below. Detailed descriptions of all results files are provided. Of particular importance are the following key results files:**
-	CompiledMetrics/Proj_*_ALL_Metrics.pdf		[summary of QC metrics]
-	Variants/annotated_exonic_variants.txt		[all variants and allele counts]
-	CopyNumber/Proj_*_copynumber_segclusp.pdf	[copy number plots]
-	SampleInfo/Proj_*_title_file.txt				[sample IDs and limited metadata]



----------------
QC Metrics Files
----------------

Proj_*_All_Metrics.pdf
----------------------
PDF file with graphical representation of all calculated quality and performance metrics for all samples in project. Numerical values are contained in the text files below.

Proj_*_ALL_basequalities.txt
----------------------------
Illumina base quality score by cycle (following GATK base quality recalibration), displayed for each sample.

Proj_*_ALL_Canonical_exoncoverage.txt
-------------------------------------
Mean coverage for each target interval corresponding to exons in canonical transcripts, displayed for each sample. Coverage is computed for reads with mapping quality > 20.

Proj_*_ALL_exoncoverage.txt
---------------------------
Mean coverage for each target interval corresponding to all protein-coding exons in all transcripts (including flanking splice sites), displayed for each sample. Coverage is computed for reads with mapping quality > 20.

Proj_*_ALL_exonnomapqcoverage.txt
---------------------------------
Mean coverage for each target interval corresponding to all protein-coding exons in all transcripts (including flanking splice sites), displayed for each sample. Coverage is computed for reads regardless of mapping quality (i.e., including reads mapping to multiple locations).

Proj_*_ALL_FPavgHom.txt
------------------------
Average minor allele fraction across all tiling SNPs that are homozygous in the given sample. This is an indication of the degree of contamination from unrelated DNA. (The IMPACT and HemePACT panels contain >1,000 “tiling SNPs” where both alleles are common in the population.)

Proj_*_ALL_FPCResultsUnMatch.txt
--------------------------------
Pairs of samples from different individuals with concordant fingerprint genotypes. (The number represents the fraction of tiling SNPs homozygous in one sample which are homozygous for the alternate allele in the other sample.)

Proj_*_ALL_FPCResultsUnMismatch.txt
-----------------------------------
Pairs of samples from the same individual with discordant fingerprint genotypes. (The number represents the fraction of tiling SNPs homozygous in one sample which are homozygous for the alternate allele in the other sample.)

Proj_*_ALL_FPCsummary.txt
-------------------------
Matrix of all pairwise concordance values based on fingerprint genotypes. (The number represents the fraction of tiling SNPs homozygous in one sample which are homozygous for the alternate allele in the other sample.)

Proj_*_ALL_FPhet.txt
--------------------
Fraction of all tiling SNPs that are heterozygous in the given sample. Samples with >0.50 heterozygous SNPs may be contaminated with unrelated DNA. (The IMPACT and HemePACT panels contain >1,000 “tiling SNPs” where both alleles are common in the population.)

Proj_*_ALL_FPsummary.txt
------------------------
For each tiling SNP in each sample: the observed allele counts of each base, the inferred genotype, and the minor allele fraction. (The IMPACT and HemePACT panels contain >1,000 “tiling SNPs” where both alleles are common in the population.)

Proj_*_ALL_gcbias.txt
---------------------
For each sample, the average coverage for target intervals in different bins of GC content (from 25-30% to 80-85%).

Proj_*_ALL_genecoverage.txt
---------------------------
Mean coverage for each target gene (across all target exons), displayed for each sample. Coverage is computed for reads with mapping quality > 20.

Proj_*_ALL_genotypehotspotnormals.txt
-------------------------------------
Sites of known somatic mutation hotspots (COSMIC) with variants detected in a normal sample. Allele counts are also shown for the matched tumor sample as well as any other tumor where the corresponding variant was detected.

Proj_*_ALL_HSmetrics.txt
------------------------
Quality and performance metrics for hybrid selection, calculated by Picard.

Proj_*_ALL_insertsizemetrics.txt
--------------------------------
Numerical values for histogram of library insert size distribution, displayed for each sample.

Proj_*_ALL_intervalnomapqcoverage_loess.txt
-------------------------------------------
Mean coverage for each target interval, displayed for each sample. Coverage is computed for reads regardless of mapping quality (i.e., including reads mapping to multiple locations). Coverage is then normalized according to a Loess normalization to adjust for biases in GC content of target intervals.

Proj_*_ALL_intervalnomapqcoverage.txt
-------------------------------------
Mean coverage for each target interval, displayed for each sample. Coverage is computed for reads regardless of mapping quality (i.e., including reads mapping to multiple locations).

Proj_*_ALL_orgbasequalities.txt
-------------------------------
Illumina base quality score by cycle (prior to GATK base quality recalibration), displayed for each sample.


-----------------
Copy Number Files
-----------------

Proj_*_ALL_copynumber.seg
-------------------------
Segmented copy number (following %GC-based loess normalization and CBS-based segmentation) for all samples. This file can be loaded in IGV.

Proj_*_copynumber_segclusp.genes.txt
------------------------------------
Gene-level copy number (following %GC-based loess normalization and CBS-based segmentation) for all samples. Each gene is assigned the copy number ratio (fold-change) for the segment on which it falls. If a gene spans multiple segments, it is assigned the average copy number ratio. The normal sample selected for copy number normalization is shown for each sample. Both fold-change and a p-value (indicative of the signal-to-noise ratio) are given.

Proj_*_copynumber_segclusp.intragenic.txt
-----------------------------------------
List of genes with putative intragenic copy number deletions. Individual targets (exons) are listed down the right-most column. The copy number ratios for the exons in each gene are grouped into two clusters. (Cluster membership is labeled as “1” or “2”.) Genes where exons with lower copy number are consecutive are candidates for intragenic deletions, though this analysis tends to overcall events.

Proj_*_copynumber_segclusp.pdf
------------------------------
Plots of copy number profiles for all tumors. Each data point represents a target interval (blue = exon; red = tiling SNP). Each tumor is normalized against a diploid normal, either from the same project or from a historical panel of normals; the specific normal used for normalization is shown. Listed below each plot are the 20 highest-level amplifications (left) and deletions (right), though not all are significant. Significantly altered genes (fold-change > 2 and p-value < 0.05) are marked with an asterisk.

Proj_*_copynumber_segclusp.probes.txt
-------------------------------------
Target-level copy number (following %GC-based loess normalization) for all samples. The normal sample selected for copy number normalization is shown for each sample. “fc” = fold change; “lr” = log ratio.

Proj_*_discrete_CNA.txt
-----------------------
Matrix of significant copy number alterations by gene (rows are genes; columns are tumors). 2 = amplification; 0 = neutral; -2 = deletion.

Proj_*_loessnorm.pdf
--------------------
Plots showing best loess normalization fit to adjust for coverage bias related to %GC content of target intervals.

------------------------
Structural Variant Files
------------------------

Proj_*_AllAnnotatedSVs.txt
--------------------------
Annotated output from DELLY rearrangement detection algorithm. Genomic coordinates and gene annotations are provided for both breakpoints. The distance between breakpoints (SV_LENGTH) and orientation of fragments that are joined (Conection_Type) are shown. Also provided are the number of supporting paired reads and split reads and the inferred derived sequence at the breakpoint (when possible). This analysis tends to overcall events—as the vast majority of rearrangements are false positives, manual review and further filtering are necessary. Events listed (in Site2Description) as “Protein fusion: in frame” should be prioritized.

Proj_*_AllAnnotatedSVs.xlsx
----------------------------
Same information as Proj_*_AllAnnotatedCVs.txt (above) in Excel file format.

----------------
Per Sample Files
----------------

s_*_Proj_*_copynumber.seg
-------------------------
Segmented copy number for each individual sample. Each segment is listed with the genomic coordinates, number of target intervals, and copy number log ratio.
	
s_*_Proj_*_mrg_cl_aln_srt_MD_IR_BR.canonical.exon.covg.sample_interval_summary
------------------------------------------------------------------------------
Coverage statistics (total, mean, quartiles) for each exon of a canonical transcript in each individual sample. Coverage is calculated for reads with mapping quality > 20.

s_*_Proj_*_mrg_cl_aln_srt_MD_IR_BR.gene_nomapq.covg.sample_interval_summary
---------------------------------------------------------------------------
Coverage statistics (total, mean, quartiles) for each exon of any target transcript in each individual sample. Coverage is computed for reads regardless of mapping quality (i.e., including reads mapping to multiple locations).

s_*_Proj_*_mrg_cl_aln_srt_MD_IR_BR.gene.covg.sample_interval_summary
--------------------------------------------------------------------
Coverage statistics (total, mean, quartiles) for each exon of any target transcript in each individual sample. Coverage is calculated for reads with mapping quality > 20.

s_*_Proj_*_mrg_cl_aln_srt_MD_IR_BR.target.covg
----------------------------------------------
Coverage statistics (mean covearage and %GC content) for all target intervals (coding exons and tiling SNPs) in each individual sample. Coverage is calculated for reads with mapping quality > 20.

s_*_Proj_*_mrg_cl_aln_srt_MD_IR_BR.tiling_nomapq.covg.sample_interval_summary
-----------------------------------------------------------------------------
Coverage statistics (mean covearage and %GC content) for just tiling intervals (tiling SNPs and fingerprint SNPs) in each individual sample. Coverage is computed for reads regardless of mapping quality (i.e., including reads mapping to multiple locations).

s_*_Proj_*_mrg_cl_aln_srt_MD_IR_BR.tiling.covg.sample_interval_summary
----------------------------------------------------------------------
Coverage statistics (mean covearage and %GC content) for just tiling intervals (tiling SNPs and fingerprint SNPs) in each individual sample. Coverage is calculated for reads with mapping quality > 20.

----------------
Sample Info File
----------------

Proj_*_title_file.txt
---------------------
Sample IDs (provided by CMO and Investigators) and limited metadata (including sample type, input DNA amount, library yield, and capture platform) for all samples in the project.


-------------
Variant Files
-------------

annotated_exonic_variants.txt
-----------------------------
Mutations and indels (and genomic annotations and allele counts) called in all tumors in the project. The normal sample used for mutation calling is displayed in the second column (usually either the matched normal or a pool of unmatched normal). Presence in dbSNP, COSMIC, and the 1000 genomes project (depicted as the minor allele fraction in the population) is shown. For tumors called against an unmatched pool of normal, variants present in >0.01 of the 1000 genomes project (and not in COSMIC) are automatically filtered out. Allele counts across a panel of historical normal are shown to indicate potential systematic artifacts. The rightmost columns display the stats for every variant in every sample, including a panel of historical normals (usually labeled with a capital M or S). DP = depth of coverage at the variant site; RD = counts of reference allele; AD = counts of alternate (mutant) allele; VF = variant frequency.

Proj_*_AllSomaticMutIndel_withAlleleDepth_annovarAnnotated.txt
--------------------------------------------------------------
Intermediate variants file displaying all mutations and indels called prior to any filtering. Stats are included for every variant in every sample in the project (as described above).

Proj_*_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedExonic.txt
--------------------------------------------------------------------
Intermediate variants files displaying all non-silent mutations and indels called in exonic (plus splice site) regions of target genes. 

Proj_*_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedExonic.Dropped.txt
----------------------------------------------------------------------------
Intermediate variants files displaying all non-silent mutations and indels called in exonic (plus splice site) regions of target genes. The file labeled “Dropped” includes candidate mutations that were subsequently rejected based on empirical filters (e.g., present in historical normal or insufficient coverage, read support, or variant frequency). 

Proj_*_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedExonic.Filtered.txt
-----------------------------------------------------------------------------
Intermediate variants files displaying all non-silent mutations and indels called in exonic (plus splice site) regions of target genes. The file labeled “Filtered” includes candidate mutations that survived all filters.

Proj_*_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedSilent.txt
--------------------------------------------------------------------
Intermediate variants files displaying all silent mutations called in exonic (plus splice site) regions of target genes. 

Proj_*_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedSilent.Dropped.txt
----------------------------------------------------------------------------
Intermediate variants files displaying all silent mutations called in exonic (plus splice site) regions of target genes. The file labeled “Dropped” includes candidate mutations that were subsequently rejected based on empirical filters (e.g., present in historical normal or insufficient coverage, read support, or variant frequency). 

Proj_*_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedSilent.Filtered.txt
-----------------------------------------------------------------------------
Intermediate variants files displaying all silent mutations called in exonic (plus splice site) regions of target genes. The file labeled “Filtered” includes candidate mutations that survived all filters.

Proj_*_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedNonPanelExonic.txt
----------------------------------------------------------------------------
Intermediate variants files displaying all non-silent mutations and indels called in exonic (plus splice site) regions of off-target genes. 

Proj_*_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedNonPanelExonic.Dropped.txt
------------------------------------------------------------------------------------
Intermediate variants files displaying all non-silent mutations and indels called in exonic (plus splice site) regions of off-target genes. The file labeled “Dropped” includes candidate mutations that were subsequently rejected based on empirical filters (e.g., present in historical normal or insufficient coverage, read support, or variant frequency).

Proj_*_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedNonPanelExonic.Filtered.txt
-------------------------------------------------------------------------------------
Intermediate variants files displaying all non-silent mutations and indels called in exonic (plus splice site) regions of off-target genes. The file labeled “Filtered” includes candidate mutations that survived all filters.

Proj_*_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedNonPanelSilent.txt
----------------------------------------------------------------------------
Intermediate variants files displaying all silent mutations and indels called in off-target regions, including introns, intergenic regions, and off-target genes. 

Proj_*_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedNonPanelSilent.Dropped.txt
------------------------------------------------------------------------------------
Intermediate variants files displaying all silent mutations and indels called in off-target regions, including introns, intergenic regions, and off-target genes. The file labeled “Dropped” includes candidate mutations that were subsequently rejected based on empirical filters (e.g., present in historical normal or insufficient coverage, read support, or variant frequency). 

Proj_*_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedNonPanelSilent.Filtered.txt
-------------------------------------------------------------------------------------
Intermediate variants files displaying all silent mutations and indels called in off-target regions, including introns, intergenic regions, and off-target genes. The file labeled “Filtered” includes candidate mutations that survived all filters.

Proj_*_AllSomaticMutIndel_withAlleleDepth_annovarAnnotated_All_Filtered.txt
---------------------------------------------------------------------------
Aggregation of all mutations and indels called (on-target and off-target) that survived all empirical filters (in the “Filtered” files above).

Proj_*_AllSomaticMutIndel_withAlleleDepth.txt
---------------------------------------------
Candidate mutations called in all samples, prior to genomic annotation.
	
Proj_*_AllSomaticMutIndel_withAlleleDepth_mergedDNP.txt
-------------------------------------------------------
Candidate mutations called, prior to genomic annotation, with adjacent SNVs with similar coverage and allele frequencies merged into dinucleotide substitutions.

	

