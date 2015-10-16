============
Requirements
============

.. sidebar:: Note:

	 This will only run on any Linux system that has an SGE or LSF cluster.

Please see template.conf file in the configuration folder.

:perl: `v5.20.2 <http://perl5.git.perl.org/perl.git/tag/2c93aff028f866699beb26e5e7504e531c31b284>`_
:python: `v2.7.8 <https://www.python.org/download/releases/2.7.8/>`_
:R: `v3.1.2 <http://cran.r-project.org/src/base/R-3/R-3.1.2.tar.gz>`_


Purpose and Tools Used
======================

.. sidebar:: Note:

	For Timmomatic we are using a custom old version which uses `cutadapt v1.1 <https://cutadapt.readthedocs.org>`_ internally. To get this please contact us.
	
:Trimming: `Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_
:Alignment: `BWA v0.7.5a <https://github.com/lh3/bwa/tree/0.7.5a>`_
:Somatic SNV calling: `MuTect v1.1.4 <https://github.com/broadinstitute/mutect/tree/1.1.4>`_
:Somatic INDEL calling: `SomaticIndelDetector in GATK v2.3-9 <http://www.broadinstitute.org/gatk/download>`_
:Somatic INDEL calling: `PINDEL v0.2.5a7 <https://github.com/genome/pindel/tree/v0.2.5a7>`_
:Indel Realignment: `ABRA v0.92 <https://github.com/mozack/abra/tree/v0.92>`_
:Mark Duplicates and Various Statistics: `Picard Tools v1.96 <https://github.com/broadinstitute/picard/tree/1.96>`_
:Base Quality Recalibration and Find Covered Intervals: `GATK v3.3.0 <http://www.broadinstitute.org/gatk/download>`_
:Genotyping Position: `SAMTOOLS v0.1.19 <https://github.com/samtools/samtools/tree/0.1.19>`_
:Somatic Structural Variant Framework: `IMPACT-SV v1.0.1 <https://github.com/rhshah/IMPACT-SV/tree/1.0.1>`_


Inside the config file
======================

There are three sections:

+-----------+-----------+-----------+
| Section 1 | Section 2 | Section 3 |
+===========+===========+===========+
| Locations | Parameters| Versions  |
+-----------+-----------+-----------+

All of this section start with ``>`` sign.


Inside each of the section here are the things that need to be set:

Locations
---------

:ZCAT: Location of the ``zcat`` program on linux 
:TMPDIR: Set the temporary directory for all tools please set somthing other then ``/tmp``
:JAVA_1_6: Set JAVA version 1.6
:JAVA_1_7: Set JAVA version 1.7
:GATK_SomaticIndel: Path to GATK somatic indel detector (GATK version 2.3-9)
:GATK: Path to GATK (GATK version 3.3.0)
:Reference: Path to fasta referece file to be used (GRCh37)
:Refseq: Path to refgene file to be used
:PICARD: Path to picard tools (Picard version 1.19)
:Mutect: Path to MuTect (MuTect version 1.1.4)
:BaitInterval: Bait file to be used for analysis. Please make the interval file based on Picard's HSmetrics tool format. 
:TargetInterval: Target file to be used for analysis. Please make the interval file based on Picard's HSmetrics tool format. 
:ABRA: Path to ABRA (ABRA version 0.92) 
:TargetRegionLIST: Target File in GATK's List format.
:PINDELBIN: Path to pindel installation (Pindel version 0.2.5a7)
:SVpipeline: Path to structural variant framework (only required if running the structural variant detection framework)
:CAT: Location of ``cat`` program on Linux 
:PYTHON: Location of ``python`` program on Linux (Python version 2.7.8)
:TrimGalore: Path to trimgalore tool
:PERL: Location of ``perl`` program on Linux (Perl version 5.20.2)
:BWA: Path to bwa tool
:GeneInterval: Gene interval file 
:GeneIntervalAnn: Gene interval annotated file
:GeneCoord: Path to Gene Coordinate file
:TilingInterval: = Path to tiling intervals
:TilingIntervalAnn: Path to tiling intervals - annotated for cytoband, for copy number
:FingerPrintInterval: Path to FingerPrint Interval file
:dbSNP: Path to db snp vcf file
:COSMIC: Path to cosmic vcf (version 0.68)
:Mills_1000G_Indels: Path to Mills 1000G Indels
:dbSNP_bitset: Path to dbsnp bitset file
:AnnotateAssessFilterVariants: Path to Annotate Assess and Filter variants script
:LoessNormalization: Path ot Loess Normalization for copynumber
:GCBiasFile: Path to GCbias file for copy number
:HistNormDir: Path to Historiacal Normal dir for Copy number
:BestCopyNumber: Path to Copy number script
:NormVsNormCopyNumber: Path to Normal vs. Normal Copy number script
:StdNormalLoess_TM: Standard Normals for copy number analysis - FFPE for tumor samples#
:StdNormalLoess_NVN: Standard Normals for copy number normal vs normal analysis
:AllMetrics: Path to all metrics R script 
:SAMTOOLS: Path to samtools
:BEDTOOLS: Path to bedtools
:GenotypeAllele: Path to Genotype allele script
:CosmicHotspotVcf: Path to cosmic hotspot vcf
:Annovar: Path to Annovar script
:Annovar_db: Path to Annovar DB
:Canonical_refFlat_file: Path to canonical reflat file
:IGVtools: Path to IGV tools
:TranslationFolder: Path to translation folder
:HotSpot_mutations: Path to hotspot mutations for 2 tiered filtering
:clinicalExons: ListOfClinicalExon 
:Validated_Exons: File with List Of Clinically Validated Exons
:Tumor_supressor_list: Path to list of tumor supressor genes 
:Canonical_Exon_Interval_table_with_aa: Path to exon interval table 
:Canonical_Exon_Interval_list: Path to canonical exon interval table for DoC
:NormalVariantsVCF: Path to compiled variants found in mixed normals
:QSUB: Path to qsub for SGE
:BSUB: Path to bsub for LSF
:RHOME: Path to R bin directory
:RLIBS: Path to R library directory
:RSYNC: Path to ``rsyn`` on system 
:BarcodeKey: Path to barcode key file
:AdaptorKey: Path to adaptor key file
:StandardNormalsDirectory: Directory where the standard normals are stored

Parameters
----------

Set the parameters to different file/folders/values required by the IMPACT pipeline

:StdNormalForMutationCalling: Path to standard normal to be used for mutation calling
:ListOfFiles: File of Files(FOF) for different steps for the pipeline (only required when the process dont start from merging fastq)
:Process: Which process to run the pipeline on ( can be 1,2,3,4,5,6,7 independently or continuous combination in ascending order )
:FastqSource: Where are the fastq file from (can be ``GCL`` or ``DMP``)
:MAPQ: Mapping Quality Threshold (Used by DMP-IMPACT:0.2)0
:BASQ: Base Quality Threshold (Used by DMP-IMPACT:0.2)
:MergeDinucleotide: Flag to Merge di-nucleotide mutation(can be 1(True) or 2(False))
:MoveFiles: Flag to Move file in folders (can be 1(True) or 2(False))
:DeleteIntermediateFiles: Flag ti Delete Intermediate Files (can be 1(True) or 2(False))
:TNfreqRatio_MutectStdFilter: TN freq Ratio for mutect std filter (Used by DMP-IMPACT:5)
:TNfreqRatio_SomIndelStdFilter: TN freq Ratio for SID std filter (Used by DMP-IMPACT:5)
:VF_threshold_hotspot: Variant Frequency threshold for SNV hotspot (Used by DMP-IMPACT:0.01)
:AD_SomIndelSTDFilter: Allele Depth Threshold for SID standard filter (Used by DMP-IMPACT:5)
:DP_SomIndelSTDFilter: Total Depth Threshold for SID standard filter (Used by DMP-IMPACT:0)
:VF_SomIndelSTDilter: Variant Frequency Threshold for SID standard filter (Used by DMP-IMPACT:0.01)
:AD_MutectSTDFilter: Allele Depth Threshold for Mutect standard filter (Used by DMP-IMPACT:5)
:DP_MutectSTDFilter: Total Depth Threshold for Mutect standard filter (Used by DMP-IMPACT:0)
:VF_MutectSTDFilter: Variant Frequency Threshold for Mutect standard filter (Used by DMP-IMPACT:0.01)
:TNfreqRatio_AnnotationFilter: Tumor to Normal frequency ratio therehold for Annotation (Used by DMP-IMPACT:5)
:PON_AD_Threshold: Panel of Normal Allele Depth Threshold (Used by DMP-IMPACT:3)
:PON_TPVF_Threshold: Panel of Normal TPVF Threshold (Used by DMP-IMPACT:10)
:Pindel_Min_Indel_Len: Minimum Length of INDEL called by PINDEL(Used by DMP-IMPACT:25)
:Pindel_Max_Indel_Len: Maximum Length of INDEL called by PINDEL (Used by DMP-IMPACT:2000)
:MAFthreshold_AnnotationFilter: Maf threshold for Annotation (Used by DMP-IMPACT:0.01)
:minimumDPforSNV: Minimum Total Depth for Novel SNVs  (Used by DMP-IMPACT:20)
:minimumADforSNV: Minimum Allele Depth for Novel SNVs (Used by DMP-IMPACT:10)
:minimumVFforSNV: Minimum Variant Frequency for Novel SNVs (Used by DMP-IMPACT:0.05)
:minimumDPforSNVhs: Minimum Total Depth for Hotspot SNVs (Used by DMP-IMPACT:20)
:minimumADforSNVhs: Minimum Allele Depth for Hotspot SNVs (Used by DMP-IMPACT:8)
:minimumVFforSNVhs: Minimum Variant Frequency for Hotspot SNVs (Used by DMP-IMPACT:0.02)
:minimumDPforINDEL: Minimum Total Depth for Novel INDELs (Used by DMP-IMPACT:20)
:minimumADforINDEL: Minimum Allele Depth for Novel INDELs (Used by DMP-IMPACT:10)
:minimumVFforINDEL: Minimum Variant Frequency for Novel INDELs (Used by DMP-IMPACT:0.05)
:minimumDPforINDELhs: Minimum Total Depth for Hotspot INDELs (Used by DMP-IMPACT:20)
:minimumADforINDELhs: Minimum Allele Depth for Hotspot INDELs (Used by DMP-IMPACT:8)
:minimumVFforINDELhs: Minimum Variant Frequnecy for Hotspot INDELs (Used by DMP-IMPACT:0.02)
:occurrencePercent: Minimum Percentage For Occurrence In Other Normals (Used by DMP-IMPACT:0.2)
:Coverage_threshold_darwin_report: Coverage threshold for darwin reports(good coverage vs bad coverage) (Used by DMP-IMPACT:100)
:QUEUE_NAME: Name of the queue on the SGE or LSF
:CLUSTER: Flag for what cluster to be used (can ``SGE`` or ``LSF``)
:runABRA: Flag to whether use ABRA or GATK indel realignment(can be 1(True) or 2(False))

Versions
--------

.. sidebar:: Note: 

	This section is just to print what version of things you are using so you can have all the dependencies with the respective versions listed here.

Inside the version there are version that are being used for each tool. This is just for consistency in reports. 


Description for title_file.txt
==============================

Headers for this tab-delimited file should be exactly with this names:

:Barcode: Has to start with bc and end with any number [for example: bc01 or bc101 should match the **adaptor & barcode** file mentioned in configuration file
:Pool:	Can be any string **joined by ``-``** and **not ``_``** and all entries should be from same pool
:Sample_ID:	Can be any string **joined by ``-``** and **not ``_``** 
:Collab_ID: Can be any string or ``-``

.. sidebar:: Note: 

	Patient with multiple samples should have **same Patient_ID**
	
:Patient_ID: Can be any string **joined by ``-``** and **not ``_``** 
:Class: Can be Tumor or Normal.
:Sample_type: Can be any string or ``-``
:Input_ng: Can be any float or ``-``
:Library_yield:	Can be float or ``-``
:Pool_input: Can be float or ``-``
:Bait_version: Can be any string or ``-``
:Gender: Can be any Male/Female or ``-``
:PatientName: Can be any string or ``-``
:MAccession: Can be any string or ``-``
:Extracted_DNA_Yield: Can be a float or ``-``

For analysis to start the **outputDirectory** will be required to have this file with ``title_file.txt`` as the name or this file needs to be present in the **configuration** file with either ``title_file.txt`` as then name or ``Pool_title.txt`` as the name where **Pool** is the string used above for that category.

Description for SampleSheet.csv
===============================

This is a comma separated file is created by the illumina sequencer and it is used to merge the fastq files. 

Headers for this tab-delimited file should be exactly with this names:

:FCID: Flowcell ID (required)
:Lane: Lane Number, this is used to merge the fastq files across lanes (required)
:SampleID: Sample ID, this is used to merge the files (required)
:SampleRef: Sample Reference is from [example:HUMAN]
:Index: Index used to sequence the sample (require)
:Description: Description of the samples
:Control: Can be any string or ``-``
:Recipe: Can be any string or ``-``
:Operator: Can be any string or ``-``
:SampleProject: Can be any string or ``-``

For analysis to start the **outputDirectory** will be required to have this file with ``SampleSheet.csv`` as the name or this file needs to be present in the **configuration** file with ``SampleSheet.csv`` as the name.


Description for adaptor file in the configuration file
======================================================

The adaptor file is the tab-delimited file with two columns:

1. Barcode Key to which the adaptor belongs which should always start with ``bc``

2. Adaptor sequence itself

There is **no header** in this file.

For Example:

	+-------+-----------------------------------------------------------------------+
	| bc01  |     GATCGGAAGAGCACACGTCTGAACTCCAGTCACAACGTGATATCTCGTATGCCGTCTTCTGCTTG |
	+-------+-----------------------------------------------------------------------+
	
	
Description for barcode file in the configuration file
======================================================

The barcode file is the tab-delimited file with two columns:

1. Barcode Sequece

2. Barcode Number that sequence represent.

There is **a header** in this file.

For Example:
	+---------+--------------+
	|Sequence | TruSeqBarcode|
	+=========+==============+
	|AACGTGAT |       bc01   |
	+---------+--------------+
