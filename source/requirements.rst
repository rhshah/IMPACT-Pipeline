============
Requirements
============

Please See template.conf file in the configuration folder.

:perl: `v5.20.2 <http://perl5.git.perl.org/perl.git/tag/2c93aff028f866699beb26e5e7504e531c31b284>`_
:python: `v2.7.8 <https://www.python.org/download/releases/2.7.8/>`_
:R: `v3.1.2 <http://cran.r-project.org/src/base/R-3/R-3.1.2.tar.gz>`_

Tools Used
----------
:Somatic SNV calling: `MuTect v1.1.4 <https://github.com/broadinstitute/mutect/tree/1.1.4>`_
:Somatic INDEL calling: `SomaticIndelDetector GATK v2.3-9 <http://www.broadinstitute.org/gatk/download>`_
:Somatic INDEL calling: `PINDEL v0.2.5a7 <https://github.com/genome/pindel/tree/v0.2.5a7>`_
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

Versions
--------
Inside the version there are version that are being used for each tool. This is just for consistency in reports. 
Note that this section is just to print what version of things you are using so you can have all the dependencies with the respective versions listed here.

   
   