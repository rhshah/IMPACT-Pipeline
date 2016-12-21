#
# MSKCC_DMP_Constants.pm
#
# This module contains the constants used by the scripts in the MSKCC DMP
# development project.  Constants will be packaged by groups.
#
# Initial Version:     4/30/2013
# Latest Revision:     4/30/2013
# Author:              Aijaz Syed
# Contact:             syeda1@mskcc.org
#




#####################################################################
# Module Header                                                     #
#####################################################################

#
# Some boilerplate statements.
#
package MSKCC_DMP_Constants;

use Exporter;
@ISA = qw(Exporter);

%EXPORT_TAGS = (
        Return_values => [qw( 
                      SUCCESS
                      FAILURE
                      UNDEFINED
                      FLAG_ON
                      FLAG_OFF
                      TRUE
                      FALSE
                      )],
        Log           => [qw( 
                      LOG_LEVEL_ERROR
                      LOG_LEVEL_WARNING
                      LOG_LEVEL_INFO
                      LOG_LEVEL_DEBUG
                      )],
        Paths         => [qw(
                      BWA_INDEX_FS_LOC
                      PICARD_FS_LOC
                      PINDEL_FS_LOC
                      ANNOVAR_FS_LOC
                      VARSCAN_FS_LOC
                      PERL_SCRIPTS_FS_LOC
                      GATK_FS_LOC
                      HG19_FA_FS_LOC
                      INTERVAL_LIST_HEADER_FS_LOC
                      INTERVALS_HG19_GENE_FS_LOC
                      INTERVALS_V4_TARGETS_FS_LOC
                      INTERVALS_RD_TARGETS_FS_LOC
                      INTERVALS_RD_BAITS_FS_LOC
                      BED_RD_AMPLICONS_FS_LOC
                      VCF_DBSNP_FS_LOC
                      KNOWN_CANONICAL_ENSEMBL_FS_LOC
                      KNOWN_CANONICAL_ENSEMBL_FS_LOC
                      
                      VARSCAN_MIN_DP
                      VARSCAN_MIN_AD
                      VARSCAN_MIN_VF
                      
                      DMP_SW_PROD_FS_LOC
                      DMP_RD_PROD_FS_LOC
                      
                      DBSNP_BITSET_FS_LOC
                      )],
        Misc          => [qw(
                        HG19_FA_NAME
                        HG19_FA_DATE
                        
                        RUN_TYPE_LOCAL
                        RUN_TYPE_SGE
                        
                        ROOT_DATA
                        ROOT_SW
                      )],
        StatusFiles   => [qw(
                      ALIGNMENT_COMPLETE_FLAG
                      ALIGNMENT_IN_PROGRESS_FLAG
                      BWA_ALN_COMPLETE_FLAG
                      BWA_SAMPE_COMPLETE_FLAG
                      BWA_ALN_IN_PROGRESS_FLAG
                      BWA_SAMPE_IN_PROGRESS_FLAG
                      SORT_SAM_COMPLETE_FLAG
                      SORT_SAM_IN_PROGRESS_FLAG
                      BAM_INDEX_COMPLETE_FLAG
                      BAM_INDEX_IN_PROGRESS_FLAG
                      BAM_STATS_COMPLETE_FLAG
                      BAM_STATS_IN_PROGRESS_FLAG
                      REALIGN_COMPLETE_FLAG
                      REALIGN_IN_PROGRESS_FLAG
                      RECAL_COMPLETE_FLAG
                      RECAL_IN_PROGRESS_FLAG
                      CTPCR_METRICS_COMPLETE_FLAG
                      CTPCR_METRICS_IN_PROGRESS_FLAG
                      UNIGENO_COMPLETE_FLAG
                      UNIGENO_IN_PROGRESS_FLAG
                      PILEUP_COMPLETE_FLAG
                      PILEUP_IN_PROGRESS_FLAG
                      VARSCAN_COMPLETE_FLAG
                      VARSCAN_IN_PROGRESS_FLAG
                      PINDEL_COMPLETE_FLAG
                      PINDEL_IN_PROGRESS_FLAG
                      DINDEL_COMPLETE_FLAG
                      DINDEL_IN_PROGRESS_FLAG
                      MERGE_VCFS_IN_PROGRESS_FLAG
                      MERGE_VCFS_COMPLETE_FLAG
                      CLINICAL_REPORT_IN_PROGRESS_FLAG
                      CLINICAL_REPORT_COMPLETE_FLAG
                      GENOTYPE_MPILEUP_R_FS_LOC
                      MUTECT_IN_PROGRESS_FLAG
                      MUTECT_COMPLETE_FLAG
                      SOMATIC_IN_PROGRESS_FLAG
                      SOMATIC_COMPLETE_FLAG

                      TUMOR_ONLY
                      TUMOR_HAPMAP

                      INTERVALS_CL_RD_BAITS_FS_LOC
                      RAINDANCE_GENES_FS_LOC
                      RAINDANCE_HOTSPOT_FS_LOC

                      )],
        );


Exporter::export_tags('Return_values','Log', 'Paths', 'Misc','StatusFiles');



#####################################################################
# Included Modules                                                  #
#####################################################################

#
# Enforce good programming style.
#
use strict;
    
#####################################################################
# Exported Constant Definitions                                     #
#####################################################################


# RETURN VALUE CONSTANTS
#
# Extracted from General.pm
use constant SUCCESS => 0;
use constant FAILURE => 1;
use constant UNDEFINED => 2;

#
# Constants to represent the allowed states of flag variables.
#
use constant FLAG_OFF => 0;
use constant FLAG_ON => 1;
use constant TRUE => 1;
use constant FALSE => 0;

#
# LOG CONSTANTS
use constant LOG_LEVEL_ERROR => "ERROR";
use constant LOG_LEVEL_WARNING => "WARNING";
use constant LOG_LEVEL_INFO => "INFO";
use constant LOG_LEVEL_DEBUG => "DEBUG";

%MSKCC_DMP_Constants::ALLOWED_USERS = 
    (
        "syeda1" => (FLAG_ON),
        "bergerm1" => (FLAG_ON),
        "chengd1" => (FLAG_ON),
        "zehira" => (FLAG_ON),
        "shahr2" => (FLAG_ON),
    );
    
############################
######## TO CHANGE #########
############################
use constant ROOT_DATA => "/analysis/dmplab/prod/data/";
use constant ROOT_SW => "/analysis/dmplab/prod/tools/";


use constant PUBLIC_DATA => (ROOT_DATA)."/pubdata/";
use constant MSK_DATA => (ROOT_DATA)."/mskdata/";
use constant BWA_INDEX_FS_LOC => (PUBLIC_DATA)."/hg19/bwa/hg19.bwa.idx";
use constant PICARD_FS_LOC => (ROOT_SW)."/picard";
use constant PERL_SCRIPTS_FS_LOC => (ROOT_SW)."/perl_scripts";
use constant ANNOVAR_FS_LOC => (ROOT_SW)."/annovar";
use constant GATK_FS_LOC => (ROOT_SW);
#use constant PINDEL_FS_LOC => "/home/chengd1/Analysis/Software/pindel_0.2.4t/pindel024t";
#use constant VARSCAN_FS_LOC => "/home/chengd1/Analysis/Software/VarScan";
#use constant MUTECT_FS_LOC => "/home/shahr2/Software/Mutect/muTect-1.1.4/muTect-1.1.4.jar";

#pubdata
use constant HG19_FA_FS_LOC => (PUBLIC_DATA)."/hg19/fasta/hg19.reordered.fasta";
use constant INTERVAL_LIST_HEADER_FS_LOC => (MSK_DATA)."/interval_list.header";
use constant INTERVALS_HG19_GENE_FS_LOC => (MSK_DATA)."/v4_hg19_gene_intervals_modified.list";
use constant INTERVALS_V4_TARGETS_FS_LOC => (MSK_DATA)."/v4_hg19_picard_baits.interval_list";
use constant INTERVALS_RD_TARGETS_FS_LOC => (MSK_DATA)."/raindance_regions.canonical.sorted.interval_list";
use constant INTERVALS_RD_BAITS_FS_LOC => (MSK_DATA)."/raindance_amplicons.canonical.sorted.interval_list";
use constant BED_RD_AMPLICONS_FS_LOC => (MSK_DATA)."/raindance_amplicons.canonical.sorted.bed";
use constant VCF_DBSNP_FS_LOC => (PUBLIC_DATA)."/dbSNP/dbsnp_137.replaced.vcf";
use constant KNOWN_CANONICAL_ENSEMBL_FS_LOC => (PUBLIC_DATA)."/hg19/hg19_knownCanonical_ensembl";
#use constant EXON_NUMBER_TABLE_FS_LOC => "/home/chengd1/Analysis/Software/snpEff/snpEff_2_0_5/data/GRCh37.64/exon.number.table";
use constant DBSNP_BITSET_FS_LOC => (PUBLIC_DATA)."/hg19/hg19.dbsnp130.nocosmic47.bitset";
use constant INTERVALS_CL_RD_BAITS_FS_LOC => (MSK_DATA)."/raindanceBidirectionallyCoveredPrimersExcluded.sorted.interval_list";
use constant RAINDANCE_GENES_FS_LOC => (MSK_DATA)."/raindance_clinical_and_investigational_genes.txt";
use constant RAINDANCE_HOTSPOT_FS_LOC => (MSK_DATA)."/raindance_clinical_hotspots.txt";

use constant GENOTYPE_MPILEUP_R_FS_LOC => (ROOT_SW)."/R_scripts"."/genotype_mpileup.dmp.R ";

# varscan options
use constant VARSCAN_MIN_DP => "100";
use constant VARSCAN_MIN_AD => "5";
use constant VARSCAN_MIN_VF => "0.01";

use constant HG19_FA_NAME => "hg19";
use constant HG19_FA_DATE => "Feb2009";

use constant RUN_TYPE_LOCAL => "local";
use constant RUN_TYPE_SGE => "sge";

use constant DMP_SW_PROD_FS_LOC => "/analysis/dmplab/prod/pipeline";
use constant DMP_RD_PROD_FS_LOC => "/analysis/dmplab/prod/pipeline/rd/dmp-rd-pipeline";

#use constant DMP_SW_PROD_FS_LOC => "/home/syeda1/prod/";
#use constant DMP_RD_PROD_FS_LOC => "/home/syeda1/prod/pipeline/dmp-rd-pipeline";

use constant ALIGNMENT_COMPLETE_FLAG => ".ALIGNMENT_COMPLETE";
use constant ALIGNMENT_IN_PROGRESS_FLAG => ".ALIGNMENT_IN_PROGRESS";
use constant BWA_ALN_COMPLETE_FLAG => ".BWA_ALN_COMPLETE";
use constant BWA_SAMPE_COMPLETE_FLAG => ".BWA_SAMPE_COMPLETE";
use constant BWA_ALN_IN_PROGRESS_FLAG => ".BWA_ALN_IN_PROGRESS";
use constant BWA_SAMPE_IN_PROGRESS_FLAG => ".BWA_SAMPE_IN_PROGRESS";
use constant SORT_SAM_COMPLETE_FLAG => ".SORT_SAM_COMPLETE";
use constant SORT_SAM_IN_PROGRESS_FLAG => ".SORT_SAM_IN_PROGRESS";
use constant BAM_INDEX_COMPLETE_FLAG => ".BAM_INDEX_COMPLETE";
use constant BAM_INDEX_IN_PROGRESS_FLAG => ".BAM_INDEX_IN_PROGRESS";
use constant BAM_STATS_COMPLETE_FLAG => ".BAM_STATS_COMPLETE";
use constant BAM_STATS_IN_PROGRESS_FLAG => ".BAM_STATS_IN_PROGRESS";
use constant REALIGN_COMPLETE_FLAG => ".REALIGN_COMPLETE";
use constant REALIGN_IN_PROGRESS_FLAG => ".REALIGN_IN_PROGRESS";
use constant RECAL_COMPLETE_FLAG => ".RECAL_COMPLETE";
use constant RECAL_IN_PROGRESS_FLAG => ".RECAL_IN_PROGRESS";
use constant CTPCR_METRICS_COMPLETE_FLAG => ".CTPCR_METRICS_COMPLETE";
use constant CTPCR_METRICS_IN_PROGRESS_FLAG => ".CTPCR_METRICS_IN_PROGRESS";
use constant UNIGENO_COMPLETE_FLAG => ".UNIGENO_COMPLETE";
use constant UNIGENO_IN_PROGRESS_FLAG => ".UNIGENO_IN_PROGRESS";
use constant PILEUP_COMPLETE_FLAG => ".PILEUP_COMPLETE";
use constant PILEUP_IN_PROGRESS_FLAG => ".PILEUP_IN_PROGRESS";
use constant VARSCAN_COMPLETE_FLAG => ".VARSCAN_COMPLETE";
use constant VARSCAN_IN_PROGRESS_FLAG => ".VARSCAN_IN_PROGRESS";
use constant PINDEL_COMPLETE_FLAG => ".PINDEL_COMPLETE";
use constant PINDEL_IN_PROGRESS_FLAG => ".PINDEL_IN_PROGRESS";
use constant DINDEL_COMPLETE_FLAG => ".DINDEL_COMPLETE";
use constant DINDEL_IN_PROGRESS_FLAG => ".DINDEL_IN_PROGRESS";
use constant MUTECT_IN_PROGRESS_FLAG => ".MUTECT_IN_PROGRESS";
use constant MUTECT_COMPLETE_FLAG => ".MUTECT_COMPLETE";
use constant SOMATIC_IN_PROGRESS_FLAG => ".SOMATIC_IN_PROGRESS";
use constant SOMATIC_COMPLETE_FLAG => ".SOMATIC_COMPLETE";
use constant MERGE_VCFS_IN_PROGRESS_FLAG => ".MERGE_FILTER_VCFS_IN_PROGRESS";
use constant MERGE_VCFS_COMPLETE_FLAG => ".MERGE_FILTER_VCFS_COMPLETE";
use constant CLINICAL_REPORT_IN_PROGRESS_FLAG => ".CLINICAL_REPORT_IN_PROGRESS";
use constant CLINICAL_REPORT_COMPLETE_FLAG => ".CLINICAL_REPORT_COMPLETE";

use constant TUMOR_ONLY => "tumor-only";
use constant TUMOR_HAPMAP => "tumor-hapmap"; 

#
# A seemingly gratuitous statement, included so that the interpreter
# doesn't complain.
1;
