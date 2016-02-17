#!/usr/bin/perl -w
##########RunIlluminaProcess.pl########
#Author: Ronak Shah(RS),Donavan Cheng(DC),Ahmet Zehir(AZ), Aijazuddin Syed(AS)
#Date: 11/08/2012
#LastModified: 08/27/2015
#Version: 3.0
#Description: Get In the data and run the mapping and mutation calling process.
#############
###Log:
##01/23/2013
#v1.1
#RS:Seprate Calculate & Compile metrics
##01/27/2013
#v1.2
#RS:Activate Standard Normal Option
##01/29/2013
#v1.2
#RS:Added sorting option for all places where and entry is a list of BAM file
#Updated the House Keeping module
#Added to change the permissions of the files at end.
#Added Deleting option after each process
#Made HouseKeeping properly active
##01/30/2013
#v1.3
#RS:Activate the option of annotation and running mutation assesor
##02/04/2013
#v1.3
#RS:Seprate First Filter of SNPs and Indel and Making of merged Annotation file.
##02/05/2013
#v.1.3
#RS:Commented Annotation Calling. Problem of SGE.
##03/07/2013
#v1.3
#RS:Removed the bug in individual bed files
##03/25/2013
#1.4
#RS:Updating to accomadate vcf filtering
#Merging of VCF & calling Allele Depth using Unified Genotyper
#Merging results for all samples
##04/04/2013
#v1.5
#DC:Added DMP support
#RS:Using mpileup to genotype the variants.
##04/09/2013
#v1.5
#RS:Merging of genotyed records.
#Updated the deletion of files.Changed it to deletion in their specific process
##04/12/2013
#v1.5
#RS:Updated the ouputs for Annotation.
##04/22/2013
#v1.5
#RS:Updated the ouputs for Annotation,Genotype,And criteria for merging
#dinucleotide.
##04/30/2013
#v1.5
#RS:Updated the total normal frequency based on median
#RS:Updated linking of stdNormals now it copies them.
#RS:Removed the merge fastq data.
##07/24/2013
#v1.6
#RS:Updated Barcode for adaptor sequences
#DC: Added funtionality for genotyping hotspot mutation.
##07/29/2013
#v1.6
#RS:Updated the qsub commands for logs.
#AZ:Adding annovar functionality and removing oncotator.
##08/06/2013
#v1.6
#RS:Updated the process retrival & multiple other fixes.
##08/09/2013
#v1.7
#RS:Updated SGE command for SABA2.
#v2.0.3
#Added Pindel
#RS:Changed the PVF filters
#v2.0.4
#RS:Added Call method to output
#Added Merging of all Filter Variants
#v2.0.5
#RS:Added support to skip indels when we dont have VCF files created.
#v2.0.6
#RS:Added intervals back for SID running.
#v2.0.7
#RS:Added Find Covered Intervals from GATK for ABRA & SID running.
#v2.0.8
#RS:Modified the SID to take new intervals. Removed TargeRegionBed TargetRegionList Variables.
#RS:Made changes the way ABRA will be ran.
#RS:Added code that will run FindCoveredIntervals from GATK.
#RS:Modified configuration files to reflect ABRA version 0.92. GATK version 3.3-0
#RS:Modified Parameters for Genotype the normal.
#v2.09
#RS:Added option to keep the symlinks by default; Use -sf 1 to remove them.
#RS:Made SID compatible with older samples as we need to have a target file in SID.
#RS:Modified the Configuration files to have a traget file for SID.
#v2.10
#RS: #RS: Modified the parameters for GenotypeHotspotNormal to vf>=0.2 or ad>=5 previously it was only for Normals now added for Tumors also
##08/25/2015
#v3.0
#RS: Adding funtionality to run on LSF
#############################
use strict;
use Getopt::Long;
use IO::File;
use List::Util qw(max sum);
use Tie::IxHash;
use Vcf;
use File::Basename;
use MSKCC_DMP_Logger;
use Data::Dumper;
use FindBin::Real qw(Bin);
use version;
my $logger = MSKCC_DMP_Logger->get_logger('IMPACT_Pipeline_Logger');
$logger->start_local();
my ( $config_file, $sv_config_file, $datadir, $outdir, $svTag,
	$removeSymLinksFlag );
undef $!;
undef $@;

#--This variable holds the current time
my $now = time;

# This is variable is the path to the bin folder
my $binPath = Bin();
my ($root) = $binPath =~ /(.*\/).*/;    #get the path above bin
my ($revision) = `svnversion $root`;
$revision =~ s/\s//g;
my $supportScriptPath = $root
  . "support-scripts/"
  ;    #path to support scripts, independent of installation location
if (
	@ARGV < 1
	or !GetOptions(
		'config|c:s'          => \$config_file,
		'svConfig|sc:s'       => \$sv_config_file,
		'dataDirectory|d:s'   => \$datadir,
		'outputDirectory|o:s' => \$outdir,
		'symLinkFlag|sf:i'    => \$removeSymLinksFlag,
		'version|v:s'         => \&VersionMessage
	)
  )
{
	Usage();
}
if ($config_file) {
	$logger->info("The configration file in use is $config_file");
}
else {
	$logger->fatal("The configration file in use is $config_file");
	exit(1);
}
if ($sv_config_file) {
	$logger->info("The sv configration file in use is $config_file");
	$svTag = 2;
}
else {
	$logger->fatal("The sv configration file in use is $config_file");
	$svTag = 1;
}
if ($removeSymLinksFlag) {
	if ( $removeSymLinksFlag == 1 ) {
		$logger->info("SYMLINKS:The symLinks will be removed");
	}
	else {
		$logger->info("SYMLINKS:The symLinks will be not be removed");
	}
}
else {
	$removeSymLinksFlag = 2;
	$logger->info("SYMLINKS:The symLinks will not be removed");
}

#Define Global Variables
my (
	$sampleFile,                  $stdNormal,
	$titleFile,                   $fof,
	$runBQSR,                     $mailID,
	$list,                        $poolName,
	$projectName,                 $barcodeFile,
	$process,                     $standardNormalList,
	$mvFiles,                     $mergeDinucleotide,
	$fastqSource,                 $adaptorFile,
	$TMPDIR,                      $JAVA_1_6,
	$JAVA_1_7,                    $ExonToGenCov,
	$FPGenotypesScript,           $GenerateMetricsFilesScript,
	$FP_genotypes,                $GATK_SomaticIndel,
	$GATK,                        $Reference,
	$Refseq,                      $PICARD,
	$Mutect,                      $filter_Mutect,
	$filter_SomaticIndel,         $BaitInterval,
	$TargetInterval,              $CompileMetrics,
	$CAT,                         $PYTHON,
	$TrimGalore,                  $PERL,
	$BWA,                         $GeneInterval,
	$GeneIntervalAnn,             $GeneCoord,
	$TilingInterval,              $TilingIntervalAnn,
	$FingerPrintInterval,         $canonicalExonIntervalsFile,
	$dbSNP,                       $COSMIC,
	$Mills_1000G_Indels,          $dbSNP_bitset,
	$dbProperties,                $Oncotator,
	$Mutation_Assessor,           $AnnotateAssessFilterVariants,
	$LoessNormalization,          $BestCopyNumber,
	$NormVsNormCopyNumber,        $StdNormalLoess_TM,
	$StdNormalLoess_NVN,          $AllMetrics,
	$SAMTOOLS,                    $BEDTOOLS,
	$GenotypeAllele,              $cosmicHotspotsVcf,
	$GCBiasFile,                  $HistNormDir,
	$TNfreqRatio_MutectStdFilter, $TNfreqRatio_SomIndelStdFilter,
	$TNfreqRatioThreshold,        $ad_SomIndelStdFilter,
	$dp_SomIndelStdFilter,        $vf_SomIndelStdFilter,
	$ad_MutectStdFilter,          $dp_MutectStdFilter,
	$vf_MutectStdFilter,          $dp_snv,
	$ad_snv,                      $vf_snv,
	$dp_snvHS,                    $ad_snvHS,
	$vf_snvHS,                    $dp_indel,
	$ad_indel,                    $vf_indel,
	$dp_indelHS,                  $ad_indelHS,
	$vf_indelHS,                  $MAFthreshold,
	$occurrencePercent,           $NormalVariantsVCF,
	$queue,                       $deleteIntermediateFiles,
	$MAPQ,                        $BASQ,
	$QSUB,                        $RHOME,
	$RLIBS,                       $RSYNC,
	$filterGenotypedVariants,     $createPatientVCFfile,
	$createCoverageFiles,         $ClinicalExons,
	$HotspotMutations,            $IGVtools,
	$translationFolder,           $exonIntervalsFile,
	$validatedExons,              $coverageThreshold,
	$SVpipeline,                  $ABRA,
	$ABRAwrapper,                 $StandardNormalsDirectory,
	$runABRA,                     $CountErrors,
	$PindelWrapper,               $PINDELBIN,
	$filter_Pindel,               $PON_AD_threshold,
	$PON_TPVF_threshold,          $Pindel_Min_Indel_Len,
	$Pindel_Max_Indel_Len,        $triallelic_fix,
	$FindCoveredIntervals,        $TargetRegionLIST,
	$CLUSTER,                     $BSUB,
	$ZCAT
);

sub VersionMessage {
	our $VERSION = version->declare("v3.0.0");
	print "You are running $0 version : $VERSION\n";
	exit();
}

#Get Configration File details
my ($Version) = &GetConfiguration($config_file);

#Check the input parameters
#Check for Sample Sheet & Title File
if ( !$sampleFile ) {
	$logger->warn(
"Sample sheet is not given in the configuration file. Looking into the raw data directory."
	);
	$sampleFile = `ls $datadir/SampleSheet.csv`;
	chomp($sampleFile);
	if ( !$sampleFile ) {
		$logger->fatal(
"Sample sheet could not be located in the raw data directory, pipeline requires a sample sheet to run. Please see usage."
		);
		Usage();
		exit(1);
	}
}
if ( !$titleFile ) {
	$logger->warn(
"Title file is not given in the configuration file. Looking into the raw data directory."
	);
	$titleFile = `ls $datadir/title_file.txt`;
	chomp($titleFile);
	if ( !$titleFile ) {
		$logger->fatal(
"Title file could not be located in the raw data directory, pipeline requires the title file to run. Please see usage."
		);
		Usage();
		exit(1);
	}
}

#Read Sample File
my (
	$fcId,        $lane,    $sampleId, $sampleRef, $index,
	$description, $control, $recipe,   $operator,  $sampleProject
) = &ReadSampleFile( $sampleFile, $projectName, $outdir );

#Read Title File
my (
	$barcode,      $pool,      $titleSampleId, $collabId,
	$patientId,    $class,     $sampleType,    $inputNg,
	$libraryYeild, $poolInput, $baitVersion
) = &ReadTitleFile( $titleFile, $outdir );
my @sampleProject = @$sampleProject;
my @pool          = @$pool;

#Check for Project Name
if ( !$projectName ) {
	$logger->warn(
"Project name is not given in the configuration file, pipeline will use the project name stated in the SampleSheet.csv file."
	);
	$projectName = $sampleProject[0];
	if ( !$projectName ) {
		$logger->fatal("Can not read project name from the sample sheet.");
		exit(1);
	}
}

#Check for Pool Name
if ( !$poolName ) {
	$logger->warn(
"Pool name is not given in the configuration file, pipeline will use the pool name stated in the title file."
	);
	$poolName = $pool[0];
	if ( !$poolName ) {
		$logger->fatal("Can not read pool name from the title file.");
		exit(1);
	}
}
my $outdir_final;

#Check output directory
if ( !$outdir ) {
	$logger->fatal(
"Please enter the directory that can be used as output directory.See Usage"
	);
	Usage();
	exit(1);
}

#Check Cluster
if ( $CLUSTER eq "SGE" ) {
	if ( !$QSUB ) {
		$logger->fatal(
"Please enter the path of QSUB executable for SGE job submission. See Usage"
		);
		Usage();
		exit(1);
	}
	else {
		$logger->info("Pipeline will try to run on the given SGE cluster");
	}
}
else {
	if ( !$BSUB ) {
		$logger->fatal(
"Please enter the path of BSUB executable for LSF job submission. See Usage"
		);
		Usage();
		exit(1);
	}
	else {
		$logger->info("Pipeline will try to run on the given LSF cluster");
	}
}

#Check for runBQSR flag
if ( !$runBQSR ) {
	$logger->info("BQSR would be ran on the project but at sample level.");
	$runBQSR = 1;
}
else {
	if ( $runBQSR == 0 ) {
		$logger->info("BQSR would not be ran on the project.");
	}
	if ( $runBQSR == 1 ) {
		$logger->info("BQSR would be ran on the project but at sample level.");
	}
	if ( $runBQSR == 2 ) {
		$logger->info("BQSR would be ran on the project but at lane level.");
	}
}

#Check for Merge Dinucleotide flag
if ( !$mergeDinucleotide ) {
	$logger->info("Dinucleotide would be merged by default.");
	$mergeDinucleotide = 1;
}
else {
	if ( $mergeDinucleotide == 1 ) {
		$logger->info("Dinucleotide would be merged by default.");
	}
	if ( $mergeDinucleotide == 2 ) {
		$logger->info("Dinucleotide would not be merged.");
	}
}

#Check for move files flag
if ( !$mvFiles ) {
	$logger->info("Folders will be created and Files will be moved.");
	$mvFiles = 1;
}
else {
	if ( $mvFiles == 1 ) {
		$logger->info("Folders will be created and Files will be moved.");
	}
	if ( $mvFiles == 2 ) {
		$logger->info(
			"Folders will not be created and Files will not be moved.");
	}
}

#Check for delete files flag
if ( !$deleteIntermediateFiles ) {
	$logger->info(
"DeleteIntermediateFiles: Flag not given. Default will delete all intermediate files"
	);
	$deleteIntermediateFiles = 1;
}
else {
	if ( $deleteIntermediateFiles == 1 ) {
		$logger->info(
			"DeleteIntermediateFiles: Flag given to delete intermediate files");
	}
	if ( $deleteIntermediateFiles == 2 ) {
		$logger->info(
"DeleteIntermediateFiles: Flag given so that intermediate files are not deleted."
		);
	}
}

#Check for runABRA flag
if ( !$runABRA ) {
	$logger->info("GATK Re-alignment will be ran");
	$runABRA = 2;
}
else {
	if ( $runABRA == 1 ) {
		$logger->info("ABRA Re-alignment will be ran");
	}
	if ( $runABRA == 2 ) {
		$logger->info("GATK Re-alignment will be ran");
	}
}

#Check for fastq source
if ( !$fastqSource ) {
	$logger->info("Assuming that Fastq files are from GCL");
	$fastqSource = "GCL";
}
else {
	if ( $fastqSource eq "DMP" ) {
		$logger->info("Fastq files are from DMP");
	}
	elsif ( $fastqSource eq "GCL" ) {
		$logger->info("Fastq files are from GCL");
	}
	else {
		$logger->fatal("Please indicate fastqSource. See Usage");
		Usage();
		exit(1);
	}
}

#Check for process
$process =~ s/\s//g;
if ( !$process ) {
	$logger->fatal("Please enter the number of process to be ran. See Usage");
	Usage();
	exit(1);
}
else {
	$logger->info("The process selected to run are $process");
}

#Check for barcode file
if ($barcodeFile) {
	$logger->info("The barcode file in use is $barcodeFile.");
}
else {
	$barcodeFile = "/home/shahr2/Scripts/All/barcodeKey48.txt";
	$logger->info("The barcode file in use is $barcodeFile.");
}

#Check for adaptor file
if ($adaptorFile) {
	$logger->info("The adaptor file in use is $adaptorFile.");
}
else {
	$adaptorFile = "/home/shahr2/Scripts/All/adaptorKey48.txt";
	$logger->info("The adaptor file in use is $adaptorFile.");
}

#Check for Standard Normal file
if ( !$stdNormal ) {
	$logger->info(
"No standard normal is given. Thus assuming there is atleast one normal in the run"
	);
	$stdNormal = "NULL";
}
else {
	$logger->info("Starndard Normal is given: $stdNormal");
}

#Check for TNratio for mutect std filter
if ( !$TNfreqRatio_MutectStdFilter ) {
	$logger->info(
"MutectStdFilter:Minimum Tumor Normal VF ratio not given default will be used"
	);
	$TNfreqRatio_MutectStdFilter = 5;
}
else {
	$logger->info("TNfreqRatio_MutectStdFilter:$TNfreqRatio_MutectStdFilter");
}

#Check DP for Mutect std filter
if ( !$dp_MutectStdFilter ) {
	$logger->info(
		"MutectStdFilter:Minimum total depth not given default will be used");
	$dp_MutectStdFilter = 0;
}
else {
	$logger->info("DP_MutectStdFilter:$dp_MutectStdFilter");
}

#Check AD for Mutect std filter
if ( !$ad_MutectStdFilter ) {
	$logger->info(
		"MutectStdFilter:Minimum Allele depth not given default will be used");
	$ad_MutectStdFilter = 3;
}
else {
	$logger->info("AD_MutectStdFilter:$ad_MutectStdFilter");
}

#Check VF for Mutect std filter
if ( !$vf_MutectStdFilter ) {
	$logger->info(
"MutectStdFilter:Minimum variant frequenct not given default will be used"
	);
	$vf_MutectStdFilter = 0.01;
}
else {
	$logger->info("VF_MutectStdFilter:$vf_MutectStdFilter");
}

#Check TNratio for somatic indel detector std filter
if ( !$TNfreqRatio_SomIndelStdFilter ) {
	$logger->info(
"SomIndelStdFilter:Minimum Tumor Normal VF ratio not given default will be used"
	);
	$TNfreqRatio_SomIndelStdFilter = 5;
}
else {
	$logger->info(
		"TNfreqRatio_SomIndelStdFilter:$TNfreqRatio_SomIndelStdFilter");
}

#Check DP for somatic indel detector std filter
if ( !$dp_SomIndelStdFilter ) {
	$logger->info(
		"SomIndelStdFilter:Minimum total depth not given default will be used");
	$dp_SomIndelStdFilter = 0;
}
else {
	$logger->info("DP_SomIndelStdFilter:$dp_SomIndelStdFilter");
}

#Check AD for somatic indel detector std filter
if ( !$ad_SomIndelStdFilter ) {
	$logger->info(
		"SomIndelStdFilter:Minimum Allele depth not given default will be used"
	);
	$ad_SomIndelStdFilter = 3;
}
else {
	$logger->info("AD_SomIndelStdFilter:$ad_SomIndelStdFilter");
}

#Check VF for somatic indel detector std filter
if ( !$vf_SomIndelStdFilter ) {
	$logger->info(
"SomIndelStdFilter:Minimum variant frequenct not given default will be used"
	);
	$vf_SomIndelStdFilter = 0.01;
}
else {
	$logger->info("VF_SomIndelStdFilter:$vf_SomIndelStdFilter");
}

#Check Base Quality threshold
if ( !$BASQ ) {
	$logger->info("BASQ: Base Quality Threshold not given defult will be used");
	$BASQ = 5;
}
else {
	$logger->info("BASQ:$BASQ");
}

#Check Mapping Quality threshold
if ( !$MAPQ ) {
	$logger->info(
		"MAPQ: Mapping Quality Threshold not given defult will be used");
	$MAPQ = 5;
}
else {
	$logger->info("MAPQ:$MAPQ");
}

#Check  for queue
if ( !$queue ) {
	$logger->info("Name of the SGE queue not given default will be used");
	$queue = "all.q";
}
else {
	$logger->info("SGE Queue:$queue");
}

#Check PON AD
if ( !$PON_AD_threshold ) {
	$logger->info(
"PON_AD_Threshold:Panel of Normal AD threshold not given default will be used"
	);
	$PON_AD_threshold = 3;
}
else {
	$logger->info("PON_AD_Threshold:$PON_AD_threshold");
}

#Check PON TPVF
if ( !$PON_TPVF_threshold ) {
	$logger->info(
"PON_TPVF_Threshold:Panel of Normal TPVF threshold not given default will be used"
	);
	$PON_TPVF_threshold = 10;
}
else {
	$logger->info("PON_TPVF_Threshold:$PON_TPVF_threshold");
}

#Check Pindel Min Indel Length
if ( !$Pindel_Min_Indel_Len ) {
	$logger->info(
"Pindel_Min_Indel_Len:Pindel Min Indel Len not given default will be used"
	);
	$Pindel_Min_Indel_Len = 25;
}
else {
	$logger->info("Pindel_Min_Indel_Len:$Pindel_Min_Indel_Len");
}

#Check Pindel Max Indel Length
if ( !$Pindel_Max_Indel_Len ) {
	$logger->info(
"Pindel_Max_Indel_Len:Pindel Max Indel Len not given default will be used"
	);
	$Pindel_Max_Indel_Len = 2000;
}
else {
	$logger->info("Pindel_Max_Indel_Len:$Pindel_Max_Indel_Len");
}

#Check for standard normal fof for genotyping
if ( !$standardNormalList ) {
	$standardNormalList =
	  "/ifs/e63data/bergerm1/Resources/StandardNormals/StandardNormals.list";
	$logger->info(
		"Standard Normal File is not given thus the default file will be used."
	);
}
else {
	$logger->info("Standard Normal File to be used is $standardNormalList");
}

#Check if path to exon2gene coverage script is given
if ( !$ExonToGenCov ) {
	$ExonToGenCov = $supportScriptPath . "dmp_exon2gene_coverage.pl";
	$logger->info(
"Exon to gene coverage script path is not given, default location ($ExonToGenCov) will be used."
	);
}
else {
	$logger->info("Exon to gene coverage script location: $ExonToGenCov.");
}

#Check if path to FP_genotypes script is given
if ( !$FPGenotypesScript ) {
	$FPGenotypesScript = $supportScriptPath . "dmp_fp_genotypes.pl";
	$logger->info(
"Fingerprint genotyping script path is not given, default location ($FPGenotypesScript) will be used."
	);
}
else {
	$logger->info(
		"Fingerprint genotyping script location: $FPGenotypesScript.");
}

#Check if path to filter mutect script is given
if ( !$filter_Mutect ) {
	$filter_Mutect = $supportScriptPath . "dmp_filter_mutect.pl";
	$logger->info(
"Filter mutect script path is not given, default location ($filter_Mutect) will be used."
	);
}
else {
	$logger->info("Filter mutect  script location: $filter_Mutect.");
}

#Check if path to filter somatic indel script is given
if ( !$filter_SomaticIndel ) {
	$filter_SomaticIndel = $supportScriptPath . "dmp_filter_indel.pl";
	$logger->info(
"Filter indel script path is not given, default location ($filter_SomaticIndel) will be used."
	);
}
else {
	$logger->info("Filter indel script location: $filter_SomaticIndel.");
}

#Check if path to filter pindel script is given
if ( !$filter_Pindel ) {
	$filter_Pindel = $supportScriptPath . "dmp_filter_pindel.py";
	$logger->info(
"Filter pindel script path is not given, default location ($filter_Pindel) will be used."
	);
}
else {
	$logger->info("Filter pindel script location: $filter_Pindel.");
}

#Check if path to compile metrics script is given
if ( !$CompileMetrics ) {
	$CompileMetrics = $binPath . "/dmp_compile_qc_metrics.pl";
	$logger->info(
"Compile metrics script path is not given, default location ($CompileMetrics) will be used."
	);
}
else {
	$logger->info("Compile metrics script location: $CompileMetrics.");
}

#Check if path to abra wrapper script is given
if ( !$ABRAwrapper ) {
	$ABRAwrapper = $binPath . "/Run_ABRArealignment.py";
	$logger->info(
"ABRAwrapper script path is not given, default location ($ABRAwrapper) will be used."
	);
}
else {
	$logger->info("ABRA wrapper script location: $ABRAwrapper.");
}

#Check if path to pindel wrapper script is given
if ( !$PindelWrapper ) {
	$PindelWrapper = $binPath . "/Run_Pindel.py";
	$logger->info(
"PindelWrapper script path is not given, default location ($PindelWrapper) will be used."
	);
}
else {
	$logger->info("Pindel wrapper script location: $PindelWrapper.");
}

#Check if path to Annotate assess filter variants script is given
if ( !$AnnotateAssessFilterVariants ) {
	$AnnotateAssessFilterVariants =
	  $supportScriptPath . "/dmp_annotate_filter_variants.pl";
	$logger->info(
"Variant annotation and filtering script path is not given, default location ($AnnotateAssessFilterVariants) will be used."
	);
}
else {
	$logger->info(
"Variant annotation and filtering script location: $AnnotateAssessFilterVariants."
	);
}

#Check if path to genotype allele script is given
if ( !$GenotypeAllele ) {
	$GenotypeAllele = $binPath . "/dmp_genotype_allele.pl";
	$logger->info(
"Genotype allele script path is not given, default location ($GenotypeAllele) will be used."
	);
}
else {
	$logger->info("Genotype allele script location: $GenotypeAllele.");
}

#Check if path to filter genotyped variants script is given
if ( !$filterGenotypedVariants ) {
	$filterGenotypedVariants =
	  $supportScriptPath . "/dmp_filter_genotyped_variants.pl";
	$logger->info(
"Filter Genotyped Variants script path is not given, default location ($filterGenotypedVariants) will be used."
	);
}
else {
	$logger->info(
		"Filter Genotyped Variants script location: $GenotypeAllele.");
}

#Check if path to FindCoveredIntervals is given
if ( !$FindCoveredIntervals ) {
	$FindCoveredIntervals = $supportScriptPath . "/Run_FindCoveredIntervals.py";
	$logger->info(
"Run_FindCoveredIntervals.py script path is not given, default location ($filterGenotypedVariants) will be used."
	);
}
else {
	$logger->info(
		"FindCoveredIntevals script location: $FindCoveredIntervals.");
}

#Check if path to create coverage files script is given
if ( !$createCoverageFiles ) {
	$createCoverageFiles =
	  $supportScriptPath . "/dmp_create_coverage_files_and_clinical_reports.pl";
	$logger->info(
"Create coverage files script path is not given, default location ($createCoverageFiles) will be used."
	);
}
else {
	$logger->info(
		"Create coverage files  script location: $createCoverageFiles.");
}

#Check if path to create coverage files script is given
if ( !$createPatientVCFfile ) {
	$createPatientVCFfile =
	  $supportScriptPath . "/dmp_create_patient_vcf_file.pl";
	$logger->info(
"Create patient VCF files script path is not given, default location ($createPatientVCFfile) will be used."
	);
}
else {
	$logger->info(
		"Create patient VCF files  script location: $createPatientVCFfile.");
}

#Check if path to generate metrics files script is given
if ( !$GenerateMetricsFilesScript ) {
	$GenerateMetricsFilesScript =
	  $supportScriptPath . "/dmp_generate_metrics_report.pl";
	$logger->info(
"Create patient VCF files script path is not given, default location ($GenerateMetricsFilesScript) will be used."
	);
}
else {
	$logger->info(
"Create patient VCF files  script location: $GenerateMetricsFilesScript."
	);
}

#Check if path to triallelic site fix script is given
if ( !$triallelic_fix ) {
	$triallelic_fix = $supportScriptPath . "/dmp_triallelic_fix.pl";
	$logger->info(
"Add entries to raw MUTECT output files in case of a triallelic site, default location ($triallelic_fix) will be used."
	);
}
else {
	$logger->info(
"Create patient VCF files  script location: $GenerateMetricsFilesScript."
	);
}

#Prin Version of tools and Files used
$Version->{"Pipeline_Version"} = $revision;
my $PrintConfigFile = "RunConfigration.txt";
open( VERSION, ">", "$outdir/$PrintConfigFile" )
  or die( $logger->fatal("Cannot open $outdir/$PrintConfigFile. Error: $!") );
print VERSION "Tools|Files\tVersion\n";
while ( my ( $tools_files, $version ) = each(%$Version) ) {
	print VERSION "$tools_files\t$version\n";
}

#print VERSION "RawDataDir\t$datadir\n";
close(VERSION);
tie( my %classPerBarcode, 'Tie::IxHash' );
for ( my $i = 0 ; $i < scalar(@$barcode) ; $i++ ) {

	#print "$$barcode[$i] => $$class[$i]\n";
	$classPerBarcode{ $$barcode[$i] } = $$class[$i];
}

# Make a dummy wait file to run after process
&MakeCSH($outdir);

#Split porces to know how many to run
my @allProcess = split( ",", $process );

#Check for raw data directory
if ( $allProcess[0] == 1 && !$datadir ) {
	$logger->fatal(
"Please enter the directory that contains the data to be processed.See Usage"
	);
	Usage();
	exit(1);
}
elsif ( $allProcess[0] == 1 && $datadir ) {
	$logger->info("The raw data directory given is $datadir");
}
else {
	$logger->info(
"Process 1 is not selected so I am not checking to see if the raw data directory is given or not."
	);
}
my $allProcessList  = join( ",", @allProcess );
my $numberOfProcess = scalar(@allProcess);
my $processCount    = 0;
my $parseFilenames;
while ( $processCount < $numberOfProcess ) {

	#print "@allProcess\n";
	my $runProcess = shift(@allProcess);

	#print "$runProcess\n";
	($parseFilenames) = &Select( $runProcess, $parseFilenames );
	$processCount++;
}
if ($fof) {
	if ( $removeSymLinksFlag == 1 ) {
		$logger->warn("Removing Softlinked files form $outdir");
		my @softlinks = ();
		eval {
			@softlinks =
			  `find $outdir -maxdepth 1 -type l -print0 | xargs -0 ls -d`;
		};
		if ($@) {
			$logger->warn("Cannot find all softlinks in $outdir. Error:$@");
		}
		foreach my $sfLink (@softlinks) {
			eval { `unlink $sfLink`; };
			if ($@) { $logger->warn("Cannot unlink $sfLink. Error:$@"); }
		}
	}
}
if ( $mvFiles == 1 ) { &RunHouseKeeping( $outdir, $pool ) }
eval { `chmod ug+rwx -R $outdir`; };
if ($@) { $logger->warn("Cannot change permissions for $outdir. Error:$@"); }

#--Calculate total runtime (current time minus start time)
$now = time - $now;
if ( length($mailID) ) {

	# if email address exists email the user
`echo \"IMPACT pipeline for $poolName has finished running the steps: $allProcessList. The output files can be found at $outdir. Thanks!\" | mail -s \"IMPACT run has completed!\" \"$mailID\"`;
}

#--Print runtime
printf(
	"\n\nTotal Pipeline running time: %02d:%02d:%02d\n\n",
	int( $now / 3600 ),
	int( ( $now % 3600 ) / 60 ),
	int( $now % 60 )
);
$logger->info("Impact_Pipeline:Done, Thanks for using the framework");
exit(0);
#####################################
#####################################
#How to use the script.
sub Usage {
	print "Unknow option: @_\n" if (@_);
	print "\nUsage : RunIlluminaProcess.pl [options]
        [--config|c                        S Path to configration file(required)]
        [--svConfig|sc                     S Path to structural variant configration file(optional)]
        [--symLinkFlag|sf           	   I Flag for Keeping or removing the symolic links(1:Remove;2:Keep)(default:2)]
        [--dataDirectory|d                 S Path where all the files to be processed are located (required)]
        [--outputDirectory|o               S Path where all the output files will be written (required)]
        Inside Config File
        >Location
        Location of all the different files and programs required by the pipeline to run.
        >Version
        Version of each script and program used by the pipeline.
        >Parameters
        SampleFile                         S csv file describing details about the sample (required and submit with full path)
        TitleFile                          S tab-delimited title file for the samples (required and submit with full path)
        StdNormalForMutationCalling|n      S file to be used as standard normal #full path of bam file, *.bai file to be located in same folder)
        ListOfFiles                        S Name of the files with there path (fof:Paired files;one per line;one after another) (optional)
        RunBQSR                            I 0 => Skip BQSR. 1 => Run Sample Level. 2 => Run Lane Level (default:1,optional)
        ProjectName                        S Name of the project(required,e.g:Colons).
        PoolName                           S Name of the pool(required,e.g:Colon5P1).
        Process                            S 1 => Merge Fastq. 2 => Mapping. 3 => Process Bams and Metrics Calculation.  4 => Calculate and Compile Metrics. 5 => Call SNPs and Indels. 6 => Filter SNPs & Indels. 7 => Annotate & Assess SNPs and Indels. qc|QC => Quality COntrol Part of pipeline (1,2,3,4). vc|VC => Varaint Calling Part (4,5,6,7). all|ALL|All => Run all steps of the pipeline(1,2,3,4,5,6,7). \"1,2,3,4,5,6,7\" for all or in valid combination.
       
        
        FastqSource                        S Fastq Source: \"GCL\" or \"DMP\" [if blank, assume GCL, if not blank, must be \"GCL\" or \"DMP\"
        ListOfStandardNoramlsForGenotyping S File containing list of all standard Normal files to be ran(optional;default:/ifs/e63data/bergerm1/Resources/StandardNormals/StandardNormals.list)
        MergeDinucleotide                  I 2 => Do not merge Di-nucleotide 1 => Merge Dinucleotide.(default:1,optional)
        MoveFiles                          I 2 => Skip Making Folders & Moving Files. 1 => Make Folders and Move Files.(default:1,optional)
        TNfreqRatio_MutectStdFilter        I Tumor Normal varaint frequecy ratio for mutect(default:5,optional).
        AD_MutectStdFilter                 I Tumor Allele depth threshold for mutect(default:3,optional). 
        DP_MutectStdFilter                 I Tumor total depth threshold for mutect(default:0,optional).
        VF_MutectStdFilter                 F Tumor variant frequency threshold for mutect(default:0.01,optional). 
        TNfreqRatio_SomIndelStdFilter      I Tumor Normal varaint frequecy ratio for Somatic Indel Detector(default:5,optional).
        AD_SomIndelStdFilter               I Tumor Allele depth threshold for Somatic Indel Detector(default:3,optional).
        DP_SomIndelStdFilter               I Tumor total depth threshold for Somatic Indel Detector(default:0,optional).
        VF_SomIndelStdFilter               F Tumor variant frequency threshold for Somatic Indel Detector(default:0.01,optional).
        TNfreqRatio_AnnotationFilter       I Tumor Normal varaint frequecy ratio for Annotation(default:5,optional).
        Tfreq_AnnotationFilter             F Tumor variant frequency threshold for Annotation(default:0.02,optional)
        MAFthreshold_AnnotationFilter      F 1000 genomes Maf threshold for Annotation.(default:0.01,optional)
        AD_Tthreshold_indel_high_AnnotationFilter   I High Tumor Allele depth threshold for indels(default:10,optional) 
        AD_Tthreshold_indel_low_AnnotationFilter    I Low Tumor Allele depth threshold for indels(default:3,optional)
        VF_Tthreshold_indel_high_AnnotationFilter   F High Tumor Varaint Frequency threshold for indels(default:0.1,optional) 
        VF_Tthreshold_indel_high_AnnotationFilter   F Low Tumor Varaint Frequency threshold for indels(default:0.02,optional) 
        QUEUE_NAME                                 S Name of the Queue where the pipeline needs to run (default:all.q,optional)
	\n";

	print "For any question please contact Ronak Shah (shahr2\@mskcc.org)\n";
	print "!~!~! Comments, flames and suggestion are welcome !~!~!\n";
	exit();
}
###################################################
###################################################
#--GET CONFIGRATION DETAIL
sub GetConfiguration {
	my ($config_file) = @_;
	my @data = ();
	tie( my %config,     'Tie::IxHash' );
	tie( my %location,   'Tie::IxHash' );
	tie( my %version,    'Tie::IxHash' );
	tie( my %parameters, 'Tie::IxHash' );
	$logger->info("Reading the configuration file");

	# change the default Input Record Separator.
	$/ = ">";
	open( CONFIG, "$config_file" )
	  or die( $logger->fatal("Cannot open $config_file. Error: $!") );
	my $junk = <CONFIG>;
	while (<CONFIG>) {
		next if ( $_ =~ /^#/ );
		chomp($_);
		my ( $defLine, @configLines ) = split /\n/, $_;
		if ( $defLine =~ /Locations/ ) {
			foreach my $config (@configLines) {
				next if ( $config =~ /^#/ );
				@data = split( "=", $config, 2 );
				$data[0] =~ s/\s//g;
				$data[1] =~ s/\s//g;
				$location{ $data[0] } = $data[1];
			}
		}
		if ( $defLine =~ /Parameters/ ) {
			foreach my $config (@configLines) {
				next if ( $config =~ /^#/ );
				@data = split( "=", $config, 2 );
				$data[0] =~ s/\s//g;
				$data[1] =~ s/\s//g;
				$parameters{ $data[0] } = $data[1];
			}
		}
		if ( $defLine =~ /Versions/ ) {
			foreach my $config (@configLines) {
				next if ( $config =~ /^#/ );
				@data = split( "=", $config, 2 );
				$data[0] =~ s/\s//g;
				$data[1] =~ s/\s//g;
				$version{ $data[0] } = $data[1];
			}
		}
	}
	close(CONFIG);
	$logger->info("Completed reading the configuration file");

	# Change back the Input Record Separator to default.
	$/ = "\n";
	eval {
		##Set Locations
		$TMPDIR              = $location{"TMPDIR"};
		$ZCAT                = $location{"ZCAT"};
		$JAVA_1_6            = $location{"JAVA_1_6"};
		$JAVA_1_7            = $location{"JAVA_1_7"};
		$RHOME               = $location{"RHOME"};
		$RLIBS               = $location{"RLIBS"};
		$ExonToGenCov        = $location{"ExonToGenCov"};
		$FPGenotypesScript   = $location{"FPGenotypesScript"};
		$FP_genotypes        = $location{"FP_genotypes"};
		$GATK_SomaticIndel   = $location{"GATK_SomaticIndel"};
		$GATK                = $location{"GATK"};
		$Reference           = $location{"Reference"};
		$Refseq              = $location{"Refseq"};
		$PICARD              = $location{"PICARD"};
		$Mutect              = $location{"Mutect"};
		$filter_Mutect       = $location{"filter_Mutect"};
		$filter_SomaticIndel = $location{"filter_SomaticIndel"};
		$filter_Pindel       = $location{"filter_Pindel"};
		$BaitInterval        = $location{"BaitInterval"};
		$TargetInterval      = $location{"TargetInterval"};
		$TargetRegionLIST    = $location{"TargetRegionLIST"};
		$CompileMetrics      = $location{"CompileMetrics"};
		$SVpipeline          = $location{"SVpipeline"};
		$CAT                 = $location{"CAT"};
		$PYTHON              = $location{"PYTHON"};
		$TrimGalore          = $location{"TrimGalore"};
		$PERL                = $location{"PERL"};
		$BWA                 = $location{"BWA"};
		$GeneInterval        = $location{"GeneInterval"};
		$GeneIntervalAnn     = $location{"GeneIntervalAnn"};
		$GeneCoord           = $location{"GeneCoord"};
		$TilingInterval      = $location{"TilingInterval"};
		$TilingIntervalAnn   = $location{"TilingIntervalAnn"};
		$FingerPrintInterval = $location{"FingerPrintInterval"};
		$dbSNP               = $location{"dbSNP"};
		$COSMIC              = $location{"COSMIC"};
		$Mills_1000G_Indels  = $location{"Mills_1000G_Indels"};
		$dbSNP_bitset        = $location{"dbSNP_bitset"};
		$dbProperties        = $location{"dbProperties"};
		$Oncotator           = $location{"Oncotator"};
		$Mutation_Assessor   = $location{"Mutation_Assessor"};
		$AnnotateAssessFilterVariants =
		  $location{"AnnotateAssessFilterVariants"};
		$LoessNormalization   = $location{"LoessNormalization"};
		$BestCopyNumber       = $location{"BestCopyNumber"};
		$NormVsNormCopyNumber = $location{"NormVsNormCopyNumber"};
		$StdNormalLoess_TM    = $location{"StdNormalLoess_TM"};
		$StdNormalLoess_NVN   = $location{"StdNormalLoess_NVN"};
		$NormalVariantsVCF    = $location{"NormalVariantsVCF"};
		$AllMetrics           = $location{"AllMetrics"};
		$SAMTOOLS             = $location{"SAMTOOLS"};
		$BEDTOOLS             = $location{"BEDTOOLS"};
		$GenotypeAllele       = $location{"GenotypeAllele"};
		$cosmicHotspotsVcf    = $location{"CosmicHotspotVcf"};
		$GCBiasFile           = $location{"GCBiasFile"};
		$HistNormDir          = $location{"HistNormDir"};
		$barcodeFile          = $location{"BarcodeKey"};
		$adaptorFile          = $location{"AdaptorKey"};
		$QSUB                 = $location{"QSUB"};
		$BSUB                 = $location{"BSUB"};
		$RSYNC                = $location{"RSYNC"};
		$ClinicalExons        = $location{"clinicalExons"};
		$HotspotMutations     = $location{"HotSpot_mutations"};
		$IGVtools             = $location{"IGVtools"};
		$translationFolder    = $location{"TranslationFolder"};
		$exonIntervalsFile = $location{"Canonical_Exon_Interval_table_with_aa"};
		$canonicalExonIntervalsFile =
		  $location{"Canonical_Exon_Interval_list"};    # this is for DoC
		$validatedExons           = $location{"Validated_Exons"};
		$ABRA                     = $location{"ABRA"};
		$ABRAwrapper              = $location{"ABRAwrapper"};
		$PINDELBIN                = $location{"PINDELBIN"};
		$PindelWrapper            = $location{"PindelWrapper"};
		$CountErrors              = $location{"CountErrors"};
		$StandardNormalsDirectory = $location{"StandardNormalsDirectory"};
		$FindCoveredIntervals     = $location{"FindCoveredIntervals"};
		##Set Parameters
		$sampleFile     = $parameters{"SampleFile"};
		$titleFile      = $parameters{"TitleFile"};
		$fastqSource    = $parameters{"FastqSource"};
		$runABRA        = $parameters{"runABRA"};
		$triallelic_fix = $parameters{"TriallelicFix"};

		#$outdir = $parameters{"OutputDir"};
		#$datadir = $parameters{"RawDataDir"};
		$MAPQ               = $parameters{"MAPQ"};
		$BASQ               = $parameters{"BASQ"};
		$stdNormal          = $parameters{"StdNormalForMutationCalling"};
		$fof                = $parameters{"ListOfFiles"};
		$poolName           = $parameters{"PoolName"};
		$projectName        = $parameters{"ProjectName"};
		$runBQSR            = $parameters{"RunBQSR"};
		$process            = $parameters{"Process"};
		$standardNormalList = $parameters{"ListOfStandardNoramlsForGenotyping"};
		$mvFiles            = $parameters{"MoveFiles"};
		$mergeDinucleotide  = $parameters{"MergeDinucleotide"};
		$MAFthreshold       = $parameters{"MAFthreshold_AnnotationFilter"};
		$TNfreqRatio_MutectStdFilter =
		  $parameters{"TNfreqRatio_MutectStdFilter"};
		$TNfreqRatio_SomIndelStdFilter =
		  $parameters{"TNfreqRatio_SomIndelStdFilter"};
		$TNfreqRatioThreshold = $parameters{"TNfreqRatio_AnnotationFilter"};
		$ad_SomIndelStdFilter = $parameters{"AD_SomIndelSTDfilter"};
		$dp_SomIndelStdFilter = $parameters{"DP_SomIndelSTDfilter"};
		$vf_SomIndelStdFilter = $parameters{"VF_SomIndelSTDfilter"};
		$ad_MutectStdFilter   = $parameters{"AD_MutectSTDfilter"};
		$dp_MutectStdFilter   = $parameters{"DP_MutectSTDfilter"};
		$vf_MutectStdFilter   = $parameters{"VF_MutectSTDfilter"};
		$dp_snv               = $parameters{"minimumDPforSNV"};
		$ad_snv               = $parameters{"minimumADforSNV"};
		$vf_snv               = $parameters{"minimumVFforSNV"};
		$dp_snvHS             = $parameters{"minimumDPforSNVhs"};
		$ad_snvHS             = $parameters{"minimumADforSNVhs"};
		$vf_snvHS             = $parameters{"minimumVFforSNVhs"};
		$dp_indel             = $parameters{"minimumDPforINDEL"};
		$ad_indel             = $parameters{"minimumADforINDEL"};
		$vf_indel             = $parameters{"minimumVFforINDEL"};
		$dp_indelHS           = $parameters{"minimumDPforINDELhs"};
		$ad_indelHS           = $parameters{"minimumADforINDELhs"};
		$vf_indelHS           = $parameters{"minimumVFforINDELhs"};
		$occurrencePercent    = $parameters{"occurrencePercent"};
		$coverageThreshold    = $parameters{"Coverage_threshold_darwin_report"};
		$mailID               = $parameters{"Email"};
		$queue                = $parameters{"QUEUE_NAME"};
		$CLUSTER              = $parameters{"CLUSTER"};
		$deleteIntermediateFiles = $parameters{"DeleteIntermediateFiles"};
		$PON_AD_threshold        = $parameters{"PON_AD_Threshold"};
		$PON_TPVF_threshold      = $parameters{"PON_TPVF_Threshold"};
		$Pindel_Min_Indel_Len    = $parameters{"Pindel_Min_Indel_Len"};
		$Pindel_Max_Indel_Len    = $parameters{"Pindel_Max_Indel_Len"};
	};
	if ($@) {
		$logger->fatal(
			"Did not find a variable in configuration file.Error: $@\n");
		exit(1);
	}
	return ( \%version );
}
###################################################
###################################################
#--Make Notification file
sub MakeCSH {
	my ($outdir) = @_;
	my $filename = $outdir . "/Notify.csh";
	my $ntmp     = new IO::File(">$filename");
	print $ntmp "#!/bin/csh\n";
	print $ntmp "#Notification File\n";
	print $ntmp "echo", " This is Done", "\n";
	$ntmp->close();
	eval { `chmod +x $filename`; };

	if ($@) {
		$logger->fatal("Cannot change permissions for $filename. Error:$@");
		exit(1);
	}
	return;
}
###################################################
###################################################
#--Waiting for the process to finish
sub WaitToFinish {
	my ( $outdir, @waitfilenames ) = @_;
	$logger->info("Waiting for the Process to finish...");
	foreach my $wfile (@waitfilenames) {
		next if ( $wfile eq "NULL" );
		sleep 10 while ( !( -e "$outdir/$wfile" ) );

		#print "$outdir/$wfile\n";
		while ( -e "$outdir/$wfile" ) {

			#print "$outdir/$wfile\n";
			open( FH, "<", "$outdir/$wfile" );
			while (<FH>) {
				if ( $_ =~ /This is Done/ig ) {

					#print "\nFinished: $wfile\n";
					last;
				}
				else {
					wait;
				}
			}
			last;
		}
		close(FH);
	}
	foreach my $wfile (@waitfilenames) {
		next if ( $wfile eq "NULL" );
		eval { `rm $outdir/$wfile`; };
		if ($@) { $logger->warn("Cannot remove $outdir/$wfile. Error:$@"); }
	}
	return;
}
###################################################
###################################################
#--Check Error Files
sub CheckErrorFiles {
	my ( $outdir, @errorfilenames ) = @_;
	$logger->info("Checking Error Files");
	foreach my $efile (@errorfilenames) {
		next if ( $efile eq "NULL" );
		if ( -e "$outdir/$efile" ) {
			if ( ( -s "$outdir/$efile" ) != 0 ) {
				$logger->warn("Please Check $efile; Something went wrong here");
			}
			else {
				eval { `rm $outdir/$efile`; };
				if ($@) {
					$logger->warn("Cannot remove $outdir/$efile. Error:$@");
				}
				next;
			}
		}
		else {
			$logger->warn(
"Please check, no error file was created named $efile; Something went wrong here."
			);
		}
	}
	return;
}
###################################################
###################################################
#--Check Output Files
sub CheckOutputFiles {
	my ( $outdir, $prog, @Outputfilenames ) = @_;
	$logger->info("Checking Output Files");
	$prog = "NULL" if ( !$prog );
	if ( ( $prog eq "STDfilter" ) or ( $prog eq "RawMutationFiles" ) ) {
		foreach my $ofile (@Outputfilenames) {
			if ( ( $ofile =~ /\// ) and ( !( $ofile =~ /,/ ) ) ) {
				if ( -e "$ofile" ) {
					next;
				}
				else {
					$logger->fatal(
"CheckFile:Please check, $ofile file was not created;Something went wrong here."
					);
					exit(1);
				}
			}
			if ( -e "$outdir/$ofile" ) {
				next;
			}
			else {
				$logger->fatal(
"CheckFile:Please check, $outdir/$ofile file was not created;Something went wrong here."
				);
				exit(1);
			}
		}
	}
	else {
		foreach my $ofile (@Outputfilenames) {
			if ( ( $ofile =~ /\// ) and ( !( $ofile =~ /,/ ) ) ) {
				if ( ( -e "$ofile" ) and ( ( -s "$ofile" ) != 0 ) ) {
					next;
				}
				else {
					$logger->fatal(
"CheckFile:Please check, $ofile file was not created or the file size is zero; Something went wrong here."
					);
					exit(1);
				}
			}
			elsif ( ( $ofile =~ /,/ ) and ( !( $ofile =~ /\// ) ) ) {
				my ( $file1, $file2 ) = split( ",", $ofile );
				if (    ( -e "$outdir/$file1" )
					and ( ( -s "$outdir/$file1" ) != 0 )
					and ( -e "$outdir/$file2" )
					and ( ( -s "$outdir/$file2" ) != 0 ) )
				{
					next;
				}
				else {
					$logger->fatal(
"CheckFile:Please check,$file1 & $file2 files was not created or the file size is zero; Something went wrong here."
					);
					exit(1);
				}
			}
			elsif ( ( $ofile =~ /,/ ) and ( $ofile =~ /\// ) ) {
				my ( $file1, $file2 ) = split( ",", $ofile );
				if (    ( -e "$file1" )
					and ( ( -s "$file1" ) != 0 )
					and ( -e "$file2" )
					and ( ( -s "$file2" ) != 0 ) )
				{
					next;
				}
				else {
					$logger->fatal(
"CheckFile:Please check,$file1 & $file2 files was not created or the file size is zero; Something went wrong here."
					);
					exit(1);
				}
			}
			else {
				if (    ( -e "$outdir/$ofile" )
					and ( ( -s "$outdir/$ofile" ) != 0 ) )
				{
					next;
				}
				else {
					$logger->fatal(
"CheckFile:Please check, $outdir/$ofile file was not created or the file size is zero;Something went wrong here."
					);
					exit(1);
				}
			}
		}
	}
	$logger->info("Finished Checking Output Files");
	return;
}
###################################################
###################################################
#--Make array of file of files list from the outdir
sub GetNames {
	my ( $fof, $outdir ) = @_;
	my (@filenames) = ();
	open( FOF, "$outdir/$fof" )
	  or die( $logger->fatal("GetNames:Cannot open $outdir/$fof Error:$!") );
	while (<FOF>) {
		chomp;
		next
		  if (/^\s$/)
		  ;   # this is to fix the extra line at the end of a list of files file
		$_ =~ s/\s//g;
		my $filename = pop @{ [ split( "/", $_ ) ] };
		push( @filenames, $filename );
	}
	return (@filenames);
}
###################################################
###################################################
#--Make Pairs of the files.
sub MAKEPAIRS {
	my ( $filenames, $outdir ) = @_;
	my @names      = @$filenames;
	my $count      = scalar @names;
	my (@newnames) = ();
	if ( $count % 2 != 0 ) {
		$logger->fatal(
			"MAKEPAIRS:Odd number of files given, please check Input file.");
		exit(1);
	}
	else {
		for ( my $i = 0 ; $i < scalar(@names) ; $i += 2 ) {
			chomp( $names[$i] );
			chomp( $names[ $i + 1 ] );
			push( @newnames, "$names[$i],$names[$i+1]" );
		}
	}
	return (@newnames);
}
###################################################
###################################################
#--Reverse Complement of the sequence.
sub RevComp {
	my ($seq) = @_;
	my $revseq = reverse($seq);
	$revseq = uc($revseq);
	$revseq =~ tr/ACGT/TGCA/;
	return ($revseq);
}
#####################################
#####################################
#Read data related to samples as well as barcodes.
sub ReadSampleFile {
	my ( $sampleFile, $projectName, $outdir ) = @_;
	my (
		@fcId,        @lane,    @sampleId, @sampleRef, @index,
		@description, @control, @recipe,   @operator,  @sampleProject
	) = ();
	my $sampleFileName = "";
	if ( $sampleFile =~ /\// ) {
		$sampleFileName = pop @{ [ split( "/", $sampleFile ) ] };
	}
	else {
		$sampleFileName = $sampleFile;
	}
	open( SAMPLEFILE, $sampleFile )
	  or
	  die( $logger->fatal("ReadSampleFile:Cannot open $sampleFile. Error:$!") );
	while (<SAMPLEFILE>) {
		next if $. == 1;
		my @oldDataCols = split( ",", $_ );
		my @dataCols =
		  grep( s/\s*$//g, @oldDataCols );    #remove whitespace if any
		if ( $dataCols[0] ) { push( @fcId,          $dataCols[0] ); }
		if ( $dataCols[1] ) { push( @lane,          $dataCols[1] ); }
		if ( $dataCols[2] ) { push( @sampleId,      $dataCols[2] ); }
		if ( $dataCols[3] ) { push( @sampleRef,     $dataCols[3] ); }
		if ( $dataCols[4] ) { push( @index,         $dataCols[4] ); }
		if ( $dataCols[5] ) { push( @description,   $dataCols[5] ); }
		if ( $dataCols[6] ) { push( @control,       $dataCols[6] ); }
		if ( $dataCols[7] ) { push( @recipe,        $dataCols[7] ); }
		if ( $dataCols[8] ) { push( @operator,      $dataCols[8] ) }
		if ( $dataCols[9] ) { push( @sampleProject, $dataCols[9] ); }
	}
	close(SAMPLEFILE);
	if ( !-e "$outdir/$sampleFileName" ) {
		eval { `cp $sampleFile $outdir/$sampleFileName`; };
		if ($@) {
			$logger->warn("Cannot copy $outdir/$sampleFileName. Error:$@");
		}
	}
	return (
		\@fcId,     \@lane,        \@sampleId, \@sampleRef,
		\@index,    \@description, \@control,  \@recipe,
		\@operator, \@sampleProject
	);
}
#####################################
#####################################
#Read data related to samples as well as barcodes from title file.
sub ReadTitleFile {
	my ( $titleFile, $outdir ) = @_;
	my @barcode      = ();
	my @pool         = ();
	my @sampleId     = ();
	my @collabId     = ();
	my @patientId    = ();
	my @class        = ();
	my @sampleType   = ();
	my @inputNg      = ();
	my @libraryYeild = ();
	my @poolInput    = ();
	my @baitVersion  = ();
	my @fof          = ();
	my @newfof       = ();
	open( TFH, $titleFile )
	  or die(
		$logger->fatal("ReadTitleFile:Cannot open file $titleFile. Error:$!") );

	while (<TFH>) {
		next if ( $. == 1 );
		my @dataCols = split( "\t", $_ );
		my @newDatacols =
		  grep( s/\s*$//g, @dataCols );    #remove whitespace if any
		push( @barcode,      $newDatacols[0] );
		push( @pool,         $newDatacols[1] );
		push( @sampleId,     $newDatacols[2] );
		push( @collabId,     $newDatacols[3] );
		push( @patientId,    $newDatacols[4] );
		push( @class,        $newDatacols[5] );
		push( @sampleType,   $newDatacols[6] );
		push( @inputNg,      $newDatacols[7] );
		push( @libraryYeild, $newDatacols[8] );
		push( @poolInput,    $newDatacols[9] );
		push( @baitVersion,  $newDatacols[10] );
	}
	close(TFH);
	my $poolName         = $pool[0];
	my $newtitleFileName = $poolName . "_title.txt";
	if ( !-e "$outdir/$newtitleFileName" ) {
		eval { `cp $titleFile $outdir/$newtitleFileName`; };
		if ($@) {
			$logger->warn("Cannot copy $outdir/$newtitleFileName. Error:$@");
		}
	}
	return (
		\@barcode,      \@pool,      \@sampleId,   \@collabId,
		\@patientId,    \@class,     \@sampleType, \@inputNg,
		\@libraryYeild, \@poolInput, \@baitVersion
	);
}
#####################################
#####################################
#Run only single process at a time.
sub Select {
	my ( $process, $parseFilenames ) = @_;
	if ( $process == 1 ) {
		$parseFilenames = &MergeDataFromDirectory();
	}
	elsif ( $process == 2 ) {
		$parseFilenames = &DoMapping($parseFilenames);
	}
	elsif ( $process == 3 ) {
		$parseFilenames = &ProcessBams($parseFilenames);
	}
	elsif ( $process == 4 ) {
		$parseFilenames = &CalculateAndCompileMetrics($parseFilenames);
		if ( $svTag == 2 ) { $parseFilenames = &CallSV($parseFilenames); }
	}
	elsif ( $process == 5 ) {
		$parseFilenames = &CallingSNPsAndIndels($parseFilenames);
	}
	elsif ( $process == 6 ) {
		$parseFilenames = &FilterSNPsAndIndels($parseFilenames);
	}
	elsif ( $process == 7 ) {
		$parseFilenames =
		  &AnnotateSNPsAndIndels( $parseFilenames, $outdir, $fof );
	}
	elsif ( ( $process eq "qc" ) or ( $process eq "QC" ) ) {
		$parseFilenames = &MergeDataFromDirectory();
		$parseFilenames = &DoMapping($parseFilenames);
		$parseFilenames = &ProcessBams($parseFilenames);
		$parseFilenames = &CalculateAndCompileMetrics($parseFilenames);
	}
	elsif ( ( $process eq "vc" ) or ( $process eq "VC" ) ) {
		$parseFilenames = &CalculateAndCompileMetrics($parseFilenames);
		if ( $svTag == 2 ) { $parseFilenames = &CallSV($parseFilenames); }
		$parseFilenames = &CallingSNPsAndIndels($parseFilenames);
		$parseFilenames = &FilterSNPsAndIndels($parseFilenames);
		$parseFilenames = &AnnotateSNPsAndIndels($parseFilenames);
	}
	elsif (( $process eq "all" )
		or ( $process eq "ALL" )
		or ( $process eq "All" ) )
	{
		$parseFilenames = &MergeDataFromDirectory();
		$parseFilenames = &DoMapping($parseFilenames);
		$parseFilenames = &ProcessBams($parseFilenames);
		$parseFilenames = &CalculateAndCompileMetrics($parseFilenames);
		if ( $svTag == 2 ) { $parseFilenames = &CallSV($parseFilenames); }
		$parseFilenames = &CallingSNPsAndIndels($parseFilenames);
		$parseFilenames = &FilterSNPsAndIndels($parseFilenames);
		$parseFilenames = &AnnotateSNPsAndIndels($parseFilenames);
	}
	else {
		$logger->fatal(
"Select:The process number entered ($process) does not exist, See Usage."
		);
		exit(1);
	}
	return ($parseFilenames);
}

sub uniq (@) {
	my %seen = ();
	grep { not $seen{$_}++ } @_;
}
#####################################
#####################################
#Merge data from reading data from the directory
sub MergeDataFromDirectory {
	my @lane            = @$lane;
	my @sampleId        = @$sampleId;
	my @index           = @$index;
	my @titleBarcode    = @$barcode;
	my @titlePool       = @$pool;
	my @titleSampleId   = @$titleSampleId;
	my %barcodes        = ();
	my %indexHash       = ();
	my %titleInfo       = ();
	my @notifyNames     = ();
	my @checkErrorFiles = ();
	my $newIndex;
	my $name;
	my $Null           = "NULL";
	my @parseFilenames = ();
	my $now            = time;
	open( BARCODEFILE, $barcodeFile )
	  or
	  die( $logger->fatal("BarcodeFile:Cannot open $barcodeFile. Error:$!") );

	while (<BARCODEFILE>) {
		next if ( $. == 1 );
		my @dataCols = split( "\t", $_ );
		$dataCols[0] =~ s/\s//g;
		$dataCols[1] =~ s/\s//g;
		$barcodes{ $dataCols[0] }  = $dataCols[1];
		$indexHash{ $dataCols[1] } = $dataCols[0];
	}
	close(BARCODEFILE);
	$logger->info("Running merge jobs on SGE/LSF");
	if ( $fastqSource eq "DMP" ) {
		foreach my $i ( 0 .. $#titleBarcode ) {
			my $newIndex = $titleBarcode[$i];
			my $name =
			    $titleSampleId[$i] . "_"
			  . $titleBarcode[$i] . "_"
			  . $titlePool[$i];
			my $read1ListName = "";
			my $read2ListName = "";
			foreach my $j ( 0 .. $#sampleId ) {
				if ( $index[$j] eq $indexHash{ $titleBarcode[$i] } ) {
					$read1ListName .=
					    $datadir . "/"
					  . $sampleId[$j] . "_"
					  . $index[$j] . "_L00"
					  . $lane[$j]
					  . "_R1_001.fastq.gz ";
					$read2ListName .=
					    $datadir . "/"
					  . $sampleId[$j] . "_"
					  . $index[$j] . "_L00"
					  . $lane[$j]
					  . "_R2_001.fastq.gz ";
				}
			}
			my $read1Name = $outdir . "/" . $name . "_L000_R1_mrg.fastq.gz";
			my $read2Name = $outdir . "/" . $name . "_L000_R2_mrg.fastq.gz";
			if (    ( -e $read1Name )
				and ( ( -s $read1Name ) != 0 )
				and ( -e $read2Name )
				and ( ( -s $read2Name ) != 0 ) )
			{
				$logger->warn(
"Files:\n$read1Name\n$read2Name\n they exists and process will not run to merge files."
				);
				push( @notifyNames,     $Null );
				push( @checkErrorFiles, $Null );
				push( @notifyNames,     $Null );
				push( @checkErrorFiles, $Null );
				push( @parseFilenames,  "$read1Name,$read2Name" );
			}
			else {
				eval {
					if ( $CLUSTER eq "SGE" )
					{
						my $cmd =
"$QSUB -q $queue -V -wd $outdir -N MergeRead1.$newIndex.$i.$$ -e MergeRead1.$newIndex.$i.$$.stderr -o /dev/null -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y '$ZCAT $read1ListName | gzip > $read1Name'";
						$logger->debug("COMMAND: $cmd");

						#Read1
`$QSUB -q $queue -V -wd $outdir -N MergeRead1.$newIndex.$i.$$ -e MergeRead1.$newIndex.$i.$$.stderr -o /dev/null -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$ZCAT $read1ListName | gzip > $read1Name"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid MergeRead1.$newIndex.$i.$$ -N NotifyMR.Read1.$i.$$ -e NotifyMR.Read1.$i.$$.stderr -o NotifyMR.Read1.$i.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;

						#Read2
						$cmd =
"$QSUB -q $queue -V -wd $outdir -N MergeRead2.$newIndex.$i.$$ -e MergeRead2.$newIndex.$i.$$.stderr -o /dev/null -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y '$ZCAT $read2ListName | gzip > $read2Name'";
						$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N MergeRead2.$newIndex.$i.$$ -e MergeRead2.$newIndex.$i.$$.stderr -o /dev/null -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$ZCAT $read2ListName | gzip > $read2Name"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid MergeRead2.$newIndex.$i.$$ -N NotifyMR.Read2.$i.$$ -e NotifyMR.Read2.$i.$$.stderr -o NotifyMR.Read2.$i.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
					}
					else {

						#Read1
						my $cmd =
"$BSUB -q $queue -cwd $outdir -J MergeRead1.$newIndex.$i.$$ -e MergeRead1.$newIndex.$i.$$.stderr -R \"rusage[mem=2]\" -M 4 -n 1 \"$ZCAT $read1ListName | gzip > $read1Name\"";
						$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J MergeRead1.$newIndex.$i.$$ -e MergeRead1.$newIndex.$i.$$.stderr -We 0:59 -R "rusage[mem=2]" -M 4 -n 1 "$ZCAT $read1ListName | gzip > $read1Name"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(MergeRead1.$newIndex.$i.$$)" -J NotifyMR.Read1.$i.$$ -e NotifyMR.Read1.$i.$$.%J.stderr -o NotifyMR.Read1.$i.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
						$cmd =
"$BSUB -q $queue -cwd $outdir -J MergeRead2.$newIndex.$i.$$ -e MergeRead2.$newIndex.$i.$$.stderr -R \"rusage[mem=2]\" -M 4 -n 1 \"$ZCAT $read2ListName | gzip > $read2Name\"";
						$logger->debug("COMMAND: $cmd");

						#Read2
`$BSUB -q $queue -cwd $outdir -J MergeRead2.$newIndex.$i.$$ -e MergeRead2.$newIndex.$i.$$.stderr -We 0:59 -R "rusage[mem=2]" -M 4 -n 1 "$ZCAT $read2ListName | gzip > $read2Name"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(MergeRead2.$newIndex.$i.$$)" -J NotifyMR.Read2.$i.$$ -e NotifyMR.Read2.$i.$$.%J.stderr -o NotifyMR.Read2.$i.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
					}
				};
				if ($@) {
					$logger->fatal("MergeFastq:Job Sumission failed");
					exit(1);
				}
				push( @notifyNames,     "NotifyMR.Read1.$i.$$.stat" );
				push( @checkErrorFiles, "MergeRead1.$newIndex.$i.$$.stderr" );
				push( @notifyNames,     "NotifyMR.Read2.$i.$$.stat" );
				push( @checkErrorFiles, "MergeRead2.$newIndex.$i.$$.stderr" );
				push( @parseFilenames,  "$read1Name,$read2Name" );
			}
		}
		&WaitToFinish( $outdir, @notifyNames );
		&CheckErrorFiles( $outdir, @checkErrorFiles );
		$now = time - $now;
		$logger->info("Finished running merge jobs on SGE");
		printf(
			"Total Merge Fastq run time: %02d:%02d:%02d\n\n",
			int( $now / 3600 ),
			int( ( $now % 3600 ) / 60 ),
			int( $now % 60 )
		);
		my (@sortedparseFilenames) =
		  sort { lowestNumber($a) <=> lowestNumber($b) } @parseFilenames;

		#Check fastq files
		my $prog;
		&CheckOutputFiles( $outdir, $prog, @sortedparseFilenames );
		return ( \@sortedparseFilenames );
	}
	elsif ( $fastqSource == "GCL" ) {
		for ( my $i = 0 ; $i < scalar(@titleBarcode) ; $i++ ) {
			$titleInfo{ $titleBarcode[$i] } =
			    $titleSampleId[$i] . "_"
			  . $titleBarcode[$i] . "_"
			  . $titlePool[$i];
		}
		@sampleId = uniq(@sampleId);
		@index    = uniq(@index);
		for ( my $sampleNum = 0 ;
			$sampleNum < scalar(@sampleId) ; $sampleNum++ )
		{
			my $read1ListName =
			    $datadir . "/"
			  . $sampleId[$sampleNum] . "_"
			  . $index[$sampleNum] . "_L00*"
			  . "_R1_*.fastq.gz";
			my $read2ListName =
			    $datadir . "/"
			  . $sampleId[$sampleNum] . "_"
			  . $index[$sampleNum] . "_L00*"
			  . "_R2_*.fastq.gz";
			if ( exists $barcodes{ $index[$sampleNum] } ) {
				$newIndex = $barcodes{ $index[$sampleNum] };
				if ( exists $titleInfo{$newIndex} ) {
					$name = $titleInfo{$newIndex};
				}
				else {
					$logger->fatal(
"MergeFastq:The barcode $newIndex doesnot exists in the title file. Cannot move ahead. Please check and rerun."
					);
					exit(1);
				}
			}
			else {
				$logger->fatal(
"MergeFastq:The barcode sequence $barcodes{$index[$sampleNum]} does not exists in barcode file. Cannot move ahead. Please check and rerun."
				);
				exit(1);
			}

			my $read1Name =
			  $outdir . "/" . $name . "_L000" . "_R1_mrg.fastq.gz";
			my $read2Name =
			  $outdir . "/" . $name . "_L000" . "_R2_mrg.fastq.gz";

			#Run the qsub command to merge the files.
			#Read1
			if (    ( -e $read1Name )
				and ( ( -s $read1Name ) != 0 )
				and ( -e $read2Name )
				and ( ( -s $read2Name ) != 0 ) )
			{
				$logger->warn(
"Files:\n$read1Name\n$read2Name\n they exists and process will not run to merge files."
				);
				push( @notifyNames,     $Null );
				push( @checkErrorFiles, $Null );
				push( @notifyNames,     $Null );
				push( @checkErrorFiles, $Null );
				push( @parseFilenames,  "$read1Name,$read2Name" );
				next;
			}
			else {
				eval {
					if ( $CLUSTER eq "SGE" )
					{
						my $cmd =
"$QSUB -q $queue -V -wd $outdir -N MergeRead1.$newIndex.$sampleNum.$$ -e MergeRead1.$newIndex.$sampleNum.$$.stderr -o /dev/null -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y \"$ZCAT $read1ListName | gzip > $read1Name\"";
						$logger->debug("COMMAND: $cmd");

						#Read1
`$QSUB -q $queue -V -wd $outdir -N MergeRead1.$newIndex.$sampleNum.$$ -e MergeRead1.$newIndex.$sampleNum.$$.stderr -o /dev/null -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$ZCAT $read1ListName | gzip > $read1Name"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid MergeRead1.$newIndex.$sampleNum.$$ -N NotifyMR.Read1.$sampleNum.$$ -e NotifyMR.Read1.$sampleNum.$$.stderr -o NotifyMR.Read1.$sampleNum.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
						$cmd =
"$QSUB -q $queue -V -wd $outdir -N MergeRead2.$newIndex.$sampleNum.$$ -e MergeRead2.$newIndex.$sampleNum.$$.stderr -o /dev/null -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y \"$ZCAT $read2ListName | gzip > $read2Name\"";
						$logger->debug("COMMAND: $cmd");

						#Read2
`$QSUB -q $queue -V -wd $outdir -N MergeRead2.$newIndex.$sampleNum.$$ -e MergeRead2.$newIndex.$sampleNum.$$.stderr -o /dev/null -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$ZCAT $read2ListName | gzip > $read2Name"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid MergeRead2.$newIndex.$sampleNum.$$ -N NotifyMR.Read2.$sampleNum.$$ -e NotifyMR.Read2.$sampleNum.$$.stderr -o NotifyMR.Read2.$sampleNum.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
					}
					else {

						#Read1
						my $cmd =
"$BSUB -q $queue -cwd $outdir -J MergeRead1.$newIndex.$sampleNum.$$ -e MergeRead1.$newIndex.$sampleNum.$$.stderr -R \"rusage[mem=2]\" -M 4 -n 1 \"$ZCAT $read1ListName | gzip > $read1Name\"";
						$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J MergeRead1.$newIndex.$sampleNum.$$ -e MergeRead1.$newIndex.$sampleNum.$$.stderr -We 0:59 -R "rusage[mem=2]" -M 4 -n 1 "$ZCAT $read1ListName | gzip > $read1Name"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(MergeRead1.$newIndex.$sampleNum.$$)" -J MergeRead1.$newIndex.$sampleNum.$$ -e NotifyMR.Read1.$sampleNum.$$.%J.stderr -o NotifyMR.Read1.$sampleNum.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;

						#Read2
						$cmd =
"$BSUB -q $queue -cwd $outdir -J MergeRead2.$newIndex.$sampleNum.$$ -e MergeRead2.$newIndex.$sampleNum.$$.stderr -R \"rusage[mem=2]\" -M 4 -n 1 \"$ZCAT $read2ListName | gzip > $read2Name\"";
						$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J MergeRead2.$newIndex.$sampleNum.$$ -e MergeRead2.$newIndex.$sampleNum.$$.stderr -We 0:59 -R "rusage[mem=2]" -M 4 -n 1 "$ZCAT $read2ListName | gzip > $read2Name"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(MergeRead2.$newIndex.$sampleNum.$$)" -J MergeRead2.$newIndex.$sampleNum.$$ -e NotifyMR.Read2.$sampleNum.$$.%J.stderr -o NotifyMR.Read2.$sampleNum.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
					}
				};
				if ($@) {
					$logger->fatal("MergeFastq:Job Sumission failed");
					exit(1);
				}
				push( @notifyNames, "NotifyMR.Read1.$sampleNum.$$.stat" );
				push( @checkErrorFiles,
					"MergeRead1.$newIndex.$sampleNum.$$.err" );
				push( @notifyNames, "NotifyMR.Read2.$sampleNum.$$.stat" );
				push( @checkErrorFiles,
					"MergeRead2.$newIndex.$sampleNum.$$.err" );
				push( @parseFilenames, "$read1Name,$read2Name" );
			}
		}
		&WaitToFinish( $outdir, @notifyNames );
		&CheckErrorFiles( $outdir, @checkErrorFiles );
		$now = time - $now;
		$logger->info("Finished running merge jobs on SGE");
		printf(
			"Total Merge Fastq run time: %02d:%02d:%02d\n\n",
			int( $now / 3600 ),
			int( ( $now % 3600 ) / 60 ),
			int( $now % 60 )
		);
		my (@sortedparseFilenames) =
		  sort { lowestNumber($a) <=> lowestNumber($b) } @parseFilenames;

		#Check fastq files
		my $prog;
		&CheckOutputFiles( $outdir, $prog, @sortedparseFilenames );
		return ( \@sortedparseFilenames );
	}
	else {
		for ( my $sampleNum = 0 ; $sampleNum < scalar(@titleBarcode) ; $sampleNum++ ) {
			my $name =
			    $titleSampleId[$sampleNum] . "_"
			  . $titleBarcode[$sampleNum] . "_"
			  . $titlePool[$sampleNum];

			my $read1ListName =
			    $datadir . "/"
			  . $titleSampleId[$sampleNum] . "_" . "*" . "_L00*"
			  . "_R1_*.fastq.gz";
			my $read2ListName =
			    $datadir . "/"
			  . $titleSampleId[$sampleNum] . "_" . "*" . "_L00*"
			  . "_R2_*.fastq.gz";

			my $read1Name =
			  $outdir . "/" . $name . "_L000" . "_R1_mrg.fastq.gz";
			my $read2Name =
			  $outdir . "/" . $name . "_L000" . "_R2_mrg.fastq.gz";

			#Run the qsub command to merge the files.
			#Read1
			if (    ( -e $read1Name )
				and ( ( -s $read1Name ) != 0 )
				and ( -e $read2Name )
				and ( ( -s $read2Name ) != 0 ) )
			{
				$logger->warn(
"Files:\n$read1Name\n$read2Name\n they exists and process will not run to merge files."
				);
				push( @notifyNames,     $Null );
				push( @checkErrorFiles, $Null );
				push( @notifyNames,     $Null );
				push( @checkErrorFiles, $Null );
				push( @parseFilenames,  "$read1Name,$read2Name" );
				next;
			}
			else {
				eval {
					if ( $CLUSTER eq "SGE" )
					{
						my $cmd =
"$QSUB -q $queue -V -wd $outdir -N MergeRead1.$newIndex.$sampleNum.$$ -e MergeRead1.$newIndex.$sampleNum.$$.stderr -o /dev/null -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y \"$ZCAT $read1ListName | gzip > $read1Name\"";
						$logger->debug("COMMAND: $cmd");

						#Read1
`$QSUB -q $queue -V -wd $outdir -N MergeRead1.$newIndex.$sampleNum.$$ -e MergeRead1.$newIndex.$sampleNum.$$.stderr -o /dev/null -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$ZCAT $read1ListName | gzip > $read1Name"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid MergeRead1.$newIndex.$sampleNum.$$ -N NotifyMR.Read1.$sampleNum.$$ -e NotifyMR.Read1.$sampleNum.$$.stderr -o NotifyMR.Read1.$sampleNum.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
						$cmd =
"$QSUB -q $queue -V -wd $outdir -N MergeRead2.$newIndex.$sampleNum.$$ -e MergeRead2.$newIndex.$sampleNum.$$.stderr -o /dev/null -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y \"$ZCAT $read2ListName | gzip > $read2Name\"";
						$logger->debug("COMMAND: $cmd");

						#Read2
`$QSUB -q $queue -V -wd $outdir -N MergeRead2.$newIndex.$sampleNum.$$ -e MergeRead2.$newIndex.$sampleNum.$$.stderr -o /dev/null -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$ZCAT $read2ListName | gzip > $read2Name"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid MergeRead2.$newIndex.$sampleNum.$$ -N NotifyMR.Read2.$sampleNum.$$ -e NotifyMR.Read2.$sampleNum.$$.stderr -o NotifyMR.Read2.$sampleNum.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
					}
					else {

						#Read1
						my $cmd =
"$BSUB -q $queue -cwd $outdir -J MergeRead1.$newIndex.$sampleNum.$$ -e MergeRead1.$newIndex.$sampleNum.$$.stderr -R \"rusage[mem=2]\" -M 4 -n 1 \"$ZCAT $read1ListName | gzip > $read1Name\"";
						$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J MergeRead1.$newIndex.$sampleNum.$$ -e MergeRead1.$newIndex.$sampleNum.$$.stderr -We 0:59 -R "rusage[mem=2]" -M 4 -n 1 "$ZCAT $read1ListName | gzip > $read1Name"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(MergeRead1.$newIndex.$sampleNum.$$)" -J MergeRead1.$newIndex.$sampleNum.$$ -e NotifyMR.Read1.$sampleNum.$$.%J.stderr -o NotifyMR.Read1.$sampleNum.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;

						#Read2
						$cmd =
"$BSUB -q $queue -cwd $outdir -J MergeRead2.$newIndex.$sampleNum.$$ -e MergeRead2.$newIndex.$sampleNum.$$.stderr -R \"rusage[mem=2]\" -M 4 -n 1 \"$ZCAT $read2ListName | gzip > $read2Name\"";
						$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J MergeRead2.$newIndex.$sampleNum.$$ -e MergeRead2.$newIndex.$sampleNum.$$.stderr -We 0:59 -R "rusage[mem=2]" -M 4 -n 1 "$ZCAT $read2ListName | gzip > $read2Name"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(MergeRead2.$newIndex.$sampleNum.$$)" -J MergeRead2.$newIndex.$sampleNum.$$ -e NotifyMR.Read2.$sampleNum.$$.%J.stderr -o NotifyMR.Read2.$sampleNum.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
					}
				};
				if ($@) {
					$logger->fatal("MergeFastq:Job Sumission failed");
					exit(1);
				}
				push( @notifyNames, "NotifyMR.Read1.$sampleNum.$$.stat" );
				push( @checkErrorFiles,
					"MergeRead1.$newIndex.$sampleNum.$$.err" );
				push( @notifyNames, "NotifyMR.Read2.$sampleNum.$$.stat" );
				push( @checkErrorFiles,
					"MergeRead2.$newIndex.$sampleNum.$$.err" );
				push( @parseFilenames, "$read1Name,$read2Name" );
			}
		}
		&WaitToFinish( $outdir, @notifyNames );
		&CheckErrorFiles( $outdir, @checkErrorFiles );
		$now = time - $now;
		$logger->info("Finished running merge jobs on SGE");
		printf(
			"Total Merge Fastq run time: %02d:%02d:%02d\n\n",
			int( $now / 3600 ),
			int( ( $now % 3600 ) / 60 ),
			int( $now % 60 )
		);
		my (@sortedparseFilenames) =
		  sort { lowestNumber($a) <=> lowestNumber($b) } @parseFilenames;

		#Check fastq files
		my $prog;
		&CheckOutputFiles( $outdir, $prog, @sortedparseFilenames );
		return ( \@sortedparseFilenames );
	}
}
#####################################
#####################################
#sort by barcode name:
sub lowestNumber {
	my $files = shift;
	my @filenames = split( ",", $files );
	my ($number) = $filenames[0] =~ m/.*_bc(\d+)_.*/g;
	return $number;
}
#####################################
#####################################
#Do Mapping which includes:
#Adapter Clipping;Matching PE fastq file
#Run BWA
#Filter ProperPairs
#Mark Duplicates
sub DoMapping {
	my ($filenames) = @_;
	my @names = ();
	if ($filenames) { (@names) = @$filenames; }
	if ( ( scalar(@names) == 0 ) and ($fof) ) {
		my @fnames = &GetNames( $fof, $outdir );
		@names = &MAKEPAIRS( \@fnames, $outdir );
	}
	my (@sortedparseFilenames) =
	  sort { lowestNumber($a) <=> lowestNumber($b) } @names;
	@names = @sortedparseFilenames;
	my @clippedFilenames = ();

	#my @SAFilenames = ();
	my @SamFilenames       = ();
	my @sortedBamFilenames = ();
	my @notifyNames        = ();
	my $now                = time;
	my %adaptorList        = ();
	open( ADAPTORFILE, $adaptorFile )
	  or die( $logger->fatal("DoMapping:Cannot open $adaptorFile,Error:$!") );
	while (<ADAPTORFILE>) {
		my @dataCols = split( "\t", $_ );
		$dataCols[0] =~ s/\s//g;
		$dataCols[1] =~ s/\s//g;
		$adaptorList{ $dataCols[0] } = $dataCols[1];
	}
	close(ADAPTORFILE);

	#Running Cutadapt through Trim Galore
	$logger->info("Started runing clipping jobs on SGE");
	for ( my $i = 0 ; $i < scalar @names ; $i++ ) {
		my ( $file1, $file2 ) = split( ",", $names[$i] );

		#print "$file1\n$file2\n";
		my ( $read1clipped, $read2clipped, $notifyname ) =
		  &RunTrimGalore( $file1, $file2, $outdir, \%adaptorList, $i );
		push( @notifyNames,      $notifyname );
		push( @clippedFilenames, "$read1clipped,$read2clipped" );
	}

	#waiting for adapter trimming to finish
	&WaitToFinish( $outdir, @notifyNames );
	$logger->info("Finished running clipping jobs on SGE");
	$now = time - $now;
	printf(
		"Total Adaptor Clipping run time: %02d:%02d:%02d\n\n",
		int( $now / 3600 ),
		int( ( $now % 3600 ) / 60 ),
		int( $now % 60 )
	);

	#Check clipped read files
	my $prog;
	&CheckOutputFiles( $outdir, $prog, @clippedFilenames );

	#Running BwaMem
	$now = time;
	$logger->info("Started runing bwa mem jobs on SGE");
	@notifyNames = ();
	for ( my $i = 0 ; $i < scalar(@clippedFilenames) ; $i++ ) {
		my ( $file1, $file2 ) = split( ",", $clippedFilenames[$i] );

		#print "$file1\n$file2\n";
		my ( $samFilename, $notifyname ) =
		  &RunBwaMem( $file1, $file2, $outdir, $i );
		push( @notifyNames, $notifyname );
		$samFilename = basename($samFilename);
		push( @SamFilenames, "$samFilename" );
	}

	#waiting for bwa aln to finish
	&WaitToFinish( $outdir, @notifyNames );
	$logger->info("Finished running bwa mem jobs on SGE");
	$now = time - $now;
	printf(
		"Total running time: %02d:%02d:%02d\n\n",
		int( $now / 3600 ),
		int( ( $now % 3600 ) / 60 ),
		int( $now % 60 )
	);

	#Check bwa sampe files
	&CheckOutputFiles( $outdir, $prog, @SamFilenames );

	#Run Sort Sam
	$now = time;
	$logger->info("Started running Sort Sam jobs on SGE");
	@notifyNames = ();
	for ( my $i = 0 ; $i < scalar(@SamFilenames) ; $i++ ) {
		my ( $sortedBamFile, $notifyname ) =
		  &RunSortSam( $SamFilenames[$i], $outdir, $i );
		push( @notifyNames,        $notifyname );
		push( @sortedBamFilenames, $sortedBamFile );
	}

	#waiting for sort sam to finish
	&WaitToFinish( $outdir, @notifyNames );
	$now = time - $now;
	$logger->info("Finished running Sort Sam jobs on SGE");

	#Check sort sam files
	&CheckOutputFiles( $outdir, $prog, @sortedBamFilenames );
	#################
	###Files to delete
	my ( @filesToDelete, @srtSam, @clfastq, @mrgfastq ) = ();
	eval {
		(@srtSam)   = `ls $outdir/*.sam`;
		(@clfastq)  = `ls $outdir/*_cl.fastq.gz`;
		(@mrgfastq) = `ls $outdir/*_mrg.fastq.gz`;
	};
	if ($@) {
		$logger->warn(
			"Cannot list files for deletion at mapping step.Error:$@");
	}
	my $flag = 0;
	foreach my $file (@sortedBamFilenames) {
		if ( ( -e "$outdir/$file" ) and ( ( -s "$outdir/$file" ) != 0 ) ) {
			next;
		}
		else {
			$flag = 1;
		}
	}
	@filesToDelete = ( @srtSam, @clfastq, @mrgfastq );
	if ( ( $flag == 0 ) and ( $deleteIntermediateFiles == 1 ) ) {
		$logger->info("Deleting files from mapping step.");
		&DeleteFiles( \@filesToDelete );
	}
	elsif (
		( $flag == 1 )
		and (  ( $deleteIntermediateFiles == 1 )
			or ( $deleteIntermediateFiles == 2 ) )
	  )
	{
		$logger->fatal(
"Sorry, some of the bam files are not present from Mapping step, thus files will not be deleted please review and then rerun it"
		);
		exit(1);
	}
	elsif ( ( $flag == 0 ) and ( $deleteIntermediateFiles == 2 ) ) {
		$logger->info("Files from mapping step will not be deleted.");
	}
	else {
		$logger->info("Deleting files from mapping step.");
		&DeleteFiles( \@filesToDelete );
	}
	printf(
		"Total Mapping step run time: %02d:%02d:%02d\n\n",
		int( $now / 3600 ),
		int( ( $now % 3600 ) / 60 ),
		int( $now % 60 )
	);
	return ( \@sortedBamFilenames );
}
#####################################
#####################################
#Do Mark Duplicates Realignment and Recalibration which includes:
#Target Creator
#Indel Realignment
#Mark Duplicates
#Table Recalibration
#Count Covariates
#BAM File Metrics Calculation
sub ProcessBams {
	my ($filenames) = @_;
	my @names = ();
	if ($filenames) { (@names) = @$filenames; }
	if ( ( scalar @names == 0 ) and ($fof) ) {
		@names = &GetNames( $fof, $outdir );
	}
	my (@sortedparseFilenames) =
	  sort { lowestNumber($a) <=> lowestNumber($b) } @names;
	@names = @sortedparseFilenames;
	my @notifyNames = ();
	tie( my %groupedFilenames,        'Tie::IxHash' );
	tie( my %indelRealignerIntervals, 'Tie::IxHash' );
	my @MarkDuplicatesBamFilenames  = ();
	my @RTCintervalFiles            = ();
	my @realignedBams               = ();
	my @sortedRealignedBamFilenames = ();
	my @recalTableFiles             = ();
	my @recalibratedBams            = ();
	my $now                         = time;
	##################
	#Run Mark Duplicates
	$logger->info("Started running Mark Duplicates jobs on SGE");
	for ( my $i = 0 ; $i < scalar(@names) ; $i++ ) {
		my ( $MarkDuplicatesBamFile, $notifyname ) =
		  &RunMarkDuplicates( $names[$i], $outdir, $i );
		push( @notifyNames,                $notifyname );
		push( @MarkDuplicatesBamFilenames, $MarkDuplicatesBamFile );
	}

	#waiting for Mark Duplicates to finish
	&WaitToFinish( $outdir, @notifyNames );
	$now = time - $now;
	$logger->info("Finished running Mark Duplicates jobs on SGE");
	printf(
		"Total MarkDuplicates run time: %02d:%02d:%02d\n\n",
		int( $now / 3600 ),
		int( ( $now % 3600 ) / 60 ),
		int( $now % 60 )
	);

	#Check Mark duplicate files
	my $prog;
	&CheckOutputFiles( $outdir, $prog, @MarkDuplicatesBamFilenames );

	#Run ABRA or GATK
	if ( $runABRA == 1 ) {

		#Group files of Indel Realignment
		my $fCount = 0;
		foreach my $file (@MarkDuplicatesBamFilenames) {
			my $irFile = $file;
			$irFile =~ s/\.bam/_IR\.bam/;
			push( @realignedBams, $irFile );
			my $indelRealignerList =
			  @$patientId[$fCount] . "_IndelRealignerInput.list";
			if ( exists $groupedFilenames{$indelRealignerList} ) {
				my $files = $groupedFilenames{$indelRealignerList};
				$files = "$files" . ",$file";
				$groupedFilenames{$indelRealignerList} = "$files";
			}
			else {
				$groupedFilenames{$indelRealignerList} = "$file";
			}
			$indelRealignerIntervals{$file} = $indelRealignerList;
			$fCount++;
		}
		####################
		#ABRA Realignment
		@notifyNames = ();
		$now         = time;
		$logger->info("Started running ABRA Realignment jobs on SGE");
		my $pcount = 0;
		while ( ( my $outputListFile, my $bamFiles ) = each(%groupedFilenames) )
		{

			#print "$outputListFile=>$bamFiles\n\n";
			my $inputFilename = $outputListFile;
			my ($patient_ID) =
			  $inputFilename =~ /(.*)_IndelRealignerInput\.list/;

			#$inputFilename =~ s/\.intervals/Input\.list/;
			#print "$inputFilename\n";
			open( FH, ">", "$outdir/$inputFilename" )
			  or die(
				$logger->fatal(
					"ABRA_OpenFile:Cannot open $inputFilename, ERROR:$!")
			  );
			my @data = split( ",", $bamFiles );
			foreach my $file (@data) {
				print FH "$outdir/$file\n";
			}

			#print FH "$bamFiles";
			close(FH);
			my ($notifyname) =
			  &RunABRA( $inputFilename, $patient_ID, $outdir, $pcount );
			push( @notifyNames, $notifyname );
			$pcount++;
			if ( $pcount >= 8 ) {
				&WaitToFinish( $outdir, @notifyNames );
				@notifyNames = ();
				$pcount      = 0;
			}
		}
		&WaitToFinish( $outdir, @notifyNames );
		$now = time - $now;
		$logger->info("Finished running ABRA Realignment jobs on SGE");
		printf(
			"Total ABRA run time: %02d:%02d:%02d\n\n",
			int( $now / 3600 ),
			int( ( $now % 3600 ) / 60 ),
			int( $now % 60 )
		);

		#Run FixMate
		$now = time;
		$logger->info("Started running Fix Mate jobs on SGE");
		@notifyNames = ();
		for ( my $i = 0 ; $i < scalar(@realignedBams) ; $i++ ) {
			my ( $sortedRealignedBamFile, $notifyname ) =
			  &RunFixMateInformation( $realignedBams[$i], $outdir, $i );
			push( @notifyNames,                 $notifyname );
			push( @sortedRealignedBamFilenames, $sortedRealignedBamFile );
		}

		#waiting for sort sam to finish
		&WaitToFinish( $outdir, @notifyNames );
		$now = time - $now;
		$logger->info("Finished running Sort Sam jobs on SGE");
	}
	else {

		#Group files of Indel Realignment
		my $fCount = 0;
		foreach my $file (@MarkDuplicatesBamFilenames) {
			my $indelIntervals =
			  @$patientId[$fCount] . "_IndelRealigner.intervals";
			if ( exists $groupedFilenames{$indelIntervals} ) {
				my $files = $groupedFilenames{$indelIntervals};
				$files = "$files" . ",$file";
				$groupedFilenames{$indelIntervals} = "$files";
			}
			else {
				$groupedFilenames{$indelIntervals} = "$file";
			}
			$indelRealignerIntervals{$file} = $indelIntervals;
			$fCount++;
		}
		####################
		#Realign Target Creator
		@notifyNames = ();
		$now         = time;
		$logger->info("Started running Realigner Target Creator jobs on SGE");
		while ( ( my $outputIntervalFile, my $bamFiles ) =
			each(%groupedFilenames) )
		{

			#print "$outputIntervalFile=>$bamFiles\n\n";
			my $inputFilename = $outputIntervalFile;
			$inputFilename =~ s/\.intervals/Input\.list/;

			#print "$inputFilename\n";
			open( FH, ">", "$outdir/$inputFilename" )
			  or die(
				$logger->fatal(
					"RTC_OpenFile:Cannot open $inputFilename, ERROR:$!")
			  );
			my @data = split( ",", $bamFiles );
			foreach my $file (@data) {
				print FH "$outdir/$file\n";
			}

			#print FH "$bamFiles";
			close(FH);
			my ($notifyname) =
			  &RunRealignerTargetCreator( $outputIntervalFile, $inputFilename,
				$outdir );
			push( @notifyNames,      $notifyname );
			push( @RTCintervalFiles, $outputIntervalFile );
		}
		&WaitToFinish( $outdir, @notifyNames );
		$now = time - $now;
		$logger->info("Finished running Realigner Target Creator jobs on SGE");
		printf(
			"Total Realign Target Creator run time: %02d:%02d:%02d\n\n",
			int( $now / 3600 ),
			int( ( $now % 3600 ) / 60 ),
			int( $now % 60 )
		);

		#Check Realign Target Creator files
		&CheckOutputFiles( $outdir, $prog, @RTCintervalFiles );
		##################
		#Indel Realignment
		@notifyNames = ();
		$now         = time;
		$logger->info("Started running Indel Realignment jobs on SGE");
		for ( my $i = 0 ; $i < scalar(@MarkDuplicatesBamFilenames) ; $i++ ) {
			my $indelRealignInterval =
			  $indelRealignerIntervals{ $MarkDuplicatesBamFilenames[$i] };
			my ( $bamFile, $notifyname ) =
			  &RunIndelRealigner( $MarkDuplicatesBamFilenames[$i],
				$indelRealignInterval, $outdir, $i );
			push( @notifyNames,   $notifyname );
			push( @realignedBams, $bamFile );
		}
		&WaitToFinish( $outdir, @notifyNames );
		$now = time - $now;
		$logger->info("Finished running Indel Realignment jobs on SGE");
		printf(
			"Total Indel Realignment run time: %02d:%02d:%02d\n\n",
			int( $now / 3600 ),
			int( ( $now % 3600 ) / 60 ),
			int( $now % 60 )
		);
		@sortedRealignedBamFilenames = @realignedBams;

		#Check Realigned bam files
		&CheckOutputFiles( $outdir, $prog, @sortedRealignedBamFilenames );
	}
	################
	#Base Quality Recalibration (BQSR)
	@notifyNames = ();

	#Sample level Processing
	if ( $runBQSR == 1 ) {
		$now = time;
		$logger->info(
			"Started running Base Quality Score Recalibration jobs on SGE");
		for ( my $i = 0 ; $i < scalar(@sortedRealignedBamFilenames) ; $i++ ) {
			my ( $recalTable, $notifyname ) =
			  &RunBaseQualityRecalibration( $sortedRealignedBamFilenames[$i],
				$outdir, $i );
			push( @notifyNames,     $notifyname );
			push( @recalTableFiles, $recalTable );
		}
		&WaitToFinish( $outdir, @notifyNames );
		$now = time - $now;
		$logger->info(
			"Finished running Base Quality Score Recalibration jobs on SGE");
		printf(
"\n\nTotal Base Quality Score Recalibration run time: %02d:%02d:%02d\n\n",
			int( $now / 3600 ),
			int( ( $now % 3600 ) / 60 ),
			int( $now % 60 )
		);

		#Check Recalibrated table files
		&CheckOutputFiles( $outdir, $prog, @recalTableFiles );
		#################
		#Print BQSR reads
		@notifyNames = ();
		$now         = time;
		$logger->info("Started running print BQSR reads jobs on SGE");
		for ( my $i = 0 ; $i < scalar(@sortedRealignedBamFilenames) ; $i++ ) {
			my ( $recalBamFile, $notifyname ) =
			  &PrintBQSRreads( $sortedRealignedBamFilenames[$i],
				$recalTableFiles[$i], $outdir, $i );
			push( @notifyNames,      $notifyname );
			push( @recalibratedBams, $recalBamFile );
		}
		&WaitToFinish( $outdir, @notifyNames );
		$now = time - $now;
		$logger->info("Finished running print BQSR reads jobs on SGE");
		printf(
			"Total Print BQSR Reads run time: %02d:%02d:%02d\n\n",
			int( $now / 3600 ),
			int( ( $now % 3600 ) / 60 ),
			int( $now % 60 )
		);

		#Check Recalibrated bam files
		&CheckOutputFiles( $outdir, $prog, @recalibratedBams );
	}

	#Lane Level Processing
	elsif ( $runBQSR == 2 ) {
		$now = time;
		$logger->info(
			"Started running Base Quality Score Recalibration jobs on SGE");
		my $allBams = "$outdir/RecalibrationInputBams.list";
		if ( -e $allBams ) {
			$logger->info(
				"File $allBams exists and the program will overwrite it.");
		}
		open( FH, ">", $allBams )
		  or
		  die( $logger->fatal("BQSR:Cannot open AllBAMs:$allBams, Error:$!") );
		foreach my $file (@sortedRealignedBamFilenames) {
			print FH "$outdir/$file\n";
		}
		close(FH);
		my ( $BQSRtable, $notifyname ) =
		  &RunBaseQualityRecalibration( $allBams, $outdir, "LaneLevel" );
		push( @notifyNames,     $notifyname );
		push( @recalTableFiles, $BQSRtable );
		&WaitToFinish( $outdir, @notifyNames );
		$now = time - $now;
		$logger->info(
			"Finished running Base Quality Score Recalibration jobs on SGE");
		printf(
"Total Base Quality Score Recalibration run time: %02d:%02d:%02d\n\n",
			int( $now / 3600 ),
			int( ( $now % 3600 ) / 60 ),
			int( $now % 60 )
		);

		#Check Recalibrated table files
		&CheckOutputFiles( $outdir, $prog, @recalTableFiles );
		#################
		#Print BQSR reads
		@notifyNames = ();
		$now         = time;
		$logger->info("Started running print BQSR reads jobs on SGE");
		for ( my $i = 0 ; $i < scalar(@sortedRealignedBamFilenames) ; $i++ ) {
			my ( $bamFile, $notifyname ) =
			  &PrintBQSRreads( $sortedRealignedBamFilenames[$i],
				$BQSRtable, $outdir, $i );
			push( @notifyNames,      $notifyname );
			push( @recalibratedBams, $bamFile );
		}
		&WaitToFinish( $outdir, @notifyNames );
		$now = time - $now;
		$logger->info("Finished running Print BQSR Reads jobs on SGE");
		printf(
			"Total Print BQSR Reads run time: %02d:%02d:%02d\n\n",
			int( $now / 3600 ),
			int( ( $now % 3600 ) / 60 ),
			int( $now % 60 )
		);

		#Check Recalibrated bam files
		&CheckOutputFiles( $outdir, $prog, @recalibratedBams );
	}

	#No BQSR
	else {
		$logger->warn("BQSR is being skipped");
		@recalibratedBams = @realignedBams;
	}
	#################
	###Files to delete
	my ( @filesToDelete, @srtBam, @MDFiles, @IRFiles, @FXFiles, @InterValList,
		@recalGrp )
	  = ();
	my $flag = 0;
	eval {
		(@srtBam)       = `ls $outdir/*_srt.ba*`;
		(@MDFiles)      = `ls $outdir/*_MD.ba*`;
		(@IRFiles)      = `ls $outdir/*_IR.ba*`;
		(@FXFiles)      = `ls $outdir/*_FX.ba*`;
		(@InterValList) = `ls $outdir/*_IndelRealigner*`;
		(@recalGrp)     = `ls $outdir/*_recalReport.grp`;
	};
	if ($@) {
		$logger->warn(
			"Cannot list files for deletion at processing bam step.Error:$@");
	}
	foreach my $file (@recalibratedBams) {
		if ( ( -e "$outdir/$file" ) and ( ( -s "$outdir/$file" ) != 0 ) ) {
			next;
		}
		else {
			$flag = 1;
		}
	}
	@filesToDelete =
	  ( @srtBam, @MDFiles, @IRFiles, @FXFiles, @InterValList, @recalGrp );
	if ( ( $flag == 0 ) and ( $deleteIntermediateFiles == 1 ) ) {
		$logger->info("Deleting files from processing bams step.");
		&DeleteFiles( \@filesToDelete );
	}
	elsif (
		( $flag == 1 )
		and (  ( $deleteIntermediateFiles == 1 )
			or ( $deleteIntermediateFiles == 2 ) )
	  )
	{
		$logger->fatal(
"Sorry, some of the bam files are not present from processing bams step, thus files will not be deleted please review and then rerun it"
		);
		exit(1);
	}
	elsif ( ( $flag == 0 ) and ( $deleteIntermediateFiles == 2 ) ) {
		$logger->info("Files from processing bams step will not be deleted.");
	}
	else {
		$logger->info("Deleting files from processing bams step.");
		&DeleteFiles( \@filesToDelete );
	}
	$now = time;
	printf(
		"Total run time: %02d:%02d:%02d\n\n",
		int( $now / 3600 ),
		int( ( $now % 3600 ) / 60 ),
		int( $now % 60 )
	);
	return ( \@recalibratedBams );
}
#####################################
#####################################
#This will calculate and compile metrics for BAM files:
sub CalculateAndCompileMetrics {
	my ($filenames) = @_;
	my @names = ();
	if ($filenames) { (@names) = @$filenames; }

	#print "F:$fof\n";
	if ( ( scalar(@names) == 0 ) and ($fof) ) {
		@names = &GetNames( $fof, $outdir );
	}
	my ( @notifyNames, @metricsOutNames ) = ();
	my (@sortedparseFilenames) =
	  sort { lowestNumber($a) <=> lowestNumber($b) } @names;
	@names = @sortedparseFilenames;
	##################
	#Calculate Metrics
	$now = time;
	$logger->info("Started running metrics calculation jobs on SGE");
	for ( my $i = 0 ; $i < scalar(@names) ; $i++ ) {
		my ( $metricsOutFiles, $waitFileNames ) =
		  &RunMetricsCalculations( $names[$i], $outdir, $i );
		foreach my $waitName (@$waitFileNames) {
			push( @notifyNames, $waitName );
		}
		foreach my $outName (@$metricsOutFiles) {
			push( @metricsOutNames, $outName );
		}
	}

	#waiting for metrics calculations to finish
	&WaitToFinish( $outdir, @notifyNames );
	$now = time - $now;
	$logger->info("Finished running metrics calculation jobs on SGE");
	printf(
		"Total Metrics Calculation run time: %02d:%02d:%02d\n\n",
		int( $now / 3600 ),
		int( ( $now % 3600 ) / 60 ),
		int( $now % 60 )
	);

	#Check metrics out files
	my $prog;
	&CheckOutputFiles( $outdir, $prog, @metricsOutNames );
	##################
	#Genotype COSMIC variants on normal samples
	$now = time;
	$logger->info(
		"Started genotyping COSMIC variants on normal samples on SGE");
	my @normalSamples = ();
	my @tumorSamples  = ();
	my %ptNormHash    = ();
	my %tmPtHash      = ();
	my $fCount        = 0;

	foreach my $name (@names) {
		my ($fileBarcode) = $name =~ /.*_(bc\d+)_.*/;
		my $fileClass = $classPerBarcode{$fileBarcode};
		if ( $fileClass =~ m/Normal/i ) {
			push( @normalSamples, $name );
			$ptNormHash{$name} = @$patientId[$fCount];
		}
		if ( $fileClass =~ m/Tumor/i ) {
			push( @tumorSamples, $name );
			$tmPtHash{ @$patientId[$fCount] } = $name;
		}
		$fCount++;
	}
	my @matchedTm = ();
	foreach my $name (@normalSamples) {
		if (   exists $ptNormHash{$name}
			&& exists $tmPtHash{ $ptNormHash{$name} } )
		{
			push @matchedTm, $tmPtHash{ $ptNormHash{$name} };
		}
		else {
			push @matchedTm, "";
		}
	}
	&GenotypeHotspotsOnNormals( \@normalSamples, \@tumorSamples, \@matchedTm,
		$cosmicHotspotsVcf );
	$now = time - $now;
	$logger->info("Finished genotyping hotspot mutations on SGE");
	printf(
		"Total Genotype Hotspots on Normals run time: %02d:%02d:%02d\n\n",
		int( $now / 3600 ),
		int( ( $now % 3600 ) / 60 ),
		int( $now % 60 )
	);
	##################
	# Gentoype expected normal variants on pooled normal sample
	$now = time;
	$logger->info(
		"Started genotyping expected variants on the pooled normal sample");
	my $pooledNormal;
	foreach my $name (@names) {
		my ($fileBarcode) = $name =~ /.*_(bc\d+)_.*/;
		my $fileClass = $classPerBarcode{$fileBarcode};
		if ( $fileClass =~ m/PoolNormal/i ) {
			$pooledNormal = $name;
		}
	}
	if ($pooledNormal) {
		&GenotypePooledNormal( $pooledNormal, $NormalVariantsVCF );
		$now = time - $now;
		$logger->info("Finished genotyping pooled normal sample");
	}
	else {
		$logger->info(
"Pooled normal sample is not present in this run, skipping this step"
		);
	}
	##################
	#Complie Metrics
	$now = time;
	$logger->info("Started running compile metrics calculation jobs on SGE");
	@notifyNames = ();
	my ($notifyname) = &CompileMetrics( \@names, $titleFile, $outdir );
	push( @notifyNames, $notifyname );

	#waiting for compile metrics calculations to finish
	&WaitToFinish( $outdir, @notifyNames );
	$now = time - $now;
	$logger->info("Finished running compile metrics calculation jobs on SGE");
	printf(
		"Total Compile Metrics run time: %02d:%02d:%02d\n\n",
		int( $now / 3600 ),
		int( ( $now % 3600 ) / 60 ),
		int( $now % 60 )
	);
	return ( \@names );
}
#####################################
#####################################
#Genotype COSMIC variants on normal samples
sub GenotypeHotspotsOnNormals {
	my ( $normalSamples, $tumorSamples, $matchedTm, $cosmicHotspotVcf ) = @_;
	my @bamFilesToGenotypeAgainst = ();
	my %descHash                  = ();
	my %matchedTmHash             = ();
	foreach my $i ( 0 .. $#$normalSamples ) {
		my $normalSample = @$normalSamples[$i];
		push @bamFilesToGenotypeAgainst, $normalSample;
		$descHash{$normalSample}      = "normal";
		$matchedTmHash{$normalSample} = @$matchedTm[$i];
	}
	foreach my $tumorSample (@$tumorSamples) {
		push @bamFilesToGenotypeAgainst, $tumorSample;
		$descHash{$tumorSample} = "tumor";
	}
	my @notifyNames    = ();
	my @GenotypedFiles = ();
	my %cosmicNameHash = ();
	my ($GHNout) = $bamFilesToGenotypeAgainst[0] =~ /.*_bc\d+_(.*)_L\d{1,3}.*/;
	$GHNout = $GHNout . "_ALL_genotypehotspotnormals.txt";
	if ( -e "$outdir/$GHNout" ) {
		$logger->warn(
"GenotypeHotSpots:$GHNout file exists and genotyping hotspots will not be run."
		);
		return;
	}
	open IN, "<${cosmicHotspotVcf}"
	  or die(
		$logger->fatal(
			"GenotypeHotSpots:Cannot open file $cosmicHotspotVcf. Error:$!")
	  );
	while (<IN>) {
		next if ( $_ =~ /^\#/ );
		chomp($_);
		my @f = split( '\t', $_ );
		my $key = $f[0] . ":" . $f[1] . ":" . $f[3] . ":" . $f[4];
		$cosmicNameHash{$key} = $f[2] . "|" . $f[7];
	}
	close IN;

	#Traverse through all Mutation Files
	$logger->info("Genotyping COSMIC hotspot mutations");
	my $count    = 0;
	my $totcount = 0;
	my $total    = scalar(@bamFilesToGenotypeAgainst);
	foreach my $bamFile (@bamFilesToGenotypeAgainst) {
		my ($bamID) = $bamFile =~ /(.*)_bc\d+/;
		my $mpileUpMutationOut = $bamFile;
		$mpileUpMutationOut =~ s/\.bam/.mpileup/;
		my $alleleDepthMutationOut = $bamFile;
		$alleleDepthMutationOut =~ s/\.bam/_mpileup\.alleledepth/;
		$descHash{$alleleDepthMutationOut} =
		  $descHash{$bamFile};    # tumor or normal
		if ( exists $matchedTmHash{$bamFile} ) {
			$matchedTmHash{$alleleDepthMutationOut} = $matchedTmHash{$bamFile};
			$matchedTmHash{$alleleDepthMutationOut} =~
			  s/\.bam/_mpileup\.alleledepth/;
		}

		#Mutation File
		if (    ( -e "$outdir/$alleleDepthMutationOut" )
			and ( ( -s "$outdir/$alleleDepthMutationOut" ) != 0 ) )
		{
			push( @GenotypedFiles, $alleleDepthMutationOut );
			$totcount++;
			next;
		}
		else {
			eval {
				if ( $CLUSTER eq "SGE" )
				{
					my $cmd =
"$QSUB -q $queue -V -N AD.$bamID.$$ -wd $outdir -e AD.$bamID.$$.stderr -o AD.$bamID.$$.stdout -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y $PERL $GenotypeAllele -fmv $cosmicHotspotVcf -bam $outdir/$bamFile -rf $Reference -s $SAMTOOLS -o $outdir -of $alleleDepthMutationOut -mof $mpileUpMutationOut -bi $bamID -q $queue";
					$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -N AD.$bamID.$$ -wd $outdir -e AD.$bamID.$$.stderr -o AD.$bamID.$$.stdout -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y $PERL $GenotypeAllele -fmv $cosmicHotspotVcf -bam $outdir/$bamFile -rf $Reference -s $SAMTOOLS -o $outdir -of $alleleDepthMutationOut -mof $mpileUpMutationOut -bi $bamID -mmq $MAPQ -mbq $BASQ -q $queue -b $BEDTOOLS -qsub $QSUB`;
`$QSUB -q $queue -V -wd $outdir -hold_jid AD.$bamID.$$ -N NotityAD.$bamID.$$ -e NotifyAD.$bamID.$$.stderr -o NotifyAD.$bamID.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
				}
				else {
					my $cmd =
"$BSUB -q $queue -J AD.$bamID.$$ -cwd $outdir -e AD.$bamID.$$.%J.stderr -o AD.$bamID.$$.%J.stdout -R \"rusage[mem=2]\" -R \"rusage[iounits=0]\" -M 4 -n 1 \"$PERL $GenotypeAllele -fmv $cosmicHotspotVcf -bam $outdir/$bamFile -rf $Reference -s $SAMTOOLS -o $outdir -of $alleleDepthMutationOut -mof $mpileUpMutationOut -bi $bamID -mmq $MAPQ -mbq $BASQ -q $queue -b $BEDTOOLS -bsub $BSUB\"";
					$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -J AD.$bamID.$$ -cwd $outdir -e AD.$bamID.$$.%J.stderr -o AD.$bamID.$$.%J.stdout -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$PERL $GenotypeAllele -fmv $cosmicHotspotVcf -bam $outdir/$bamFile -rf $Reference -s $SAMTOOLS -o $outdir -of $alleleDepthMutationOut -mof $mpileUpMutationOut -bi $bamID -mmq $MAPQ -mbq $BASQ -q $queue -b $BEDTOOLS -bsub $BSUB"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(AD.$bamID.$$)" -J NotityAD.$bamID.$$ -e NotifyAD.$bamID.$$.%J.stderr -o NotifyAD.$bamID.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
				}
			};
			if ($@) {
				$logger->fatal(
					"GenotypeHotSpot:Job Submission Failed. Error:$!");
				exit(1);
			}
			push( @notifyNames,    "NotifyAD.$bamID.$$.stat" );
			push( @GenotypedFiles, $alleleDepthMutationOut );
			$count++;
			$totcount++;
		}
		if ( ( $count >= 10 ) or ( $totcount == $total ) ) {
			&WaitToFinish( $outdir, @notifyNames );
			@notifyNames = ();
			$count       = 0;
		}
	}

	#my @pstr = ();
	my $normposcnt  = 0;
	my %normposHash = ();
	foreach my $GenotypedFile (@GenotypedFiles) {
		if ( $descHash{$GenotypedFile} eq "tumor" ) { next; }
		my ($samplename) = $GenotypedFile =~ /(.*)_bc\d+_/;
		open IN, "<$outdir/${GenotypedFile}"
		  or die(
			$logger->fatal(
"GenotypeHotSpots:Cannot open file $outdir/$GenotypedFile. Error:$!"
			)
		  );
		while (<IN>) {
			next if ( $_ =~ /^Ref_BAM/ );
			chomp($_);
			my @f      = split( '\t', $_ );
			my $sample = $f[0];
			my $chr    = $f[2];
			my $pos    = $f[3];
			my $ref    = $f[4];
			my $alt    = $f[5];
			my $key    = $chr . ":" . $pos . ":" . $ref . ":" . $alt;
			my $dp     = $f[6];
			my $ad     = $f[8];
			my $vf     = $f[9];

			if (   $dp > 0
				&& ( $vf >= 0.02 or $ad >= 5 )
				&& exists $cosmicNameHash{$key} )
			{

#push @pstr, "$sample\t$cosmicNameHash{$key}\t$chr\t$pos\t$ref\t$alt\t$dp\t$ad\t$vf";
				push @{ $normposHash{$key} },
"$GenotypedFile\t$sample\t$cosmicNameHash{$key}\t$chr\t$pos\t$ref\t$alt\t$dp\t$ad\t$vf";
				$normposcnt++;
			}
		}
		close IN;
		eval { `rm $outdir/${GenotypedFile}`; };
		if ($@) {
			$logger->fatal(
				"GenotypeHotSpot:Cannot remove $outdir/$GenotypedFile. Error:$@"
			);
			exit(1);
		}
	}
	my %tmposHash = ();
	foreach my $GenotypedFile (@GenotypedFiles) {
		if ( $descHash{$GenotypedFile} ne "tumor" ) { next; }
		my ($samplename) = $GenotypedFile =~ /(.*)_bc\d+_/;
		open IN, "<$outdir/${GenotypedFile}"
		  or die(
			$logger->fatal(
"GenotypeHotSpots:Cannot open file $outdir/$GenotypedFile. Error:$!"
			)
		  );
		while (<IN>) {
			next if ( $_ =~ /^Ref_BAM/ );
			chomp($_);
			my @f      = split( '\t', $_ );
			my $sample = $f[0];
			my $chr    = $f[2];
			my $pos    = $f[3];
			my $ref    = $f[4];
			my $alt    = $f[5];
			my $key    = $chr . ":" . $pos . ":" . $ref . ":" . $alt;
			my $dp     = $f[6];
			my $ad     = $f[8];
			my $vf     = $f[9];

			if ( exists $normposHash{$key} ) {
				foreach my $norm_i ( @{ $normposHash{$key} } ) {
					my @g = split( '\t', $norm_i );
					if ( $GenotypedFile eq $matchedTmHash{ $g[0] } ) {
						push @{ $tmposHash{$key} },
"MATCHED_TM\t$sample\t$cosmicNameHash{$key}\t$chr\t$pos\t$ref\t$alt\t$dp\t$ad\t$vf";
					}
					elsif ( $dp > 0 && ( $vf >= 0.02 or $ad >= 5 ) ) {
						push @{ $tmposHash{$key} },
"TM\t$sample\t$cosmicNameHash{$key}\t$chr\t$pos\t$ref\t$alt\t$dp\t$ad\t$vf";
					}
				}
			}
		}
		close IN;
		eval { `rm $outdir/${GenotypedFile}`; };
		if ($@) {
			$logger->fatal(
				"GenotypeHotSpot:Cannot remove $outdir/$GenotypedFile. Error:$@"
			);
			exit(1);
		}
	}
	open GHNOUT, ">$outdir/${GHNout}"
	  or die(
		$logger->fatal(
			"GenotypeHotSpots:Cannot open file $outdir/$GHNout. Error:$!")
	  );

	#if(scalar (@pstr) > 0)
	if ( $normposcnt == 0 ) {
		print GHNOUT "No hotspots detected in any normals\n";
	}
	else {
		foreach my $key ( sort keys %normposHash ) {
			my %toPrint;
			foreach my $norm_i ( @{ $normposHash{$key} } ) {
				my @g = split( '\t', $norm_i );
				$g[0] = "NORMAL";
				print GHNOUT join( "\t", @g ) . "\n";
			}
			foreach my $tm_i ( @{ $tmposHash{$key} } ) {
				$toPrint{$tm_i} = 1;
			}
			foreach my $key ( sort keys %toPrint ) {
				print GHNOUT "$key\n";
			}
			print GHNOUT "\n";
		}
	}
	close GHNOUT;
	return;
	exit;
}
#####################################
#####################################
# Genotype Pooled Normal
sub GenotypePooledNormal {
	my ( $pooledNormal, $NormalVariantsVCF ) = @_;
	my ( %expectedVFs, %observedVFs, @notifyNames );
	my ($GPNout) = $pooledNormal =~ /.*_bc\d+_(.*)_L\d{1,3}.*/;
	$GPNout .= "_ALL_GenotypePooledNormal.txt";
	my $outfile = $pooledNormal;
	$outfile =~ s/.bam/.mpileup.alleledepth/;

	#print "$outdir/$outfile\n";
	if ( !-e $GPNout ) {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -N GenotypePooledNormal.$$ -wd $outdir -e GenotypePooledNormal.$$.stderr -o GenotypePooledNormal.$$.stdout -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y $PERL $GenotypeAllele -fmv $NormalVariantsVCF -bam $outdir/$pooledNormal -rf $Reference -s $SAMTOOLS -o $outdir -of $outfile -q $queue";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -N GenotypePooledNormal.$$ -wd $outdir -e GenotypePooledNormal.$$.stderr -o GenotypePooledNormal.$$.stdout -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y $PERL $GenotypeAllele -fmv $NormalVariantsVCF -bam $outdir/$pooledNormal -rf $Reference -s $SAMTOOLS -o $outdir -of $outfile -mmq $MAPQ -mbq $BASQ -q $queue -b $BEDTOOLS -qsub $QSUB`;
`$QSUB -q $queue -V -wd $outdir -hold_jid GenotypePooledNormal.$$ -N NotifyGenotypePooledNormal.$$ -e NotifyGenotypePooledNormal.$$.stderr -o NotifyGenotypePooledNormal.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -J GenotypePooledNormal.$$ -cwd $outdir -e GenotypePooledNormal.$$.%J.stderr -o GenotypePooledNormal.$$.%J.stdout -R \"rusage[mem=2]\" -R \"rusage[iounits=0]\" -M 4 -n 1 \"$PERL $GenotypeAllele -fmv $NormalVariantsVCF -bam $outdir/$pooledNormal -rf $Reference -s $SAMTOOLS -o $outdir -of $outfile -mmq $MAPQ -mbq $BASQ -q $queue -b $BEDTOOLS -bsub $BSUB\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -J GenotypePooledNormal.$$ -cwd $outdir -e GenotypePooledNormal.$$.%J.stderr -o GenotypePooledNormal.$$.%J.stdout -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$PERL $GenotypeAllele -fmv $NormalVariantsVCF -bam $outdir/$pooledNormal -rf $Reference -s $SAMTOOLS -o $outdir -of $outfile -mmq $MAPQ -mbq $BASQ -q $queue -b $BEDTOOLS -bsub $BSUB"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(GenotypePooledNormal.$$)" -J NotifyGenotypePooledNormal.$$ -e NotifyGenotypePooledNormal.%J.$$.stderr -o NotifyGenotypePooledNormal.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal(
				"GenotypePooledNormal:Job Submission Failed. Error:$!");
			exit(1);
		}
		push( @notifyNames, "NotifyGenotypePooledNormal.$$.stat" );
		&WaitToFinish( $outdir, @notifyNames );
		@notifyNames = ();
		open NORMALVCF, "<", $NormalVariantsVCF
		  || die $logger->fatal(
			" Can not open $NormalVariantsVCF file. Error: $!");
		while (<NORMALVCF>) {
			chomp;
			next if (/^#/);
			my @line  = split("\t");
			my $chr   = $line[0];
			my $pos   = $line[1];
			my $ref   = $line[3];
			my $alt   = $line[4];
			my @annot = split( ";", $line[7] );
			my ($VFexp) = $annot[0] =~ /VF_EXPECTED=(.*)/;
			$expectedVFs{"$chr:$pos:$ref:$alt"} = $VFexp;
		}
		close NORMALVCF;
		open POOLEDALLELEDEPTH, "<", "$outdir/$outfile"
		  || die $logger->fatal(" Can not open $outfile file. Error: $!");
		while (<POOLEDALLELEDEPTH>) {
			chomp;
			next if (/^Ref_BAM\t/);
			my @line  = split("\t");
			my $chr   = $line[2];
			my $pos   = $line[3];
			my $ref   = $line[4];
			my $alt   = $line[5];
			my $VFobs = $line[9];
			$observedVFs{"$chr:$pos:$ref:$alt"} = $VFobs;
		}
		close POOLEDALLELEDEPTH;
		open OUT, ">", "$outdir/$GPNout"
		  || die $logger->info("Can not open $GPNout file: Error: $!");
		print OUT
"Chromosome\tPosition\tReference\tAlternative\tExpectedVF\tObservedVF\n";
		foreach my $key ( sort keys %observedVFs ) {
			my ( $chr, $pos, $ref, $alt ) = split( ":", $key );
			my $VFexp = $expectedVFs{$key};
			my $VFobs = $observedVFs{$key};
			print OUT "$chr\t$pos\t$ref\t$alt\t$VFexp\t$VFobs\n";
		}
		close OUT;
	}
}
#####################################
#####################################
#This will help to call:
#GermLine SNPs: Unified Genotyper
#GermLine Indels: Unified Genotyper
#Somatic SNPs: Mutect
#Somatics Indels: SomaticIndelCaller
sub CallingSNPsAndIndels {
	my ($filenames) = @_;
	my @names = ();
	if ($filenames) { (@names) = @$filenames; }

	#print "F:$fof\n";
	if ( ( scalar(@names) == 0 ) and ($fof) ) {
		@names = &GetNames( $fof, $outdir );
	}
	my @notifyNames  = ();
	my @germLineVcfs = ();
	tie( my %groupedFilenames, 'Tie::IxHash' );
	my @somaticMutationFiles    = ();
	my @somaticIndelFiles       = ();
	my @somaticPindelFiles      = ();
	my @NTCsomaticMutationFiles = ();
	my @NTCsomaticIndelFiles    = ();
	my @NTCsomaticPindelFiles   = ();
	my @CoveragePerSample       = ();
	my @meanCoverageValues      = ();
	tie( my %coverageForNormal, 'Tie::IxHash' );
	tie( my %NormalPerFile,     'Tie::IxHash' );
	my $standardNormal;
	my $now = time;
	my (@sortedparseFilenames) =
	  sort { lowestNumber($a) <=> lowestNumber($b) } @names;
	@names = @sortedparseFilenames;
	my ($poolName) = $names[0] =~ /.*_bc\d+_(.*)_L\d{1,3}.*/;
	my $NormalUsed = $poolName . "_NormalUsedInMutationCalling.txt";

	#Call Somatic SNPs and Indels
	$logger->info("Started running Somatics Variant jobs on SGE");
	my $fCount = 0;

	#Group the files
	foreach my $file (@names) {
		if ( exists $groupedFilenames{ @$patientId[$fCount] } ) {
			my $files = $groupedFilenames{ @$patientId[$fCount] };
			$files = "$files" . ",$file";
			$groupedFilenames{ @$patientId[$fCount] } = "$files";
		}
		else {
			$groupedFilenames{ @$patientId[$fCount] } = "$file";
		}
		$fCount++;
	}

	#Get Mean Coverage from HSmetrics file.
	my $poolNormalMeanCov;
	my $poolNormal;
	foreach my $file (@names) {
		my ($fileBarcode) = $file =~ /.*_(bc\d+)_.*/;
		my $fileClass = $classPerBarcode{$fileBarcode};
		if ( $fileClass =~ m/Normal/i ) {
			if ( $fileClass =~ m/PoolN/i ) {
				$poolNormal = $file;
				my $HSmetricsFile = $file;
				$HSmetricsFile =~ s/\.bam/\.HSmetrics\.txt/g;

				#print "HS:$HSmetricsFile\n";
				open( FH, "$outdir/$HSmetricsFile" )
				  or die(
					$logger->fatal(
"CallSomaticSNP:Cannot Open $outdir/$HSmetricsFile, Error:$!"
					)
				  );
				while (<FH>) {
					next until ( $_ =~ /^BAIT_SET/ );
					while (<FH>) {
						next if ( ( $_ =~ /^BAIT_SET/ ) or ( $_ =~ /^\s$/ ) );
						my (@values) = split( "\t", $_ );
						$poolNormalMeanCov = $values[21];
					}
				}
				close(FH);
				next;
			}
			else {
				my $HSmetricsFile = $file;
				$HSmetricsFile =~ s/\.bam/\.HSmetrics\.txt/g;

				#print "HS:$HSmetricsFile\n";
				open( FH, "$outdir/$HSmetricsFile" )
				  or die(
					$logger->fatal(
"CallSomaticSNP:Cannot Open $outdir/$HSmetricsFile, Error:$!"
					)
				  );
				my $meanCov;
				my $CovForFile;
				while (<FH>) {
					next until ( $_ =~ /^BAIT_SET/ );
					while (<FH>) {
						next if ( ( $_ =~ /^BAIT_SET/ ) or ( $_ =~ /^\s$/ ) );
						my (@values) = split( "\t", $_ );

						#print "MeanCOv:$values[21]\n";
						$CovForFile = $values[21];
						$meanCov    = $values[21];
					}
				}
				close(FH);
				$coverageForNormal{$file} = $CovForFile;
				push( @CoveragePerSample, $meanCov );
			}
		}
		else {
			next;
		}
	}

	#Get file that will be used as standard normal
	my $maxCoverage = max @CoveragePerSample;

	#print "MAX:$maxCoverage\n";
	while ( ( my $key, my $value ) = each(%coverageForNormal) ) {
		if ( $value == $maxCoverage ) {
			$standardNormal = $key;
		}
		else {
			next;
		}
	}
	if ( !$standardNormal ) {
		$standardNormal = $stdNormal;
	}

	#print "SN:$standardNormal:$maxCoverage\n";
	#print "PN:$poolNormal:$poolNormalMeanCov\n";
	#Running Mutect and Somatic Indel Caller
	my $count = 0;
	while ( ( my $key, my $value ) = each(%groupedFilenames) ) {
		my @files = split( ",", $value );

		# Section of Normal
		my @normalSamples = ();
		tie( my %coverageForSampleNormals, 'Tie::IxHash' );
		my @CoverageForMultipleNormal = ();
		my $normal;

		#	my $poolNormal;
		foreach my $file (@files) {
			my ($fileBarcode) = $file =~ /.*_(bc\d+)_.*/;
			my $fileClass = $classPerBarcode{$fileBarcode};
			if ( $fileClass =~ m/Normal/i ) {
				push( @normalSamples, $file );
			}
			else {
				next;
			}
		}
		if ( scalar @normalSamples != 0 ) {
			foreach my $file (@normalSamples) {
				my $HSmetricsFile = $file;
				$HSmetricsFile =~ s/\.bam/\.HSmetrics\.txt/g;

				#print "HS:$HSmetricsFile\n";
				open( FH, "$outdir/$HSmetricsFile" )
				  or die(
					$logger->fatal(
"CallSomaticSNP:Cannot Open $outdir/$HSmetricsFile, Error:$!"
					)
				  );
				my $meanCov;
				my $CovForFile;
				while (<FH>) {
					next until ( $_ =~ /^BAIT_SET/ );
					while (<FH>) {
						next if ( ( $_ =~ /^BAIT_SET/ ) or ( $_ =~ /^\s$/ ) );
						my (@values) = split( "\t", $_ );

						#print "MeanCOv:$values[21]\n";
						$CovForFile = $values[21];
						$meanCov    = $values[21];
					}
					$coverageForSampleNormals{$file} = $CovForFile;
					push( @CoverageForMultipleNormal, $meanCov );
				}
				close(FH);
			}

			#Get file that will be used as normal
			my $maxCoverage = max @CoverageForMultipleNormal;

			#print "MAX:$maxCoverage\n";
			if ( scalar @normalSamples > 1 ) {
				while ( ( my $key, my $value ) =
					each(%coverageForSampleNormals) )
				{
					if ( ( $value == $maxCoverage ) and ( $value >= 50 ) ) {
						$normal = $key;
					}
					else {
						if ( $poolNormalMeanCov <= 50 ) {
							$normal = $standardNormal;
						}
						else {
							$normal = $poolNormal;
						}
					}
				}
			}
			else {
				if ( scalar @normalSamples == 1 ) {
					my $coverage =
					  $coverageForSampleNormals{ $normalSamples[0] };
					if ( $coverage >= 50 ) {
						$normal = $normalSamples[0];
					}
					else {
						if ( $poolNormalMeanCov <= 50 ) {
							$normal = $standardNormal;
						}
						else {
							$normal = $poolNormal;
						}
					}
				}
				else {
					$normal = $poolNormal;
				}
			}
		}
		else {
			if ( $poolNormalMeanCov <= 50 ) {
				$normal = $standardNormal;
			}
			else {
				$normal = $poolNormal;
			}
		}

		#Check if the normal file is with full path
		if ( $normal =~ /\// ) {
			$normal = pop @{ [ split( "/", $normal ) ] };
		}
		else {
			$normal = $normal;
		}

		#RUN Mutation Calling Jobs
		my $pcount = 0;
		foreach my $file (@files) {
			my ($fileBarcode) = $file =~ /.*_(bc\d+)_.*/;
			my $fileClass = $classPerBarcode{$fileBarcode};
			next if ( $fileClass =~ m/Normal/i );
			print "Final2:$file\t$normal\n";
			my ($tFileId) = $file   =~ /(.*)_bc\d+_/;
			my ($nFileId) = $normal =~ /(.*)_bc\d+_/;
			$NormalPerFile{$tFileId} = $nFileId;
			##PindelCofig
			my $pconfig = "$outdir/$tFileId" . "_pindel_config.cfg";
			open( PFH, ">", "$pconfig" )
			  or die(
				$logger->fatal(
					"CallSomaticSNP:Cannot Open $pconfig, Error:$!"
				)
			  );
			print PFH
			  "$outdir/$file\t150\t$tFileId\n$outdir/$normal\t150\t$nFileId\n";
			close(PFH);
			my (
				$somaticMutationFile, $somaticIndelFile,
				$somaticPindelFile,   $waitFileNames
			  )
			  = &RunMutect_SomaticIndelDetector( $key, $normal, $file, $pconfig,
				$outdir, $count );
			foreach my $waitName (@$waitFileNames) {
				push( @notifyNames, $waitName );
			}

			if ( $fileClass =~ m/NTC/i ) {
				push( @NTCsomaticMutationFiles,
					"$somaticMutationFile" . ":" . "$count" );
				push( @NTCsomaticIndelFiles,
					"$somaticIndelFile" . ":" . "$count" );
			}
			else {
				push( @somaticMutationFiles, $somaticMutationFile );
				push( @somaticIndelFiles,    $somaticIndelFile );
				push( @somaticPindelFiles,   $somaticPindelFile );
			}
			$count++;
			$pcount++;
			if ( $pcount >= 8 ) {
				&WaitToFinish( $outdir, @notifyNames );
				@notifyNames = ();
				$pcount      = 0;
			}
		}
	}
	&WaitToFinish( $outdir, @notifyNames );
	my $prog = "RawMutationFiles";

	#Check UG vcf files
	#&CheckOutputFiles($outdir,$prog, @germLineVcfs);
	#Check Mutation out files
	&CheckOutputFiles( $outdir, $prog, @somaticMutationFiles );

	#Check Indel out files
	#&CheckOutputFiles($outdir,$prog, @somaticIndelFiles);
	#push the NTC files back to the place they came from.
	foreach (@NTCsomaticMutationFiles) {
		my ( $name, $num ) = split( ":", $_ );
		splice @somaticMutationFiles, $num, 0, $name;
	}
	foreach (@NTCsomaticIndelFiles) {
		my ( $name, $num ) = split( ":", $_ );
		splice @somaticIndelFiles, $num, 0, $name;
	}
	foreach (@NTCsomaticPindelFiles) {
		my ( $name, $num ) = split( ":", $_ );
		splice @somaticPindelFiles, $num, 0, $name;
	}

	#Output what normal file is used
	open( NFH, ">", "$outdir/$NormalUsed" )
	  or die(
		$logger->fatal("NormalUsed:Cannot open $outdir/$NormalUsed. Error:$!")
	  );
	while ( my ( $key, $value ) = each(%NormalPerFile) ) {
		print NFH "$key\t$value\n";
	}
	close(NFH);
	$now = time - $now;
	$logger->info("Finished running Germline and Somatic Variant jobs on SGE");
	printf(
"Total Germline and Somatic Variant Calling run time: %02d:%02d:%02d\n\n",
		int( $now / 3600 ),
		int( ( $now % 3600 ) / 60 ),
		int( $now % 60 )
	);
	return ( \@names );
}
#####################################
#####################################
#This will help to Filter SNPS:
#Somatic SNPs: Mutect Filter
#Somatics Indels: SomaticIndelCaller Filter
sub FilterSNPsAndIndels {
	my ($filenames) = @_;
	my @names = ();
	if ($filenames) { (@names) = @$filenames; }
	if ( ( scalar(@names) == 0 ) and ($fof) ) {
		@names = &GetNames( $fof, $outdir );
	}
	my (@sortedparseFilenames) =
	  sort { lowestNumber($a) <=> lowestNumber($b) } @names;
	@names = @sortedparseFilenames;
	my @notifyNames                = ();
	my @bedFileNames               = ();
	my @somaticMutationVcfFiles    = ();
	my @somaticIndelVcfFiles       = ();
	my @somaticPindelVcfFiles      = ();
	my @somaticMutationTxtFiles    = ();
	my @somaticIndelTxtFiles       = ();
	my @somaticPindelTxtFiles      = ();
	my @NTCbedFileNames            = ();
	my @NTCsomaticMutationVcfFiles = ();
	my @NTCsomaticIndelVcfFiles    = ();
	my @NTCsomaticPindelVcfFiles   = ();
	my @NTCsomaticMutationTxtFiles = ();
	my @NTCsomaticIndelTxtFiles    = ();
	my @NTCsomaticPindelTxtFiles   = ();
	tie( my %NormalPerFile, 'Tie::IxHash' );
	my $now = time;
	$logger->info("Started running Filter SNPs and Indels jobs on SGE");
	my $count = 0;

	for ( my $i = 0 ; $i < scalar(@names) ; $i++ ) {
		my ($fileBarcode) = $names[$i] =~ /.*_(bc\d+)_.*/;
		my $fileClass = $classPerBarcode{$fileBarcode};
		next if ( $fileClass =~ m/Normal/i );
		my (
			$somaticMutationVcfFile, $somaticMutationTxtFile,
			$somaticIndelVcfFile,    $somaticIndelTxtFile,
			$somaticPindelVcfFile,   $somaticPindelTxtFile,
			$waitFileNames
		) = &RunSomaticMutIndelFilter( $names[$i], $outdir, $i );
		foreach my $waitName (@$waitFileNames) {
			push( @notifyNames, $waitName );
		}
		my ($filename) = $names[$i] =~ /(.*)\.bam/;
		$filename = $filename . ".bed";
		if ( $fileClass =~ m/NTC/i ) {
			push( @NTCbedFileNames, "$filename" . ":" . "$count" );
			push( @NTCsomaticMutationVcfFiles,
				"$somaticMutationVcfFile" . ":" . "$count" );
			push( @NTCsomaticIndelVcfFiles,
				"$somaticIndelVcfFile" . ":" . "$count" );
			push( @NTCsomaticPindelVcfFiles,
				"$somaticPindelVcfFile" . ":" . "$count" );
			push( @NTCsomaticMutationTxtFiles,
				"$somaticMutationTxtFile" . ":" . "$count" );
			push( @NTCsomaticIndelTxtFiles,
				"$somaticIndelTxtFile" . ":" . "$count" );
			push( @NTCsomaticPindelTxtFiles,
				"$somaticPindelTxtFile" . ":" . "$count" );
		}
		else {
			push( @bedFileNames,            $filename );
			push( @somaticMutationVcfFiles, $somaticMutationVcfFile );
			push( @somaticIndelVcfFiles,    $somaticIndelVcfFile );
			push( @somaticPindelVcfFiles,   $somaticPindelVcfFile );
			push( @somaticMutationTxtFiles, $somaticMutationTxtFile );
			push( @somaticIndelTxtFiles,    $somaticIndelTxtFile );
			push( @somaticPindelTxtFiles,   $somaticPindelTxtFile );
		}
		$count++;
	}
	&WaitToFinish( $outdir, @notifyNames );
	$logger->info("Completed Filtering Mutations & Indels");
	my $prog = "STDfilter";

	#Check Mutation out files
	&CheckOutputFiles( $outdir, $prog, @somaticMutationVcfFiles );
	&CheckOutputFiles( $outdir, $prog, @somaticMutationTxtFiles );

	#Check Indel out files
	&CheckOutputFiles( $outdir, $prog, @somaticIndelVcfFiles );
	&CheckOutputFiles( $outdir, $prog, @somaticIndelTxtFiles );

	#Check Indel out files
	&CheckOutputFiles( $outdir, $prog, @somaticPindelVcfFiles );
	&CheckOutputFiles( $outdir, $prog, @somaticPindelTxtFiles );

	#push the NTC files back to the place they came from.
	foreach (@NTCbedFileNames) {
		my ( $name, $num ) = split( ":", $_ );
		splice @bedFileNames, $num, 0, $name;
	}
	foreach (@NTCsomaticMutationVcfFiles) {
		my ( $name, $num ) = split( ":", $_ );
		splice @somaticMutationVcfFiles, $num, 0, $name;
	}
	foreach (@NTCsomaticIndelVcfFiles) {
		my ( $name, $num ) = split( ":", $_ );
		splice @somaticIndelVcfFiles, $num, 0, $name;
	}
	foreach (@NTCsomaticPindelVcfFiles) {
		my ( $name, $num ) = split( ":", $_ );
		splice @somaticPindelVcfFiles, $num, 0, $name;
	}
	foreach (@NTCsomaticMutationTxtFiles) {
		my ( $name, $num ) = split( ":", $_ );
		splice @somaticMutationTxtFiles, $num, 0, $name;
	}
	foreach (@NTCsomaticIndelTxtFiles) {
		my ( $name, $num ) = split( ":", $_ );
		splice @somaticIndelTxtFiles, $num, 0, $name;
	}
	foreach (@NTCsomaticPindelTxtFiles) {
		my ( $name, $num ) = split( ":", $_ );
		splice @somaticPindelTxtFiles, $num, 0, $name;
	}

	#Concatenate all Coding Entries
	my ($poolName) = $names[0] =~ /.*_bc\d+_(.*)_L\d{1,3}.*/;
	my $AllStdFilterEntriesTxt = $poolName . "_AllSomaticMutIndel.txt";
	my $AllStdFilterEntriesVcf = $poolName . "_AllSomaticMutIndel.vcf";
	my $NormalUsed             = $poolName . "_NormalUsedInMutationCalling.txt";

	#store normal used file
	open( NFH, "<", "$outdir/$NormalUsed" )
	  or die(
		$logger->fatal("FilterSNPs:Cannot open $outdir/$NormalUsed. Error:$!")
	  );
	while (<NFH>) {
		chomp($_);
		my @data = split( "\t", $_ );
		$NormalPerFile{ $data[0] } = $data[1];
	}
	close(NFH);
	open( TFH, ">", "$outdir/$AllStdFilterEntriesTxt" )
	  or die(
		$logger->fatal(
			"FilterSNPs:Cannot open $outdir/$AllStdFilterEntriesTxt, Error:$!")
	  );
	open( VFH, ">", "$outdir/$AllStdFilterEntriesVcf" )
	  or die(
		$logger->fatal(
			"FilterSNPs:Cannot open $outdir/$AllStdFilterEntriesTxt, Error:$!")
	  );

	#Get Header Lines of both files
	my $headerMutationVcf = $somaticMutationVcfFiles[0];
	my $headerIndelVcf    = $somaticIndelVcfFiles[0];
	my $headerPindelVcf   = $somaticPindelVcfFiles[0];
	my @headerM           = ();
	my @headerI           = ();
	my @filter            = ();
	my @format            = ();
	my @info              = ();
	my @contigs           = ();
	my @headerLine        = ();
	my $vcfVersion;
	my $customInfo =
"##INFO=<ID=Sample,Number=1,Type=String,Description=\"Name of the Tumor Sample\">";
	my $customInfoMethod =
"##INFO=<ID=CallMethod,Number=1,Type=String,Description=\"Name of the tool used to make the call\">";
	my %SIDrecords = ();
	open( HMV, "$outdir/$headerMutationVcf" )
	  or die(
		$logger->fatal(
			"FilterSNPs:Cannot Open $outdir/$headerMutationVcf, Error:$!")
	  );
	open( HIV, "$outdir/$headerIndelVcf" )
	  or die(
		$logger->fatal(
			"FilterSNPs:Cannot Open $outdir/$headerIndelVcf, Error:$!")
	  );
	open( HPV, "$outdir/$headerPindelVcf" )
	  or die(
		$logger->fatal(
			"FilterSNPs:Cannot Open $outdir/$headerPindelVcf, Error:$!")
	  );
	print "before here\n";
	$vcfVersion = <HMV>;
	chomp($vcfVersion);

	while (<HMV>) {
		if ( $_ =~ /^#/ ) {
			chomp($_);
			if ( $_ =~ /^##FILTER/ ) {
				push( @filter, $_ );
			}
			if ( $_ =~ /^##FORMAT/ ) {
				push( @format, $_ );
			}
			if ( $_ =~ /^##INFO/ ) {
				push( @info, $_ );
			}
			if ( $_ =~ /^##contig/ ) {
				push( @contigs, $_ );
			}
			if ( $_ =~ /^#CHROM/ ) {
				my @splitHeaderLine = split( "\t", $_ );
				my $elementsNum = scalar(@splitHeaderLine);
				splice @splitHeaderLine, $elementsNum - 2, 2;
				push( @splitHeaderLine, $poolName );
				my $joinLine = join( "\t", @splitHeaderLine );
				push( @headerLine, $joinLine );
			}
		}
		else {
			next;
		}
	}
	close(HMV);
	while (<HIV>) {
		if ( $_ =~ /^#/ ) {
			chomp($_);
			if ( ( $_ =~ /ID=AD/ ) or ( $_ =~ /ID=DP/ ) ) {
				next;
			}
			if ( $_ =~ /^##FILTER/ ) {
				push( @filter, $_ );
			}
			if ( $_ =~ /^##FORMAT/ ) {
				push( @format, $_ );
			}
			if ( $_ =~ /^##INFO/ ) {
				push( @info, $_ );
				push( @info, $customInfo );
				push( @info, $customInfoMethod );
			}
		}
		else {
			next;
		}
	}
	close(HIV);
	while (<HPV>) {
		if ( $_ =~ /^#/ ) {
			chomp($_);
			if ( ( $_ =~ /ID=AD/ ) or ( $_ =~ /ID=DP/ ) ) {
				next;
			}
			if ( $_ =~ /^##FILTER/ ) {
				push( @filter, $_ );
			}
			if ( $_ =~ /^##FORMAT/ ) {
				push( @format, $_ );
			}
			if ( $_ =~ /^##INFO/ ) {
				push( @info, $_ );
				push( @info, $customInfo );
				push( @info, $customInfoMethod );
			}
		}
		else {
			next;
		}
	}
	close(HPV);

	#VCF header
	my (@VCFheader) = ();
	push( @VCFheader, $vcfVersion );
	@VCFheader = ( @VCFheader, @filter, @format, @info, @contigs, @headerLine );
	@VCFheader = uniq(@VCFheader);

	#Txt Header
	print TFH
	  "Sample\tNormalUsed\tChrom\tStart\tRef\tAlt\tFailureReason\tCallMethod\n";

	#
	#populate vcf file
	foreach my $line (@VCFheader) {
		print VFH "$line\n";
	}
	my $sampleId = "";
	for (
		my $fileNum = 0 ;
		$fileNum < scalar(@somaticMutationVcfFiles) ;
		$fileNum++
	  )
	{
		if ( -s "$outdir/$somaticMutationVcfFiles[$fileNum]" ) {

			#Mutation VCF
			my $Mvcf =
			  Vcf->new( file => "$outdir/$somaticMutationVcfFiles[$fileNum]" );
			($sampleId) = $somaticMutationVcfFiles[$fileNum] =~ /(.*)_bc\d+/;
			$Mvcf->parse_header();

			#Get Sample ID
			my $Mname1 = $Mvcf->get_column_name(9);
			if ( $Mname1 =~ /\// ) {
				$Mname1 = basename($Mname1);
			}
			my $Mname2 = $Mvcf->get_column_name(10);
			if ( $Mname2 =~ /\// ) {
				$Mname2 = basename($Mname2);
			}

			#Match Sample Id to column
			my ($mCol);
			if ( $sampleId eq $Mname1 ) {
				$mCol = 9;
			}
			if ( $sampleId eq $Mname2 ) {
				$mCol = 10;
			}
			while ( my $rec = $Mvcf->next_data_array() ) {
				$$rec[7] =
				  "Sample=" . $sampleId . ";CallMethod=MuTect;" . $$rec[7];
				print VFH
"$$rec[0]\t$$rec[1]\t$$rec[2]\t$$rec[3]\t$$rec[4]\t$$rec[5]\t$$rec[6]\t$$rec[7]\t$$rec[8]\t$$rec[$mCol]\n";
			}
			$Mvcf->close();
		}
		else {
			$logger->fatal(
"FILTER SNPS: $sampleId: MuTect vcf was not created. You wont have any SNPs Results"
			);
		}

		#Indels VCF
		#Store record for pindel
		if ( -s "$outdir/$somaticIndelVcfFiles[$fileNum]" ) {
			%SIDrecords = ();
			my $Ivcf =
			  Vcf->new( file => "$outdir/$somaticIndelVcfFiles[$fileNum]" );
			($sampleId) = $somaticIndelVcfFiles[$fileNum] =~ /(.*)_bc\d+/;
			$Ivcf->parse_header();

			#Get Sample ID
			my $Iname1 = $Ivcf->get_column_name(9);
			if ( $Iname1 =~ /\// ) {
				$Iname1 = basename($Iname1);
			}
			my $Iname2 = $Ivcf->get_column_name(10);
			if ( $Iname2 =~ /\// ) {
				$Iname2 = basename($Iname2);
			}

			#Match Sample Id to Column
			my ($iCol);
			if ( $sampleId eq $Iname1 ) {
				$iCol = 9;
			}
			if ( $sampleId eq $Iname2 ) {
				$iCol = 10;
			}
			while ( my $rec = $Ivcf->next_data_array() ) {
				$SIDrecords{"$$rec[0]\t$$rec[1]\t$$rec[2]\t$$rec[3]\t$$rec[4]"}
				  = 0;
				$$rec[7] =
				    "Sample="
				  . $sampleId
				  . ";CallMethod=SomaticIndelDetector;"
				  . $$rec[7];
				print VFH
"$$rec[0]\t$$rec[1]\t$$rec[2]\t$$rec[3]\t$$rec[4]\t$$rec[5]\t$$rec[6]\t$$rec[7]\t$$rec[8]\t$$rec[$iCol]\n";
			}
			$Ivcf->close();

			#Pindel VCF
			my $Pvcf =
			  Vcf->new( file => "$outdir/$somaticPindelVcfFiles[$fileNum]" );
			($sampleId) = $somaticIndelVcfFiles[$fileNum] =~ /(.*)_bc\d+/;
			$Pvcf->parse_header();

			#Get Sample ID
			my $Pname1 = $Pvcf->get_column_name(9);
			if ( $Pname1 =~ /\// ) {
				$Pname1 = basename($Pname1);
			}
			my $Pname2 = $Pvcf->get_column_name(10);
			if ( $Pname2 =~ /\// ) {
				$Pname2 = basename($Pname2);
			}

			#Match Sample Id to Column
			my ($pCol);
			if ( $sampleId eq $Pname1 ) {
				$pCol = 9;
			}
			if ( $sampleId eq $Pname2 ) {
				$pCol = 10;
			}
			while ( my $rec = $Pvcf->next_data_array() ) {
				my $pindelKey =
				  "$$rec[0]\t$$rec[1]\t$$rec[2]\t$$rec[3]\t$$rec[4]";
				$$rec[7] =
				  "Sample=" . $sampleId . ";CallMethod=Pindel;" . $$rec[7];
				if ( exists $SIDrecords{$pindelKey} ) {
					$logger->info(
						"$sampleId: $pindelKey was also called by Pindel");
				}
				else {
					print VFH
"$$rec[0]\t$$rec[1]\t$$rec[2]\t$$rec[3]\t$$rec[4]\t$$rec[5]\t$$rec[6]\t$$rec[7]\t$$rec[8]\t$$rec[$iCol]\n";
				}
			}
			$Pvcf->close();
		}
		else {
			$logger->fatal(
"FILTER INDEL: $sampleId: SID vcf was not created. You wont have any Indel Results"
			);
		}

		#Populate Txt File
		open( SMTFH, "$outdir/$somaticMutationTxtFiles[$fileNum]" )
		  or die(
			$logger->fatal(
"FilterSNPs:Cannot open $outdir/$somaticMutationTxtFiles[$fileNum], Error:$!"
			)
		  );
		open( SITFH, "$outdir/$somaticIndelTxtFiles[$fileNum]" )
		  or die(
			$logger->fatal(
"FilterSNPs:Cannot open $outdir/$somaticIndelTxtFiles[$fileNum], Error:$!"
			)
		  );
		open( PITFH, "$outdir/$somaticPindelTxtFiles[$fileNum]" )
		  or die(
			$logger->fatal(
"FilterSNPs:Cannot open $outdir/$somaticPindelTxtFiles[$fileNum], Error:$!"
			)
		  );
		open( BFH, ">", "$outdir/$bedFileNames[$fileNum]" )
		  or die(
			$logger->fatal(
"FilterSNPs:Cannot open $outdir/$bedFileNames[$fileNum], Error:$!"
			)
		  );
		my (@dataCols)    = ();
		my (@newDataCols) = ();
		my $normalUsed;
		my $mafStart;
		my $mafEnd;

		while (<SMTFH>) {
			chomp($_);
			@newDataCols = ();
			(@newDataCols) = split( "\t", $_ );
			$normalUsed = "";
			if ( exists $NormalPerFile{ $newDataCols[0] } ) {
				$normalUsed = $NormalPerFile{ $newDataCols[0] };
			}
			else {
				$normalUsed = "NULL";
			}
			$mafStart = $newDataCols[2];
			$mafEnd   = $newDataCols[2];
			print TFH
"$newDataCols[0]\t$normalUsed\t$newDataCols[1]\t$newDataCols[2]\t$newDataCols[3]\t$newDataCols[4]\t$newDataCols[5]\tMuTect\n";
			print BFH
"$newDataCols[1]\t$mafStart\t$mafEnd\t$newDataCols[0]:$normalUsed\n";
		}
		close(SMTFH);
		%SIDrecords = ();
		while (<SITFH>) {
			chomp($_);
			@dataCols = ();
			(@dataCols) = split( "\t", $_ );
			@newDataCols = ();
			@newDataCols = grep( s/\s*$//g, @dataCols );
			$normalUsed  = "";
			if ( exists $NormalPerFile{ $newDataCols[0] } ) {
				$normalUsed = $NormalPerFile{ $newDataCols[0] };
			}
			else {
				$normalUsed = "NULL";
			}
			if ( ( length( $newDataCols[3] ) ) > ( length( $newDataCols[4] ) ) )
			{
				my $deletionLength = ( length( $newDataCols[3] ) ) - 1;
				$mafStart = $newDataCols[2];
				$mafEnd   = $newDataCols[2] + $deletionLength;
			}
			else {
				$mafStart = $newDataCols[2];
				$mafEnd   = $newDataCols[2];
			}
			$SIDrecords{
"$newDataCols[0]\t$normalUsed\t$newDataCols[1]\t$newDataCols[2]\t$newDataCols[3]\t$newDataCols[4]"
			  } = 0;

#print TFH "$newDataCols[0]\t$normalUsed\t$newDataCols[1]\t$newDataCols[2]\t$newDataCols[3]\t$newDataCols[4]\tSID\n";
			print BFH
"$newDataCols[1]\t$mafStart\t$mafEnd\t$newDataCols[0]:$normalUsed\n";
		}
		close(SITFH);
		while (<PITFH>) {
			chomp($_);
			@dataCols = ();
			(@dataCols) = split( "\t", $_ );
			@newDataCols = ();
			@newDataCols = grep( s/\s*$//g, @dataCols );
			$normalUsed  = "";
			if ( exists $NormalPerFile{ $newDataCols[0] } ) {
				$normalUsed = $NormalPerFile{ $newDataCols[0] };
			}
			else {
				$normalUsed = "NULL";
			}
			if ( ( length( $newDataCols[3] ) ) > ( length( $newDataCols[4] ) ) )
			{
				my $deletionLength = ( length( $newDataCols[3] ) ) - 1;
				$mafStart = $newDataCols[2];
				$mafEnd   = $newDataCols[2] + $deletionLength;
			}
			else {
				$mafStart = $newDataCols[2];
				$mafEnd   = $newDataCols[2];
			}
			my $pindelKey =
"$newDataCols[0]\t$normalUsed\t$newDataCols[1]\t$newDataCols[2]\t$newDataCols[3]\t$newDataCols[4]";
			if ( exists $SIDrecords{$pindelKey} ) {
				$logger->info(
					"$newDataCols[0]:$pindelKey was called by Pindel also");
				print TFH
"$newDataCols[0]\t$normalUsed\t$newDataCols[1]\t$newDataCols[2]\t$newDataCols[3]\t$newDataCols[4]\t.\tSID,PINDEL\n";
				delete $SIDrecords{$pindelKey};
			}
			else {
				print TFH
"$newDataCols[0]\t$normalUsed\t$newDataCols[1]\t$newDataCols[2]\t$newDataCols[3]\t$newDataCols[4]\t.\tPINDEL\n";
				print BFH
"$newDataCols[1]\t$mafStart\t$mafEnd\t$newDataCols[0]:$normalUsed\n";
			}
		}
		close(PITFH);
		close(BFH);
		while ( my ( $sidKey, $sidVal ) = each(%SIDrecords) ) {
			print TFH "$sidKey\t.\tSID\n";
		}
	}
	close(TFH);
	close(VFH);

	#Genotype these variants.
	my ($mergeOutFile) =
	  &GenotypeAllele_And_MergeFiles( \@names, $poolName,
		$AllStdFilterEntriesVcf, $AllStdFilterEntriesTxt, $NormalUsed );
	$now = time - $now;
	$logger->info("Finished running Somatic Variant Filter jobs on SGE");
	printf(
		"Total Somatic Variant Filter run time: %02d:%02d:%02d\n\n",
		int( $now / 3600 ),
		int( ( $now % 3600 ) / 60 ),
		int( $now % 60 )
	);
	return ( \@names );
}
#####################################
#####################################
#Genotype varinats from each filtered vcf file
sub GenotypeAllele_And_MergeFiles {
	my ( $bamFiles, $poolName, $AllStdFilterEntriesVcf, $AllStdFilterEntriesTxt,
		$NormalUsed )
	  = @_;
	my @bamFilesToGenotypeAgainst = @$bamFiles;
	my @notifyNames               = ();
	my @GenotypedFiles            = ();
	my @stdNormals                = ();
	my @stdNormalIds              = ();
	my @typeOfSample              = ();
	foreach my $file (@bamFilesToGenotypeAgainst) {
		my ($fileBarcode) = $file =~ /.*_(bc\d+)_.*/;
		my $fileClass = $classPerBarcode{$fileBarcode};
		push( @typeOfSample, $fileClass );
	}

	#get standard Normal List
	$logger->info(
"GenotypeMutations:Globbing standard normals location ($StandardNormalsDirectory)"
	);
	my @stdNormalBams = glob("$StandardNormalsDirectory/*bam");

	#Added if mutation calling normal exists
	if ( $stdNormal ne "NULL" ) {
		push( @stdNormalBams, $stdNormal );
	}
	foreach (@stdNormalBams) {
		chomp($_);
		my $bamfilename = basename($_);
		my ($stdBamId) = $bamfilename =~ /(.*)_bc\d+/;
		push( @bamFilesToGenotypeAgainst, $bamfilename );
		push( @typeOfSample,              "Normal" );
		push( @stdNormals,                $bamfilename );
		push( @stdNormalIds,              $stdBamId );
		if ( -e "$outdir/$bamfilename" ) {
			next;
		}
		else {
			chop($_);
			eval {
				if ( $CLUSTER eq "SGE" )
				{
					my $cmd =
"$QSUB -q $queue -V -N RSYNC.$stdBamId.$$ -wd $outdir -e RSYNC.$stdBamId.$$.stderr -o RSYNC.$stdBamId.$$.stdout -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y $RSYNC -a -L $_\* $outdir/.";
					$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -N RSYNC.$stdBamId.$$ -wd $outdir -e RSYNC.$stdBamId.$$.stderr -o RSYNC.$stdBamId.$$.stdout -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y $RSYNC -a -L $_\* $outdir/.`;

#`$QSUB -q $queue -V -N RSYNC.$stdBamId.$$ -wd $outdir -e RSYNC.$stdBamId.$$.stderr -o RSYNC.$stdBamId.$$.stdout -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y echo "Run Rsync"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid RSYNC.$stdBamId.$$ -N NotityRsync.$stdBamId.$$ -e NotifyRsync.$stdBamId.$$.stderr -o NotifyRsync.$stdBamId.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
				}
				else {
					my $cmd =
"$BSUB -q $queue -J RSYNC.$stdBamId.$$ -cwd $outdir -e RSYNC.$stdBamId.$$.%J.stderr -o RSYNC.$stdBamId.$$.%J.stdout -R \"rusage[mem=2]\" -M 4 -n 1 \"$RSYNC -a -L $_\* $outdir/.\"";
					$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -J RSYNC.$stdBamId.$$ -cwd $outdir -e RSYNC.$stdBamId.$$.%J.stderr -o RSYNC.$stdBamId.$$.%J.stdout -We 0:59 -R "rusage[mem=2]" -M 4 -n 1 "$RSYNC -a -L $_\* $outdir/."`;
`$BSUB -q $queue -cwd $outdir -w "post_done(RSYNC.$stdBamId.$$)" -J NotityRsync.$stdBamId.$$ -e NotifyRsync.$stdBamId.$$.%J.stderr -o NotifyRsync.$stdBamId.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
				}
			};
			if ($@) {
				$logger->fatal(
"GenotypeMutations:Cannot rsync standard normal list. Error:$@"
				);
				exit(1);
			}
			push( @notifyNames, "NotifyRsync.$stdBamId.$$.stat" );
		}
	}

	#close(SNFH);
	&WaitToFinish( $outdir, @notifyNames );
	@notifyNames = ();

	#Traverse through all Mutation Files
	$logger->info("Genotyping Found Mutations");
	my $count    = 0;
	my $totcount = 0;
	my $index    = 0;
	my $total    = scalar(@bamFilesToGenotypeAgainst);
	foreach my $bamFile (@bamFilesToGenotypeAgainst) {

		#print "$bamFile\n";
		my ($bamID) = $bamFile =~ /(.*)_bc\d+/;
		my $mpileUpMutationOut = $bamFile;
		$mpileUpMutationOut =~ s/\.bam/.mpileup/;
		my $alleleDepthMutationOut = $bamFile;
		$alleleDepthMutationOut =~ s/\.bam/_mpileup\.alleledepth/;

		#Mutation File
		if (    ( -e "$outdir/$alleleDepthMutationOut" )
			and ( ( -s "$outdir/$alleleDepthMutationOut" ) != 0 ) )
		{
			push( @GenotypedFiles, $alleleDepthMutationOut );
			$totcount++;
			next;
		}
		else {
			eval {
				if ( $CLUSTER eq "SGE" )
				{
					my $cmd =
"$QSUB -q $queue -V -N AD.$bamID.$$ -wd $outdir -e AD.$bamID.$$.stderr -o AD.$bamID.$$.stdout -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y $PERL $GenotypeAllele -fmv $outdir/$AllStdFilterEntriesVcf -bam $outdir/$bamFile -rf $Reference -s $SAMTOOLS -o $outdir -of $alleleDepthMutationOut -mof $mpileUpMutationOut -bi $bamID -q $queue -tos $typeOfSample[$index]";
					$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -N AD.$bamID.$$ -wd $outdir -e AD.$bamID.$$.stderr -o AD.$bamID.$$.stdout -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y $PERL $GenotypeAllele -fmv $outdir/$AllStdFilterEntriesVcf -bam $outdir/$bamFile -rf $Reference -s $SAMTOOLS -o $outdir -of $alleleDepthMutationOut -mof $mpileUpMutationOut -bi $bamID -mmq $MAPQ -mbq $BASQ -q $queue -b $BEDTOOLS -qsub $QSUB -tos $typeOfSample[$index]`;
`$QSUB -q $queue -V -wd $outdir -hold_jid AD.$bamID.$$ -N NotifyAD.$bamID.$$ -e NotifyAD.$bamID.$$.stderr -o NotifyAD.$bamID.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
				}
				else {
					my $cmd =
"$BSUB -q $queue -J AD.$bamID.$$ -cwd $outdir -e AD.$bamID.$$.%J.stderr -o AD.$bamID.$$.%J.stdout -R \"rusage[mem=2]\" -R \"rusage[iounits=0]\" -M 4 -n 1 \"$PERL $GenotypeAllele -fmv $outdir/$AllStdFilterEntriesVcf -bam $outdir/$bamFile -rf $Reference -s $SAMTOOLS -o $outdir -of $alleleDepthMutationOut -mof $mpileUpMutationOut -bi $bamID -mmq $MAPQ -mbq $BASQ -q $queue -b $BEDTOOLS -bsub $BSUB -tos $typeOfSample[$index]\"";
					$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -J AD.$bamID.$$ -cwd $outdir -e AD.$bamID.$$.%J.stderr -o AD.$bamID.$$.%J.stdout -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$PERL $GenotypeAllele -fmv $outdir/$AllStdFilterEntriesVcf -bam $outdir/$bamFile -rf $Reference -s $SAMTOOLS -o $outdir -of $alleleDepthMutationOut -mof $mpileUpMutationOut -bi $bamID -mmq $MAPQ -mbq $BASQ -q $queue -b $BEDTOOLS -bsub $BSUB -tos $typeOfSample[$index]"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(AD.$bamID.$$)" -J NotifyAD.$bamID.$$ -e NotifyAD.$bamID.$$.%J.stderr -o NotifyAD.$bamID.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
				}
			};
			if ($@) {
				$logger->fatal("Genotype_Merge:Job Sumission Failed. Error:$@");
				exit(1);
			}
			push( @notifyNames,    "NotifyAD.$bamID.$$.stat" );
			push( @GenotypedFiles, $alleleDepthMutationOut );
			$count++;
			$totcount++;
		}
		if ( ( $count >= 10 ) or ( $totcount == $total ) ) {
			&WaitToFinish( $outdir, @notifyNames );
			@notifyNames = ();
			$count       = 0;
		}
		$index++;
	}
	print "Completed Running Allele Depth on the files\n";

	#Check Genotyped out files
	my $prog;
	&CheckOutputFiles( $outdir, $prog, @GenotypedFiles );
	print "Deleting Standard Normal Files\n";
	foreach my $stdFiles (@stdNormals) {
		chop($stdFiles);
		eval { `rm $outdir/$stdFiles\*`; };
		if ($@) {
			$logger->warn(
"GenotypeMutations:Cannot remove standard normal list files. Error:$@"
			);
		}
	}

	#Merging of all files
	$logger->info("Merging Mutations with genotyped files");

	#Hash storing mpileup data for Normal Samples
	tie( my %NormalDataHash, 'Tie::IxHash' );

	#Hash storing mpileup data for Tumor Samples
	tie( my %TumorDataHash, 'Tie::IxHash' );

	#Title file used to decide Normal and tumor samples and makking a tied hash
	my $titleFile = "$outdir/$poolName" . "_title.txt";
	open( TFH, "<", $titleFile )
	  or die(
		$logger->fatal("GenotypeMutations:Cannot open $titleFile, Error:$!") );
	while (<TFH>) {
		next if ( $. == 1 );
		my @ndataCols = split( "\t", $_ );
		my @dataCols = grep( s/\s*$//g, @ndataCols );

		#print "$dataCols[2]:$dataCols[5]\n";
		if ( $dataCols[5] =~ /Normal/i ) {
			tie %{ $NormalDataHash{ $dataCols[2] } }, 'Tie::IxHash';
		}
		else {
			tie %{ $TumorDataHash{ $dataCols[2] } }, 'Tie::IxHash';
		}
	}
	foreach my $stdNormal (@stdNormalIds) {
		tie %{ $NormalDataHash{$stdNormal} }, 'Tie::IxHash';
	}

	#populating Normal and Tumor Hash from the genotyped files
	for ( my $i = 0 ; $i < scalar(@GenotypedFiles) ; $i++ ) {
		open( FH, "$outdir/$GenotypedFiles[$i]" )
		  or die(
			$logger->fatal(
"GenotypeMutations:Cannot open $outdir/$GenotypedFiles[$i], Error:$!"
			)
		  );
		my $header = <FH>;
		while (<FH>) {
			next if ( $_ =~ /^\s*$/ );
			chomp($_);
			my @ndataCols = split( "\t", $_ );
			my @dataCols = grep( s/\s*$//g, @ndataCols );
			my (
				$refBamId, $sampleId, $chr,        $pos,
				$alt,      $ref,      $totalDepth, $totalRefCount,
				$altCount, $altFreq,  $refFwd,     $refRev,
				$altFwd,   $altRev
			);
			if   ( $dataCols[0] ) { $refBamId = $dataCols[0]; }
			else                  { $refBamId = ""; }
			if   ( $dataCols[1] ) { $sampleId = $dataCols[1]; }
			else                  { $sampleId = ""; }
			if   ( $dataCols[2] ) { $chr = $dataCols[2]; }
			else                  { $chr = ""; }
			if   ( $dataCols[3] ) { $pos = $dataCols[3]; }
			else                  { $pos = 0; }
			if   ( $dataCols[4] ) { $ref = $dataCols[4]; }
			else                  { $ref = ""; }
			if   ( $dataCols[5] ) { $alt = $dataCols[5]; }
			else                  { $alt = ""; }
			if   ( $dataCols[6] ) { $totalDepth = $dataCols[6]; }
			else                  { $totalDepth = 0; }
			if   ( $dataCols[7] ) { $totalRefCount = $dataCols[7]; }
			else                  { $totalRefCount = 0; }
			if   ( $dataCols[8] ) { $altCount = $dataCols[8]; }
			else                  { $altCount = 0; }
			if   ( $dataCols[9] ) { $altFreq = $dataCols[9]; }
			else                  { $altFreq = 0.0; }
			if   ( $dataCols[10] ) { $refFwd = $dataCols[10]; }
			else                   { $refFwd = 0; }
			if   ( $dataCols[11] ) { $refRev = $dataCols[11]; }
			else                   { $refRev = 0; }
			if   ( $dataCols[12] ) { $altFwd = $dataCols[12]; }
			else                   { $altFwd = 0; }
			if   ( $dataCols[13] ) { $altRev = $dataCols[13]; }
			else                   { $altRev = 0; }

			if ( exists $NormalDataHash{$refBamId} ) {
				$NormalDataHash{$refBamId}{ $sampleId . ":" 
					  . $chr . ":" 
					  . $pos . ":" 
					  . $ref . ":"
					  . $alt } =
				  "$totalDepth\t$totalRefCount\t$altCount\t$altFreq";
			}
			if ( exists $TumorDataHash{$refBamId} ) {
				$TumorDataHash{$refBamId}{ $sampleId . ":" 
					  . $chr . ":" 
					  . $pos . ":" 
					  . $ref . ":"
					  . $alt } =
"$totalDepth\t$totalRefCount\t$altCount\t$altFreq\t$refFwd\t$refRev\t$altFwd\t$altRev";
			}
		}
		close(FH);
	}

	#Merging Records
	my $mergedOutFile = $AllStdFilterEntriesTxt;
	$mergedOutFile =~ s/\.txt/_withAlleleDepth\.txt/;
	open( OMFMH, ">", "$outdir/$mergedOutFile" )
	  or die(
		$logger->fatal(
			"GenotypeMutations:Cannot open $outdir/$mergedOutFile, Error:$!")
	  );
	open( FMH, "<", "$outdir/$AllStdFilterEntriesTxt" )
	  or die(
		$logger->fatal(
"GenotypeMutations:Cannot open $outdir/$AllStdFilterEntriesTxt, Error:$!"
		)
	  );
	my $Orgheader = <FMH>;
	chomp($Orgheader);
	my ( $sampleId, $NsampleId, $chr, $pos, $ref, $alt, $failureReason,
		$callMethod )
	  = split( "\t", $Orgheader );
	my $OncotatorHeader =
"Hugo_Symbol\tCall_Confidence\tComments\t$failureReason\t$callMethod\tVariant_Class\tdbSNP_RS\tdbSNP_Val\tProtein_Change\tCOSMIC_site";
	my $sampleTumorHeader =
"T_TotalDepth\tT_RefCount\tT_AltCount\tT_AltFreq\tT_Ref+\tT_Ref-\tT_Alt+\tT_Alt-";
	my $sampleNormalHeader = "N_TotalDepth\tN_RefCount\tN_AltCount\tN_AltFreq";
	my @keyNormalHash      = keys(%NormalDataHash);
	my $allNormalHeader    = join( "\t", @keyNormalHash );
	my @keyTumorHash       = keys(%TumorDataHash);
	my $allTumorHeader     = join( "\t", @keyTumorHash );
	my $AllnormalAggregate = "All_N_Aggregate_AlleleDepth";
	my $AllnormalFreq      = "All_N_Median_AlleleFreq";
	my $AlltumornormalFreq = "T_freq/All_N_Freq";
	my $OccurenceInNormals = "Occurence_in_Normals";

#my $MutationAssessorHeader = "MA:FImpact\tMA:link.MSA\tMA:link.PDB\tMA:link.var";
	print OMFMH
"$sampleId\t$NsampleId\t$chr\t$pos\t$ref\t$alt\t$OncotatorHeader\t$sampleNormalHeader\t$sampleTumorHeader\t$AllnormalAggregate\t$AllnormalFreq\t$AlltumornormalFreq\t$OccurenceInNormals\t$allNormalHeader\t$allTumorHeader\n";
	while (<FMH>) {
		chomp($_);
		my ( $sampleId, $NsampleId, $chr, $pos, $ref, $alt, $failureReason,
			$callMethod )
		  = split( "\t", $_ );
		my $orgKey =
		  $sampleId . ":" . $chr . ":" . $pos . ":" . $ref . ":" . $alt;
		my $SampleTumorMpileUpvalues  = $TumorDataHash{$sampleId}{$orgKey};
		my $SampleNormalMpileUpvalues = $NormalDataHash{$NsampleId}{$orgKey};

		#print "Tumor:$sampleId:$orgKey:$SampleTumorMpileUpvalues\n";
		#print "Normal:$NsampleId:$orgKey:$SampleNormalMpileUpvalues\n";
		my @addNcols  = ();
		my @addTcols  = ();
		my @AllRefAD  = ();
		my @AllAltAD  = ();
		my @AllFreqAD = ();

		#Processing for Normal Samples
		while ( my ( $key, $value ) = each(%NormalDataHash) ) {
			if ( exists $NormalDataHash{$key}{$orgKey} ) {
				my $eachNormalSampleMpileup = $NormalDataHash{$key}{$orgKey};

				#print "Each:$eachNormalSampleMpileup\n";
				my ( $totalDepth, $totalRefCount, $altCount, $altFreq ) =
				  split( "\t", $eachNormalSampleMpileup );
				push( @AllRefAD,  $totalDepth );
				push( @AllAltAD,  $altCount );
				push( @AllFreqAD, $altFreq );
				$altFreq = sprintf( "%.3f", $altFreq );
				my $altLine =
				    "DP=$totalDepth"
				  . ";RD=$totalRefCount"
				  . ";AD=$altCount"
				  . ";VF=$altFreq";
				push( @addNcols, $altLine );
			}
		}

		#Adding Total Counts
		my $totnRefAD  = sum(@AllRefAD);
		my $totnAltAD  = sum(@AllAltAD);
		my $totnFreqAD = 0;

		#All normal median freq
		$totnFreqAD = &median( \@AllFreqAD );
		$totnFreqAD = sprintf( "%.3f", $totnFreqAD );

		#Sample Tumor divided my All normal median frequency
		my ( $T_totalDepth, $T_totalRefCount, $T_altCount, $T_altFreq,
			$T_refFwd, $T_refRev, $T_altFwd, $T_altRev )
		  = split( "\t", $SampleTumorMpileUpvalues );
		my $Tfreq_AllNfreq = 0;
		$Tfreq_AllNfreq = sprintf( "%.5f", ( $T_altFreq / $totnFreqAD ) )
		  if ( $totnFreqAD != 0 );
		my $allADtogather = $totnRefAD . "/" . $totnAltAD;

# AZ Aug 29th, changed the dash to pipe. with small numbers, excel changes dash to dates
		my $Tpvf = sprintf( "%.5f", ( $T_altFreq / $PON_TPVF_threshold ) );
		if ( $Tpvf < 0.01 ) {
			$Tpvf = 0.01;
		}
		my $OccurenceCount      = 0;
		my $OccurencePercentage = 0;

		#calculating Occurence
		for ( my $sCount = 0 ; $sCount < scalar(@AllFreqAD) ; $sCount++ ) {
			if (    ( $AllAltAD[$sCount] >= $PON_AD_threshold )
				and ( $AllFreqAD[$sCount] >= $Tpvf ) )
			{

				#print "$AllAltAD[$sCount] $AllFreqAD[$sCount]\n";
				#print "$_";
				$OccurenceCount++;
			}
		}
		if ( ($OccurenceCount) and ( $OccurenceCount != 0 ) ) {
			$OccurencePercentage =
			  sprintf( "%.3f", $OccurenceCount / keys(%NormalDataHash) );
		}
		else {
			$OccurencePercentage = 0.0;
		}
		my $allNormalCols = join( "\t", @addNcols );

		#Processing for Tumor Samples
		while ( my ( $key, $value ) = each(%TumorDataHash) ) {
			if ( exists $TumorDataHash{$key}{$orgKey} ) {
				my $eachTumorSampleMpileup = $TumorDataHash{$key}{$orgKey};

				#print "$key\t$orgKey\t$eachTumorSampleMpileup\n";
				my ( $totalDepth, $totalRefCount, $altCount, $altFreq, $refFwd,
					$refRev, $altFwd, $altRev )
				  = split( "\t", $eachTumorSampleMpileup );
				$altFreq = sprintf( "%.3f", $altFreq );
				my $altLine =
				    "DP=$totalDepth"
				  . ";RD=$totalRefCount"
				  . ";AD=$altCount"
				  . ";VF=$altFreq";
				push( @addTcols, $altLine );
			}
		}
		my $allTumorCols = join( "\t", @addTcols );
		print OMFMH
"$sampleId\t$NsampleId\t$chr\t$pos\t$ref\t$alt\t\t\t\t$failureReason\t$callMethod\t\t\t\t\t\t$SampleNormalMpileUpvalues\t$SampleTumorMpileUpvalues\t$allADtogather\t$totnFreqAD\t$Tfreq_AllNfreq\t$OccurenceCount;$OccurencePercentage\t$allNormalCols\t$allTumorCols\n";
	}
	close(OMFMH);
	close(FMH);
	$logger->info("Competed Merging of Mutation & Genotyped Files");
	if ( $mergeDinucleotide == 1 ) {
		&MergeDinucleotides( $mergedOutFile, $poolName, $outdir );
	}
	else {
		$logger->warn("Not merging the Dinucleotides");
	}
	return ($mergedOutFile);
}
#####################################
#####################################
#calculate median from the array
sub median {
	@_ == 1 or die('Sub usage: $median = median(\@array);');
	my ($array_ref) = @_;
	my $count = scalar @$array_ref;

	# Sort a COPY of the array, leaving the original untouched
	my @array = sort { $a <=> $b } @$array_ref;
	if ( $count % 2 ) {
		return $array[ int( $count / 2 ) ];
	}
	else {
		return ( $array[ $count / 2 ] + $array[ $count / 2 - 1 ] ) / 2;
	}
}
#####################################
#####################################
#Merge Dinucleotide
sub MergeDinucleotides {
	my ( $genotypeFile, $poolName, $outdir ) = @_;
	my $mdnFile = $genotypeFile;
	$mdnFile =~ s/\.txt/_mergedDNP\.txt/;
	my $InputToOncotator = $poolName . "_InputToOncotator.maf";
	my $MergedStats      = $poolName . "_DNPmerged.stats";
	my $UnMergedStats    = $poolName . "_DNPunmerged.stats";

	#Hash storing sample data for Samples
	tie( my %DataHash, 'Tie::IxHash' );

	#Hash storing sample-coordinate information
	tie( my %SCHash, 'Tie::IxHash' );

	#Hash storing sample-coordinate-Allele information
	tie( my %SCAHash, 'Tie::IxHash' );

	#Open Merge Dinucleotide File
	open( OFH, ">", "$outdir/$mdnFile" )
	  or die(
		$logger->fatal(
			"MergeDinucleotide:Cannot open $outdir/$mdnFile, Error:$!")
	  );

	#Open Input to oncotator file
	open( IFH, ">", "$outdir/$InputToOncotator" )
	  or die(
		$logger->fatal(
			"MergeDinucleotide:Cannot open $outdir/$InputToOncotator, Error:$!")
	  );

	#Open Merge Dinucleotide stats File
	open( MFH, ">", "$outdir/$MergedStats" )
	  or die(
		$logger->fatal(
			"MergeDinucleotide:Cannot open $outdir/$MergedStats, Error:$!")
	  );

	#Open UnMerge Dinucleotide stats File
	open( UMFH, ">", "$outdir/$UnMergedStats" )
	  or die(
		$logger->fatal(
			"MergeDinucleotide:Cannot open $outdir/$UnMergedStats, Error:$!")
	  );

	#Open Genotype File
	open( FH, "$outdir/$genotypeFile" )
	  or die(
		$logger->fatal(
			"MergeDinucleotide:Cannot open $outdir/$genotypeFile, Error:$!")
	  );
	my $header = <FH>;
	while (<FH>) {
		chomp($_);
		my (@dataCols) = split( "\t", $_ );
		if ( ( length( $dataCols[4] ) ) == ( length( $dataCols[5] ) ) ) {
			if ( exists $SCHash{ $dataCols[0] } ) {
				my $val = $SCHash{ $dataCols[0] };
				$val = $val . ";" . $dataCols[2] . ":" . $dataCols[3];
				$SCHash{ $dataCols[0] } = $val;
			}
			else {
				$SCHash{ $dataCols[0] } = $dataCols[2] . ":" . $dataCols[3];
			}
		}
		$SCAHash{ $dataCols[0] }{ $dataCols[2] . ":" . $dataCols[3] } =
		    $dataCols[4] . "#"
		  . $dataCols[5] . "#"
		  . $dataCols[24] . "#"
		  . $dataCols[26] . "#"
		  . $dataCols[27];
	}
	close(FH);
	print MFH
"Sample\tChr1:Pos1\tChr2:Pos2\tRef1\tRef2\tAlt1\tAlt2\tDP1\tDP2\tAD1\tAD2\tVF1\tVF2\n";
	while ( my ( $key, $val ) = each(%SCHash) ) {
		my (@coords) = ();
		@coords = split( ";", $val );
		my %searchHash = ();
		foreach my $coord (@coords) { $searchHash{$coord} = () }
		foreach my $coord (@coords) {
			my ( $chr, $adCoord ) = split( ":", $coord );
			$adCoord = $adCoord + 1;
			$adCoord = $chr . ":" . $adCoord;
			if ( exists $searchHash{$adCoord} ) {
				my $dataval = $SCAHash{$key}{$coord};
				my ( $ref1, $alt1, $dp1, $ad1, $vf1 ) = split( "#", $dataval );
				$dataval = $SCAHash{$key}{$adCoord};
				my ( $ref2, $alt2, $dp2, $ad2, $vf2 ) = split( "#", $dataval );
				my ( $ref, $alt );
				if (   ( abs( $dp1 - $dp2 ) <= 10 )
					&& ( abs( $ad1 - $ad2 ) <= 5 ) )
				{
					$ref = $ref1 . $ref2;
					$alt = $alt1 . $alt2;
					$SCAHash{$key}{$coord} =
					  $ref . "#" . $alt . "#" . $dp1 . "#" . $ad1 . "#" . $vf1;
					print MFH
"$key\t$coord\t$adCoord\t$ref1\t$ref2\t$alt1\t$alt2\t$dp1\t$dp2\t$ad1\t$ad2\t$vf1\t$vf2\tMerged\n";
					delete $SCAHash{$key}{$adCoord};
				}
				else {
					print UMFH
"$key\t$coord\t$adCoord\t$ref1\t$ref2\t$alt1\t$alt2\t$dp1\t$dp2\t$ad1\t$ad2\t$vf1\t$vf2\tUnMerged\n";
					next;
				}
			}
		}
	}
	close(MFH);
	close(UMFH);
	print IFH
"Chromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tNCBI_Build\n";
	open( FH, "$outdir/$genotypeFile" )
	  or die(
		$logger->fatal(
			"MergeDinucleotide:Cannot open $outdir/$genotypeFile, Error:$!")
	  );
	$header = <FH>;
	print OFH "$header";
	while (<FH>) {
		chomp($_);
		my (@dataCols) = split( "\t", $_ );
		if (
			exists $SCAHash{ $dataCols[0] }
			{ $dataCols[2] . ":" . $dataCols[3] } )
		{
			my ( $ref, $alt, $dp, $ad, $vf ) =
			  split( "#",
				$SCAHash{ $dataCols[0] }{ $dataCols[2] . ":" . $dataCols[3] } );
			$dataCols[4]  = $ref;
			$dataCols[5]  = $alt;
			$dataCols[24] = $dp;
			$dataCols[26] = $ad;
			$dataCols[27] = $vf;
			my $data = join( "\t", @dataCols );
			my ( $mafStart, $mafEnd );

			if ( ( length($ref) ) > ( length($alt) ) ) {
				my $deletionLength = ( length($ref) ) - 1;
				$mafStart = $dataCols[3];
				$mafEnd   = $dataCols[3] + $deletionLength;
			}
			elsif ( ( length($ref) == 2 ) && ( length($alt) == 2 ) ) {
				$mafStart = $dataCols[3];
				$mafEnd   = $dataCols[3] + 1;
			}
			else {
				$mafStart = $dataCols[3];
				$mafEnd   = $dataCols[3];
			}
			print IFH "$dataCols[2]\t$mafStart\t$mafEnd\t$ref\t$alt\t\t37\n";
			print OFH "$data\n";
		}
	}
	close(IFH);
	close(OFH);
	close(FH);
	return;
}
#####################################
#####################################
#Annotate as well as Assess the
#SNPs & Indels
sub AnnotateSNPsAndIndels {
	my ( $filenames, $outdir, $fof ) = @_;

	#print "$filenames\t$outdir\t$fof\n";
	my @names = ();
	my @notifyNames;
	if ($filenames) { (@names) = @$filenames; }
	if ( ( scalar(@names) == 0 ) and ($fof) ) {
		@names = &GetNames( $fof, $outdir );
	}
	my (@sortedparseFilenames) =
	  sort { lowestNumber($a) <=> lowestNumber($b) } @names;
	@names = @sortedparseFilenames;
	my $somaticMutIndelFile =
	  $poolName . "_AllSomaticMutIndel_withAlleleDepth_mergedDNP.txt";

	#my $oncotatorInput = $poolName . "_InputToOncotator.maf";
	my $exonCoverage          = $poolName . "_ALL_exoncoverage.txt";
	my $geneCoverage          = $poolName . "_ALL_genecoverage.txt";
	my $canonicalExonCoverage = $poolName . "_ALL_Canonical_exoncoverage.txt";

# if(-e "$outdir/$somaticMutIndelFile")
# {
#     `$QSUB -q $queue -V -wd $outdir -N AnnotateVariants.$$ -o AnnotateVariants.$$.stdout -e AnnotateVariants.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1  -b y $PERL $AnnotateAssessFilterVariants -si $outdir/$somaticMutIndelFile -c $config_file -t $titleFile -o $outdir -ec $outdir/$exonCoverage -gc $outdir/$geneCoverage -d $deleteIntermediateFiles`;
#     `$QSUB -q $queue -V -wd $outdir -hold_jid AnnotateVariants.$$ -N NotifyAV.$$ -e NotifyAV.$$.stderr -o NotifyAV.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
#     push(@notifyNames, "NotifyAV.$$.stat");
# }
# else
# {
# Annotate variants:
	$somaticMutIndelFile =
	  $poolName . "_AllSomaticMutIndel_withAlleleDepth.txt";
	if ( $CLUSTER eq "SGE" ) {
		my $cmd =
"$QSUB -q $queue -V -wd $outdir -N AnnotateVariants.$$ -o AnnotateVariants.$$.stdout -e AnnotateVariants.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1  -b y $PERL $AnnotateAssessFilterVariants -si $outdir/$somaticMutIndelFile -c $config_file -t $titleFile -o $outdir -d $deleteIntermediateFiles";
		$logger->debug("COMMAND : $cmd");
`$QSUB -q $queue -V -wd $outdir -N AnnotateVariants.$$ -o AnnotateVariants.$$.stdout -e AnnotateVariants.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1  -b y $PERL $AnnotateAssessFilterVariants -si $outdir/$somaticMutIndelFile -c $config_file -t $titleFile -o $outdir -d $deleteIntermediateFiles`;
`$QSUB -q $queue -V -wd $outdir -hold_jid AnnotateVariants.$$ -N NotifyAV.$$ -e NotifyAV.$$.stderr -o NotifyAV.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
	}
	else {
		my $cmd =
"$BSUB -q $queue -cwd $outdir -J AnnotateVariants.$$ -o AnnotateVariants.$$.%J.stdout -e AnnotateVariants.$$.%J.stderr -R \"rusage[mem=2]\" -M 4 -n 1 \"$PERL $AnnotateAssessFilterVariants -si $outdir/$somaticMutIndelFile -c $config_file -t $titleFile -o $outdir -d $deleteIntermediateFiles\"";
		$logger->debug("COMMAND : $cmd");
`$BSUB -q $queue -cwd $outdir -J AnnotateVariants.$$ -o AnnotateVariants.$$.%J.stdout -e AnnotateVariants.$$.%J.stderr -We 0:59 -R "rusage[mem=2]" -M 4 -n 1 "$PERL $AnnotateAssessFilterVariants -si $outdir/$somaticMutIndelFile -c $config_file -t $titleFile -o $outdir -d $deleteIntermediateFiles"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(AnnotateVariants.$$)" -J NotifyAV.$$ -e NotifyAV.$$.%J.stderr -o NotifyAV.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
	}
	push( @notifyNames, "NotifyAV.$$.stat" );
	&WaitToFinish( $outdir, @notifyNames );

	# Filter annotated and gentoyped variants
	$somaticMutIndelFile = $poolName
	  . "_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedExonic.txt";
	if ( $CLUSTER eq "SGE" ) {
		my $cmd =
"$QSUB -q $queue -V -wd $outdir -N FilterExonicVariants.$$ -o FilterVariants.$$.stdout -e FilterVariants.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1  -b y $PERL $filterGenotypedVariants -i $outdir/$somaticMutIndelFile -h $HotspotMutations -ce $ClinicalExons -t $titleFile -dp_snv $dp_snv -ad_snv $ad_snv -vf_snv $vf_snv -dp_snvHS $dp_snvHS -ad_snvHS $ad_snvHS -vf_snvHS $vf_snvHS -dp_indel $dp_indel -ad_indel $ad_indel -vf_indel $vf_indel -dp_indelHS $dp_indelHS -ad_indelHS $ad_indelHS -vf_indelHS $vf_indelHS -occurrence $occurrencePercent -tn_ratio $TNfreqRatioThreshold -mt $MAFthreshold ";
		$logger->debug("COMMAND : $cmd");
		@notifyNames = ();
`$QSUB -q $queue -V -wd $outdir -N FilterExonicVariants.$$ -o FilterVariants.$$.stdout -e FilterVariants.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1  -b y $PERL $filterGenotypedVariants -i $outdir/$somaticMutIndelFile -h $HotspotMutations -ce $ClinicalExons -t $titleFile -dp_snv $dp_snv -ad_snv $ad_snv -vf_snv $vf_snv -dp_snvHS $dp_snvHS -ad_snvHS $ad_snvHS -vf_snvHS $vf_snvHS -dp_indel $dp_indel -ad_indel $ad_indel -vf_indel $vf_indel -dp_indelHS $dp_indelHS -ad_indelHS $ad_indelHS -vf_indelHS $vf_indelHS -occurrence $occurrencePercent -tn_ratio $TNfreqRatioThreshold -mt $MAFthreshold`;
`$QSUB -q $queue -V -wd $outdir -hold_jid FilterExonicVariants.$$ -N NotifyFV.$$ -e NotifyFV.$$.stderr -o NotifyFV.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
	}
	else {
		my $cmd =
"$BSUB -q $queue -cwd $outdir -J FilterExonicVariants.$$ -o FilterVariants.$$.%J.stdout -e FilterVariants.$$.%J.stderr -R \"rusage[mem=2]\" -M 4 -n 1 \"$PERL $filterGenotypedVariants -i $outdir/$somaticMutIndelFile -h $HotspotMutations -ce $ClinicalExons -t $titleFile -dp_snv $dp_snv -ad_snv $ad_snv -vf_snv $vf_snv -dp_snvHS $dp_snvHS -ad_snvHS $ad_snvHS -vf_snvHS $vf_snvHS -dp_indel $dp_indel -ad_indel $ad_indel -vf_indel $vf_indel -dp_indelHS $dp_indelHS -ad_indelHS $ad_indelHS -vf_indelHS $vf_indelHS -occurrence $occurrencePercent -tn_ratio $TNfreqRatioThreshold -mt $MAFthreshold\"";
		$logger->debug("COMMAND : $cmd");
		@notifyNames = ();
`$BSUB -q $queue -cwd $outdir -J FilterExonicVariants.$$ -o FilterVariants.$$.%J.stdout -e FilterVariants.$$.%J.stderr -We 0:59 -R "rusage[mem=2]" -M 4 -n 1 "$PERL $filterGenotypedVariants -i $outdir/$somaticMutIndelFile -h $HotspotMutations -ce $ClinicalExons -t $titleFile -dp_snv $dp_snv -ad_snv $ad_snv -vf_snv $vf_snv -dp_snvHS $dp_snvHS -ad_snvHS $ad_snvHS -vf_snvHS $vf_snvHS -dp_indel $dp_indel -ad_indel $ad_indel -vf_indel $vf_indel -dp_indelHS $dp_indelHS -ad_indelHS $ad_indelHS -vf_indelHS $vf_indelHS -occurrence $occurrencePercent -tn_ratio $TNfreqRatioThreshold -mt $MAFthreshold"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(FilterExonicVariants.$$)" -J NotifyFV.$$ -e NotifyFV.$$.%J.stderr -o NotifyFV.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
	}
	push( @notifyNames, "NotifyFV.$$.stat" );
	&WaitToFinish( $outdir, @notifyNames );
	my $filteredExonicVariants = $poolName
	  . "_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedExonic.Filtered.txt";
	`cp $filteredExonicVariants annotated_exonic_variants.txt`;
	$somaticMutIndelFile = $poolName
	  . "_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedNonPanelExonic.txt";
	if ( $CLUSTER eq "SGE" ) {
		my $cmd =
"$QSUB -q $queue -V -wd $outdir -N FilterNpVariants.$$ -o FilterVariants.$$.stdout -e FilterVariants.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1  -b y $PERL $filterGenotypedVariants -i $outdir/$somaticMutIndelFile -h $HotspotMutations -ce $ClinicalExons -t $titleFile -dp_snv $dp_snv -ad_snv $ad_snv -vf_snv $vf_snv -dp_snvHS $dp_snvHS -ad_snvHS $ad_snvHS -vf_snvHS $vf_snvHS -dp_indel $dp_indel -ad_indel $ad_indel -vf_indel $vf_indel -dp_indelHS $dp_indelHS -ad_indelHS $ad_indelHS -vf_indelHS $vf_indelHS -occurrence $occurrencePercent -tn_ratio $TNfreqRatioThreshold -mt $MAFthreshold ";
		$logger->debug("COMMAND : $cmd");
		@notifyNames = ();
`$QSUB -q $queue -V -wd $outdir -N FilterNpVariants.$$ -o FilterVariants.$$.stdout -e FilterVariants.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1  -b y $PERL $filterGenotypedVariants -i $outdir/$somaticMutIndelFile -h $HotspotMutations -ce $ClinicalExons -t $titleFile -dp_snv $dp_snv -ad_snv $ad_snv -vf_snv $vf_snv -dp_snvHS $dp_snvHS -ad_snvHS $ad_snvHS -vf_snvHS $vf_snvHS -dp_indel $dp_indel -ad_indel $ad_indel -vf_indel $vf_indel -dp_indelHS $dp_indelHS -ad_indelHS $ad_indelHS -vf_indelHS $vf_indelHS -occurrence $occurrencePercent -tn_ratio $TNfreqRatioThreshold -mt $MAFthreshold`;
`$QSUB -q $queue -V -wd $outdir -hold_jid FilterNpVariants.$$ -N NotifyNpFV.$$ -e NotifyFV.$$.stderr -o NotifyNpFV.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
	}
	else {
		my $cmd =
"$BSUB -q $queue -cwd $outdir -J FilterNpVariants.$$ -o FilterVariants.$$.%J.stdout -e FilterVariants.$$.%J.stderr -R \"rusage[mem=2]\" -M 4 -n 1 \"$PERL $filterGenotypedVariants -i $outdir/$somaticMutIndelFile -h $HotspotMutations -ce $ClinicalExons -t $titleFile -dp_snv $dp_snv -ad_snv $ad_snv -vf_snv $vf_snv -dp_snvHS $dp_snvHS -ad_snvHS $ad_snvHS -vf_snvHS $vf_snvHS -dp_indel $dp_indel -ad_indel $ad_indel -vf_indel $vf_indel -dp_indelHS $dp_indelHS -ad_indelHS $ad_indelHS -vf_indelHS $vf_indelHS -occurrence $occurrencePercent -tn_ratio $TNfreqRatioThreshold -mt $MAFthreshold\"";
		$logger->debug("COMMAND : $cmd");
		@notifyNames = ();
`$BSUB -q $queue -cwd $outdir -J FilterNpVariants.$$ -o FilterVariants.$$.%J.stdout -e FilterVariants.$$.%J.stderr -We 0:59 -R "rusage[mem=2]" -M 4 -n 1 "$PERL $filterGenotypedVariants -i $outdir/$somaticMutIndelFile -h $HotspotMutations -ce $ClinicalExons -t $titleFile -dp_snv $dp_snv -ad_snv $ad_snv -vf_snv $vf_snv -dp_snvHS $dp_snvHS -ad_snvHS $ad_snvHS -vf_snvHS $vf_snvHS -dp_indel $dp_indel -ad_indel $ad_indel -vf_indel $vf_indel -dp_indelHS $dp_indelHS -ad_indelHS $ad_indelHS -vf_indelHS $vf_indelHS -occurrence $occurrencePercent -tn_ratio $TNfreqRatioThreshold -mt $MAFthreshold"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(FilterNpVariants.$$)" -J NotifyFV.$$ -e NotifyFV.$$.stderr -o NotifyNpFV.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
	}
	push( @notifyNames, "NotifyNpFV.$$.stat" );
	&WaitToFinish( $outdir, @notifyNames );
	$filteredExonicVariants = $poolName
	  . "_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedNonPanelExonic.Filtered.txt";
	`cp $filteredExonicVariants annotated_nonpanel_exonic_variants.txt`;

	# Filter annotated and gentoyped SILENT variants
	$somaticMutIndelFile = $poolName
	  . "_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedSilent.txt";
	if ( $CLUSTER eq "SGE" ) {
		my $cmd =
"$QSUB -q $queue -V -wd $outdir -N FilterSilentVariants.$$ -o FilterVariants.$$.stdout -e FilterVariants.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1  -b y $PERL $filterGenotypedVariants -i $outdir/$somaticMutIndelFile -h $HotspotMutations -ce $ClinicalExons -t $titleFile -dp_snv $dp_snv -ad_snv $ad_snv -vf_snv $vf_snv -dp_snvHS $dp_snvHS -ad_snvHS $ad_snvHS -vf_snvHS $vf_snvHS -dp_indel $dp_indel -ad_indel $ad_indel -vf_indel $vf_indel -dp_indelHS $dp_indelHS -ad_indelHS $ad_indelHS -vf_indelHS $vf_indelHS -occurrence $occurrencePercent -tn_ratio $TNfreqRatioThreshold -mt $MAFthreshold ";
		$logger->debug("COMMAND : $cmd");
		@notifyNames = ();
`$QSUB -q $queue -V -wd $outdir -N FilterSilentVariants.$$ -o FilterSilentVariants.$$.stdout -e FilterSilentVariants.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1  -b y $PERL $filterGenotypedVariants -i $outdir/$somaticMutIndelFile -h $HotspotMutations -ce $ClinicalExons -t $titleFile -dp_snv $dp_snv -ad_snv $ad_snv -vf_snv $vf_snv -dp_snvHS $dp_snvHS -ad_snvHS $ad_snvHS -vf_snvHS $vf_snvHS -dp_indel $dp_indel -ad_indel $ad_indel -vf_indel $vf_indel -dp_indelHS $dp_indelHS -ad_indelHS $ad_indelHS -vf_indelHS $vf_indelHS -occurrence $occurrencePercent -tn_ratio $TNfreqRatioThreshold -mt $MAFthreshold`;
`$QSUB -q $queue -V -wd $outdir -hold_jid FilterSilentVariants.$$ -N NotifySilentFV.$$ -e NotifySilentFV.$$.stderr -o NotifySilentFV.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
	}
	else {
		my $cmd =
"$BSUB -q $queue -cwd $outdir -J FilterSilentVariants.$$ -o FilterSilentVariants.$$.%J.stdout -e FilterSilentVariants.$$.%J.stderr -R \"rusage[mem=2]\" -M 4 -n 1 \"$PERL $filterGenotypedVariants -i $outdir/$somaticMutIndelFile -h $HotspotMutations -ce $ClinicalExons -t $titleFile -dp_snv $dp_snv -ad_snv $ad_snv -vf_snv $vf_snv -dp_snvHS $dp_snvHS -ad_snvHS $ad_snvHS -vf_snvHS $vf_snvHS -dp_indel $dp_indel -ad_indel $ad_indel -vf_indel $vf_indel -dp_indelHS $dp_indelHS -ad_indelHS $ad_indelHS -vf_indelHS $vf_indelHS -occurrence $occurrencePercent -tn_ratio $TNfreqRatioThreshold -mt $MAFthreshold\"";
		$logger->debug("COMMAND : $cmd");
		@notifyNames = ();
`$BSUB -q $queue -cwd $outdir -J FilterSilentVariants.$$ -o FilterSilentVariants.$$.%J.stdout -e FilterSilentVariants.$$.%J.stderr -We 0:59 -R "rusage[mem=2]" -M 4 -n 1 "$PERL $filterGenotypedVariants -i $outdir/$somaticMutIndelFile -h $HotspotMutations -ce $ClinicalExons -t $titleFile -dp_snv $dp_snv -ad_snv $ad_snv -vf_snv $vf_snv -dp_snvHS $dp_snvHS -ad_snvHS $ad_snvHS -vf_snvHS $vf_snvHS -dp_indel $dp_indel -ad_indel $ad_indel -vf_indel $vf_indel -dp_indelHS $dp_indelHS -ad_indelHS $ad_indelHS -vf_indelHS $vf_indelHS -occurrence $occurrencePercent -tn_ratio $TNfreqRatioThreshold -mt $MAFthreshold"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(FilterSilentVariants.$$)" -J NotifySilentFV.$$ -e NotifySilentFV.$$.%J.stderr -o NotifySilentFV.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
	}
	push( @notifyNames, "NotifySilentFV.$$.stat" );
	&WaitToFinish( $outdir, @notifyNames );
	$filteredExonicVariants = $poolName
	  . "_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedSilent.Filtered.txt";
	`cp $filteredExonicVariants annotated_silent_variants.txt`;
	$somaticMutIndelFile = $poolName
	  . "_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedNonPanelSilent.txt";
	if ( $CLUSTER eq "SGE" ) {
		my $cmd =
"$QSUB -q $queue -V -wd $outdir -N FilterNpsVariants.$$ -o FilterVariants.$$.stdout -e FilterVariants.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1  -b y $PERL $filterGenotypedVariants -i $outdir/$somaticMutIndelFile -h $HotspotMutations -ce $ClinicalExons -t $titleFile -dp_snv $dp_snv -ad_snv $ad_snv -vf_snv $vf_snv -dp_snvHS $dp_snvHS -ad_snvHS $ad_snvHS -vf_snvHS $vf_snvHS -dp_indel $dp_indel -ad_indel $ad_indel -vf_indel $vf_indel -dp_indelHS $dp_indelHS -ad_indelHS $ad_indelHS -vf_indelHS $vf_indelHS -occurrence $occurrencePercent -tn_ratio $TNfreqRatioThreshold -mt $MAFthreshold ";
		$logger->debug("COMMAND : $cmd");
		@notifyNames = ();
`$QSUB -q $queue -V -wd $outdir -N FilterNpsVariants.$$ -o FilterSilentVariants.$$.stdout -e FilterSilentVariants.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1  -b y $PERL $filterGenotypedVariants -i $outdir/$somaticMutIndelFile -h $HotspotMutations -ce $ClinicalExons -t $titleFile -dp_snv $dp_snv -ad_snv $ad_snv -vf_snv $vf_snv -dp_snvHS $dp_snvHS -ad_snvHS $ad_snvHS -vf_snvHS $vf_snvHS -dp_indel $dp_indel -ad_indel $ad_indel -vf_indel $vf_indel -dp_indelHS $dp_indelHS -ad_indelHS $ad_indelHS -vf_indelHS $vf_indelHS -occurrence $occurrencePercent -tn_ratio $TNfreqRatioThreshold -mt $MAFthreshold`;
`$QSUB -q $queue -V -wd $outdir -hold_jid FilterNpsVariants.$$ -N NotifyNpsFV.$$ -e NotifySilentFV.$$.stderr -o NotifyNpsFV.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
	}
	else {
		my $cmd =
"$BSUB -q $queue -cwd $outdir -J FilterNpsVariants.$$ -o FilterSilentVariants.$$.%J.stdout -e FilterSilentVariants.$$.%J.stderr -R \"rusage[mem=2]\" -M 4 -n 1 \"$PERL $filterGenotypedVariants -i $outdir/$somaticMutIndelFile -h $HotspotMutations -ce $ClinicalExons -t $titleFile -dp_snv $dp_snv -ad_snv $ad_snv -vf_snv $vf_snv -dp_snvHS $dp_snvHS -ad_snvHS $ad_snvHS -vf_snvHS $vf_snvHS -dp_indel $dp_indel -ad_indel $ad_indel -vf_indel $vf_indel -dp_indelHS $dp_indelHS -ad_indelHS $ad_indelHS -vf_indelHS $vf_indelHS -occurrence $occurrencePercent -tn_ratio $TNfreqRatioThreshold -mt $MAFthreshold\"";
		$logger->debug("COMMAND : $cmd");
		@notifyNames = ();
`$BSUB -q $queue -cwd $outdir -J FilterNpsVariants.$$ -o FilterSilentVariants.$$.%J.stdout -e FilterSilentVariants.$$.%J.stderr -We 0:59 -R "rusage[mem=2]" -M 4 -n 1 "$PERL $filterGenotypedVariants -i $outdir/$somaticMutIndelFile -h $HotspotMutations -ce $ClinicalExons -t $titleFile -dp_snv $dp_snv -ad_snv $ad_snv -vf_snv $vf_snv -dp_snvHS $dp_snvHS -ad_snvHS $ad_snvHS -vf_snvHS $vf_snvHS -dp_indel $dp_indel -ad_indel $ad_indel -vf_indel $vf_indel -dp_indelHS $dp_indelHS -ad_indelHS $ad_indelHS -vf_indelHS $vf_indelHS -occurrence $occurrencePercent -tn_ratio $TNfreqRatioThreshold -mt $MAFthreshold"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(FilterNpsVariants.$$)" -J NotifyNpsFV.$$ -e NotifySilentFV.$$.%J.stderr -o NotifyNpsFV.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
	}
	push( @notifyNames, "NotifyNpsFV.$$.stat" );
	&WaitToFinish( $outdir, @notifyNames );
	$filteredExonicVariants = $poolName
	  . "_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedNonPanelSilent.Filtered.txt";
	`cp $filteredExonicVariants annotated_nonpanel_silent_variants.txt`;

	# Generate VCF files
	$somaticMutIndelFile = $poolName
	  . "_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedExonic.Filtered.txt";
	if ( $CLUSTER eq "SGE" ) {
		my $cmd =
"$QSUB -q $queue -V -wd $outdir -N GenerateVCFfiles.$$ -o GenerateVCFfiles.$$.stdout -e GenerateVCFfiles.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1  -b y $PERL $createPatientVCFfile -v $somaticMutIndelFile -o $outdir -t $titleFile -j $JAVA_1_7 -it $IGVtools";
		$logger->debug("COMMAND : $cmd");
		@notifyNames = ();
`$QSUB -q $queue -V -wd $outdir -N GenerateVCFfiles.$$ -o GenerateVCFfiles.$$.stdout -e GenerateVCFfiles.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1  -b y $PERL $createPatientVCFfile -v $somaticMutIndelFile -o $outdir  -t $titleFile -j $JAVA_1_7 -it $IGVtools`;
`$QSUB -q $queue -V -wd $outdir -hold_jid GenerateVCFfiles.$$ -N NotifyGenVCF.$$ -e NotifyGenVCF.$$.stderr -o NotifyGenVCF.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
	}
	else {
		my $cmd =
"$BSUB -q $queue -cwd $outdir -J GenerateVCFfiles.$$ -o GenerateVCFfiles.$$.%J.stdout -e GenerateVCFfiles.$$.%J.stderr -R \"rusage[mem=2]\" -M 4 -n 1 \"$PERL $createPatientVCFfile -v $somaticMutIndelFile -o $outdir  -t $titleFile -j $JAVA_1_7 -it $IGVtools\"";
		$logger->debug("COMMAND : $cmd");
		@notifyNames = ();
`$BSUB -q $queue -cwd $outdir -J GenerateVCFfiles.$$ -o GenerateVCFfiles.$$.%J.stdout -e GenerateVCFfiles.$$.%J.stderr -We 0:59 -R "rusage[mem=2]" -M 4 -n 1 "$PERL $createPatientVCFfile -v $somaticMutIndelFile -o $outdir  -t $titleFile -j $JAVA_1_7 -it $IGVtools"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(GenerateVCFfiles.$$)" -J NotifyGenVCF.$$ -e NotifyGenVCF.$$.%J.stderr -o NotifyGenVCF.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
	}
	push( @notifyNames, "NotifyGenVCF.$$.stat" );
	&WaitToFinish( $outdir, @notifyNames );

	# Generate Coverage files and clinical reports
	$somaticMutIndelFile = $poolName
	  . "_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedExonic.Filtered.txt";
	my $geneLevelCNfile = $poolName . "_copynumber_segclusp.genes.txt";
	my $intragenicLevelCNfile =
	  $poolName . "_copynumber_segclusp.intragenic.txt";
	my $normalsUSed = $poolName . "_NormalUsedInMutationCalling.txt";
	if ( $CLUSTER eq "SGE" ) {
		my $cmd =
"$QSUB -q $queue -V -wd $outdir -N GenerateVCFfiles.$$ -o GenerateVCFfiles.$$.stdout -e GenerateVCFfiles.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1  -b y $PERL $createCoverageFiles -v $somaticMutIndelFile -o $outdir -nu $normalsUSed -t $titleFile -ss $sampleFile -id $intragenicLevelCNfile -gc $outdir/$geneCoverage -cn $geneLevelCNfile -gi $GeneIntervalAnn -ec $outdir/$exonCoverage -cc $outdir/$canonicalExonCoverage -ve $validatedExons -ce $ClinicalExons -ei $exonIntervalsFile -tf $translationFolder  -ct $coverageThreshold";
		$logger->debug("COMMAND : $cmd");
		@notifyNames = ();
`$QSUB -q $queue -V -wd $outdir -N GenerateCoverageFiles.$$ -o GenerateCoverageFiles.$$.stdout -e GenerateCoverageFiles.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1  -b y $PERL $createCoverageFiles -v $somaticMutIndelFile -o $outdir -nu $normalsUSed -t $titleFile -ss $sampleFile -id $intragenicLevelCNfile -gc $outdir/$geneCoverage -cn $geneLevelCNfile -gi $GeneIntervalAnn -ec $outdir/$exonCoverage -cc $outdir/$canonicalExonCoverage -ve $validatedExons -ce $ClinicalExons -ei $exonIntervalsFile -tf $translationFolder  -ct $coverageThreshold`;
`$QSUB -q $queue -V -wd $outdir -hold_jid GenerateCoverageFiles.$$ -N NotifyCovFile.$$ -e NotifyCovFile.$$.stderr -o NotifyCovFile.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
	}
	else {
		my $cmd =
"$BSUB -q $queue -cwd $outdir -J GenerateCoverageFiles.$$ -o GenerateCoverageFiles.$$.%J.stdout -e GenerateCoverageFiles.$$.%J.stderr -R \"rusage[mem=2]\" -M 4 -n 1 \"$PERL $createCoverageFiles -v $somaticMutIndelFile -o $outdir -nu $normalsUSed -t $titleFile -ss $sampleFile -id $intragenicLevelCNfile -gc $outdir/$geneCoverage -cn $geneLevelCNfile -gi $GeneIntervalAnn -ec $outdir/$exonCoverage -cc $outdir/$canonicalExonCoverage -ve $validatedExons -ce $ClinicalExons -ei $exonIntervalsFile -tf $translationFolder  -ct $coverageThreshold\"";
		$logger->debug("COMMAND : $cmd");
		@notifyNames = ();
`$BSUB -q $queue -cwd $outdir -J GenerateCoverageFiles.$$ -o GenerateCoverageFiles.$$.%J.stdout -e GenerateCoverageFiles.$$.%J.stderr -We 0:59 -R "rusage[mem=2]" -M 4 -n 1 "$PERL $createCoverageFiles -v $somaticMutIndelFile -o $outdir -nu $normalsUSed -t $titleFile -ss $sampleFile -id $intragenicLevelCNfile -gc $outdir/$geneCoverage -cn $geneLevelCNfile -gi $GeneIntervalAnn -ec $outdir/$exonCoverage -cc $outdir/$canonicalExonCoverage -ve $validatedExons -ce $ClinicalExons -ei $exonIntervalsFile -tf $translationFolder  -ct $coverageThreshold"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(GenerateCoverageFiles.$$)" -J NotifyCovFile.$$ -e NotifyCovFile.$$.%J.stderr -o NotifyCovFile.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
	}
	push( @notifyNames, "NotifyCovFile.$$.stat" );
	&WaitToFinish( $outdir, @notifyNames );

	#Merge Filter Variants
	my $filteredVariantsFile = $poolName
	  . "_AllSomaticMutIndel_withAlleleDepth_annovarAnnotated_All_Filtered.txt";
	open( OFH, ">", "$outdir/$filteredVariantsFile" )
	  or die(
		$logger->fatal(
"MergeFilteredVariants:Cannot open $outdir/$filteredVariantsFile, Error:$!"
		)
	  );

	#my @FilterVariantFiles = glob("$outdir/*annovarAnnotated*.Filtered.txt");
	my $filecount = 0;
	foreach my $filterFile ( glob "$outdir/*annovarAnnotated*.Filtered.txt" ) {
		open( FH, "$filterFile" )
		  or die(
			$logger->fatal(
				"MergeFilterVariants:Cannot open $outdir/$filterFile, Error:$!")
		  );
		if ( $filecount == 0 ) {
			my $header = <FH>;
			print OFH "$header";
			while (<FH>) {
				print OFH "$_";
			}
		}
		else {
			while (<FH>) {
				next if ( $. == 1 );
				print OFH "$_";
			}
		}
		close(FH);
		$filecount++;
	}
	close(OFH);
	return ( \@names );
}
#####################################
#####################################
#Clip adapter sequences.
sub RunTrimGalore {
	my ( $file1, $file2, $outdir, $adaptorList, $id ) = @_;
	my %barcodeList = %$adaptorList;
	my ($barcode) = $file1 =~ /.*_(bc\d+)_.*/;
	my $adapter1 = $barcodeList{$barcode};
	my $adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT";
	my ($basename1) = $file1 =~ /(.*)\.fastq.gz/;
	my $outFilename1 = "$basename1" . "_cl.fastq.gz";
	my ($basename2) = $file2 =~ /(.*)\.fastq.gz/;
	my $outFilename2 = "$basename2" . "_cl.fastq.gz";

	if (    ( -e "$outFilename1" )
		and ( ( -s "$outFilename1" ) != 0 )
		and ( -e "$outFilename2" )
		and ( ( -s "$outFilename2" ) != 0 ) )
	{
		$logger->info(
"Files:\n$outFilename1\n$outFilename2\n they exists and process will not run to clip adapters in them."
		);
		return ( "$outFilename1", "$outFilename2", 'NULL' );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N Clipping.$id.$$ -o Clipping.$id.$$.stdout -e Clipping.$id.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1  -b y \"$PERL $TrimGalore --paired --gzip -q 1 --suppress_warn --stringency 3 -length 25 -o $outdir -a $adapter1 -a2 $adapter2 $file1 $file2\"";
				$logger->debug("COMMAND : $cmd");
`$QSUB -q $queue -V -wd $outdir -N Clipping.$id.$$ -o Clipping.$id.$$.stdout -e Clipping.$id.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1  -b y "$PERL $TrimGalore --paired --gzip -q 1 --suppress_warn --stringency 3 -length 25 -o $outdir -a $adapter1 -a2 $adapter2 $file1 $file2"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid Clipping.$id.$$ -N NotifyCR.$id.$$ -e NotifyCR.$id.$$.stderr -o NotifyCR.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"`$BSUB -q $queue -cwd $outdir -J Clipping.$id.$$ -o Clipping.$id.$$.%J.stdout -e Clipping.$id.$$.%J.stderr -R \"rusage[mem=2]\" -M 4 -n 1 \"$PERL $TrimGalore --paired --gzip -q 1 --suppress_warn --stringency 3 -length 25 -o $outdir -a $adapter1 -a2 $adapter2 $file1 $file2\"";
				$logger->debug("COMMAND : $cmd");
`$BSUB -q $queue -cwd $outdir -J Clipping.$id.$$ -o Clipping.$id.$$.%J.stdout -e Clipping.$id.$$.%J.stderr -We 0:59 -R "rusage[mem=2]" -M 4 -n 1 "$PERL $TrimGalore --paired --gzip -q 1 --suppress_warn --stringency 3 -length 25 -o $outdir -a $adapter1 -a2 $adapter2 $file1 $file2"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(Clipping.$id.$$)" -J NotifyCR.$id.$$ -e NotifyCR.$id.$$.%J.stderr -o NotifyCR.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal("Clipping:Job Submission Failed, Error:$@");
			exit(1);
		}
	}
	return ( "$outFilename1", "$outFilename2", "NotifyCR.$id.$$.stat" );
}
#####################################
#####################################
#BWA Find Suffix Array(SA) Co-ordinates
sub RunBwaAln {
	my ( $file, $outdir, $id, $readType ) = @_;
	my ($basename) = $file =~ /(.*)\.fastq.gz/;
	my ($outFilename) = "$basename" . ".sai";
	if ( ( -e "$outFilename" ) and ( ( -s "$outFilename" ) != 0 ) ) {
		$logger->info(
"Files:\n$outFilename\n they exists and process will not run to make \".sai\" file."
		);
		return ( "$outFilename", 'NULL' );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N bwaAln.$readType.$id.$$ -o bwaAln.$readType.$id.$$.stdout -e bwaAln.$readType.$id.$$.stderr -l h_vmem=5G,virtual_free=5G -pe smp 1 -b y '$BWA aln -f $outFilename $Reference $file'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N bwaAln.$readType.$id.$$ -o bwaAln.$readType.$id.$$.stdout -e bwaAln.$readType.$id.$$.stderr -l h_vmem=5G,virtual_free=5G -pe smp 1 -b y "$BWA aln -f $outFilename $Reference $file"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid bwaAln.$readType.$id.$$ -N NotifyBwaAln.$id.$$ -e NotifyBwaAln.$readType.$id.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1 -o NotifyBwaAln.$readType.$id.$$.stat -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J bwaAln.$readType.$id.$$ -o bwaAln.$readType.$id.$$.%J.stdout -e bwaAln.$readType.$id.$$.%J.stderr -R \"rusage[mem=5]\" -M 8 -n 1 \"$BWA aln -f $outFilename $Reference $file\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J bwaAln.$readType.$id.$$ -o bwaAln.$readType.$id.$$.%J.stdout -e bwaAln.$readType.$id.$$.%J.stderr -We 24:00 -R "rusage[mem=5]" -M 8 -n 1 "$BWA aln -f $outFilename $Reference $file"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(bwaAln.$readType.$id.$$)" -J NotifyBwaAln.$id.$$ -e NotifyBwaAln.$readType.$id.$$.%J.stderr -o NotifyBwaAln.$readType.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal("BWA ALN:Job Submission Failed, Error:$@");
			exit(1);
		}
	}
	return ( "$outFilename", "NotifyBwaAln.$readType.$id.$$.stat" );
}
#####################################
#####################################
#BWA SAMPE align SA.
sub RunBwaSampe {
	my ( $clippedfile1, $clippedfile2, $SAfile1, $SAfile2, $outdir, $id ) = @_;
	my ($basename) = $clippedfile1 =~ /(.*)_R1.*\.fastq.gz/;
	my $outFilename = "$basename" . "_mrg_cl_aln.sam";
	if ( $basename =~ /\// ) {
		$basename = basename($basename);
	}
	my @sampleDetails = split( "_bc", $basename );
	my $sampleId = $sampleDetails[0];
	my ($barcode) = $basename =~ /.*_(bc\d+)_.*/;
	my ($pool)    = $basename =~ /.*_bc\d+_(.*)_L\d{1,3}_.*/;
	if ( ( -e "$outFilename" ) and ( ( -s "$outFilename" ) != 0 ) ) {
		$logger->info(
"Files:\n$outFilename\n they exists and process will not run to make \"_aln.sam\" file."
		);
		return ( "$outFilename", 'NULL' );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N bwaSampe.$id.$$ -o bwaSampe.$id.$$.stdout -e bwaSampe.$id.$$.stderr -l h_vmem=6G,virtual_free=6G -pe smp 1 -b y '$BWA sampe -f $outFilename $Reference $SAfile1 $SAfile2 $clippedfile1 $clippedfile2'";
				$logger->debug("COMMAND: $cmd");

#`$QSUB -q $queue -V -wd $outdir -N bwaSampe.$id.$$ -o bwaSampe.$id.$$.stdout -e bwaSampe.$id.$$.stderr -b y "$BWA sampe -r \'\@RG\tID:$basename\tLB:$id\tSM:$sampleId\tPL:Illumina\tPU:$barcode\tCN:BergerLab_MSKCC\' -f $outFilename $Reference $SAfile1 $SAfile2 $clippedfile1 $clippedfile2"`;
`$QSUB -q $queue -V -wd $outdir -N bwaSampe.$id.$$ -o bwaSampe.$id.$$.stdout -e bwaSampe.$id.$$.stderr -l h_vmem=6G,virtual_free=6G -pe smp 1 -b y "$BWA sampe -f $outFilename $Reference $SAfile1 $SAfile2 $clippedfile1 $clippedfile2"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid bwaSampe.$id.$$ -N NotifyBwaSampe.$id.$$ -e NotifyBwaSampe.$id.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1 -o NotifyBwaSampe.$id.$$.stat -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J bwaSampe.$id.$$ -o bwaSampe.$id.$$.%J.stdout -e bwaSampe.$id.$$.%J.stderr -R \"rusage[mem=6]\" -M 8 -n 1 \"$BWA sampe -f $outFilename $Reference $SAfile1 $SAfile2 $clippedfile1 $clippedfile2\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J bwaSampe.$id.$$ -o bwaSampe.$id.$$.%J.stdout -e bwaSampe.$id.$$.%J.stderr -We 24:00 -R "rusage[mem=6]" -M 8 -n 1 "$BWA sampe -f $outFilename $Reference $SAfile1 $SAfile2 $clippedfile1 $clippedfile2"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(bwaSampe.$id.$$)" -J NotifyBwaSampe.$id.$$ -e NotifyBwaSampe.$id.$$.%J.stderr -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 -o NotifyBwaSampe.$id.$$.stat "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal("BWA SAMPE:Job Submission Failed, Error:$@");
			exit(1);
		}
	}
	return ( "$outFilename", "NotifyBwaSampe.$id.$$.stat" );
}
#####################################
#####################################
#BWA MEM to align fastq.
sub RunBwaMem {
	my ( $fastq1, $fastq2, $outdir, $id ) = @_;
	my ($basename) = $fastq1 =~ /(.*)_R1.*\.fastq.gz/;
	my $outFilename = "$basename" . "_mrg_cl_aln.sam";
	if ( $basename =~ /\// ) {
		$basename = basename($basename);
	}
	my @sampleDetails = split( "_bc", $basename );
	my $sampleId = $sampleDetails[0];
	my ($barcode) = $basename =~ /.*_(bc\d+)_.*/;
	my ($pool)    = $basename =~ /.*bc\d+_(.*)_L\d{1,3}_.*/;
	if ( ( -e "$outFilename" ) and ( ( -s "$outFilename" ) != 0 ) ) {
		$logger->info(
"Files:\n$outFilename\n they exists and process will not run to make \"_aln.sam\" file."
		);
		return ( "$outFilename", 'NULL' );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -wd $outdir -N bwaMem.$id.$$ -l h_vmem=10G,virtual_free=10G -pe smp 4 -o $outFilename -e /dev/null -b y '$BWA mem -t 4 -PM -R \'\@RG\tID:$basename\tLB:$id\tSM:$sampleId\tPL:Illumina\tPU:$barcode\tCN:BergerLab_MSKCC\' $Reference $fastq1 $fastq2'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -wd $outdir -N bwaMem.$id.$$ -l h_vmem=10G,virtual_free=10G -pe smp 4 -o $outFilename -e /dev/null -b y "$BWA mem -t 4 -PM -R \'\@RG\tID:$basename\tLB:$id\tSM:$sampleId\tPL:Illumina\tPU:$barcode\tCN:BergerLab_MSKCC\' $Reference $fastq1 $fastq2"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid bwaMem.$id.$$ -N NotifyBwaMem.$id.$$ -l h_vmem=2G,virtual_free=2G -pe smp 1 -e /dev/null -o NotifyBwaMem.$id.$$.stat -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J bwaMem.$id.$$ -R \"rusage[mem=10]\" -M 40 -n 4 -o bwaMem.$id.$$.stdout -e bwaMem.$id.$$.stderr \"$BWA mem -t 4 -PM -R \'\@RG\tID:$basename\tLB:$id\tSM:$sampleId\tPL:Illumina\tPU:$barcode\tCN:BergerLab_MSKCC\' $Reference $fastq1 $fastq2 > $outFilename\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J bwaMem.$id.$$ -We 24:00 -R "rusage[mem=10]" -M 40 -n 4 -o bwaMem.$id.$$.stdout -e bwaMem.$id.$$.stderr "$BWA mem -t 4 -PM -R \'\@RG\tID:$basename\tLB:$id\tSM:$sampleId\tPL:Illumina\tPU:$barcode\tCN:BergerLab_MSKCC\' $Reference $fastq1 $fastq2 > $outFilename"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(bwaMem.$id.$$)" -J NotifyBwaMem.$id.$$ -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 2 -n 1 -o NotifyBwaMem.$id.$$.stat "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal("BEA MEM:Job Submission Failed, Error:$@");
			exit(1);
		}
	}
	return ( "$outFilename", "NotifyBwaMem.$id.$$.stat" );
}
#####################################
#####################################
#Sort Sam file
sub RunSortSam {
	my ( $samFile, $outdir, $id ) = @_;
	my ($basename) = $samFile =~ /(.*)_mrg_cl_aln.sam/;
	my $outFilename = $samFile;
	$outFilename =~ s/\.sam/_srt\.bam/g;
	my @sampleDetails = split( "_bc", $basename );
	my $sampleId = $sampleDetails[0];
	my ($barcode) = $basename =~ /.*_(bc\d+)_.*/;
	my ($pool)    = $basename =~ /.*_bc\d+_(.*)_L\d{1,3}_.*/;
	my $platform  = "Illumina";

	if ( ( -e "$outFilename" ) and ( ( -s "$outFilename" ) != 0 ) ) {
		$logger->info(
"Files:\n$outFilename\n they exists and process will not run to make \"_srt.bam\" file."
		);
		return ( "$outFilename", 'NULL' );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N SortSam.$id.$$ -o SortSam.$id.$$.stdout -e SortSam.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y '$JAVA_1_7 -Xmx4g -jar $PICARD/AddOrReplaceReadGroups.jar I=$samFile O=$outFilename SO=coordinate RGID=$basename RGLB=$id RGPL=$platform RGPU=$barcode RGSM=$sampleId RGCN=MSKCC TMP_DIR=$TMPDIR COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N SortSam.$id.$$ -o SortSam.$id.$$.stdout -e SortSam.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y "$JAVA_1_7 -Xmx4g -jar $PICARD/AddOrReplaceReadGroups.jar I=$samFile O=$outFilename SO=coordinate RGID=$basename RGLB=$id RGPL=$platform RGPU=$barcode RGSM=$sampleId RGCN=MSKCC TMP_DIR=$TMPDIR COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid SortSam.$id.$$ -N NotifySortSam.$id.$$ -e NotifySortSam.$id.$$.stderr -o NotifySortSam.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J SortSam.$id.$$ -o SortSam.$id.$$.%J.stdout -e SortSam.$id.$$.%J.stderr -R \"rusage[mem=8]\" -M 12 -n 1 \"$JAVA_1_7 -Xmx4g -jar $PICARD/AddOrReplaceReadGroups.jar I=$samFile O=$outFilename SO=coordinate RGID=$basename RGLB=$id RGPL=$platform RGPU=$barcode RGSM=$sampleId RGCN=MSKCC TMP_DIR=$TMPDIR COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J SortSam.$id.$$ -o SortSam.$id.$$.%J.stdout -e SortSam.$id.$$.%J.stderr -We 6:00 -R "rusage[mem=8]" -M 12 -n 1 "$JAVA_1_7 -Xmx4g -jar $PICARD/AddOrReplaceReadGroups.jar I=$samFile O=$outFilename SO=coordinate RGID=$basename RGLB=$id RGPL=$platform RGPU=$barcode RGSM=$sampleId RGCN=MSKCC TMP_DIR=$TMPDIR COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(SortSam.$id.$$)" -J NotifySortSam.$id.$$ -e NotifySortSam.$id.$$.%J.stderr -o NotifySortSam.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal("SortSam:Job Submission Failed, Error:$@");
			exit(1);
		}
	}
	return ( "$outFilename", "NotifySortSam.$id.$$.stat" );
}
#####################################
#####################################
#FxMate information & sort file
sub RunFixMateInformation {
	my ( $bamFile, $outdir, $id ) = @_;
	my $outFilename = $bamFile;
	$outFilename =~ s/\.bam/_FX\.bam/g;
	if ( ( -e "$outFilename" ) and ( ( -s "$outFilename" ) != 0 ) ) {
		$logger->info(
"Files:\n$outFilename\n they exists and process will not run to make \"_srt.bam\" file."
		);
		return ( "$outFilename", 'NULL' );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N FixMate.$id.$$ -o FixMate.$id.$$.stdout -e FixMate.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 4 -b y '$JAVA_1_7 -Xmx24g -jar $PICARD/FixMateInformation.jar I=$bamFile O=$outFilename SO=coordinate TMP_DIR=$TMPDIR COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N FixMate.$id.$$ -o FixMate.$id.$$.stdout -e FixMate.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 4 -b y "$JAVA_1_7 -Xmx24g -jar $PICARD/FixMateInformation.jar I=$bamFile O=$outFilename SO=coordinate TMP_DIR=$TMPDIR COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid FixMate.$id.$$ -N NotifyFixMate.$id.$$ -e NotifyFixMate.$id.$$.stderr -o NotifyFixMate.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J FixMate.$id.$$ -o FixMate.$id.$$.%J.stdout -e FixMate.$id.$$.%J.stderr -R \"rusage[mem=32]\" -M 40 -n 4 \"$JAVA_1_7 -Xmx24g -jar $PICARD/FixMateInformation.jar I=$bamFile O=$outFilename SO=coordinate TMP_DIR=$TMPDIR COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J FixMate.$id.$$ -o FixMate.$id.$$.%J.stdout -e FixMate.$id.$$.%J.stderr -We 6:00 -R "rusage[mem=32]" -M 40 -n 4 "$JAVA_1_7 -Xmx24g -jar $PICARD/FixMateInformation.jar I=$bamFile O=$outFilename SO=coordinate TMP_DIR=$TMPDIR COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(FixMate.$id.$$)" -J NotifyFixMate.$id.$$ -e NotifyFixMate.$id.$$.%J.stderr -o NotifyFixMate.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal("FixMate:Job Submission Failed, Error:$@");
			exit(1);
		}
	}
	return ( "$outFilename", "NotifyFixMate.$id.$$.stat" );
}
#####################################
#####################################
#Mark Duplicates in Bam
sub RunMarkDuplicates {
	my ( $bamFile, $outdir, $id ) = @_;
	my $outFilename     = $bamFile;
	my $metricsFilename = $bamFile;
	$outFilename     =~ s/\.bam/_MD\.bam/g;
	$metricsFilename =~ s/\.bam/_MD\.metrics/g;
	if (    ( -e "$outdir/$outFilename" )
		and ( ( -s "$outdir/$outFilename" ) != 0 ) )
	{
		$logger->info(
"Files:\n$outFilename\n they exists and process will not run to make \"_MD.bam\" file."
		);
		return ( "$outFilename", 'NULL' );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N MD.$id.$$ -o MD.$id.$$.stdout -e MD.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y '$JAVA_1_7 -Xmx4g -jar $PICARD/MarkDuplicates.jar I=$bamFile O=$outFilename ASSUME_SORTED=true METRICS_FILE=$metricsFilename TMP_DIR=$TMPDIR COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N MD.$id.$$ -o MD.$id.$$.stdout -e MD.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y "$JAVA_1_7 -Xmx4g -jar $PICARD/MarkDuplicates.jar I=$bamFile O=$outFilename ASSUME_SORTED=true METRICS_FILE=$metricsFilename TMP_DIR=$TMPDIR COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid MD.$id.$$ -N NotifyMD.$id.$$ -e NotifyMD.$id.$$.stderr -o NotifyMD.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J MD.$id.$$ -o MD.$id.$$.%J.stdout -e MD.$id.$$.%J.stderr -R \"rusage[mem=8]\" -M 12 -n 1 \"$JAVA_1_7 -Xmx4g -jar $PICARD/MarkDuplicates.jar I=$bamFile O=$outFilename ASSUME_SORTED=true METRICS_FILE=$metricsFilename TMP_DIR=$TMPDIR COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J MD.$id.$$ -o MD.$id.$$.%J.stdout -e MD.$id.$$.%J.stderr -We 6:00 -R "rusage[mem=8]" -M 12 -n 1 "$JAVA_1_7 -Xmx4g -jar $PICARD/MarkDuplicates.jar I=$bamFile O=$outFilename ASSUME_SORTED=true METRICS_FILE=$metricsFilename TMP_DIR=$TMPDIR COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(MD.$id.$$)" -J NotifyMD.$id.$$ -e NotifyMD.$id.$$.%J.stderr -o NotifyMD.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal("MarkDuplicates:Job Submission Failed, Error:$@");
			exit(1);
		}
	}
	return ( "$outFilename", "NotifyMD.$id.$$.stat" );
}
#####################################
#####################################
#Running Realign Target Creator from GATK
sub RunRealignerTargetCreator {
	my ( $outFilename, $inputFiles, $outdir ) = @_;
	$inputFiles =~ s/,/ /g;
	my ($jobName) = $outFilename =~ /(.*)_IndelRealigner.intervals/;
	if (    ( -e "$outdir/$outFilename" )
		and ( ( -s "$outdir/$outFilename" ) != 0 ) )
	{
		$logger->info(
"Files:\n$outFilename\n they exists and process will not run to make \"_IndelRealigner.intervals\" file."
		);
		return ("NULL");
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N RTC.$jobName.$$ -o RTC.$jobName.$$.stdout -e RTC.$jobName.$$.stderr -l h_vmem=12G,virtual_free=12G -pe smp 2 -b y '$JAVA_1_7 -Xmx20g -jar $GATK -T RealignerTargetCreator -I $outdir/$inputFiles -R $Reference -L $GeneInterval -o $outFilename --known $Mills_1000G_Indels'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N RTC.$jobName.$$ -o RTC.$jobName.$$.stdout -e RTC.$jobName.$$.stderr -l h_vmem=12G,virtual_free=12G -pe smp 2 -b y "$JAVA_1_7 -Xmx20g -jar $GATK -T RealignerTargetCreator -I $outdir/$inputFiles -R $Reference -L $GeneInterval -o $outFilename --known $Mills_1000G_Indels"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid RTC.$jobName.$$ -N NotifyRTC.$jobName.$$ -e NotifyRTC.$jobName.$$.stderr -o NotifyRTC.$jobName.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"`$BSUB -q $queue -cwd $outdir -J RTC.$jobName.$$ -o RTC.$jobName.$$.%J.stdout -e RTC.$jobName.$$.%J.stderr -R \"rusage[mem=12]\" -M 24 -n 2 \"$JAVA_1_7 -Xmx20g -jar $GATK -T RealignerTargetCreator -I $outdir/$inputFiles -R $Reference -L $GeneInterval -o $outFilename --known $Mills_1000G_Indels\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J RTC.$jobName.$$ -o RTC.$jobName.$$.%J.stdout -e RTC.$jobName.$$.%J.stderr -We 24:00 -R "rusage[mem=12]" -M 24 -n 2 "$JAVA_1_7 -Xmx20g -jar $GATK -T RealignerTargetCreator -I $outdir/$inputFiles -R $Reference -L $GeneInterval -o $outFilename --known $Mills_1000G_Indels"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(RTC.$jobName.$$)" -J NotifyRTC.$jobName.$$ -e NotifyRTC.$jobName.$$.%J.stderr -o NotifyRTC.$jobName.$$.%J.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal(
				"RealignTargetCreator:Job Submission Failed, Error:$@");
			exit(1);
		}
	}
	return ("NotifyRTC.$jobName.$$.stat");
}
##########################################################################
#Running Indel Realignment from GATK
sub RunIndelRealigner {
	my ( $file, $interval, $outdir, $id ) = @_;
	my ($basename) = $file =~ /(.*)\.bam/;
	my $outFilename = $basename . "_IR.bam";
	if (    ( -e "$outdir/$outFilename" )
		and ( ( -s "$outdir/$outFilename" ) != 0 ) )
	{
		$logger->info(
"Files:\n$outFilename\n they exists and process will not run to make \"_IR.bam\" file."
		);
		return ( "$outFilename", 'NULL' );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N IR.$id.$$ -o IR.$id.$$.stdout -e IR.$id.$$.stderr -l h_vmem=12G,virtual_free=12G -pe smp 2 -b y '$JAVA_1_7 -Xmx20g -Djava.io.tmpdir=$TMPDIR -jar $GATK -T IndelRealigner -I $file -R $Reference -targetIntervals $interval -o $outFilename -baq RECALCULATE -known $Mills_1000G_Indels'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N IR.$id.$$ -o IR.$id.$$.stdout -e IR.$id.$$.stderr -l h_vmem=12G,virtual_free=12G -pe smp 2 -b y "$JAVA_1_7 -Xmx20g -Djava.io.tmpdir=$TMPDIR -jar $GATK -T IndelRealigner -I $file -R $Reference -targetIntervals $interval -o $outFilename -baq RECALCULATE -known $Mills_1000G_Indels"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid IR.$id.$$ -N NotifyIR.$id.$$ -e NotifyIR.$id.$$.stderr -o NotifyIR.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J IR.$id.$$ -o IR.$id.$$.%J.stdout -e IR.$id.$$.%J.stderr -R \"rusage[mem=12]\" -M 30 -n 2 \"$JAVA_1_7 -Xmx20g -Djava.io.tmpdir=$TMPDIR -jar $GATK -T IndelRealigner -I $file -R $Reference -targetIntervals $interval -o $outFilename -baq RECALCULATE -known $Mills_1000G_Indels\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J IR.$id.$$ -o IR.$id.$$.%J.stdout -e IR.$id.$$.%J.stderr -We 24:00 -R "rusage[mem=12]" -M 30 -n 2 "$JAVA_1_7 -Xmx20g -Djava.io.tmpdir=$TMPDIR -jar $GATK -T IndelRealigner -I $file -R $Reference -targetIntervals $interval -o $outFilename -baq RECALCULATE -known $Mills_1000G_Indels"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(IR.$id.$$)" -J NotifyIR.$id.$$ -e NotifyIR.$id.$$.%J.stderr -o NotifyIR.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal("IndelRealigner:Job Submission Failed, Error:$@");
			exit(1);
		}
	}
	return ( "$outFilename", "NotifyIR.$id.$$.stat" );
}
#####################################
#####################################
#Running Base Quality Recalibration from GATK
sub RunBaseQualityRecalibration {
	my ( $files, $outdir, $id ) = @_;
	my $basename;
	my $outFilename1;

	# my $outFilename2;
	#my $plot_pdf;
	if ( $id eq "LaneLevel" ) {
		$outFilename1 = "recalibration-report.grp";

		#$outFilename2 = "recalibration-report.tmp";
		#$plot_pdf = "recalibration-report-plot.pdf";
		if (    ( -e "$outdir/$outFilename1" )
			and ( ( -s "$outdir/$outFilename1" ) != 0 ) )
		{
			$logger->info(
				"Files:\n$outFilename1\n they exists and process will not run."
			);
			return ( "$outFilename1", 'NULL' );
		}
		else {
			eval {
				if ( $CLUSTER eq "SGE" )
				{
					my $cmd =
"$QSUB -q $queue -V -wd $outdir -o BQSR1.$id.$$.stdout -e BQSR1.$id.$$.stderr -N BQSR1.$id.$$ -l h_vmem=3G,virtual_free=3G -pe smp 10 -b y '$JAVA_1_7 -Xmx20g -Djava.io.tmpdir=$TMPDIR -jar $GATK -T BaseRecalibrator -nct 10 -I $files -plots -R $Reference -knownSites $dbSNP -knownSites $Mills_1000G_Indels -rf BadCigar -o $outFilename1'";
					$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -o BQSR1.$id.$$.stdout -e BQSR1.$id.$$.stderr -N BQSR1.$id.$$ -l h_vmem=3G,virtual_free=3G -pe smp 10 -b y "$JAVA_1_7 -Xmx20g -Djava.io.tmpdir=$TMPDIR -jar $GATK -T BaseRecalibrator -nct 10 -I $files -plots -R $Reference -knownSites $dbSNP -knownSites $Mills_1000G_Indels -rf BadCigar -o $outFilename1"`;

#`$QSUB -q $queue -V -wd $outdir -l mem_free=10G -pe alloc 8 -o /dev/null -e /dev/null -hold_jid BQSR1.$id.$$ -N BQSR2.$id.$$ -b y "$JAVA_1_7 -Xmx20g -Djava.io.tmpdir=$TMPDIR -jar $GATK -T BaseRecalibrator -nct 8 -I $files -plots $plot_pdf -R $Reference -knownSites $dbSNP -knownSites $Mills_1000G_Indels -o $outFilename2 -BQSR $outdir/$outFilename1"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid BQSR1.$id.$$ -N NotifyBQSR.$id.$$ -e NotifyBQSR.$id.$$.stderr -o NotifyBQSR.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
				}
				else {
					my $cmd =
"$BSUB -q $queue -cwd $outdir -o BQSR1.$id.$$.%J.stdout -e BQSR1.$id.$$.%J.stderr -J BQSR1.$id.$$ -R \"rusage[mem=3]\" -M 36 -n 10 \"$JAVA_1_7 -Xmx20g -Djava.io.tmpdir=$TMPDIR -jar $GATK -T BaseRecalibrator -nct 10 -I $files -plots -R $Reference -knownSites $dbSNP -knownSites $Mills_1000G_Indels -rf BadCigar -o $outFilename1\"";
					$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -o BQSR1.$id.$$.%J.stdout -e BQSR1.$id.$$.%J.stderr -J BQSR1.$id.$$ -We 24:00 -R "rusage[mem=3]" -M 36 -n 10 "$JAVA_1_7 -Xmx20g -Djava.io.tmpdir=$TMPDIR -jar $GATK -T BaseRecalibrator -nct 10 -I $files -plots -R $Reference -knownSites $dbSNP -knownSites $Mills_1000G_Indels -rf BadCigar -o $outFilename1"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(BQSR1.$id.$$)" -J NotifyBQSR.$id.$$ -e NotifyBQSR.$id.$$.%J.stderr -o NotifyBQSR.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
				}
			};
			if ($@) {
				$logger->fatal(
					"RunBQSRLaneLevel:Job Submission Failed, Error:$@");
				exit(1);
			}
		}
	}
	else {
		($basename) = $files =~ /(.*)\.bam/;
		$outFilename1 = $basename . "_recalReport.grp";

		#$outFilename2 = $basename . "_recalReport.tmp";
		#$plot_pdf  = $basename . "_recalPlot.pdf";
		if (    ( -e "$outdir/$outFilename1" )
			and ( ( -s "$outdir/$outFilename1" ) != 0 ) )
		{
			$logger->info(
"Files:\n$outFilename1\n they exists and process will not run to make \"-recal.grp\" file."
			);
			return ( "$outFilename1", 'NULL' );
		}
		else {
			eval {
				if ( $CLUSTER eq "SGE" )
				{
					my $cmd =
"$QSUB -q $queue -V -wd $outdir -l h_vmem=8G,virtual_free=8G -pe smp 5 -o BQSR1.$id.$$.stdout -e BQSR1.$id.$$.stderr -N BQSR1.$id.$$ -b y '$JAVA_1_7 -Xmx30g -Djava.io.tmpdir=$TMPDIR -jar $GATK -T BaseRecalibrator -I $files -R $Reference -knownSites $dbSNP -knownSites $Mills_1000G_Indels -rf BadCigar -o $outFilename1 -nct 3'";
					$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -l h_vmem=8G,virtual_free=8G -pe smp 5 -o BQSR1.$id.$$.stdout -e BQSR1.$id.$$.stderr -N BQSR1.$id.$$ -b y "$JAVA_1_7 -Xmx30g -Djava.io.tmpdir=$TMPDIR -jar $GATK -T BaseRecalibrator -I $files -R $Reference -knownSites $dbSNP -knownSites $Mills_1000G_Indels -rf BadCigar -o $outFilename1 -nct 3"`;

#`$QSUB -q $queue -V -wd $outdir -l mem_free=10G -o /dev/null -e /dev/null -hold_jid BQSR1.$id.$$ -N BQSR2.$id.$$ -b y "$JAVA_1_7 -Xmx20g -Djava.io.tmpdir=$TMPDIR -jar $GATK -T BaseRecalibrator -I $files -plots $plot_pdf -R $Reference -knownSites $dbSNP -knownSites $Mills_1000G_Indels -o $outFilename2 -BQSR $outdir/$outFilename1"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid BQSR1.$id.$$ -N NotifyBQSR.$id.$$ -e NotifyBQSR.$id.$$.stderr -o NotifyBQSR.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
				}
				else {
					my $cmd =
"$BSUB -q $queue -cwd $outdir -R \"rusage[mem=8]\" -M 44 -n 5 -o BQSR1.$id.$$.%J.stdout -e BQSR1.$id.$$.%J.stderr -J BQSR1.$id.$$ \"$JAVA_1_7 -Xmx30g -Djava.io.tmpdir=$TMPDIR -jar $GATK -T BaseRecalibrator -I $files -R $Reference -knownSites $dbSNP -knownSites $Mills_1000G_Indels -rf BadCigar -o $outFilename1 -nct 3\"";
					$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -We 24:00 -R "rusage[mem=8]" -M 44 -n 5 -o BQSR1.$id.$$.%J.stdout -e BQSR1.$id.$$.%J.stderr -J BQSR1.$id.$$ "$JAVA_1_7 -Xmx30g -Djava.io.tmpdir=$TMPDIR -jar $GATK -T BaseRecalibrator -I $files -R $Reference -knownSites $dbSNP -knownSites $Mills_1000G_Indels -rf BadCigar -o $outFilename1 -nct 3"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(BQSR1.$id.$$)" -J NotifyBQSR.$id.$$ -e NotifyBQSR.$id.$$.%J.stderr -o NotifyBQSR.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
				}
			};
			if ($@) {
				$logger->fatal(
					"RunBQSRSampleLevel:Job Submission Failed, Error:$@");
				exit(1);
			}
		}
	}
	return ( "$outFilename1", "NotifyBQSR.$id.$$.stat" );
}
#####################################
#####################################
#Running PrintReads from GATK
sub PrintBQSRreads {
	my ( $file, $BQSRtable, $outdir, $id ) = @_;
	my ($basename) = $file =~ /(.*)\.bam/;
	my $outFilename = $basename . "_BR.bam";
	if (    ( -e "$outdir/$outFilename" )
		and ( ( -s "$outdir/$outFilename" ) != 0 ) )
	{
		$logger->info(
"Files:\n$outFilename\n they exists and process will not run to make \"-BR.bam\" file."
		);
		return ( "$outFilename", 'NULL' );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N PrintBQSR.$id.$$ -o PrintBQSR.$id.$$.stdout -e PrintBQSR.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y '$JAVA_1_7 -Xmx4g -Djava.io.tmpdir=$TMPDIR -jar $GATK -T PrintReads -I $file -R $Reference -baq RECALCULATE -BQSR $BQSRtable -EOQ -o $outFilename'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N PrintBQSR.$id.$$ -o PrintBQSR.$id.$$.stdout -e PrintBQSR.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 3 -b y "$JAVA_1_7 -Xmx4g -Djava.io.tmpdir=$TMPDIR -jar $GATK -T PrintReads -I $file -R $Reference -baq RECALCULATE -BQSR $BQSRtable -nct 8 -EOQ -o $outFilename"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid PrintBQSR.$id.$$ -N NotifyPrintBQSR.$id.$$ -e NotifyPrintBQSR.$id.$$.stderr -o NotifyPrintBQSR.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J PrintBQSR.$id.$$ -o PrintBQSR.$id.$$.%J.stdout -e PrintBQSR.$id.$$.%J.stderr -R \"rusage[mem=8]\" -M 30 -n 3 \"$JAVA_1_7 -Xmx4g -Djava.io.tmpdir=$TMPDIR -jar $GATK -T PrintReads -I $file -R $Reference -baq RECALCULATE -BQSR $BQSRtable -nct 8 -EOQ -o $outFilename\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J PrintBQSR.$id.$$ -o PrintBQSR.$id.$$.%J.stdout -e PrintBQSR.$id.$$.%J.stderr -We 24:00 -R "rusage[mem=8]" -M 30 -n 3 "$JAVA_1_7 -Xmx4g -Djava.io.tmpdir=$TMPDIR -jar $GATK -T PrintReads -I $file -R $Reference -baq RECALCULATE -BQSR $BQSRtable -nct 8 -EOQ -o $outFilename"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(PrintBQSR.$id.$$)" -J NotifyPrintBQSR.$id.$$ -e NotifyPrintBQSR.$id.$$.stderr -o NotifyPrintBQSR.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal("PrintReads:Job Submission Failed, Error:$@");
			exit(1);
		}
	}
	return ( "$outFilename", "NotifyPrintBQSR.$id.$$.stat" );
}
#####################################
#####################################
#Metrics Calculations for Bam
sub RunMetricsCalculations {
	my ( $bamFile, $outdir, $id ) = @_;
	my ($basename) = $bamFile =~ /(.*)\.bam/;
	my $HSmetricsFilename         = $basename . ".HSmetrics.txt";
	my $perTargetCoverageFilename = $basename . ".target.covg";
	my $InsertSizeMetricsFilename = $basename . ".insert_size_metrics";
	my $BaseQualitiesFilename     = $basename . ".quality_by_cycle_metrics";
	my $geneCoverage              = $basename . ".gene.covg";
	my $gene_nomapqCoverage       = $basename . ".gene_nomapq.covg";
	my $geneCount                 = $basename . ".gene.cn.txt";
	my $canonicalExonCoverage     = $basename . ".canonical.exon.covg";
	my $canonicalExonCoverageOutput =
	  $basename . ".canonical.exon.covg.sample_interval_summary";
	my $tilingCoverage        = $basename . ".tiling.covg";
	my $tiling_nomapqCoverage = $basename . ".tiling_nomapq.covg";
	my $tilingCoverageOutput =
	  $basename . ".tiling.covg.sample_interval_summary";
	my $fingerprintCounts  = $basename . ".FP.counts";
	my $fingerprintSummary = $basename . ".FP.summary.txt";

	#my $countErrors        = $basename . ".CountErrors.txt";
	my @notifynames   = ();
	my @metricsOutput = ();

	#Calculate Hybrid Selection specific metrics
	if (    ( -e "$outdir/$HSmetricsFilename" )
		and ( ( -s "$outdir/$HSmetricsFilename" ) != 0 ) )
	{
		$logger->info(
"Files:\n$HSmetricsFilename\n they exists and process will not run to make \"-BR.HSmetrics.txt\" file."
		);
		push( @notifynames,   "NULL" );
		push( @metricsOutput, $HSmetricsFilename );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N HSmetrics.$id.$$ -o /dev/null -e /dev/null -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y '$JAVA_1_7 -Xmx4g -jar $PICARD/CalculateHsMetrics.jar I=$bamFile O=$HSmetricsFilename BI=$BaitInterval TI=$TargetInterval REFERENCE_SEQUENCE=$Reference PER_TARGET_COVERAGE=$perTargetCoverageFilename TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=LENIENT'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N HSmetrics.$id.$$ -o /dev/null -e /dev/null -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y "$JAVA_1_7 -Xmx4g -jar $PICARD/CalculateHsMetrics.jar I=$bamFile O=$HSmetricsFilename BI=$BaitInterval TI=$TargetInterval REFERENCE_SEQUENCE=$Reference PER_TARGET_COVERAGE=$perTargetCoverageFilename TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=LENIENT"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid HSmetrics.$id.$$ -N NotifyHSmetrics.$id.$$ -e /dev/null -o NotifyHSmetrics.$id.$$.stat -l h_vmem=5G,virtual_free=5G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J HSmetrics.$id.$$ -R \"rusage[mem=8]\" -M 12 -n 1 \"$JAVA_1_7 -Xmx4g -jar $PICARD/CalculateHsMetrics.jar I=$bamFile O=$HSmetricsFilename BI=$BaitInterval TI=$TargetInterval REFERENCE_SEQUENCE=$Reference PER_TARGET_COVERAGE=$perTargetCoverageFilename TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=LENIENT\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J HSmetrics.$id.$$ -o HSmetrics.$id.$$.stdout -e HSmetrics.$id.$$.stderr -We 0:59 -R "rusage[mem=8]" -M 12 -n 1 "$JAVA_1_7 -Xmx4g -jar $PICARD/CalculateHsMetrics.jar I=$bamFile O=$HSmetricsFilename BI=$BaitInterval TI=$TargetInterval REFERENCE_SEQUENCE=$Reference PER_TARGET_COVERAGE=$perTargetCoverageFilename TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=LENIENT"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(HSmetrics.$id.$$)" -J NotifyHSmetrics.$id.$$ -o NotifyHSmetrics.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal("HSmetrics:Job Submission Failed, Error:$@");
			exit(1);
		}
		push( @notifynames,   "NotifyHSmetrics.$id.$$.stat" );
		push( @metricsOutput, $HSmetricsFilename );
	}

	#Collect Insert Size and Mean Quality by Cycle metrics
	if (    ( -e "$outdir/$InsertSizeMetricsFilename" )
		and ( ( -s "$outdir/$InsertSizeMetricsFilename" ) != 0 )
		and ( -e "$outdir/$BaseQualitiesFilename" )
		and ( ( -s "$outdir/$BaseQualitiesFilename" ) != 0 ) )
	{
		$logger->info(
"Files:\n$InsertSizeMetricsFilename\n$BaseQualitiesFilename\n they exists and process will not run to make \"-BR.insert_size_metrics\" and \"-BR.quality_by_cycle_metrics\" files."
		);
		push( @notifynames,   "NULL" );
		push( @metricsOutput, $InsertSizeMetricsFilename );
		push( @metricsOutput, $BaseQualitiesFilename );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N InsertSizeMetrics.$id.$$ -o InsertSizeMetrics.$id.$$.stdout -e InsertSizeMetrics.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y '$JAVA_1_7 -Xmx4g -jar $PICARD/CollectMultipleMetrics.jar I=$bamFile O=$basename PROGRAM=null PROGRAM=CollectInsertSizeMetrics ASSUME_SORTED=true TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=LENIENT'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N InsertSizeMetrics.$id.$$ -o InsertSizeMetrics.$id.$$.stdout -e InsertSizeMetrics.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y "$JAVA_1_7 -Xmx4g -jar $PICARD/CollectMultipleMetrics.jar I=$bamFile O=$basename PROGRAM=null PROGRAM=CollectInsertSizeMetrics ASSUME_SORTED=true TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=LENIENT"`;
				$cmd =
"$QSUB -q $queue -V -wd $outdir -N MeanQualityByCycle.$id.$$ -o MeanQualityByCycle.$id.$$.stdout -e MeanQualityByCycle.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y '$JAVA_1_7 -Xmx4g -jar $PICARD/CollectMultipleMetrics.jar I=$bamFile O=$basename PROGRAM=null PROGRAM=MeanQualityByCycle ASSUME_SORTED=true TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=LENIENT'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N MeanQualityByCycle.$id.$$ -o MeanQualityByCycle.$id.$$.stdout -e MeanQualityByCycle.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y "$JAVA_1_7 -Xmx4g -jar $PICARD/CollectMultipleMetrics.jar I=$bamFile O=$basename PROGRAM=null PROGRAM=MeanQualityByCycle ASSUME_SORTED=true TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=LENIENT"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid "InsertSizeMetrics.$id.$$,MeanQualityByCycle.$id.$$" -N NotifyMultipleMetrics.$id.$$ -e NotifyMultipleMetrics.$id.$$.stderr -o NotifyMultipleMetrics.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J InsertSizeMetrics.$id.$$ -o InsertSizeMetrics.$id.$$.%J.stdout -e InsertSizeMetrics.$id.$$.%J.stderr -R \"rusage[mem=8]\" -M 12 -n 1 \"$JAVA_1_7 -Xmx4g -jar $PICARD/CollectMultipleMetrics.jar I=$bamFile O=$basename PROGRAM=null PROGRAM=CollectInsertSizeMetrics ASSUME_SORTED=true TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=LENIENT\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J InsertSizeMetrics.$id.$$ -o InsertSizeMetrics.$id.$$.%J.stdout -e InsertSizeMetrics.$id.$$.%J.stderr -We 0:59 -R "rusage[mem=8]" -M 12 -n 1 "$JAVA_1_7 -Xmx4g -jar $PICARD/CollectMultipleMetrics.jar I=$bamFile O=$basename PROGRAM=null PROGRAM=CollectInsertSizeMetrics ASSUME_SORTED=true TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=LENIENT"`;
				$cmd =
"$BSUB -q $queue -cwd $outdir -J MeanQualityByCycle.$id.$$ -o MeanQualityByCycle.$id.$$.%J.stdout -e MeanQualityByCycle.$id.$$.%J.stderr -R \"rusage[mem=8]\" -M 12 -n 1 \"$JAVA_1_7 -Xmx4g -jar $PICARD/CollectMultipleMetrics.jar I=$bamFile O=$basename PROGRAM=null PROGRAM=MeanQualityByCycle ASSUME_SORTED=true TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=LENIENT\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J MeanQualityByCycle.$id.$$ -o MeanQualityByCycle.$id.$$.%J.stdout -e MeanQualityByCycle.$id.$$.%J.stderr -We 0:59 -R "rusage[mem=8]" -M 12 -n 1 "$JAVA_1_7 -Xmx4g -jar $PICARD/CollectMultipleMetrics.jar I=$bamFile O=$basename PROGRAM=null PROGRAM=MeanQualityByCycle ASSUME_SORTED=true TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=LENIENT"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(InsertSizeMetrics.$id.$$) && post_done(MeanQualityByCycle.$id.$$)" -J NotifyMultipleMetrics.$id.$$ -e NotifyMultipleMetrics.$id.$$.%J.stderr -o NotifyMultipleMetrics.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal(
				"InsertSize_MeanQuality:Job Submission Failed, Error:$@");
			exit(1);
		}
		push( @notifynames,   "NotifyMultipleMetrics.$id.$$.stat" );
		push( @metricsOutput, $InsertSizeMetricsFilename );
		push( @metricsOutput, $BaseQualitiesFilename );
	}

	#Collect depth of coverage information: gene/exon, tiling, FP
	#Gene Coverage
	if ( ( -e "$outdir/$geneCount" ) and ( ( -s "$outdir/$geneCount" ) != 0 ) )
	{
		$logger->info(
"File:\n$geneCount\n they exists and process will not run to make \"-BR.gene.cn.txt\" file."
		);
		push( @notifynames,   "NULL" );
		push( @metricsOutput, $geneCount );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N GeneCoverage.$id.$$ -o GeneCoverage.$id.$$.stdout -e GeneCoverage.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y '$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $geneCoverage -L $GeneInterval  -rf BadCigar -mmq $MAPQ -mbq $BASQ -omitLocusTable -omitSampleSummary -omitBaseOutput'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N GeneCoverage.$id.$$ -o GeneCoverage.$id.$$.stdout -e GeneCoverage.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y "$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $geneCoverage -L $GeneInterval  -rf BadCigar -mmq $MAPQ -mbq $BASQ -omitLocusTable -omitSampleSummary -omitBaseOutput"`;
				$logger->debug(
"$QSUB -q $queue -V -wd $outdir -hold_jid GeneCoverage.$id.$$ -N ExonToGeneCov.$id.$$ -o $geneCount -e ExonToGeneCov.$id.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y '$PERL $ExonToGenCov $geneCoverage.sample_interval_summary $GeneCoord'"
				);
`$QSUB -q $queue -V -wd $outdir -hold_jid GeneCoverage.$id.$$ -N ExonToGeneCov.$id.$$ -o $geneCount -e ExonToGeneCov.$id.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$PERL $ExonToGenCov $geneCoverage.sample_interval_summary $GeneCoord"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid ExonToGeneCov.$id.$$ -N NotifyGeneCoverage.$id.$$ -e NotifyGeneCoverage.$id.$$.stderr -o NotifyGeneCoverage.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J GeneCoverage.$id.$$ -o GeneCoverage.$id.$$.%J.stdout -e GeneCoverage.$id.$$.%J.stderr -R \"rusage[mem=8]\" -M 12 -n 1 \"$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $geneCoverage -L $GeneInterval  -rf BadCigar -mmq $MAPQ -mbq $BASQ -omitLocusTable -omitSampleSummary -omitBaseOutput\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J GeneCoverage.$id.$$ -o GeneCoverage.$id.$$.%J.stdout -e GeneCoverage.$id.$$.%J.stderr -We 0:59 -R "rusage[mem=8]" -M 12 -n 1 "$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $geneCoverage -L $GeneInterval  -rf BadCigar -mmq $MAPQ -mbq $BASQ -omitLocusTable -omitSampleSummary -omitBaseOutput"`;
				$cmd =
"$BSUB -q $queue -cwd $outdir -w \"done(GeneCoverage.$id.$$)\" -J ExonToGeneCov.$id.$$ -o ExonToGeneCov.$id.$$.%J.stdout -e ExonToGeneCov.$id.$$.%J.stderr -R \"rusage[mem=2]\" -R \"rusage[iounits=0]\" -M 4 -n 1 \"$PERL $ExonToGenCov $geneCoverage.sample_interval_summary $GeneCoord > $geneCount\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -w "post_done(GeneCoverage.$id.$$)" -J ExonToGeneCov.$id.$$ -o ExonToGeneCov.$id.$$.%J.stdout -e ExonToGeneCov.$id.$$.%J.stderr -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$PERL $ExonToGenCov $geneCoverage.sample_interval_summary $GeneCoord > $geneCount"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(ExonToGeneCov.$id.$$)" -J NotifyGeneCoverage.$id.$$ -e NotifyGeneCoverage.$id.$$.stderr -o NotifyGeneCoverage.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal("GeneCoverage:Job Submission Failed, Error:$@");
			exit(1);
		}
		push( @notifynames,   "NotifyGeneCoverage.$id.$$.stat" );
		push( @metricsOutput, $geneCount );
	}

	#Gene_Nomapq Coverage
	if (
		    ( -e "${outdir}/${gene_nomapqCoverage}.sample_interval_summary" )
		and
		( ( -s "${outdir}/${gene_nomapqCoverage}.sample_interval_summary" ) !=
			0 )
	  )
	{
		$logger->info(
"File:\n${gene_nomapqCoverage}.sample_interval_summary\n they exists and process will not run to make \".gene_nomapq.covg\.sample_interval_summary\" file."
		);
		push( @notifynames, "NULL" );
		push( @metricsOutput,
			"${gene_nomapqCoverage}.sample_interval_summary" );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N GeneNomapqCoverage.$id.$$ -o GeneNomapqCoverage.$id.$$.stdout -e GeneNomapqCoverage.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y '$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $gene_nomapqCoverage -L $GeneInterval  -rf BadCigar -mmq 0 -mbq $BASQ -omitLocusTable -omitSampleSummary -omitBaseOutput'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N GeneNomapqCoverage.$id.$$ -o GeneNomapqCoverage.$id.$$.stdout -e GeneNomapqCoverage.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y "$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $gene_nomapqCoverage -L $GeneInterval  -rf BadCigar -mmq 0 -mbq $BASQ -omitLocusTable -omitSampleSummary -omitBaseOutput"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid GeneNomapqCoverage.$id.$$ -N NotifyGeneNomapqCoverage.$id.$$ -e NotifyGeneNomapqCoverage.$id.$$.stderr -o NotifyGeneNomapqCoverage.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J GeneNomapqCoverage.$id.$$ -o GeneNomapqCoverage.$id.$$.%J.stdout -e GeneNomapqCoverage.$id.$$.%J.stderr -R \"rusage[mem=8]\" -M 12 -n 1 \"$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $gene_nomapqCoverage -L $GeneInterval  -rf BadCigar -mmq 0 -mbq $BASQ -omitLocusTable -omitSampleSummary -omitBaseOutput\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J GeneNomapqCoverage.$id.$$ -o GeneNomapqCoverage.$id.$$.%J.stdout -e GeneNomapqCoverage.$id.$$.%J.stderr -We 0:59 -R "rusage[mem=8]" -M 12 -n 1 "$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $gene_nomapqCoverage -L $GeneInterval  -rf BadCigar -mmq 0 -mbq $BASQ -omitLocusTable -omitSampleSummary -omitBaseOutput"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(GeneNomapqCoverage.$id.$$)" -J NotifyGeneNomapqCoverage.$id.$$ -e NotifyGeneNomapqCoverage.$id.$$.%J.stderr -o NotifyGeneNomapqCoverage.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal(
				"GeneNomapqCoverage:Job Submission Failed, Error:$@");
			exit(1);
		}
		push( @notifynames, "NotifyGeneNomapqCoverage.$id.$$.stat" );
		push( @metricsOutput,
			"${gene_nomapqCoverage}.sample_interval_summary" );
	}

	#Canonical Exon Coverage
	if (    ( -e "$outdir/$canonicalExonCoverageOutput" )
		and ( ( -s "$outdir/$canonicalExonCoverageOutput" ) != 0 ) )
	{
		$logger->info(
"File:\$canonicalExonCoverage\n they exists and process will not run to make \"-BR.canonical.exon.covg\" file."
		);
		push( @notifynames,   "NULL" );
		push( @metricsOutput, $canonicalExonCoverageOutput );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N CanonicalExonCoverage.$id.$$ -o CanonicalExonCoverage.$id.$$.stdout -e CanonicalExonCoverage.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y '$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $canonicalExonCoverage -L $canonicalExonIntervalsFile  -rf BadCigar -mmq $MAPQ -mbq $BASQ -omitLocusTable -omitSampleSummary -omitBaseOutput'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N CanonicalExonCoverage.$id.$$ -o CanonicalExonCoverage.$id.$$.stdout -e CanonicalExonCoverage.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y "$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $canonicalExonCoverage -L $canonicalExonIntervalsFile -rf BadCigar -mmq $MAPQ -mbq $BASQ -omitLocusTable -omitSampleSummary -omitBaseOutput"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid CanonicalExonCoverage.$id.$$ -N NotifyCanonicalExonCoverage.$id.$$ -e NotifyCanonicalExonCoverage.$id.$$.stderr -o NotifyCanonicalExonCoverage.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J CanonicalExonCoverage.$id.$$ -o CanonicalExonCoverage.$id.$$.%J.stdout -e CanonicalExonCoverage.$id.$$.%J.stderr -R \"rusage[mem=8]\" -M 12 -n 1 \"$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $canonicalExonCoverage -L $canonicalExonIntervalsFile -rf BadCigar -mmq $MAPQ -mbq $BASQ -omitLocusTable -omitSampleSummary -omitBaseOutput\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J CanonicalExonCoverage.$id.$$ -o CanonicalExonCoverage.$id.$$.%J.stdout -e CanonicalExonCoverage.$id.$$.%J.stderr -We 0:59 -R "rusage[mem=8]" -M 12 -n 1 "$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $canonicalExonCoverage -L $canonicalExonIntervalsFile -rf BadCigar -mmq $MAPQ -mbq $BASQ -omitLocusTable -omitSampleSummary -omitBaseOutput"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(CanonicalExonCoverage.$id.$$)" -J NotifyCanonicalExonCoverage.$id.$$ -e NotifyCanonicalExonCoverage.$id.$$.%J.stderr -o NotifyCanonicalExonCoverage.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal(
				"CanonicalExonCoverage:Job Submission Failed, Error:$@");
			exit(1);
		}
		push( @notifynames,   "NotifyCanonicalExonCoverage.$id.$$.stat" );
		push( @metricsOutput, $canonicalExonCoverageOutput );
	}

	#Tiling Coverage
	if (    ( -e "$outdir/$tilingCoverageOutput" )
		and ( ( -s "$outdir/$tilingCoverageOutput" ) != 0 ) )
	{
		$logger->info(
"File:\n$tilingCoverage\nthey exists and process will not run to make \"-BR.tiling.covg\" file."
		);
		push( @notifynames,   "NULL" );
		push( @metricsOutput, $tilingCoverageOutput );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N TilingCoverage.$id.$$ -o TilingCoverage.$id.$$.stdout -e TilingCoverage.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y '$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $tilingCoverage -L $TilingInterval -rf BadCigar -mmq $MAPQ -mbq $BASQ -omitLocusTable -omitSampleSummary -omitBaseOutput'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N TilingCoverage.$id.$$ -o TilingCoverage.$id.$$.stdout -e TilingCoverage.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y "$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $tilingCoverage -L $TilingInterval -rf BadCigar -mmq $MAPQ -mbq $BASQ -omitLocusTable -omitSampleSummary -omitBaseOutput"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid TilingCoverage.$id.$$ -N NotifyTilingCoverage.$id.$$ -e NotifyTilingCoverage.$id.$$.stderr -o NotifyTilingCoverage.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J TilingCoverage.$id.$$ -o TilingCoverage.$id.$$.%J.stdout -e TilingCoverage.$id.$$.%J.stderr -R \"rusage[mem=8]\" -M 12 -n 1 \"$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $tilingCoverage -L $TilingInterval -rf BadCigar -mmq $MAPQ -mbq $BASQ -omitLocusTable -omitSampleSummary -omitBaseOutput\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J TilingCoverage.$id.$$ -o TilingCoverage.$id.$$.%J.stdout -e TilingCoverage.$id.$$.%J.stderr -We 0:59 -R "rusage[mem=8]" -M 12 -n 1 "$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $tilingCoverage -L $TilingInterval -rf BadCigar -mmq $MAPQ -mbq $BASQ -omitLocusTable -omitSampleSummary -omitBaseOutput"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(TilingCoverage.$id.$$)" -J NotifyTilingCoverage.$id.$$ -e NotifyTilingCoverage.$id.$$.stderr -o NotifyTilingCoverage.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal("TitlingCoverage:Job Submission Failed, Error:$@");
			exit(1);
		}
		push( @notifynames,   "NotifyTilingCoverage.$id.$$.stat" );
		push( @metricsOutput, $tilingCoverageOutput );
	}

	#Tiling_Nomapq Coverage
	if (
		    ( -e "${outdir}/${tiling_nomapqCoverage}.sample_interval_summary" )
		and
		( ( -s "${outdir}/${tiling_nomapqCoverage}.sample_interval_summary" )
			!= 0 )
	  )
	{
		$logger->info(
"File:\n${tiling_nomapqCoverage}.sample_interval_summary\n they exists and process will not run to make \".tiling_nomapq.covg\.sample_interval_summary\" file."
		);
		push( @notifynames, "NULL" );
		push( @metricsOutput,
			"${tiling_nomapqCoverage}.sample_interval_summary" );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N TilingNomapqCoverage.$id.$$ -o TilingNomapqCoverage.$id.$$.stdout -e TilingNomapqCoverage.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y '$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $tiling_nomapqCoverage -L $TilingInterval  -rf BadCigar -mmq 0 -mbq $BASQ -omitLocusTable -omitSampleSummary -omitBaseOutput'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N TilingNomapqCoverage.$id.$$ -o TilingNomapqCoverage.$id.$$.stdout -e TilingNomapqCoverage.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y "$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $tiling_nomapqCoverage -L $TilingInterval  -rf BadCigar -mmq 0 -mbq $BASQ -omitLocusTable -omitSampleSummary -omitBaseOutput"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid TilingNomapqCoverage.$id.$$ -N NotifyTilingNomapqCoverage.$id.$$ -e NotifyTilingNomapqCoverage.$id.$$.stderr -o NotifyTilingNomapqCoverage.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J TilingNomapqCoverage.$id.$$ -o TilingNomapqCoverage.$id.$$.%J.stdout -e TilingNomapqCoverage.$id.$$.%J.stderr -R \"rusage[mem=8]\" -M 12 -n 1 \"$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $tiling_nomapqCoverage -L $TilingInterval  -rf BadCigar -mmq 0 -mbq $BASQ -omitLocusTable -omitSampleSummary -omitBaseOutput\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J TilingNomapqCoverage.$id.$$ -o TilingNomapqCoverage.$id.$$.%J.stdout -e TilingNomapqCoverage.$id.$$.%J.stderr -We 0:59 -R "rusage[mem=8]" -M 12 -n 1 "$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $tiling_nomapqCoverage -L $TilingInterval  -rf BadCigar -mmq 0 -mbq $BASQ -omitLocusTable -omitSampleSummary -omitBaseOutput"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(TilingNomapqCoverage.$id.$$)" -J NotifyTilingNomapqCoverage.$id.$$ -e NotifyTilingNomapqCoverage.$id.$$.stderr -o NotifyTilingNomapqCoverage.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal(
				"TilingNomapqCoverage:Job Submission Failed, Error:$@");
			exit(1);
		}
		push( @notifynames, "NotifyTilingNomapqCoverage.$id.$$.stat" );
		push( @metricsOutput,
			"${tiling_nomapqCoverage}.sample_interval_summary" );
	}

	#FingerPrint Counts
	if (    ( -e "$outdir/$fingerprintSummary" )
		and ( ( -s "$outdir/$fingerprintSummary" ) != 0 ) )
	{
		$logger->info(
"File:\n$fingerprintSummary\n they exists and process will not run to make \"-BR.FP.summary.txt\" file."
		);
		push( @notifynames,   "NULL" );
		push( @metricsOutput, $fingerprintSummary );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N FingerPrint.$id.$$ -o FingerPrint.$id.$$.stdout -e FingerPrint.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y '$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $fingerprintCounts -L $FingerPrintInterval -rf BadCigar -mmq $MAPQ -mbq $BASQ -omitLocusTable -omitSampleSummary -omitIntervals -baseCounts'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N FingerPrint.$id.$$ -o FingerPrint.$id.$$.stdout -e FingerPrint.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y "$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $fingerprintCounts -L $FingerPrintInterval -rf BadCigar -mmq $MAPQ -mbq $BASQ -omitLocusTable -omitSampleSummary -omitIntervals -baseCounts"`;
				$logger->debug(
"$QSUB -q $queue -V -wd $outdir -hold_jid FingerPrint.$id.$$ -N FingerPrintSummary.$id.$$ -o $fingerprintSummary -e FingerPrint.$id.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y '$PERL $FPGenotypesScript $fingerprintCounts $FP_genotypes'"
				);
`$QSUB -q $queue -V -wd $outdir -hold_jid FingerPrint.$id.$$ -N FingerPrintSummary.$id.$$ -o $fingerprintSummary -e FingerPrint.$id.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$PERL $FPGenotypesScript $fingerprintCounts $FP_genotypes"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid FingerPrintSummary.$id.$$ -N NotifyFingerPrintSummary.$id.$$ -e NotifyFingerPrintSummary.$id.$$.stderr -o NotifyFingerPrintSummary.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J FingerPrint.$id.$$ -o FingerPrint.$id.$$.%J.stdout -e FingerPrint.$id.$$.%J.stderr -R \"rusage[mem=8]\" -M 12 -n 1 \"$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $fingerprintCounts -L $FingerPrintInterval -rf BadCigar -mmq $MAPQ -mbq $BASQ -omitLocusTable -omitSampleSummary -omitIntervals -baseCounts\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J FingerPrint.$id.$$ -o FingerPrint.$id.$$.%J.stdout -e FingerPrint.$id.$$.%J.stderr -We 0:59 -R "rusage[mem=8]" -M 12 -n 1 "$JAVA_1_7 -Xmx4g -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $fingerprintCounts -L $FingerPrintInterval -rf BadCigar -mmq $MAPQ -mbq $BASQ -omitLocusTable -omitSampleSummary -omitIntervals -baseCounts"`;
				$cmd =
"$BSUB -q $queue -cwd $outdir -w \"done(FingerPrint.$id.$$)\" -J FingerPrintSummary.$id.$$ -o FingerPrint.$id.$$.%J.stdout -e FingerPrint.$id.$$.%J.stderr -R \"rusage[mem=2]\" -M 4 -n 1 \"$PERL $FPGenotypesScript $fingerprintCounts $FP_genotypes > $fingerprintSummary\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -w "post_done(FingerPrint.$id.$$)" -J FingerPrintSummary.$id.$$ -o FingerPrint.$id.$$.%J.stdout -e FingerPrint.$id.$$.%J.stderr -We 0:59 -R "rusage[mem=2]" -M 4 -n 1 "$PERL $FPGenotypesScript $fingerprintCounts $FP_genotypes > $fingerprintSummary"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(FingerPrintSummary.$id.$$)" -J NotifyFingerPrintSummary.$id.$$ -e NotifyFingerPrintSummary.$id.$$.%J.stderr -o NotifyFingerPrintSummary.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal("FingerPrint:Job Submission Failed, Error:$@");
			exit(1);
		}
		push( @notifynames,   "NotifyFingerPrintSummary.$id.$$.stat" );
		push( @metricsOutput, $fingerprintSummary );
	}

=begin
	#Calculate G>T errors
	if (     ( -e "$outdir/$countErrors" )
		 and ( ( -s "$outdir/$countErrors" ) != 0 ) )
	{
		$logger->info(
"Files:\n$countErrors\n they exists and process will not run to make \"-BR.zCounts.txt\" file."
		);
		push( @notifynames,   "NULL" );
		push( @metricsOutput, $countErrors );
	} else
	{
		eval {
			my $cmd =
"$QSUB -q $queue -V -wd $outdir -N CountErrors.$id.$$ -o /dev/null -e /dev/null -l h_vmem=5G,virtual_free=5G -pe smp 6 -b y '$CountErrors --fasta $Reference --bam $bamFile --target $TargetRegionBED --output $countErrors'";
			$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N CountErrors.$id.$$ -o /dev/null -e /dev/null -l h_vmem=5G,virtual_free=5G -pe smp 6 -b y '$CountErrors --fasta $Reference --bam $bamFile --target $TargetRegionBED --output $countErrors'`;
`$QSUB -q $queue -V -wd $outdir -hold_jid CountErrors.$id.$$ -N NotifyCountErrors.$id.$$ -e /dev/null -o NotifyCountErrors.$id.$$.stat -l h_vmem=5G,virtual_free=5G -pe smp 1 -b y "$outdir/Notify.csh"`;
		};
		if ($@)
		{
			$logger->fatal("CountErrors:Job Submission Failed, Error:$@");
			exit(1);
		}
		push( @notifynames,   "NotifyCountErrors.$id.$$.stat" );
		push( @metricsOutput, $countErrors );
	}
=cut

	return ( \@metricsOutput, \@notifynames );
}
#####################################
#####################################
#Compiling BAM metrics
sub CompileMetrics {
	my ( $recalibratedBams, $titleFile, $outdir ) = @_;
	my @bams    = @$recalibratedBams;
	my $allBams = "$outdir/RecalibrationInputBams.list";
	if ( -e "$outdir/$allBams" ) {
		$logger->info(
			"File $allBams exists and the program will overwrite it.");
		open( FH, ">", $allBams )
		  or die(
			$logger->fatal("CompileMetrics:Cannot open $allBams, Error:$!") );
		foreach my $file (@bams) {
			print FH "$outdir/$file\n";
		}
		close(FH);
	}
	else {
		open( FH, ">", $allBams )
		  or die(
			$logger->fatal("CompileMetrics:Cannot open $allBams, Error:$!") );
		foreach my $file (@bams) {
			print FH "$outdir/$file\n";
		}
		close(FH);
	}
	eval {
		if ( $CLUSTER eq "SGE" )
		{
			my $cmd =
"$QSUB -q $queue -V -wd $outdir -N CompileMetrics.$$ -o CompileMetrics.$$.stdout -e CompileMetrics.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1  -b y '$PERL $CompileMetrics -i $allBams -t $titleFile -am $AllMetrics -cn $BestCopyNumber -ncn $NormVsNormCopyNumber -ce $exonIntervalsFile -gcb $GCBiasFile -snlo $StdNormalLoess_TM -snlon $StdNormalLoess_NVN -gia $GeneIntervalAnn -tia $TilingIntervalAnn -ln $LoessNormalization -o $outdir -q $queue -qsub $QSUB -p $PERL -rh $RHOME -rl $RLIBS -ms $GenerateMetricsFilesScript'";
			$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N CompileMetrics.$$ -o CompileMetrics.$$.stdout -e CompileMetrics.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1  -b y "$PERL $CompileMetrics -i $allBams -t $titleFile -am $AllMetrics -cn $BestCopyNumber -ncn $NormVsNormCopyNumber -ce $exonIntervalsFile -gcb $GCBiasFile -snlo $StdNormalLoess_TM -snlon $StdNormalLoess_NVN -gia $GeneIntervalAnn -tia $TilingIntervalAnn -ln $LoessNormalization -o $outdir -q $queue -qsub $QSUB -p $PERL -rh $RHOME -rl $RLIBS -ms $GenerateMetricsFilesScript"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid CompileMetrics.$$ -N NotifyCM.$$ -e NotifyCM.$$.stderr -o NotifyCM.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
		}
		else {
			my $cmd =
"$BSUB -q $queue -cwd $outdir -J CompileMetrics.$$ -o CompileMetrics.$$.%J.stdout -e CompileMetrics.$$.%J.stderr -We 0:59 -R \"rusage[mem=2]\" -R \"rusage[iounits=0]\" -M 4 -n 1 \"$PERL $CompileMetrics -i $allBams -t $titleFile -am $AllMetrics -cn $BestCopyNumber -ncn $NormVsNormCopyNumber -ce $exonIntervalsFile -gcb $GCBiasFile -snlo $StdNormalLoess_TM -snlon $StdNormalLoess_NVN -gia $GeneIntervalAnn -tia $TilingIntervalAnn -ln $LoessNormalization -o $outdir -q $queue -bsub $BSUB -p $PERL -rh $RHOME -rl $RLIBS -ms $GenerateMetricsFilesScript\"";
			$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J CompileMetrics.$$ -o CompileMetrics.$$.%J.stdout -e CompileMetrics.$$.%J.stderr -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$PERL $CompileMetrics -i $allBams -t $titleFile -am $AllMetrics -cn $BestCopyNumber -ncn $NormVsNormCopyNumber -ce $exonIntervalsFile -gcb $GCBiasFile -snlo $StdNormalLoess_TM -snlon $StdNormalLoess_NVN -gia $GeneIntervalAnn -tia $TilingIntervalAnn -ln $LoessNormalization -o $outdir -q $queue -bsub $BSUB -p $PERL -rh $RHOME -rl $RLIBS -ms $GenerateMetricsFilesScript"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(CompileMetrics.$$)" -J NotifyCM.$$ -e NotifyCM.$$.%J.stderr -o NotifyCM.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
		}
	};
	if ($@) {
		$logger->fatal("CompileMetrics:Job Submission Failed, Error:$@");
		exit(1);
	}
	return ("NotifyCM.$$.stat");
}
#####################################
#####################################
#Call Germline mutations and indels
sub RunUnifiedGenotyper {
	my ( $bamFile, $outdir, $id ) = @_;
	my $outFilename1 = $bamFile;
	my $outFilename2 = $bamFile;
	$outFilename1 =~ s/\.bam/_GLM\.vcf/g;
	$outFilename2 =~ s/\.bam/_GLI\.vcf/g;
	if (    ( -e "$outdir/$outFilename1" )
		and ( ( -s "$outdir/$outFilename1" ) != 0 )
		and ( -e "$outdir/$outFilename2" )
		and ( ( -s "$outdir/$outFilename2" ) != 0 ) )
	{
		$logger->info(
"Files:\n$outFilename1\n$outFilename2\n they exists and process will not run to make \"_GLM.vcf\" and \"_GLI.vcf\"file."
		);
		return ( "$outFilename1,$outFilename2", 'NULL' );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N UGM.$id.$$ -o UGM.$id.$$.stdout -e UGM.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y '$JAVA_1_7 -Xmx4g -jar $GATK -T UnifiedGenotyper -R $Reference -I $bamFile -o $outFilename1  -rf BadCigar --dbsnp $dbSNP -glm SNP -dcov 50000 -L $GeneInterval'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N UGM.$id.$$ -o UGM.$id.$$.stdout -e UGM.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y "$JAVA_1_7 -Xmx4g -jar $GATK -T UnifiedGenotyper -R $Reference -I $bamFile -o $outFilename1  -rf BadCigar --dbsnp $dbSNP -glm SNP -dcov 50000 -L $GeneInterval"`;
				$cmd =
"$QSUB -q $queue -V -wd $outdir -N UGI.$id.$$ -o UGI.$id.$$.stdout -e UGI.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y '$JAVA_1_7 -Xmx4g -jar $GATK -T UnifiedGenotyper -R $Reference -I $bamFile -o $outFilename2  -rf BadCigar --dbsnp $dbSNP -glm INDEL -dcov 50000 -L $GeneInterval'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N UGI.$id.$$ -o UGI.$id.$$.stdout -e UGI.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y "$JAVA_1_7 -Xmx4g -jar $GATK -T UnifiedGenotyper -R $Reference -I $bamFile -o $outFilename2  -rf BadCigar --dbsnp $dbSNP -glm INDEL -dcov 50000 -L $GeneInterval"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid "UGM.$id.$$,UGI.$id.$$" -N NotifyUG.$id.$$ -e NotifyUG.$id.$$.stderr -o NotifyUG.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J UGM.$id.$$ -o UGM.$id.$$.%J.stdout -e UGM.$id.$$.%J.stderr -R \"rusage[mem=8]\" -M 12 -n 1 \"$JAVA_1_7 -Xmx4g -jar $GATK -T UnifiedGenotyper -R $Reference -I $bamFile -o $outFilename1  -rf BadCigar --dbsnp $dbSNP -glm SNP -dcov 50000 -L $GeneInterval\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J UGM.$id.$$ -o UGM.$id.$$.%J.stdout -e UGM.$id.$$.%J.stderr -We 24:00 -R "rusage[mem=8]" -M 12 -n 1 "$JAVA_1_7 -Xmx4g -jar $GATK -T UnifiedGenotyper -R $Reference -I $bamFile -o $outFilename1  -rf BadCigar --dbsnp $dbSNP -glm SNP -dcov 50000 -L $GeneInterval"`;
				$cmd =
"$BSUB -q $queue -cwd $outdir -J UGI.$id.$$ -o UGI.$id.$$.%J.stdout -e UGI.$id.$$.%J.stderr -R \"rusage[mem=8]\" -M 12 -n 1 \"$JAVA_1_7 -Xmx4g -jar $GATK -T UnifiedGenotyper -R $Reference -I $bamFile -o $outFilename2  -rf BadCigar --dbsnp $dbSNP -glm INDEL -dcov 50000 -L $GeneInterval\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J UGI.$id.$$ -o UGI.$id.$$.%J.stdout -e UGI.$id.$$.%J.stderr -We 24:00 -R "rusage[mem=8]" -M 12 -n 1 "$JAVA_1_7 -Xmx4g -jar $GATK -T UnifiedGenotyper -R $Reference -I $bamFile -o $outFilename2  -rf BadCigar --dbsnp $dbSNP -glm INDEL -dcov 50000 -L $GeneInterval"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(UGM.$id.$$) && post_done(UGI.$id.$$)" -J NotifyUG.$id.$$ -e NotifyUG.$id.$$.%J.stderr -o NotifyUG.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal("UnifiedGenotyper:Job Submission Failed, Error:$@");
			exit(1);
		}
	}
	return ( "$outFilename1,$outFilename2", "NotifyUG.$id.$$.stat" );
}
#####################################
#####################################
#Call somatic mutations and indels
sub RunMutect_SomaticIndelDetector {
	my ( $pID, $normalBamFile, $tumorBamFile, $pconfig, $outdir, $id ) = @_;

	#MutectFiles
	my $MutationVerboseOutFilename = $tumorBamFile;
	$MutationVerboseOutFilename =~ s/\.bam/\.callstats\.TN\.txt/g;
	my $MutationOutFilename = $tumorBamFile;
	$MutationOutFilename =~ s/\.bam/\.callstats\.TN\.vcf/g;

	#SID files
	my $IndelOutFilename = $tumorBamFile;
	$IndelOutFilename =~ s/\.bam/\.indel\.detailed\.TN\.matched\.vcf/g;
	my $IndelBedOutFilename = $tumorBamFile;
	$IndelBedOutFilename =~ s/\.bam/\.indel\.brief\.TN\.matched\.bed/g;
	my $IndelVerboseOutFilename = $tumorBamFile;
	$IndelVerboseOutFilename =~ s/\.bam/\.indel\.detailed\.TN\.matched\.txt/g;
	my @notifynames = ();
	my ($tFileId) = $tumorBamFile =~ /(.*)_bc\d+_/;
	my $targetFilePrefix = $pID . "_targetsToRealign";
	my $targetFile = $outdir . "/" . $targetFilePrefix . "_covered_srt.list";

	if ( -e $targetFile ) {
		$targetFile = $targetFile;
	}
	else {
		$targetFile = $TargetRegionLIST;
	}

	#bai file linking
	my ($tumorBaiFile) = $tumorBamFile =~ /(.*)\.bam/;
	$tumorBaiFile = $tumorBaiFile . ".bai";
	my $lntumorBaiFile = $tumorBamFile;
	$lntumorBaiFile = $lntumorBaiFile . ".bai";
	`ln -s $tumorBaiFile $lntumorBaiFile`;
	my ($normalBaiFile) = $normalBamFile =~ /(.*)\.bam/;
	$normalBaiFile = $normalBaiFile . ".bai";
	my $lnnormalBaiFile = $normalBamFile;
	$lnnormalBaiFile = $lnnormalBaiFile . ".bai";
	`ln -s $normalBaiFile $lnnormalBaiFile`;

	#Pindel File
	my $PindelVCFOutFilename = $tumorBamFile;
	$PindelVCFOutFilename =~ s/\.bam/\.pindel\.detailed\.TN\.matched\.vcf/g;

	#Mutect
	if (    ( -e "$outdir/$MutationOutFilename" )
		and ( ( -s "$outdir/$MutationOutFilename" ) != 0 ) )
	{
		$logger->info(
"Files:\n$MutationOutFilename they exists and process will not run to make \".callstats.TN.txt\" file."
		);
		push( @notifynames, "NULL" );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N MuTect.$id.$$ -o MuTect.$id.$$.stdout -e MuTect.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y '$JAVA_1_6 -Xmx4g -jar $Mutect -T MuTect --input_file:normal $normalBamFile --input_file:tumor $tumorBamFile --reference_sequence $Reference --dbsnp $dbSNP --cosmic $COSMIC -o $MutationVerboseOutFilename -vcf $MutationOutFilename --enable_extended_output  -dcov 50000 -rf BadCigar'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N MuTect.$id.$$ -o MuTect.$id.$$.stdout -e MuTect.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y "$JAVA_1_6 -Xmx4g -jar $Mutect -T MuTect --input_file:normal $normalBamFile --input_file:tumor $tumorBamFile --reference_sequence $Reference --dbsnp $dbSNP --cosmic $COSMIC -o $MutationVerboseOutFilename -vcf $MutationOutFilename --enable_extended_output  -dcov 50000 -rf BadCigar -rf MappingQuality -mmq $MAPQ"`;
				$cmd =
"$QSUB -q $queue -V -wd $outdir -hold_jid MuTect.$id.$$ -N Triallelic.$id.$$ -o Triallelic.$id.$$.stdout -e Triallelic.$id.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y '$PERL $triallelic_fix --output_dir $outdir --mutect_raw_file $MutationVerboseOutFilename --mutect_vcf_file $MutationOutFilename'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -hold_jid MuTect.$id.$$ -N Triallelic.$id.$$ -o Triallelic.$id.$$.stdout -e Triallelic.$id.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$PERL $triallelic_fix --output_dir $outdir --mutect_raw_file $MutationVerboseOutFilename --mutect_vcf_file $MutationOutFilename"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid MuTect.$id.$$,Triallelic.$id.$$ -N NotifyMT.$id.$$ -e NotifyMT.$id.$$.stderr -o NotifyMT.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J MuTect.$id.$$ -o MuTect.$id.$$.%J.stdout -e MuTect.$id.$$.%J.stderr -We 24:00 -R \"rusage[mem=8]\" -M 12 -n 1 \"$JAVA_1_6 -Xmx4g -jar $Mutect -T MuTect --input_file:normal $normalBamFile --input_file:tumor $tumorBamFile --reference_sequence $Reference --dbsnp $dbSNP --cosmic $COSMIC -o $MutationVerboseOutFilename -vcf $MutationOutFilename --enable_extended_output  -dcov 50000 -rf BadCigar -rf MappingQuality -mmq $MAPQ\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J MuTect.$id.$$ -o MuTect.$id.$$.%J.stdout -e MuTect.$id.$$.%J.stderr -We 24:00 -R "rusage[mem=8]" -M 12 -n 1 "$JAVA_1_6 -Xmx4g -jar $Mutect -T MuTect --input_file:normal $normalBamFile --input_file:tumor $tumorBamFile --reference_sequence $Reference --dbsnp $dbSNP --cosmic $COSMIC -o $MutationVerboseOutFilename -vcf $MutationOutFilename --enable_extended_output  -dcov 50000 -rf BadCigar -rf MappingQuality -mmq $MAPQ"`;
				$cmd =
"$BSUB -q $queue -cwd $outdir -w \"post_done(MuTect.$id.$$)\" -J Triallelic.$id.$$ -o Triallelic.$id.$$.%J.stdout -e Triallelic.$id.$$.%J.stderr -We 0:59 -R \"rusage[mem=2]\" -M 4 -n 1 \"$PERL $triallelic_fix --output_dir $outdir --mutect_raw_file $MutationVerboseOutFilename --mutect_vcf_file $MutationOutFilename\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -w "post_done(MuTect.$id.$$)" -J Triallelic.$id.$$ -o Triallelic.$id.$$.%J.stdout -e Triallelic.$id.$$.%J.stderr -We 0:59 -R "rusage[mem=2]" -M 4 -n 1 "$PERL $triallelic_fix --output_dir $outdir --mutect_raw_file $MutationVerboseOutFilename --mutect_vcf_file $MutationOutFilename"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(MuTect.$id.$$) && post_done(Triallelic.$id.$$)" -J NotifyMT.$id.$$ -e NotifyMT.$id.$$.%J.stderr -o NotifyMT.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal("Mutect:Job Submission Failed, Error:$@");
			exit(1);
		}
		push( @notifynames, "NotifyMT.$id.$$.stat" );
	}

	#Somatic Indel Detector
	if (    ( -e "$outdir/$IndelVerboseOutFilename" )
		and ( ( -s "$outdir/$IndelVerboseOutFilename" ) != 0 ) )
	{
		$logger->info(
"Files:\n$IndelVerboseOutFilename\n exists and process will not run to make \".indel.detailed.TN.matched.txt\" file."
		);
		push( @notifynames, "NULL" );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N SID.$id.$$ -o SID.$id.$$.stdout -e SID.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y '$JAVA_1_6 -Xmx4g -jar $GATK_SomaticIndel -T SomaticIndelDetector -R $Reference -I:normal $normalBamFile -I:tumor $tumorBamFile -filter 'T_COV<10||N_COV<4||T_INDEL_F<0.0001||T_INDEL_CF<0.7' -verbose $IndelVerboseOutFilename -o $IndelOutFilename -refseq $Refseq --maxNumberOfReads 100000 -rf DuplicateRead -rf FailsVendorQualityCheck -rf NotPrimaryAlignment -rf BadMate -rf MappingQualityUnavailable -rf UnmappedRead -rf BadCigar -rf MappingQuality -mmq $MAPQ -L $targetFile'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N SID.$id.$$ -o SID.$id.$$.stdout -e SID.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y "$JAVA_1_6 -Xmx4g -jar $GATK_SomaticIndel -T SomaticIndelDetector -R $Reference -I:normal $normalBamFile -I:tumor $tumorBamFile -filter 'T_COV<10||N_COV<4||T_INDEL_F<0.0001||T_INDEL_CF<0.7' -verbose $IndelVerboseOutFilename -o $IndelOutFilename -refseq $Refseq --maxNumberOfReads 100000 -rf DuplicateRead -rf FailsVendorQualityCheck -rf NotPrimaryAlignment -rf BadMate -rf MappingQualityUnavailable -rf UnmappedRead -rf BadCigar -rf MappingQuality -mmq $MAPQ -L $targetFile"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid SID.$id.$$ -N NotifySID.$id.$$ -e NotifySID.$id.$$.stderr -o NotifySID.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J SID.$id.$$ -o SID.$id.$$.%J.stdout -e SID.$id.$$.%J.stderr -R \"rusage[mem=8]\" -M 12 -n 1 \"$JAVA_1_6 -Xmx4g -jar $GATK_SomaticIndel -T SomaticIndelDetector -R $Reference -I:normal $normalBamFile -I:tumor $tumorBamFile -filter 'T_COV<10||N_COV<4||T_INDEL_F<0.0001||T_INDEL_CF<0.7' -verbose $IndelVerboseOutFilename -o $IndelOutFilename -refseq $Refseq --maxNumberOfReads 100000 -rf DuplicateRead -rf FailsVendorQualityCheck -rf NotPrimaryAlignment -rf BadMate -rf MappingQualityUnavailable -rf UnmappedRead -rf BadCigar -rf MappingQuality -mmq $MAPQ -L $targetFile\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J SID.$id.$$ -o SID.$id.$$.%J.stdout -e SID.$id.$$.%J.stderr -We 0:59 -R "rusage[mem=8]" -M 12 -n 1 "$JAVA_1_6 -Xmx4g -jar $GATK_SomaticIndel -T SomaticIndelDetector -R $Reference -I:normal $normalBamFile -I:tumor $tumorBamFile -filter 'T_COV<10||N_COV<4||T_INDEL_F<0.0001||T_INDEL_CF<0.7' -verbose $IndelVerboseOutFilename -o $IndelOutFilename -refseq $Refseq --maxNumberOfReads 100000 -rf DuplicateRead -rf FailsVendorQualityCheck -rf NotPrimaryAlignment -rf BadMate -rf MappingQualityUnavailable -rf UnmappedRead -rf BadCigar -rf MappingQuality -mmq $MAPQ -L $targetFile"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(SID.$id.$$)" -J NotifySID.$id.$$ -e NotifySID.$id.$$.%J.stderr -o NotifySID.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal(
				"SomaticIndelDetector:Job Submission Failed, Error:$@");
			exit(1);
		}
		push( @notifynames, "NotifySID.$id.$$.stat" );
	}

	#Pindel
	if (    ( -e "$outdir/$PindelVCFOutFilename" )
		and ( ( -s "$outdir/$PindelVCFOutFilename" ) != 0 ) )
	{
		$logger->info(
"Files:\n$PindelVCFOutFilename\n exists and process will not run to make \".pindel.detailed.TN.matched.txt\" file."
		);
		push( @notifynames, "NULL" );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N PID.$id.$$ -o PID.$id.$$.stdout -e PID.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y '$PYTHON $PindelWrapper -i $pconfig -pId $tFileId -t 5 -r $Reference -p $PINDELBIN -chr ALL -q $queue -qsub $QSUB -o $outdir -op $tFileId'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N PID.$id.$$ -o PID.$id.$$.stdout -e PID.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y "$PYTHON $PindelWrapper -i $pconfig -pId $tFileId -t 5 -r $Reference -p $PINDELBIN -chr ALL -q $queue -qsub $QSUB -o $outdir -op $tFileId"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid PID.$id.$$ -N NotifyPID.$id.$$ -e NotifyPID.$id.$$.stderr -o NotifyPID.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J PID.$id.$$ -o PID.$id.$$.%J.stdout -e PID.$id.$$.%J.stderr -We 24:00 -R \"rusage[mem=8]\" -R \"rusage[iounits=0]\" -M 12 -n 1 \"$PYTHON $PindelWrapper -i $pconfig -pId $tFileId -t 5 -r $Reference -p $PINDELBIN -chr ALL -q $queue -bsub $BSUB -o $outdir -op $tFileId\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J PID.$id.$$ -o PID.$id.$$.%J.stdout -e PID.$id.$$.%J.stderr -We 24:00 -R "rusage[mem=8]" -R "rusage[iounits=0]" -M 12 -n 1 "$PYTHON $PindelWrapper -i $pconfig -pId $tFileId -t 5 -r $Reference -p $PINDELBIN -chr ALL -q $queue -bsub $BSUB -o $outdir -op $tFileId"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(PID.$id.$$)" -J NotifyPID.$id.$$ -e NotifyPID.$id.$$.%J.stderr -o NotifyPID.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal("Pindel:Job Submission Failed, Error:$@");
			exit(1);
		}
		push( @notifynames, "NotifyPID.$id.$$.stat" );
	}
	return (
		"$MutationOutFilename",  "$IndelVerboseOutFilename",
		"$PindelVCFOutFilename", \@notifynames
	);
}
#####################################
#####################################
#Filter somatic mutations and indels
sub RunSomaticMutIndelFilter {
	my ( $tumorBamFile, $outdir, $id ) = @_;
	my $MutationVerboseFilename = $tumorBamFile;
	$MutationVerboseFilename =~ s/\.bam/\.callstats\.TN\.txt/g;
	my $MutationVcfFilename = $tumorBamFile;
	$MutationVcfFilename =~ s/\.bam/\.callstats\.TN\.vcf/g;
	my $MutationOutVcfFilename = $tumorBamFile;
	$MutationOutVcfFilename =~ s/\.bam/\.callstats\.TN\_STDfilter\.vcf/g;
	my $MutationOutTxtFilename = $tumorBamFile;
	$MutationOutTxtFilename =~ s/\.bam/\.callstats\.TN\_STDfilter\.txt/g;
	my $IndelVerboseFilename = $tumorBamFile;
	$IndelVerboseFilename =~ s/\.bam/\.indel\.detailed\.TN\.matched\.txt/g;
	my $IndelVcfFilename = $tumorBamFile;
	$IndelVcfFilename =~ s/\.bam/\.indel\.detailed\.TN\.matched\.vcf/g;
	my $IndelOutVcfFilename = $tumorBamFile;
	$IndelOutVcfFilename =~
	  s/\.bam/\.indel\.detailed\.TN\.matched\_STDfilter\.vcf/g;
	my $IndelOutTxtFilename = $tumorBamFile;
	$IndelOutTxtFilename =~
	  s/\.bam/\.indel\.detailed\.TN\.matched\_STDfilter\.txt/g;
	my $PindelVcfFilename = $tumorBamFile;
	$PindelVcfFilename =~ s/\.bam/\.pindel\.detailed\.TN\.matched\.vcf/g;
	my $PindelOutVcfFilename = $tumorBamFile;
	$PindelOutVcfFilename =~
	  s/\.bam/\.pindel\.detailed\.TN\.matched\_STDfilter\.vcf/g;
	my $PindelOutTxtFilename = $tumorBamFile;
	$PindelOutTxtFilename =~
	  s/\.bam/\.pindel\.detailed\.TN\.matched\_STDfilter\.txt/g;
	my ($tFileId) = $tumorBamFile =~ /(.*)_bc\d+_/;
	my @notifynames = ();

	#Filter Somatic Mutation
	if (    ( -e "$outdir/$MutationOutVcfFilename" )
		and ( ( -s "$outdir/$MutationOutVcfFilename" ) != 0 ) )
	{
		$logger->info(
"Files:\n$MutationOutVcfFilename they exists and process will not run to make \"_STDfilter.vcf\" file."
		);
		push( @notifynames, "NULL" );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N FilterSID.$id.$$ -o FilterSID.$id.$$.stdout -e FilterSID.$id.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y '$PERL $filter_SomaticIndel -t $outdir/$IndelVerboseFilename -v $outdir/$IndelVcfFilename -s $tFileId -o $outdir -dp $dp_SomIndelStdFilter -ad $ad_SomIndelStdFilter -vf $vf_SomIndelStdFilter -tnr $TNfreqRatio_SomIndelStdFilter'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N FilterMuTect.$id.$$ -o FilterMuTect.$id.$$.stdout -e FilterMuTect.$id.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$PERL $filter_Mutect -t $outdir/$MutationVerboseFilename -v $outdir/$MutationVcfFilename -s $tFileId -o $outdir -dp $dp_MutectStdFilter -ad $ad_MutectStdFilter -vf $vf_MutectStdFilter -tnr $TNfreqRatio_MutectStdFilter"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid FilterMuTect.$id.$$ -N NotifyFMT.$id.$$ -e NotifyFMT.$id.$$.stderr -o NotifyFMT.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J FilterMuTect.$id.$$ -o FilterMuTect.$id.$$.%J.stdout -e FilterMuTect.$id.$$.%J.stderr -R \"rusage[mem=2]\" -M 4 -n 1 \"$PERL $filter_Mutect -t $outdir/$MutationVerboseFilename -v $outdir/$MutationVcfFilename -s $tFileId -o $outdir -dp $dp_MutectStdFilter -ad $ad_MutectStdFilter -vf $vf_MutectStdFilter -tnr $TNfreqRatio_MutectStdFilter\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J FilterMuTect.$id.$$ -o FilterMuTect.$id.$$.%J.stdout -e FilterMuTect.$id.$$.%J.stderr -We 0:59 -R "rusage[mem=2]" -M 4 -n 1 "$PERL $filter_Mutect -t $outdir/$MutationVerboseFilename -v $outdir/$MutationVcfFilename -s $tFileId -o $outdir -dp $dp_MutectStdFilter -ad $ad_MutectStdFilter -vf $vf_MutectStdFilter -tnr $TNfreqRatio_MutectStdFilter"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(FilterMuTect.$id.$$)" -J NotifyFMT.$id.$$ -e NotifyFMT.$id.$$.%J.stderr -o NotifyFMT.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal("FilerMutect:Job Submission Failed, Error:$@");
			exit(1);
		}
		push( @notifynames, "NotifyFMT.$id.$$.stat" );
	}

	#Filter Somatic Indel
	if (    ( -e "$outdir/$IndelOutVcfFilename" )
		and ( ( -s "$outdir/$IndelOutVcfFilename" ) != 0 ) )
	{
		$logger->info(
"Files:$IndelOutVcfFilename\nthey exists and process will not run to make \"_STDfilter.vcf\" file."
		);
		push( @notifynames, "NULL" );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N FilterSID.$id.$$ -o FilterSID.$id.$$.stdout -e FilterSID.$id.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y '$PERL $filter_SomaticIndel -t $outdir/$IndelVerboseFilename -v $outdir/$IndelVcfFilename -s $tFileId -o $outdir -dp $dp_SomIndelStdFilter -ad $ad_SomIndelStdFilter -vf $vf_SomIndelStdFilter -tnr $TNfreqRatio_SomIndelStdFilter'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N FilterSID.$id.$$ -o FilterSID.$id.$$.stdout -e FilterSID.$id.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$PERL $filter_SomaticIndel -t $outdir/$IndelVerboseFilename -v $outdir/$IndelVcfFilename -s $tFileId -o $outdir -dp $dp_SomIndelStdFilter -ad $ad_SomIndelStdFilter -vf $vf_SomIndelStdFilter -tnr $TNfreqRatio_SomIndelStdFilter"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid FilterSID.$id.$$ -N NotifyFSID.$id.$$ -e NotifyFSID.$id.$$.stderr -o NotifySID.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J FilterSID.$id.$$ -o FilterSID.$id.$$.%J.stdout -e FilterSID.$id.$$.%J.stderr -R \"rusage[mem=2]\" -M 4 -n 1 \"$PERL $filter_SomaticIndel -t $outdir/$IndelVerboseFilename -v $outdir/$IndelVcfFilename -s $tFileId -o $outdir -dp $dp_SomIndelStdFilter -ad $ad_SomIndelStdFilter -vf $vf_SomIndelStdFilter -tnr $TNfreqRatio_SomIndelStdFilter\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J FilterSID.$id.$$ -o FilterSID.$id.$$.%J.stdout -e FilterSID.$id.$$.%J.stderr -We 0:59 -R "rusage[mem=2]" -M 4 -n 1 "$PERL $filter_SomaticIndel -t $outdir/$IndelVerboseFilename -v $outdir/$IndelVcfFilename -s $tFileId -o $outdir -dp $dp_SomIndelStdFilter -ad $ad_SomIndelStdFilter -vf $vf_SomIndelStdFilter -tnr $TNfreqRatio_SomIndelStdFilter"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(FilterSID.$id.$$)" -J NotifyFSID.$id.$$ -e NotifyFSID.$id.$$.%J.stderr -o NotifySID.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal("FilterIndel:Job Submission Failed, Error:$@");
			exit(1);
		}
		push( @notifynames, "NotifySID.$id.$$.stat" );
	}

	#Filter Somatic Pindel
	if (    ( -e "$outdir/$PindelOutVcfFilename" )
		and ( ( -s "$outdir/$PindelOutVcfFilename" ) != 0 ) )
	{
		$logger->info(
"Files:$PindelOutVcfFilename\nthey exists and process will not run to make \"_STDfilter.vcf\" file."
		);
		push( @notifynames, "NULL" );
	}
	else {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N FilterPID.$id.$$ -o FilterPID.$id.$$.stdout -e FilterPID.$id.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y '$PYTHON $filter_Pindel -i $outdir/$PindelVcfFilename -tsn $tFileId -o $outdir -dp $dp_SomIndelStdFilter -ad $ad_SomIndelStdFilter -vf $vf_SomIndelStdFilter -tnr $TNfreqRatio_SomIndelStdFilter -min $Pindel_Min_Indel_Len -max $Pindel_Max_Indel_Len -o $outdir'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N FilterPID.$id.$$ -o FilterPID.$id.$$.stdout -e FilterPID.$id.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$PYTHON $filter_Pindel -i $outdir/$PindelVcfFilename -tsn $tFileId -o $outdir -dp $dp_SomIndelStdFilter -ad $ad_SomIndelStdFilter -vf $vf_SomIndelStdFilter -tnr $TNfreqRatio_SomIndelStdFilter -min $Pindel_Min_Indel_Len -max $Pindel_Max_Indel_Len -o $outdir"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid FilterPID.$id.$$ -N NotifyFPID.$id.$$ -e NotifyFPID.$id.$$.stderr -o NotifyPID.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J FilterPID.$id.$$ -o FilterPID.$id.$$.%J.stdout -e FilterPID.$id.$$.%J.stderr -R \"rusage[mem=2]\" -M 4 -n 1 \"$PYTHON $filter_Pindel -i $outdir/$PindelVcfFilename -tsn $tFileId -o $outdir -dp $dp_SomIndelStdFilter -ad $ad_SomIndelStdFilter -vf $vf_SomIndelStdFilter -tnr $TNfreqRatio_SomIndelStdFilter -min $Pindel_Min_Indel_Len -max $Pindel_Max_Indel_Len -o $outdir\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J FilterPID.$id.$$ -o FilterPID.$id.$$.%J.stdout -e FilterPID.$id.$$.%J.stderr -We 0:59 -R "rusage[mem=2]" -M 4 -n 1 "$PYTHON $filter_Pindel -i $outdir/$PindelVcfFilename -tsn $tFileId -o $outdir -dp $dp_SomIndelStdFilter -ad $ad_SomIndelStdFilter -vf $vf_SomIndelStdFilter -tnr $TNfreqRatio_SomIndelStdFilter -min $Pindel_Min_Indel_Len -max $Pindel_Max_Indel_Len -o $outdir"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(FilterPID.$id.$$)" -J NotifyFPID.$id.$$ -e NotifyFPID.$id.$$.%J.stderr -o NotifyPID.$id.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal("FilterPindel:Job Submission Failed, Error:$@");
			exit(1);
		}
		push( @notifynames, "NotifyPID.$id.$$.stat" );
	}
	return (
		"$MutationOutVcfFilename", "$MutationOutTxtFilename",
		"$IndelOutVcfFilename",    "$IndelOutTxtFilename",
		"$PindelOutVcfFilename",   "$PindelOutTxtFilename",
		\@notifynames
	);
}
#####################################
#####################################
#Calll Structural Variant Pipeline
sub CallSV {
	$logger->info("Starting Structural Variant Pipeline");
	my ($filenames) = @_;
	my @names = ();
	if ($filenames) { (@names) = @$filenames; }
	if ( ( scalar @names == 0 ) and ($fof) ) {
		@names = &GetNames( $fof, $outdir );
	}
	my (@sortedparseFilenames) =
	  sort { lowestNumber($a) <=> lowestNumber($b) } @names;
	@names = @sortedparseFilenames;
	my @notifynames = ();
	eval {
		if ( $CLUSTER eq "SGE" )
		{
			my $cmd =
"$QSUB -q $queue -V -wd $outdir -N RunSV.$$ -o RunSV.$$.stdout -e RunSV.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y '$PERL $SVpipeline -c $sv_config_file -d $datadir -o $outdir'";
			$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N RunSV.$$ -o RunSV.$$.stdout -e RunSV.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$PERL $SVpipeline -c $sv_config_file -d $datadir -o $outdir"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid RunSV.$$ -N NotifyRunSV.$$ -e NotifyRunSV.$$.stderr -o NotifyRunSV.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
		}
		else {
			my $cmd =
"$BSUB -q $queue -cwd $outdir -J RunSV.$$ -o RunSV.$$.%J.stdout -e RunSV.$$.%J.stderr -R \"rusage[mem=2]\" -R \"rusage[iounits=0]\" -M 4 -n 1 \"$PERL $SVpipeline -c $sv_config_file -d $datadir -o $outdir\"";
			$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J RunSV.$$ -o RunSV.$$.%J.stdout -e RunSV.$$.%J.stderr -We 24:00 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$PERL $SVpipeline -c $sv_config_file -d $datadir -o $outdir"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(RunSV.$$)" -J NotifyRunSV.$$ -e NotifyRunSV.$$.%J.stderr -o NotifyRunSV.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
		}
	};
	if ($@) {
		$logger->fatal("RunSV:Job Submission Failed, Error:$@");
		exit(1);
	}
	push( @notifynames, "NotifyRunSV.$$.stat" );
	if ($fof) {
		&WaitToFinish( $outdir, @notifynames );
	}
	return ( \@names );
}
#####################################
#####################################
#Calll Structural Variant Pipeline
sub RunABRA {
	$logger->info("Running ABRA realignment");
	my ( $inputFiles, $patientID, $outdir, $i ) = @_;
	my @notifynames      = ();
	my @processed        = ();
	my @unprocessed      = ();
	my $BWA_PATH         = dirname($BWA);
	my $targetFilePrefix = $patientID . "_targetsToRealign";
	my $targetFile = $outdir . "/" . $targetFilePrefix . "_covered_srt.bed";
	open( FH, $inputFiles )
	  or die( $logger->fatal("RunABRA:Cannot open $inputFiles, Error:$!") );

	while (<FH>) {
		chomp($_);
		my $outFilename = $_;
		$outFilename =~ s/\.bam/_IR\.bam/;
		if ( ( -e "$outFilename" ) and ( ( -s "$outFilename" ) != 0 ) ) {
			push( @processed, $outFilename );
		}
		else {
			push( @unprocessed, $outFilename );
		}
	}
	close(FH);
	if (@unprocessed) {
		eval {
			if ( $CLUSTER eq "SGE" )
			{
				my $cmd =
"$QSUB -q $queue -V -wd $outdir -N RunFindCoveredInterval_IR.$i.$$ -o RunFindCoveredInterval_IR.$i.$$.stdout -e RunFindCoveredInterval_IR.$i.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y '$PYTHON $FindCoveredIntervals -i $inputFiles -of $targetFilePrefix -t 4 -r $Reference -g $GATK -j $JAVA_1_7 -q $queue -qsub $QSUB -o $outdir'";
				$logger->debug("COMMAND: $cmd");
				$cmd =
"$QSUB -q $queue -V -wd $outdir -N RunABRA_IR.$i.$$ -o RunABRA_IR.$i.$$.stdout -e RunABRA.$i.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y '$PYTHON $ABRAwrapper -i $inputFiles -p $patientID -t 12 -k 43 53 63 83 93 -r $Reference -tr $targetFile -j $JAVA_1_7 -b $BWA_PATH -a $ABRA -q $queue -qsub $QSUB -o $outdir -temp $TMPDIR'";
				$logger->debug("COMMAND: $cmd");
`$QSUB -q $queue -V -wd $outdir -N RunFindCoveredInterval_IR.$i.$$ -o RunFindCoveredInterval_IR.$i.$$.stdout -e RunFindCoveredInterval.$i.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y "$PYTHON $FindCoveredIntervals -i $inputFiles -of $targetFilePrefix -t 4 -r $Reference -g $GATK -j $JAVA_1_7 -q $queue -qsub $QSUB -o $outdir"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid RunFindCoveredInterval_IR.$i.$$ -N RunABRA_IR.$i.$$ -o RunABRA_IR.$i.$$.stdout -e RunABRA_IR.$i.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y "$PYTHON $ABRAwrapper -i $inputFiles -p $patientID -t 12 -k 43 53 63 83 93 -r $Reference -tr $targetFile -j $JAVA_1_7 -b $BWA_PATH -a $ABRA -q $queue -qsub $QSUB -o $outdir -temp $TMPDIR"`;
`$QSUB -q $queue -V -wd $outdir -hold_jid RunABRA_IR.$i.$$ -N NotifyRunABRA_IR.$i.$$ -e NotifyRunABRA_IR.$i.$$.stderr -o NotifyRunABRA_IR.$i.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
			}
			else {
				my $cmd =
"$BSUB -q $queue -cwd $outdir -J RunFindCoveredInterval_IR.$i.$$ -o RunFindCoveredInterval_IR.$i.$$.%J.stdout -e RunFindCoveredInterval.$i.$$.%J.stderr -R \"rusage[mem=8]\" -R \"rusage[iounits=0]\" -M 12 -n 1 \"$PYTHON $FindCoveredIntervals -i $inputFiles -of $targetFilePrefix -t 4 -r $Reference -g $GATK -j $JAVA_1_7 -q $queue -bsub $BSUB -o $outdir\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -J RunFindCoveredInterval_IR.$i.$$ -o RunFindCoveredInterval_IR.$i.$$.%J.stdout -e RunFindCoveredInterval.$i.$$.%J.stderr -We 24:00 -R "rusage[mem=8]" -R "rusage[iounits=0]" -M 12 -n 1 "$PYTHON $FindCoveredIntervals -i $inputFiles -of $targetFilePrefix -t 4 -r $Reference -g $GATK -j $JAVA_1_7 -q $queue -bsub $BSUB -o $outdir"`;
				$cmd =
"$BSUB -q $queue -cwd $outdir -w \"post_done(RunFindCoveredInterval_IR.$i.$$)\" -J RunABRA_IR.$i.$$ -o RunABRA_IR.$i.$$.%J.stdout -e RunABRA_IR.$i.$$.%J.stderr -R \"rusage[mem=8]\" -R \"rusage[iounits=0]\" -M 12 -n 1 \"$PYTHON $ABRAwrapper -i $inputFiles -p $patientID -t 12 -k 43 53 63 83 93 -r $Reference -tr $targetFile -j $JAVA_1_7 -b $BWA_PATH -a $ABRA -q $queue -bsub $BSUB -o $outdir -temp $TMPDIR\"";
				$logger->debug("COMMAND: $cmd");
`$BSUB -q $queue -cwd $outdir -w "post_done(RunFindCoveredInterval_IR.$i.$$)" -J RunABRA_IR.$i.$$ -o RunABRA_IR.$i.$$.%J.stdout -e RunABRA_IR.$i.$$.%J.stderr -We 24:00 -R "rusage[mem=8]" -R "rusage[iounits=0]" -M 12 -n 1 "$PYTHON $ABRAwrapper -i $inputFiles -p $patientID -t 12 -k 43 53 63 83 93 -r $Reference -tr $targetFile -j $JAVA_1_7 -b $BWA_PATH -a $ABRA -q $queue -bsub $BSUB -o $outdir -temp $TMPDIR"`;
`$BSUB -q $queue -cwd $outdir -w "post_done(RunABRA_IR.$i.$$)" -J NotifyRunABRA_IR.$i.$$ -e NotifyRunABRA_IR.$i.$$.%J.stderr -o NotifyRunABRA_IR.$i.$$.stat -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -M 4 -n 1 "$outdir/Notify.csh"`;
			}
		};
		if ($@) {
			$logger->fatal("RunABRA:Job Submission Failed, Error:$@");
			exit(1);
		}
		return ("NotifyRunABRA_IR.$i.$$.stat");
	}
	else {
		$logger->info(
"Files:@processed\nthey exists and process will not run to make \"_IR.bam\" file."
		);
		return ("NULL");
	}
}
#####################################
#####################################
#Delete Intermediate files
sub DeleteFiles {
	my ($filesToDelete) = @_;
	foreach my $file (@$filesToDelete) {
		chomp($file);
		if ( -e $file ) {
			eval { `rm $file`; };
			if ($@) {
				$logger->warn(
					"DeleteFile:Sorry not able to delete $file, Error:$@");
			}
		}
		else {
			next;
		}
	}
	return;
}
#####################################
#####################################
#Make Folder and Move files
sub RunHouseKeeping {
	my ( $outdir, $pool ) = @_;
	my (@poolNames)      = @$pool;
	my $poolName         = $poolNames[0];
	my $FinalBams        = "$outdir/FinalBams";
	my $Results          = "$outdir/Results";
	my $CompileMetrics   = "$outdir/Results/CompileMetrics";
	my $AllSampleResults = "$outdir/AllSampleResults";

	#my $vcfs = "$outdir/VCFs";
	my $mutations   = "$outdir/Mutations";
	my $AlleleDepth = "$outdir/AlleleDepth";
	my $StdLogFiles = "$outdir/StdLogFiles";
	eval {
		if ( -d $FinalBams )
		{
			$logger->info("DIR:$FinalBams exits and won\'t be created.");
		}
		else {
			`mkdir $FinalBams`;
		}
		if ( -d $Results ) {
			$logger->info("DIR:$Results exits and won\'t be created.");
		}
		else {
			`mkdir $Results`;
		}
		if ( -d $CompileMetrics ) {
			$logger->info("DIR:$CompileMetrics exits and won\'t be created.");
		}
		else {
			`mkdir $CompileMetrics`;
		}
		if ( -d $AlleleDepth ) {
			$logger->info("DIR:$AlleleDepth exits and won\'t be created.");
		}
		else {
			`mkdir $AlleleDepth`;
		}
		if ( -d $AllSampleResults ) {
			$logger->info("DIR:$AllSampleResults exits and won\'t be created.");
		}
		else {
			`mkdir $AllSampleResults`;
		}

		#if(-d $vcfs)
		#{
		#    $logger->info("DIR:$vcfs exits and won\'t be created.");
		#}
		#else
		#{
		#    `mkdir $vcfs`;
		#}
		if ( -d $mutations ) {
			$logger->info("DIR:$mutations exits and won\'t be created.");
		}
		else {
			`mkdir $mutations`;
		}
		if ( -d $StdLogFiles ) {
			$logger->info("DIR:$StdLogFiles exits and won\'t be created.");
		}
		else {
			`mkdir $StdLogFiles`;
		}
	};
	if ($@) {
		$logger->warn("MoveFiles:Error while making folders. Error:$@");
	}
	eval {
		my $mvTitleFiles        = $poolName . "_title.txt";
		my @mvTitleData         = <$outdir/$mvTitleFiles>;
		my $mvCMPFiles          = $poolName . "_ALL_*";
		my $mvCMPFiles_All      = $poolName . "_All_*";
		my $mvCopyNumber        = $poolName . "_copynumber_*";
		my $mvNormVsNorm        = $poolName . "_normvsnorm_*";
		my $mvLoessNorm         = $poolName . "_loessnorm*";
		my $BAMFiles            = "*_BR.ba*";
		my @BAMData             = <$outdir/$BAMFiles>;
		my $BedFiles            = "*_BR.bed";
		my @BedData             = <$outdir/$BedFiles>;
		my $mutectFiles         = "*_BR.callstats.TN*";
		my @mutectData          = <$outdir/$mutectFiles>;
		my $indelFiles          = "*_BR.*indel.detailed.TN.matched*";
		my @indelData           = <$outdir/$indelFiles>;
		my $MetricsResultsFiles = "*_BR.*";
		my @MetricsResultsData  = <$outdir/$MetricsResultsFiles>;
		my $ClippingResults     = "*_cl.stats";
		my @ClippingData        = <$outdir/$ClippingResults>;
		my $DuplicationResults  = "*_MD.metrics";
		my @DuplicationData     = <$outdir/$DuplicationResults>;
		my $SegmentationResults = "*_copynumber.seg";
		my @SegmentationData    = <$outdir/$SegmentationResults>;
		my $cnvJsons            = "*_copynumber.json";
		my @cnvJson_arr         = <$outdir/$cnvJsons>;
		my $ADfiles             = "*_mpileup.alleledepth";
		my @ADData              = <$outdir/$ADfiles>;

		#my $gMfiles = "*_GLM.*";
		#my @gMData = <$outdir/$gMfiles>;
		#my $gIfiles = "*_GLI.*";
		#my @gIData = <$outdir/$gIfiles>;
		my $stdoutFiles = "*.stdout";
		my @stdoutData  = <$outdir/$stdoutFiles>;
		my $stderrFiles = "*.stderr";
		my @stderrData  = <$outdir/$stderrFiles>;
		my $QCfigures   = "figure";
		if ( scalar(@BAMData) > 0 ) {
			`mv $outdir/$BAMFiles $FinalBams`;
		}
		if ( scalar(@BedData) > 0 ) {
			`mv $outdir/$BedFiles $FinalBams`;
		}
		if ( scalar(@mutectData) > 0 ) {
			`mv $outdir/$mutectFiles $mutations`;
		}
		if ( scalar(@indelData) > 0 ) {
			`mv $outdir/$indelFiles $mutations`;
		}

		#if(scalar(@gMData) > 0)
		#{
		#    `mv $outdir/$gMfiles $vcfs`;
		#}
		#if(scalar(@gIData) > 0)
		#{
		#    `mv $outdir/$gIfiles $vcfs`;
		#}
		if ( -d $QCfigures ) {
			`mv $outdir/$QCfigures/ Results/`;
		}
		if ( scalar(@mvTitleData) > 0 ) {
			`mv $outdir/$mvTitleFiles $Results`;
		}
		my @mvCMPData = <$outdir/$mvCMPFiles>;
		if ( scalar(@mvCMPData) > 0 ) {
			`mv $outdir/$mvCMPFiles $Results`;
		}
		my @mvCMPData_All = <$outdir/$mvCMPFiles_All>;
		if ( scalar(@mvCMPData_All) > 0 ) {
			`mv $outdir/$mvCMPFiles_All $Results`;
		}
		my @mvCopyNumberFiles = glob("$outdir/$mvCopyNumber");
		if ( scalar(@mvCopyNumberFiles) > 0 ) {
			`mv $outdir/$mvCopyNumber $Results`;
			my $src1 =
			  $Results . "/" . $poolName . "_copynumber_segclusp.probes.txt";
			my $dest1 = $Results . "/cnv_variants_probe_level.txt";
			my $src2 =
			  $Results . "/" . $poolName . "_copynumber_segclusp.genes.txt";
			my $dest2 = $Results . "/cnv_variants_gene_level.txt";
			`ln -s $src1 $dest1`;
			`ln -s $src2 $dest2`;
		}
		my @mvNormVsNormFiles = <$outdir/$mvNormVsNorm>;
		if ( scalar(@mvNormVsNormFiles) > 0 ) {
			`mv $outdir/$mvNormVsNorm $Results`;
		}
		my @mvLoessNormFiles = <$outdir/$mvLoessNorm>;
		if ( scalar(@mvLoessNormFiles) > 0 ) {
			`mv $outdir/$mvLoessNorm $Results`;
		}
		if ( scalar(@MetricsResultsData) > 0 ) {
			`mv $outdir/$MetricsResultsFiles $AllSampleResults`;
		}
		if ( scalar(@ClippingData) > 0 ) {
			`mv $outdir/$ClippingResults $AllSampleResults`;
		}
		if ( scalar(@DuplicationData) > 0 ) {
			`mv $outdir/$DuplicationResults $AllSampleResults`;
		}
		if ( scalar(@SegmentationData) > 0 ) {
			`mv $outdir/$SegmentationResults $AllSampleResults`;
		}
		if ( scalar(@cnvJson_arr) > 0 ) {
			`mv $outdir/$cnvJsons $AllSampleResults`;
		}
		if ( scalar(@ADData) > 0 ) {
			`mv $outdir/$ADfiles $AlleleDepth`;
		}
		if ( scalar(@stdoutData) > 0 ) {
			`mv $outdir/$stdoutFiles $StdLogFiles`;
		}
		if ( scalar(@stderrData) > 0 ) {
			`mv $outdir/$stderrFiles $StdLogFiles`;
		}
		`mv $outdir/SampleSheet.csv $Results`;
	};
	if ($@) {
		$logger->warn("MoveFiles:Error while moving files. Error:$@");
	}
	eval {

		#`rm $outdir/*.csh`;
	};
	if ($@) {
		$logger->warn("MoveFiles:Error while removing .csh files. Error:$@");
	}
	return;
}
