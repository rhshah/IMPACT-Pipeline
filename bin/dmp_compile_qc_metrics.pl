#!/usr/bin/perl
##########ComplieMetrics.pl########
#Author: Ronak Shah(RS),Donvan Cheng(DC)
#Date: 12/11/2012
#LastModified: 07/24/2013
#Version:1.3
#Description: Compile all the metrics files.
##04/24/2013
#v1.2
#RS:Fixed bugs in names
##05/14/2013
#v1.2
#RS:Fixed bugs for version 5 baits
##07/24/2013
#v1.3
#DC:Updated the copynumber scripts
##08/08/2013
#v1.3
#RS:Added option of SGE queue
##08/10/2013
#v1.3
#AZ: Added logger options
##08/27/2013
#RS:Changed chosing of the normal for contamination check.
##09/10/2013
# AZ: Added option for canonical exon coverage file
#############################
##08/27/2015
# RS: Added option for using bsub
#############################

use strict;
use warnings;
use Getopt::Long;
use IO::File;
use Tie::Autotie 'Tie::IxHash';
use Cwd;
use Statistics::Lite qw(:all);
use MSKCC_DMP_Logger;
my $logger = MSKCC_DMP_Logger->get_logger('COMPILE_QC_METRICS');
$logger->start_local();
#--This variable holds the current time
my $now = time;
my ($bamList,$titleFile,$allMetrics,$loessNorm,$QSUB,$BSUB,$PERL,$RHOME,$RLIBS,$bestCN,$nvnCN,$outdir,$gcbiasfile,$histnormdir, $stdnormloess_tm, $stdnormloess_nvn, $queue, $exonIntervals, $MetricsScript);
my ($geneIntervalAnn,$tilingIntervalAnn);

if (@ARGV < 1 or !GetOptions (
	    'bamFileList|i:s'   	=> \$bamList,
	    'titleFile|t:s'     	=> \$titleFile,
            'AllMetrics|am:s'	=> \$allMetrics,
	    'LoessNorm|ln:s'    	=> \$loessNorm,
	    'BestCN|cn:s'       	=> \$bestCN,
	    'NormVsNormCN|ncn:s'  	=> \$nvnCN,
	    'PERL|p:s'       		=> \$PERL,
	    'RHOME|rh:s'       		=> \$RHOME,
	    'RLIBS|rl:s'  			=> \$RLIBS,
            'GCBias|gcb:s'      => \$gcbiasfile,
            'StdNormLoess_TM|snlo:s'      => \$stdnormloess_tm,
            'StdNormLoess_NVN|snlon:s'      => \$stdnormloess_nvn,
	    'GeneIntervalAnn|gia:s' => \$geneIntervalAnn,
	    'TilingIntervalAnn|tia:s' => \$tilingIntervalAnn,
	    'HistNorm|his:s'    => \$histnormdir,
	    'qsub:s'						=> \$QSUB,
	    'bsub:s'						=> \$BSUB,
	    'queue|q:s'         => \$queue,
	    'canonicalExons|ce:s' => \$exonIntervals,
	    'metricsScript|ms:s'   => \$MetricsScript,
	    'outdir|o:s'        => \$outdir))
	{
		Usage();
	}


if((!defined $bamList) or (!defined $titleFile))
{
    $logger->warn("Please provide the file of file containing the name of Bam's as well as the Title File.See Usage\n");
    Usage();
    exit;
}
if (!defined $outdir)
{
    $logger->warn("Output directory is not specified, current working directory will be used as output directory.\n");
    $outdir = getcwd;
}
if (!defined $allMetrics)
{
    $logger->warn("All metrics script is not specified, default (v1.2) will be used instead.\n");
    $allMetrics = "/home/shahr2/Scripts/All/AllMetrics_v1.2.r";
}
if (!defined $loessNorm)
{
    $logger->warn("Loess normalization script is not supplied, default (/home/chengd1/Analysis/LabPipeline/CopyNumberWork/LoessNormalizeExonData_d1.3.R) will be used.\n");
    $loessNorm = "/home/chengd1/Analysis/LabPipeline/CopyNumberWork/LoessNormalizeExonData_d1.3.R";
}
if (!defined $bestCN)
{
    $logger->warn("BestCN script not supplied, default (/home/chengd1/Analysis/LabPipeline/CopyNumberWork/CopyNumber_d2.1.no_hist.R) will be used.\n");
    $bestCN = "/home/chengd1/Analysis/LabPipeline/CopyNumberWork/CopyNumber_d2.1.no_hist.R";
}
if (!defined $nvnCN)
{
    $logger->warn("NormVsNormCN script not supplied, default (/home/chengd1/Analysis/LabPipeline/CopyNumberWork/CopyNumber_d2.1.normal_vs_normal.R) will be used.\n");
    $nvnCN = "/home/chengd1/Analysis/LabPipeline/CopyNumberWork/CopyNumber_d2.1.normal_vs_normal.R";
}
if (!defined $PERL)
{
    $logger->warn("PERL not supplied, default(/usr/bin/perl) will be used.\n");
    $PERL = "/usr/bin/perl";
}
if (!defined $RHOME)
{
    $logger->warn("RHOME not supplied, default(/ifs/e63data/bergerm1/Resources/SupportTools/R/Source/R-3.0.1/bin) will be used.\n");
    $RHOME = "/ifs/e63data/bergerm1/Resources/SupportTools/R/Source/R-3.0.1/bin";
}
if (!defined $RLIBS)
{
    $logger->warn("RHOME not supplied, default(/home/shahr2/R/library) will be used.\n");
    $RLIBS = "/home/shahr2/R/library";
}
if (!defined $gcbiasfile)
{
    $logger->warn("GC bias file is not supplied, default(/home/chengd1/Analysis/LabPipeline/All-2013Apr16th/v4_hg19_all_GC200bp.txt) will be used.\n");
    $gcbiasfile = "/home/chengd1/Analysis/LabPipeline/All-2013Apr16th/v4_hg19_all_GC200bp.txt";
}
if (!defined $stdnormloess_tm)
{
    $logger->warn("Std Norm Loess file TM is not supplied, default(/dmp/data/pubdata/std-normal-loess/production/FFPE_CtrlPool_BAM_ALL_intervalcoverage_loess.txt) will be used.\n");
    $stdnormloess_tm = "/dmp/data/pubdata/std-normal-loess/production/FFPE_CtrlPool_BAM_ALL_intervalcoverage_loess.txt";
}
if (!defined $stdnormloess_nvn)
{
    $logger->warn("Std Norm Loess file NVN is not supplied, default(/dmp/data/pubdata/std-normal-loess/production/CombinedStdNormals_ALL_intervalcoverage_loess.txt) will be used.\n");
    $stdnormloess_nvn = "/dmp/data/pubdata/std-normal-loess/production/CombinedStdNormals_ALL_intervalcoverage_loess.txt";
}
if (!defined $geneIntervalAnn)
{
    $logger->warn("Gene intervals annotated file is not supplied, default(/dmp/data/mskdata/interval-lists/production/gene_intervals.list.annotated) will be used.\n");
    $geneIntervalAnn = "/dmp/data/mskdata/interval-lists/production/gene_intervals.list.annotated";
}
if (!defined $tilingIntervalAnn)
{
    $logger->warn("Tiling intervals annotated file is not supplied, default(/dmp/data/mskdata/interval-lists/production/tiling_intervals.list.annotated) will be used.\n");
    $tilingIntervalAnn = "/dmp/data/mskdata/interval-lists/production/tiling_intervals.list.annotated";
}
if (!defined $histnormdir)
{
    $logger->warn("Historical normal data is not provided, default (/home/chengd1/Analysis/LabPipeline/CopyNumberWork/GoodNormals) will be used.\n");
    $histnormdir = "/home/chengd1/Analysis/LabPipeline/CopyNumberWork/GoodNormals";
}
if(!defined $exonIntervals){
    $logger->warn("Canonical exon interval file is not provided, default (//home/zehira/PubData/cv1.genelist.b37.with_aa.interval_list) will be used.");
    $exonIntervals = "/home/zehira/PubData/cv1.genelist.b37.with_aa.interval_list";
}
#Check  for qsub
if(!defined $QSUB)
{
    $logger->warn("Path to QSUB command not given\n");
    $logger->warn("Path to BSUB command not given. Will assume we are running BSUB");
}
else
{
    $logger->info( "SGE QSUB:$QSUB\n");
}
#Check  for qsub
if(!defined $BSUB)
{
    $logger->warn("Path to BSUB command not given. Will assume we are running QSUB");
}
else
{
    $logger->info( "LSF BSUB:$BSUB\n");
}
if(defined $BSUB && defined $QSUB){
	$logger->fatal("Please provide path either for BSUB or QSUB not of both. See Usage");
	Usage();
	exit(1);
	
}
#Check  for queue
if(!defined $queue)
{
    $logger->warn("Name of the SGE/LSF queue not given default will be used\n");
    $queue = "all.q";
}
else
{
     $logger->info("SGE/LSF Queue:$queue will be used\n");
}

#Read Title File
my($barcode,$pool,$titleSampleId,$collabId,$patientId,$class,$sampleType,$inputNg,$libraryYeild,$poolInput,$baitVersion) = &ReadTitleFile($titleFile,$outdir);
tie (my %classPerBarcode, 'Tie::IxHash');
for(my $i = 0; $i < scalar(@$barcode); $i++)
{
    #print "$$barcode[$i] => $$class[$i]\n";
    $classPerBarcode{$$barcode[$i]} = $$class[$i];
}
&RunCompileMetrics($bamList,$patientId,$outdir);

#--Calculate total runtime (current time minus start time) 
$now = time - $now;

#--Print runtime 
$logger->info(printf("\n\nTotal running time: %02d:%02d:%02d\n\n", int($now / 3600), int(($now % 3600) / 60), int($now % 60)));

$logger->info("!!!!Done, Thanks for using the script!!!\n");
exit;

#####################################
#####################################
#How to use the script.
sub Usage
{
    print "Unknow option: @_\n" if (@_);

	print "\nUsage : RunIlluminaProcess.pl [options]
        [--bamList|i       S File of files having list of all bam files (required)]
        [--titleFile|t     S tab-delimited title file for the samples (required and submit with full path)]
        [--AllMetrics|am   S Path to AllMetrics script (optional;default:/home/shahr2/Scripts/All/AllMetrics_v1.2.r)]
        [--LoessNorm|ln    S Path to Loess Normalization Script  (optional;default:/home/chengd1/Analysis/LabPipeline/CopyNumberWork/LoessNormalizeExonData_d1.3.R)]
        [--BestCN|cn       S Path to Best Copy Number Script (optional;default:/home/chengd1/Analysis/LabPipeline/CopyNumberWork/CopyNumber_d1.2.R)] 
        [--GCBias|gcb      S Path to GC bias file (optional;default:/home/chengd1/Analysis/LabPipeline/All-2013Apr16th/v4_hg19_all_GC200bp.txt)]
        [--HistNorm|his    S Path to Directory with all historical normal files (/home/chengd1/Analysis/LabPipeline/CopyNumberWork/GoodNormals)]
        [--queue|q         S Name of the Sun Grd Engine Queue where the pipeline needs to run (default:all.q,optional)]
        [--qsub            S Path to qsub executable for SGE(default:None,optional)]
        [--bsub            S Path to bsub executable for LSF(default:None,required)]
        [--metricsScript|ms         S Name of the script used to generate .html and .pdf files]
        [--outdir|o        S Path where all the output files will be written (optional) [default:cwd]]
	\n";

	print "For any question please contact Ronak Shah (shahr2\@mskcc.org)\n";
	print "!~!~! Comments, flames and suggestion are welcome !~!~!\n";
	exit;
}
###################################################
###################################################
#--Make Notification file
sub MakeCSH
 {
        my($outdir) = @_;
        my $filename = $outdir . "/Notify.csh";
        my $ntmp = new IO::File(">$filename");
        print $ntmp "#!/bin/csh\n";
        print $ntmp "#Notification File\n";
        print $ntmp "echo"," This is Done","\n";
        $ntmp->close();
        eval{`chmod +x $filename`;};
	if($@){$logger->fatal("Cannot change permissions for $filename. Error:$@");exit(1);}
        return;
}
###################################################
###################################################
#--Waiting for the process to finish

sub WaitToFinish
{
	my($outdir,@waitfilenames) = @_;
	$logger->info("Waiting for the Process to finish...\n");
	foreach my $wfile(@waitfilenames)
	{
	    next if($wfile eq "NULL");
	    wait while(! -e "$outdir/$wfile");
	    #print "$outdir/$wfile\n";
	    while(-e "$outdir/$wfile")
	    {
		#print "$outdir/$wfile\n";
		open(FH,"<","$outdir/$wfile");
		while(<FH>)
		{
		    if($_ =~ /This is Done/ig)
		    {
			#print "\nFinished: $wfile\n";
			last;
		    }
		    else
		    {
			wait;
		    }
		}
		last;
		}
	    close(FH);
	}
	foreach my $wfile(@waitfilenames)
	{
	   next if($wfile eq "NULL");
	   eval{`rm $outdir/$wfile`;};
	   if($@){$logger->warn("Cannot remove $outdir/$wfile. Error:$@");}
	}
	return;
}
#####################################
#####################################
#Read data related to samples as well as barcodes from title file.

sub ReadTitleFile
{
    my($titleFile,$outdir) = @_;
    my @barcode = ();
    my @pool = ();
    my @sampleId = ();
    my @collabId = ();
    my @patientId = ();
    my @class = ();
    my @sampleType = ();
    my @inputNg = ();
    my @libraryYeild = ();
    my @poolInput = ();
    my @baitVersion = ();
    my @fof = ();
    my @newfof = ();

    open(TFH,$titleFile)||die"Cannot open file $titleFile, $!\n";
    while(<TFH>)
    {
	next if($. == 1);
	my @dataCols = split("\t",$_);
	my @newDatacols = grep(s/\s*$//g, @dataCols);#remove whitespace if any
	
	push(@barcode,$newDatacols[0]);
	push(@pool,$newDatacols[1]);
	push(@sampleId,$newDatacols[2]);
	push(@collabId,$newDatacols[3]);
	push(@patientId,$newDatacols[4]);
	push(@class,$newDatacols[5]);
	push(@sampleType,$newDatacols[6]);
	push(@inputNg,$newDatacols[7]);
	push(@libraryYeild,$newDatacols[8]);
	push(@poolInput,$newDatacols[9]);
	push(@baitVersion,$newDatacols[10]);
    
    }
    close(TFH);
    my $poolName = $pool[0];
    my $newtitleFileName = $poolName . "_title.txt";
    if(! -e "$outdir/$newtitleFileName")
    {
	eval{`cp $titleFile $outdir/$newtitleFileName`;};
	if($@){$logger->warn("Cannot copy $outdir/$newtitleFileName. Error:$@");}
    }
    return(\@barcode,\@pool,\@sampleId,\@collabId,\@patientId,\@class,\@sampleType,\@inputNg,\@libraryYeild,\@poolInput,\@baitVersion);

}
#####################################
#####################################
#Run Compile Metrics Calculations
sub RunCompileMetrics
{
    my($bamList,$patientId,$outdir) = @_;
    my @bams = ();
    ##getting the list of files ready.
    open(FH,$bamList) or die $logger->fatal("Cannot open $bamList. Error: $!\n");
    while(<FH>)
    {
	my $file;
	if($_ =~ /\//)
	{
	    $file = pop @{[split("/",$_)]};
	    push(@bams,$file);
	}
	else
 	{
	    $file = $_;
	    push(@bams,$file);
	}
    }
    close(FH);
    my(@sortedparseFilenames) = sort {lowestNumber($a) <=>  lowestNumber($b)} @bams;
    @bams = @sortedparseFilenames;
    my ($HSout) = &RunHSandDuplicationMetrics(\@bams,$outdir);
    my $INSout = &RunInsertSizeMetrics(\@bams,$outdir);
    my ($BQRout,$BQOout) = &RunBaseQualityMetrics(\@bams,$outdir);
    my $ECout = &RunExonCoverageMetrics(\@bams,$outdir);
    my $ECout_nomapq = &RunExonNomapqCoverageMetrics(\@bams,$outdir);
    my $CECout = &RunCanonicalExonCoverageMetrics(\@bams,$exonIntervals,$outdir);
    my $GCout = &RunGeneCoverageMetrics(\@bams,$outdir);
    my ($FPG1out,$FPG2out,$FPG3out,$FPGC1out,$FPGC2out) = &RunFPGenotypeMetrics(\@bams,$patientId,$outdir);
    my $FPTout = &RunFPtilingCoverageMetrics(\@bams,$outdir);
    my $FPTout_nomapq = &RunFPtilingNomapqCoverageMetrics(\@bams,$outdir);
    my ($HST1out,$HST2out) = &RunTargetCoverageMetrics(\@bams,$outdir);

    # DEBUG: For testing
    &PlotGraphs($FPG1out,$ECout,$outdir);
    return;
}
#####################################
#####################################
#sort by barcode name from bams:

sub lowestNumber
{
    my $files = shift;
    my $number;
    my @filenames;
    if($files =~ /,/)
    {
	@filenames = split(",",$files);
	($number) = $filenames[0] =~ m/.*_bc(\d{1,2})_.*/g;
    }
    else
    {
	($number) = $files =~ m/.*_bc(\d{1,2})_.*/g;
    }
    return $number;
}
#####################################
#####################################
#sort by barcode name from barcode:

sub lowestNumberBC
{
    my $lines = shift;
    my @filenames = split("\t",$lines);
    my ($number) = $filenames[0] =~ m/.*bc(\d{1,2}).*/g;
    return $number;
}
#####################################
#####################################
#Compile HS metrics and Duplication metrics and Clipping Metrics

sub RunHSandDuplicationMetrics
{
    my ($bamList,$outdir) = @_;
    my @bams = @$bamList;
    my @HSmetricsFiles = ();
    my @DupMetricsFiles = ();
    my @Clip1MetricsFiles = (); 
    my @Clip2MetricsFiles = ();
    my $HSheader;
    my @HSmetrics;
    my $dupheader;
    my @dupmetrics;
    my $clip1header;
    my @clip1metrics;
    my $clip2header;
    my @clip2metrics;
    my($HSout) =  $bams[0] =~ /.*_bc\d+_(.*)_L\d{1,3}.*/;
    $HSout =  $HSout . "_ALL_HSmetrics.txt";
    foreach my $file (@bams)
    {
	my $hsFile = $file;
	$hsFile =~ s/\.bam/\.HSmetrics\.txt/;
	#print "H:$hsFile\n";
	push(@HSmetricsFiles, "$outdir/$hsFile");
	my $dupFile = $file;
	if($dupFile =~ /_IR_FX_BR/){
	$dupFile =~ s/_IR_FX_BR\.bam/\.metrics/;
	}
	if($dupFile =~ /_IR_BR/){
		$dupFile =~ s/_IR_BR\.bam/\.metrics/;
	}
	#print "D:$dupFile\n";	
	push (@DupMetricsFiles,"$outdir/$dupFile");
	my $clip1File = $file;
	my $clip2File = $file;
	if($clip1File =~ /_IR_FX_BR/){
		$clip1File =~ s/_mrg_cl_aln_srt_MD_IR_FX_BR\.bam/\_R1_mrg_cl.stats/;
	
	}
	if($clip1File =~ /_IR_BR/){
		$clip1File =~ s/_mrg_cl_aln_srt_MD_IR_BR\.bam/\_R1_mrg_cl.stats/;
	
	}
	if($clip2File =~ /_IR_FX_BR/){
		$clip2File =~ s/_mrg_cl_aln_srt_MD_IR_FX_BR\.bam/\_R2_mrg_cl.stats/;
	
	}
	if($clip2File =~ /_IR_BR/){
		$clip2File =~ s/_mrg_cl_aln_srt_MD_IR_BR\.bam/\_R2_mrg_cl.stats/;
	
	}
	#print "C1:$clip1File\nC2:$clip2File\n";
	push (@Clip1MetricsFiles,"$outdir/$clip1File");
	push (@Clip2MetricsFiles,"$outdir/$clip2File");
    }

    open (OUT, ">$HSout") or die $logger->fatal("can't open output file $HSout. Error: $!\n");
    for (my $j=0; $j<=$#HSmetricsFiles; $j++) 
    {
	open (FH, "<$HSmetricsFiles[$j]") or die $logger->fatal("can't open $HSmetricsFiles[$j]. Error: $!\n");
	while (my $text = <FH>) 
	{
	    chomp $text;
	    if ($text =~ "BAIT_SET") 
	    {
		$HSheader = $text;
		my $newline = <FH>;
		chomp $newline;
		$HSmetrics[$j] = $newline;
	    }
	}
	close (FH);
	open (FH, "<$DupMetricsFiles[$j]") or die $logger->fatal("Can't open $DupMetricsFiles[$j]. Error: $!\n");
	while (my $text = <FH>)
	{
	    chomp $text;
	    if ($text =~ "LIBRARY")
	    {
		$dupheader = "$text\tboth reads align\tone read aligns\tneither read aligns";
		my $newline = <FH>;
	    chomp $newline;
		my @val = split "\t", $newline;
		my $both = $val[2];
		my $one = $val[1];
		my $none = ($val[3]-$val[1])/2;
		$dupmetrics[$j] = "$newline\t$both\t$one\t$none";
	    }
	}
	close (FH);

	open (FH, "<$Clip1MetricsFiles[$j]") or die $logger->fatal("Can't open $Clip1MetricsFiles[$j]. Error:  $!\n");
	$clip1header = "Read1TotalReads\tRead1TrimmedReads\tPerRead1Trimmed";
	my $processedReads1 = 0;
	my $trimReads1 = 0;
	my $perTrimmedReads1 = 0;
	while (my $text = <FH>)
	{
	    chomp $text;
	    
	    if ($text =~ "Processed reads:") 
	    {
		#print "I:$text\n";
		my @data = split(":",$text);
		$data[1] =~ s/\s//g;
		$processedReads1 = $data[1];
	    }
	    if ($text =~ "Trimmed reads:")
	    {
		my @data = split(":",$text);
		my @realdata = split(/\(/,$data[1]);
		$realdata[0] =~ s/\s//g;
		$realdata[1] =~ s/\s|\)|\%//g;
		$trimReads1 = $realdata[0];
		$perTrimmedReads1 = $realdata[1];
	    }
	}
	#print "T:$processedReads1\t$trimReads1\n";
	$clip1metrics[$j] = "$processedReads1\t$trimReads1\t$perTrimmedReads1";
	close (FH);
	open (FH, "<$Clip2MetricsFiles[$j]") or die $logger->fatal("Can't open $Clip2MetricsFiles[$j]. Error: $!\n");
	$clip2header = "Read2TotalReads\tRead2TrimmedReads\tPerRead2Trimmed";
	my $processedReads2 = 0;
	my $trimReads2 = 0;
	my $perTrimmedReads2 = 0;
	while (my $text = <FH>)
	{
	    chomp $text;
	    if ($text =~ "Processed reads:")
	    {
		my @data = split(":",$text);
		$data[1] =~ s/\s//g;
		$processedReads2 = $data[1];
	    }
	    if ($text =~ "Trimmed reads:")
	    {
		my @data = split(":",$text);
		my @realdata = split(/\(/,$data[1]);
		$realdata[0] =~ s/\s//g;
		$realdata[1] =~ s/\s|\)|\%//g;
		$trimReads2 = $realdata[0];
		$perTrimmedReads2 = $realdata[1];
	    }
	}
	$clip2metrics[$j] = "$processedReads2\t$trimReads2\t$perTrimmedReads2";
	close (FH);
    }
    print OUT "Sample\t$HSheader\t$dupheader\t$clip1header\t$clip2header\n";
    for (my $i=0; $i<scalar(@bams); $i++) 
    {
	#print "TT:$bams[$i]\n";
	my($barcode) = $bams[$i] =~ m/.*_(bc\d+).*/g;
	print OUT "$barcode\t$HSmetrics[$i]\t$dupmetrics[$i]\t$clip1metrics[$i]\t$clip2metrics[$i]\n";
    }

    close(OUT);
    return($HSout);
}

#####################################
#####################################
#Compile Insert Size Metrics
sub RunInsertSizeMetrics
{
    my ($bamList,$outdir) = @_;
    my @bams = @$bamList;
    my @InsMetricsFiles = ();
    my @insertmetrics;
    my @insertsum;
    my @insertpeak;
    my($INSout) =  $bams[0] =~ /.*_bc\d+_(.*)_L\d{1,3}.*/;
    $INSout =  $INSout . "_ALL_insertsizemetrics.txt";
    foreach my $file (@bams)
    {
	my $insFile = $file;
	$insFile =~ s/\.bam/\.insert_size_metrics/;
	#print "I:$insFile\n";
	push(@InsMetricsFiles, "$outdir/$insFile");
    }
    for (my $j=0; $j<=$#InsMetricsFiles; $j++)
    {
	$insertpeak[$j][1]=0;
	open (FH, "<$InsMetricsFiles[$j]") or die $logger->fatal("can't open $InsMetricsFiles[$j]. Error:  $!\n");
	while (my $text = <FH>)
	{
	    chomp $text;
	    if ($text =~ "insert_size")
	    {
		while (my $newline = <FH>)
		{
		    chomp $newline;
		    if ($newline) 
		    {
			my @count = split "\t", $newline;
			$insertmetrics[$j][$count[0]] = $count[1];
			$insertsum[$j] += $count[1];
			if ($count[1]>$insertpeak[$j][1]) 
			{
			    $insertpeak[$j][0] = $count[0]; ##argmax (insert)
			    $insertpeak[$j][1] = $count[1]; ##max (count)
			}
		    }
		}
	    }
	}
	close (FH);
    }
    open (OUT, ">$INSout") or die $logger->fatal("Can't open output file $INSout. Error: $!\n");
    print OUT "insert_size";
    for (my $i=0; $i<=$#bams; $i++)
    {
	my($barcode) = $bams[$i] =~ m/.*_(bc\d+).*/g;
	print OUT "\t$barcode";
    }
    print OUT "\n";
    for (my $ins=1; $ins<=500; $ins++)
    {
	print OUT "$ins";
	for (my $j=0; $j<=$#insertmetrics; $j++)
	{
	    if ($insertmetrics[$j][$ins])
	    {
		my $ratio = $insertmetrics[$j][$ins]/$insertsum[$j];
		print OUT "\t$ratio";
	    }
	    else
	    {
		print OUT "\t0";
	    }
	}
	print OUT "\n";
    }
    print OUT "\n";
    print OUT "Peak";
    for (my $j=0; $j<=$#bams; $j++)
    {
	print OUT "\t$insertpeak[$j][0]";
    }
    print OUT "\n";
    close(OUT);
    return($INSout);
}

#####################################
#####################################
#Compile Base quality metrics

sub RunBaseQualityMetrics
{
    my ($bamList,$outdir) = @_;
    my @bams = @$bamList;
    my @BQmetricsFiles;
    my($BQRout) =  $bams[0] =~ /.*_bc\d+_(.*)_L\d{1,3}.*/;
    $BQRout =  $BQRout . "_ALL_basequalities.txt"; 
    my($BQOout) =  $bams[0] =~ /.*_bc\d+_(.*)_L\d{1,3}.*/;
    $BQOout =  $BQOout . "_ALL_orgbasequalities.txt";
    foreach my $file (@bams)
    {
	my $bqFile = $file;
	$bqFile =~ s/\.bam/\.quality_by_cycle_metrics/;
	#print "bq:$bqFile\n";
	push(@BQmetricsFiles, "$outdir/$bqFile");
    }
    open (OUT1, ">$BQRout") or die $logger->fatal("can't open output file $BQRout, $!\n");
    open (OUT2, ">$BQOout") or die $logger->fatal("can't open output file $BQOout, $!\n");
    my $baseheader = "cycle";
    my @basequalities = ();
    my @baseoriginalqualities = ();
    for (my $j=0; $j<=$#BQmetricsFiles; $j++)
    {
	my($barcode) = $BQmetricsFiles[$j] =~ m/.*_(bc\d+).*/g;
	open (FH, "<$BQmetricsFiles[$j]") or die "can't open $BQmetricsFiles[$j]\n";
	while (my $text = <FH>)
	{
	    chomp $text;
	    if ($text =~ "CYCLE")
	    {
		#$baseheader = $baseheader."\t$barcode\_MeanQuality\t$barcode\_MeanOriginalQuality";
		$baseheader = $baseheader."\t$barcode";
		while (my $newline = <FH>)
		{
		    chomp $newline;
		    my @scores = split "\t", $newline;
		    push @{$basequalities[$j]}, $scores[1];
		    push @{$baseoriginalqualities[$j]}, $scores[2];
		}
	    }
	}
	close (FH);
    }
    print OUT1 "$baseheader\n"; 
    print OUT2 "$baseheader\n";
    for (my $cyc=0; $cyc<=$#{$basequalities[0]}; $cyc++)
    {
	my $truecycle = $cyc+1;
	print OUT1 "$truecycle";
	print OUT2 "$truecycle";
	for (my $bc=0; $bc<=$#BQmetricsFiles; $bc++)
	{
	    if ($basequalities[$bc][$cyc])
	    {
		print OUT1 "\t$basequalities[$bc][$cyc]";
	    }
	    if ($baseoriginalqualities[$bc][$cyc])
	    {
		print OUT2 "\t$baseoriginalqualities[$bc][$cyc]";
	    }
	}
	print OUT1 "\n";
	print OUT2 "\n";
    }
    close(OUT1);
    close(OUT2);
    return($BQRout,$BQOout);
}

#####################################
#####################################
#Compile Exon Coverage Metrics
sub RunExonCoverageMetrics
{
    my ($bamList,$outdir) = @_;
    my @bams = @$bamList;
    my @ECmetricsFiles = ();
    my $covgheader;
    my @covgvalues;
    my @intervals; 
    my($ECout) =  $bams[0] =~ /.*_bc\d+_(.*)_L\d{1,3}.*/;
    $ECout =  $ECout . "_ALL_exoncoverage.txt"; 
    foreach my $file (@bams)
    {
	my $ecFile = $file;
	$ecFile =~ s/\.bam/\.gene\.covg\.sample_interval_summary/;
	#print "ec:$ecFile\n";
	push(@ECmetricsFiles, "$outdir/$ecFile");
    }
    open (OUT, ">$ECout") or die $logger->fatal("can't open output file $ECout, $!\n");
    for (my $j=0; $j<=$#ECmetricsFiles; $j++)
    {
	@intervals = ();
	open (FH, "<$ECmetricsFiles[$j]") or die $logger->fatal("can't open $ECmetricsFiles[$j], $!\n");
	my $topline = <FH>;
	while (my $text = <FH>)
	{
	    chomp $text;
	    my @cov = split "\t", $text;
	    push @intervals, $cov[0];
	    push @{$covgvalues[$j]}, $cov[2];
	}
	close (FH);
    }
    print OUT "0\tTarget";
    for (my $bc=0; $bc<=$#ECmetricsFiles; $bc++)
    {
	my($barcode) = $ECmetricsFiles[$bc] =~ m/.*_(bc\d+).*/g;
	print OUT "\t$barcode";
    }
    print OUT "\n";
    for (my $i=0; $i<=$#{$covgvalues[0]}; $i++)
    {
	my $tally = $i+1;
	print OUT "$tally\t$intervals[$i]";
	for (my $bc=0; $bc<=$#ECmetricsFiles; $bc++) 
	{
	    print OUT "\t$covgvalues[$bc][$i]";
	}
	print OUT "\n";
    }
    close(OUT);
    return($ECout);
}

#####################################
#####################################
#Compile Nomapq Exon Coverage Metrics
sub RunExonNomapqCoverageMetrics
{
    my ($bamList,$outdir) = @_;
    my @bams = @$bamList;
    my @ECmetricsFiles = ();
    my $covgheader;
    my @covgvalues;
    my @intervals; 
    my($ECout) =  $bams[0] =~ /.*_bc\d+_(.*)_L\d{1,3}.*/;
    $ECout =  $ECout . "_ALL_exonnomapqcoverage.txt"; 
    foreach my $file (@bams)
    {
	my $ecFile = $file;
	$ecFile =~ s/\.bam/\.gene_nomapq\.covg\.sample_interval_summary/;
	#print "ec:$ecFile\n";
	push(@ECmetricsFiles, "$outdir/$ecFile");
    }
    open (OUT, ">$ECout") or die $logger->fatal("can't open output file $ECout, $!\n");
    for (my $j=0; $j<=$#ECmetricsFiles; $j++)
    {
	@intervals = ();
	open (FH, "<$ECmetricsFiles[$j]") or die $logger->info("can't open $ECmetricsFiles[$j], Skip $!\n");
	my $topline = <FH>;
	while (my $text = <FH>)
	{
	    chomp $text;
	    my @cov = split "\t", $text;
	    push @intervals, $cov[0];
	    push @{$covgvalues[$j]}, $cov[2];
	}
	close (FH);
    }
    print OUT "0\tTarget";
    for (my $bc=0; $bc<=$#ECmetricsFiles; $bc++)
    {
	my($barcode) = $ECmetricsFiles[$bc] =~ m/.*_(bc\d+).*/g;
	print OUT "\t$barcode";
    }
    print OUT "\n";
    for (my $i=0; $i<=$#{$covgvalues[0]}; $i++)
    {
	my $tally = $i+1;
	print OUT "$tally\t$intervals[$i]";
	for (my $bc=0; $bc<=$#ECmetricsFiles; $bc++) 
	{
	    print OUT "\t$covgvalues[$bc][$i]";
	}
	print OUT "\n";
    }
    close(OUT);
    return($ECout);
}


#####################################
#####################################
# Compile Canonical Exon Coverage Metrics
sub RunCanonicalExonCoverageMetrics {
    my ($bamList,$exonIntervals,$outdir) = @_;
    my @bams = @$bamList;
    my @ECmetricsFiles = ();
    my $covgheader;
    my @covgvalues;
    my @intervals; 
    my %canonicalExons;
    my($ECout) =  $bams[0] =~ /.*_bc\d+_(.*)_L\d{1,3}.*/;
    $ECout =  $ECout . "_ALL_Canonical_exoncoverage.txt"; 
    foreach my $file (@bams)
    {
		my $ecFile = $file;
		$ecFile =~ s/\.bam/\.canonical\.exon\.covg\.sample_interval_summary/;
		print "ec:$ecFile\n";
		push(@ECmetricsFiles, "$outdir/$ecFile");
    }
    open (OUT, ">$ECout") or die $logger->fatal("can't open output file $ECout, $!\n");
    for (my $j=0; $j<=$#ECmetricsFiles; $j++)
    {
		@intervals = ();
		open (FH, "<$ECmetricsFiles[$j]") or die $logger->fatal("can't open $ECmetricsFiles[$j], $!\n");
			my $topline = <FH>;
			while (my $text = <FH>)
			{
			    chomp $text;
			    my @cov = split "\t", $text;
			    push @intervals, $cov[0];
			    push @{$covgvalues[$j]}, $cov[2];
			}
		close (FH);
    }
    
    
    print OUT "0\tTarget";
    for (my $bc=0; $bc<=$#ECmetricsFiles; $bc++)
    {
		my($barcode) = $ECmetricsFiles[$bc] =~ m/.*_(bc\d+).*/g;
		print OUT "\t$barcode";
    }
    print OUT "\n";
    for (my $i=0; $i<=$#{$covgvalues[0]}; $i++)
    {
		my $tally = $i+1;
		print OUT "$tally\t$intervals[$i]";
		for (my $bc=0; $bc<=$#ECmetricsFiles; $bc++) 
		{
		    print OUT "\t$covgvalues[$bc][$i]";
		}
		print OUT "\n";
    }
    close(OUT);
    return($ECout);
}





#####################################
#####################################
#Compile Gene Coverage Metrics
sub RunGeneCoverageMetrics
{
    my ($bamList,$outdir) = @_;
    my @bams = @$bamList;
    my @GCmetricsFiles = ();
    my @genecovgvalues;
    my @genes;
    my @gene_intervals;
    my($GCout) =  $bams[0] =~ /.*_bc\d+_(.*)_L\d{1,3}.*/;
    $GCout =  $GCout . "_ALL_genecoverage.txt";
    foreach my $file (@bams)
    {
	my $gcFile = $file;
	$gcFile =~ s/\.bam/\.gene\.cn\.txt/;
	#print "gc:$gcFile\n";
	push(@GCmetricsFiles, "$outdir/$gcFile");
    }
    open (OUT, ">$GCout") or die $logger->fatal("can't open output file $GCout, $!\n");
    for (my $j=0; $j<=$#GCmetricsFiles; $j++)
    {
	@genes=();
	@gene_intervals=();
	open (FH, "<$GCmetricsFiles[$j]") or die $logger->fatal("can't open $GCmetricsFiles[$j], $!\n");
	while (my $text = <FH>)
	{
	    chomp $text;
	    my @cov = split "\t", $text;
	    push @genes, $cov[0];
	    push @{$genecovgvalues[$j]}, $cov[1];
	    push @gene_intervals, $cov[2];
	}
	close (FH);
    }
    print OUT "Gene";
    for (my $bc=0; $bc<=$#GCmetricsFiles; $bc++)
    {
	my($barcode) = $GCmetricsFiles[$bc] =~ m/.*_(bc\d+).*/g;
	print OUT "\t$barcode";
    }
    print OUT "\tLocus\n";
    for (my $i=0; $i<=$#{$genecovgvalues[0]}; $i++)
    {
	print OUT "$genes[$i]";
	for (my $bc=0; $bc<=$#GCmetricsFiles; $bc++)
	{
	    print OUT "\t$genecovgvalues[$bc][$i]";
	}
	print OUT "\t$gene_intervals[$i]\n";
    }
    close(OUT);
    return($GCout);
}



#####################################
#####################################
#Compile Finger Print Genotype Metrics
sub RunFPGenotypeMetrics
{
    my ($bamList,$patientId,$outdir) = @_;
    my @bams = @$bamList;
    my @FPGmetricsFiles = ();
    my $FPheader = "Locus";
    my @FPcounts;
    my @FPgenotype;
    my @FPratios;
    my @FPcoords;
    my @FPhom;
    #Calculate Homozygous position based on availability of normal sample.
    
    tie (my %groupedFilenames, 'Tie::IxHash');
    tie (my %groupedbarcode, 'Tie::IxHash');
    my($FPGout) =  $bams[0] =~ /.*_bc\d+_(.*)_L\d{1,3}.*/;
    my $FPG1out =  $FPGout . "_ALL_FPsummary.txt";
    my $FPG2out =  $FPGout . "_ALL_FPavgHom.txt";
    my $FPG3out =  $FPGout . "_ALL_FPhet.txt";
    open (OUT1, ">$FPG1out") or die $logger->fatal("can't open output file $FPG1out\n");
    #open (OUTF, ">$FPGFout") or die "can't open output file $FPGFout\n";
    #Combine Same Patient Files
    my $fCount = 0;
    foreach my $file (@bams)
    {
	#my ($commonSampleId) = $file =~ /(.*)([A-Z]\d|[A-Z]|[A-Z]-\d{1,3}ng|[A-Z]\d-\d{1,3}ng)_bc/;
        #my($commonPatient) = $commonSampleId;
	my ($barcode) = $file =~ /.*_(bc\d+)_.*/;
	#print "$file\n";
	#print "$barcode\n";
	if(exists $groupedbarcode{@$patientId[$fCount]})
	{
	    my $barcodes =  $groupedbarcode{@$patientId[$fCount]};
	    $barcodes = "$barcodes" . ",$barcode";
	    $groupedbarcode{@$patientId[$fCount]} = "$barcodes";
	}
	else
	{
	    $groupedbarcode{@$patientId[$fCount]} = "$barcode";
	}
	if(exists $groupedFilenames{@$patientId[$fCount]})
	{
	    my $files =  $groupedFilenames{@$patientId[$fCount]};
	    $files = "$files" . ",$file";
	    $groupedFilenames{@$patientId[$fCount]} = "$files";
	}
	else
	{
	    $groupedFilenames{@$patientId[$fCount]} = "$file";
	}
	$fCount++;
    }   
    my $count = 0;
    while((my $key, my $value) = each (%groupedFilenames))
    {
	tie ( my %StdFPcoords, 'Tie::IxHash');
	my @files = split(",",$value);
	# Section of Normal
	my $normal;
	my @normalSamples =();
	tie (my %coverageForSampleNormals, 'Tie::IxHash');
	my @CoverageForMultipleNormal = ();
	foreach my $file (@files)
	{
	    my($fileBarcode) = $file =~ /.*_(bc\d+)_.*/;
	    my $fileClass = $classPerBarcode{$fileBarcode};
            #print "FileClass:$fileClass\n";
            
	    if($fileClass =~ m/Normal/i)
	    {
		push (@normalSamples, $file);
	    }
	}
	#print "NS:",scalar(@normalSamples),"\n";
	
	if (scalar @normalSamples >= 1)
	{
	    foreach my $file (@normalSamples)
	    {
		my $HSmetricsFile = $file;
		$HSmetricsFile =~ s/\.bam/\.HSmetrics\.txt/g;
		#print "HS:$HSmetricsFile\n";
		open(FH,"$outdir/$HSmetricsFile") or die $logger->fatal("Cannot Open $outdir/$HSmetricsFile, $!\n");
		my $meanCov;
		my $CovForFile;
		while(<FH>)
		{
		    next until($_ =~ /^BAIT_SET/);
		    while(<FH>)
		    {
			next if(($_ =~ /^BAIT_SET/) or ($_ =~ /^\s$/));
			my(@values) = split("\t",$_);
			#print "MeanCOv:$values[21]\n";
			$CovForFile = $values[21];
			$meanCov = $values[21];
		    }
		    $coverageForSampleNormals{$file} = $CovForFile;
		    push (@CoverageForMultipleNormal, $meanCov);
		}
		close(FH);
	    }
	    #Get file that will be used as normal
            
	    my $maxCoverage = max @CoverageForMultipleNormal;
	    #print "MAX:$maxCoverage\n"; 
	    #if (scalar @normalSamples > 1)
	    #{
	    while((my $key, my $value) = each (%coverageForSampleNormals))
                {
                    if($value == $maxCoverage)
                        {
                            if ($value >= 50) {
                                $normal = $key;
                            }
                            else {
                                $normal = "NULL";
                            }
                        }
                    
		    else
                        {
                            next;
                        }
		}
        }
        else
            {
                $normal = "NULL";
            }
	#print "NOrmal:","$normal\n";
	
	if ($normal ne "NULL")
	{
	    my $normalFPGfile = $normal;
	   # print "Normal:$normalFPGfile\n";
	    $normalFPGfile =~ s/\.bam/\.FP\.summary\.txt/;
	    open (FH, "<$outdir/$normalFPGfile") or die $logger->fatal("Can't open $outdir/$normalFPGfile. Error: $!\n");
	    my $row = 0;
	    while (my $text = <FH>)
	    {
		#print "$normalFPGfile\t$text\n";
		my @line = split ("\t", $text);
		$line[3] = 0.0 if (! $line[3]);
		if($line[3] < 0.10)
		{
		    $StdFPcoords{ $line[0]}=$line[2];
		}
		$row++;
	    }
	    close(FH);
	    foreach my $file(@files)
	    {
		my $fpgFile = $file; 
		#print "All:$fpgFile\n";
		$fpgFile =~ s/\.bam/\.FP\.summary\.txt/;
		my($barcode) = $fpgFile =~ m/.*_(bc\d+).*/g;
		$FPheader = "$FPheader" . "\t$barcode" . "_Counts\t$barcode" . "_Genotypes\t$barcode" . "_MinorAlleleFreq";
		push(@FPGmetricsFiles, $fpgFile);
		open (FH, "<$outdir/$fpgFile") or die $logger->fatal("Can't open $outdir/$fpgFile. Error: $!\n");
		$row = 0;
		while (my $text = <FH>)
		{
		    chomp $text; 
		    my @line = split "\t", $text;
		    $line[3] = 0.0 if (! $line[3]);
		    $FPcoords[$row] = $line[0];
		    push @{$FPcounts[$count]}, $line[1];
		    push @{$FPgenotype[$count]}, $line[2];
		    push @{$FPratios[$count]}, $line[3];
		    if (exists $StdFPcoords{$line[0]})
		    {
			#print "1:$count:$line[3]\n";
			push @{$FPhom[$count]}, $line[3];
		    }
		    $row++;
		}
		$count++;
		close(FH);
	    }
	}
	else
	{
	    foreach my $file(@files)
	    {
		my $fpgFile = $file;	
		#print "All:$fpgFile\n";
		$fpgFile =~ s/\.bam/\.FP\.summary\.txt/;
		push(@FPGmetricsFiles,$fpgFile);
		my($barcode) = $fpgFile =~ m/.*_(bc\d+).*/g;
		$FPheader = "$FPheader" . "\t$barcode" . "_Counts\t$barcode" . "_Genotypes\t$barcode" . "_MinorAlleleFreq";
		open (FH, "<$outdir/$fpgFile") or die $logger->fatal("Can't open $outdir/$fpgFile. Error: $!\n");
		my $row = 0;
		while (my $text = <FH>)
		{
		    chomp $text; 
		    my @line = split "\t", $text;
		    $line[3] = 0.0 if (! $line[3]);
		    $FPcoords[$row] = $line[0];
		    push @{$FPcounts[$count]}, $line[1];
		    push @{$FPgenotype[$count]}, $line[2];
		    push @{$FPratios[$count]}, $line[3];
		    if ($line[3] < 0.10)
		    {
			#print "2:$count:$line[3]\n";
			push @{$FPhom[$count]}, $line[3];
		    }
		    $row++;
		}
		$count++;
		close(FH);
	    }
	}
    }
    #Calculate Average Minor Allele Freq for Homozygous positions
    my @FPcontamination;
    for (my $j=0; $j<=$#FPGmetricsFiles; $j++)
    {
	$FPcontamination[$j]=0;
	for (my $hom=0; $hom<=$#{$FPhom[$j]}; $hom++)
	{
	    $FPcontamination[$j] += $FPhom[$j][$hom];
	}
	#print "HOM:$j:$#{$FPhom[$j]}\n";
	#print "CON:$j:$FPcontamination[$j]\n";
	$FPcontamination[$j] = $FPcontamination[$j] / ($#{$FPhom[$j]} + 1);
    }
    print OUT1 "$FPheader\n";
    for (my $snp=0; $snp<=$#{$FPcounts[0]}; $snp++)
    {
	print OUT1 "$FPcoords[$snp]";
	for (my $bc=0; $bc<=$#FPGmetricsFiles; $bc++)
	{
	    print OUT1 "\t$FPcounts[$bc][$snp]\t$FPgenotype[$bc][$snp]\t$FPratios[$bc][$snp]";
	}
	print OUT1 "\n";
    }
    close(OUT1);
    
    ####Print Avg Homozygous frequency
    open (OUT2, ">$FPG2out") or die $logger->fatal("Can't open output file $FPG2out. Error: $!\n");
    print OUT2 "Sample\tAvgMinorHomFreq\n";
    my @dataforAvgMinorHomFreq = ();
    for (my $j=0; $j<=$#FPGmetricsFiles; $j++)
    {
	my($barcode) = $FPGmetricsFiles[$j] =~ m/.*_(bc\d+)_.*/g;
	my $FPcon = sprintf "%.5f",$FPcontamination[$j];
	push(@dataforAvgMinorHomFreq,"$barcode\t$FPcon");
	#print OUT2 "$barcode\t$FPcontamination[$j]\n";
    }
    my(@sorteddataforAvgMinorHomFreq) = sort {lowestNumberBC($a) <=>  lowestNumberBC($b)} @dataforAvgMinorHomFreq;
    foreach my $SampleAvgHomFreq (@sorteddataforAvgMinorHomFreq)
    {
	print OUT2 "$SampleAvgHomFreq\n";
    }
    close(OUT2);

    ####Farction of Hetrozygous Position
    open (OUT3, ">$FPG3out") or die $logger->fatal("Can't open output file $FPG3out. Error: $!\n");
    #Calculate Average hetrozygosity.
    my @dataFracHetPos = ();
    print OUT3 "Sample\tPerHetrozygosPos\n";
    foreach my $file(@FPGmetricsFiles)
    {
	my($barcode) = $file =~ m/.*_(bc\d+)_.*/g;
	open (FH, "<$file") or die "can't open $file, $!\n";
	my $hetCount = 0;
	my $totalCount = 0;
	while(<FH>)
	{
	    $totalCount++;
	    my @line = split("\t",$_);
	    if(($line[2] eq "AA") or ($line[2] eq "TT") or ($line[2] eq "GG") or ($line[2] eq "CC") or ($line[2] eq "--"))
	    {
		next
	    }
	    else
	    {
		$hetCount++;
	    }
	}
	close(FH);
	#print OUT3 "$barcode\t";
	my $hetrozygosity = ($hetCount/$totalCount);
	$hetrozygosity = sprintf "%.3f",$hetrozygosity;
	push(@dataFracHetPos,"$barcode\t$hetrozygosity");
	#printf OUT3 ("%.2f\n",$hetrozygosity );
    }
    my(@sorteddataFracHetPos) = sort {lowestNumberBC($a) <=>  lowestNumberBC($b)} @dataFracHetPos;
    foreach my $SampleFracHetPos (@sorteddataFracHetPos)
    {
	print OUT3 "$SampleFracHetPos\n";
    }
    close(OUT3);

    #############
    #Sample Mixup.
    #############
    tie (my %RefvsSampleGenotypes, 'Tie::IxHash');
    my $Cheader;
    tie (my %unexpectedmismatch, 'Tie::IxHash');
    tie (my %unexpectedmatch, 'Tie::IxHash');
    my @unexpectedmismatch = ();
    my @unexpectedmatch = ();
    my($FPGCout) =  $bams[0] =~ /.*_bc\d+_(.*)_L\d{1,3}.*/;
    my $FPGC1out =  $FPGCout . "_ALL_FPCsummary.txt";
    my $FPGC2out =  $FPGCout . "_ALL_FPCResultsUnMatch.txt";
    my $FPGC3out =  $FPGCout . "_ALL_FPCResultsUnMismatch.txt";
    open (OUTA, ">$FPGC1out") or die $logger->fatal("Can't open output file $FPGC1out. Error: $!\n");
    open (OUTC, ">$FPGC2out") or die $logger->fatal("Can't open output file $FPGC2out. Error: $!\n");
    open (OUTB, ">$FPGC3out") or die $logger->fatal("Can't open output file $FPGC3out. Error: $!\n");
	my(@sortedFPGmetricsFiles) = sort {lowestNumber($a) <=>  lowestNumber($b)} @FPGmetricsFiles;
    foreach my $refFile(@sortedFPGmetricsFiles)
    {
	chomp($refFile);
	tie %{$RefvsSampleGenotypes{$refFile}},'Tie::IxHash';
    }
    foreach my $refFile(@sortedFPGmetricsFiles)
    {
	chomp($refFile);
	my($refbarcode) = $refFile =~ m/.*_(bc\d+).*/g;
	$Cheader = $Cheader . "\t$refbarcode";
	open (FH, "<$refFile") or die $logger->fatal("Can't open $refFile. Error: $!\n");
	tie (my %refGenotypes, 'Tie::IxHash');
	my $refCount = 0;
	while(<FH>)
	{
	    my @line = split("\t",$_);
	    if(($line[2] eq "AA") or ($line[2] eq "TT") or ($line[2] eq "GG") or ($line[2] eq "CC"))
	    {
		$refGenotypes{$line[0]} = $line[2];
		$refCount++;
	    }
	    else
	    {
		next;
	    }
	}
	close(FH);
	#tie (my %sampleGenotypes, 'Tie::IxHash');
	foreach my $file(@sortedFPGmetricsFiles)
	{
	    chomp($file);
	   # my($barcode) = $file =~ m/.*-(bc\d+).*/g;
	   # push(@Crows,$barcode);
	    open (FH, "<$file") or die $logger->fatal("Can't open $file. Error: $!\n");
	    my $genotypeCount = 0;
	    while(<FH>)
	    {
		my @line = split("\t",$_);
		if(exists $refGenotypes{$line[0]})
		{
		    if(($line[2] eq "AA") or ($line[2] eq "TT") or ($line[2] eq "GG") or ($line[2] eq "CC"))
		    {
			my $refValue = $refGenotypes{$line[0]};
			my $value = $line[2];
			if($refValue ne $value)
			{
			    $genotypeCount++;
			}
		    }
		    else
		    {
			next;
		    }
		}
		else
		{
		    next;
		}
	    }
	    close(FH);
	    my $avgGenotypeCount = 0;
	    $avgGenotypeCount = $genotypeCount/$refCount if ($refCount != 0);
	    $avgGenotypeCount = sprintf "%.2f", $avgGenotypeCount;
	    $RefvsSampleGenotypes{$refFile}{$file} = $avgGenotypeCount;
	}
    }
    print OUTA "$Cheader\n";
    my $sCount = 0;
    #print "@$patientId\n";
    while ( my($refSamples, $querySamples) = each %RefvsSampleGenotypes ) 
    {
	my ($refPatientId) = @$patientId[$sCount];
	my($refbarcode) = $refSamples =~ m/.*_(bc\d+).*/g;
	#print "$refPatientId", "-$refbarcode", "-MixupAnalysis:\n";
	print OUTA "$refbarcode";
	my $PateintBarcodes = $groupedbarcode{$refPatientId};
	#print "$PateintBarcodes\n";
	while ( my($sample, $avgAltAlleleCount) = each %$querySamples ) 
	{
	    #print "$sample=$avgAltAlleleCount\n";
	    print OUTA "\t$avgAltAlleleCount";
	    my($qbarcode) = $sample =~ /.*_(bc\d+)_.*/;
	    if($PateintBarcodes =~ /$qbarcode/)
	    {
		if($avgAltAlleleCount > 0.05)
		{
		    #print OUT2 "Unexpected Mismatch:$refPatientId","-$refbarcode"," vs $qPatientId","-$qbarcode"," have high fraction of discordant alleles i.e $avgAltAlleleCount\n";
		    my ($printRefId) = $refSamples =~ /(.*_bc\d+)_.*/;
		    my ($printqueryId) = $sample =~ /(.*_bc\d+)_.*/;
		    my $unMismatch1 = "$printRefId"."\t$printqueryId";
		    my $unMismatch2 = "$printqueryId". "\t$printRefId";
		    my $avgAltAlleleCounts;
		    if((exists $unexpectedmismatch{$unMismatch1}) or (exists $unexpectedmismatch{$unMismatch2}))
		    {
			my $avgAltAlleleCounts1 = $unexpectedmismatch{$unMismatch1};
			my $avgAltAlleleCounts2 = $unexpectedmismatch{$unMismatch2};
			if($avgAltAlleleCounts1)
			{
			    $avgAltAlleleCounts = $avgAltAlleleCounts1;
			    $avgAltAlleleCounts = $avgAltAlleleCounts . ",$avgAltAlleleCount";
			    $unexpectedmismatch{$unMismatch1} = "$avgAltAlleleCounts";
			}
			if($avgAltAlleleCounts2)
			{
			    $avgAltAlleleCounts = $avgAltAlleleCounts2;
			    $avgAltAlleleCounts = $avgAltAlleleCounts . ",$avgAltAlleleCount";
			    $unexpectedmismatch{$unMismatch2} = "$avgAltAlleleCounts";
			}
		    }
		    else
		    {
			 $unexpectedmismatch{$unMismatch1} = "$avgAltAlleleCount";
		    }
		}
	    }
	    else
	    {
		if($avgAltAlleleCount < 0.05)
		{
		    #print OUT2 "Unxpected Match:$refPatientId","-$refbarcode"," vs $qPatientId","-$qbarcode"," have low fraction of discordant alleles i.e $avgAltAlleleCount\n";
		    #my $unMatch = "$refPatientId"."-$refbarcode".",$qPatientId"."-$qbarcode"; 
		    my ($printRefId) = $refSamples =~ /(.*_bc\d+)_.*/;
		    my ($printqueryId) = $sample =~ /(.*_bc\d+)_.*/;
		    my $unMatch1 = "$printRefId"."\t$printqueryId";
		    my $unMatch2 = "$printqueryId". "\t$printRefId";
		    #my $unMatch1 = "$refPatientId"."-$refbarcode"."\t$qPatientId"."-$qbarcode";
		    #my $unMatch2 = "$qPatientId"."-$qbarcode". "\t$refPatientId"."-$refbarcode";
		    my $avgAltAlleleCounts;
		    if((exists $unexpectedmatch{$unMatch1}) or (exists $unexpectedmatch{$unMatch2}))
		    {
			my $avgAltAlleleCounts1 = $unexpectedmatch{$unMatch1};
			my $avgAltAlleleCounts2 = $unexpectedmatch{$unMatch2};
			if($avgAltAlleleCounts1)
			{
			    $avgAltAlleleCounts = $avgAltAlleleCounts1;
			    $avgAltAlleleCounts = $avgAltAlleleCounts . ",$avgAltAlleleCount";
			    $unexpectedmatch{$unMatch1} = "$avgAltAlleleCounts";
			}
			if($avgAltAlleleCounts2)
			{
			    $avgAltAlleleCounts = $avgAltAlleleCounts2;
			    $avgAltAlleleCounts = $avgAltAlleleCounts . ",$avgAltAlleleCount";
			    $unexpectedmatch{$unMatch2} = "$avgAltAlleleCounts";
			}
		    }
		    else
		    {
			 $unexpectedmatch{$unMatch1} = "$avgAltAlleleCount";
		    }
		}
	    }
	}
	print OUTA "\n";
	$sCount++;
    }
    close(OUTA);
 
    print OUTB "#Unexpected Mismatch:\n";
    my (@keysUnMismatch) = keys(%unexpectedmismatch);
    my (@keysUnMatch) = keys(%unexpectedmatch);
    if((scalar @keysUnMismatch) != 0)
    {
	print OUTB "Sample1\tSample2\tFractionOfDicordantAlleles\n";
	while(my ($key,$value) = each %unexpectedmismatch)
	{
	    print OUTB "$key\t$value\n";
	}
    }
    else
    {
	 print OUTB "\n";
	 print OUTB "---There are no Unexpected Mismatches.---\n";
    }
    close(OUTB);
    print OUTC "#Unexpected Match:\n";
    if((scalar @keysUnMatch) != 0)
    {
	print OUTC "Sample1\tSample2\tFractionOfDiscordantAlleles\n";
	while(my ($key,$value) = each %unexpectedmatch)
	{
	    print OUTC "$key\t$value\n";
	}
    }
    else
    {
	print OUTC "\n";
	print OUTC "---There are no Unexpected Matches.---\n";
    }
    close(OUTC);
    return($FPG1out,$FPG2out,$FPG3out,$FPGC1out,$FPGC2out);
}

#####################################
#####################################
#Compile Finger Print Tiling Coverage Metrics
sub RunFPtilingCoverageMetrics
{
    my ($bamList,$outdir) = @_;
    my @bams = @$bamList;
    my @FPTmetricsFiles = ();
    my $tilingheader;
    my @tilingvalues;
    my @tilingintervals; 
    my($FPTout) =  $bams[0] =~ /.*_bc\d+_(.*)_L\d{1,3}.*/;
    $FPTout =  $FPTout . "_ALL_tilingcoverage.txt";
    foreach my $file (@bams)
    {
	my $fptFile = $file;
	$fptFile =~ s/\.bam/\.tiling\.covg\.sample_interval_summary/;
	#print "fpt:$fptFile\n";
	push(@FPTmetricsFiles, "$outdir/$fptFile");
    }
    open (OUT, ">$FPTout") or die $logger->fatal("Can't open output file $FPTout. Error: $!\n");
    for (my $j=0; $j<=$#FPTmetricsFiles; $j++)
    {
	@tilingintervals = ();
	open (FH, "<$FPTmetricsFiles[$j]") or die $logger->fatal("Can't open $FPTmetricsFiles[$j]. Error:  $!\n");
	my $topline = <FH>;
	while (my $text = <FH>)
	{
	    chomp $text;
	    my @cov = split "\t", $text;
	    push @tilingintervals, $cov[0];
	    push @{$tilingvalues[$j]}, $cov[2];
	}
	close (FH);
    }
    print OUT "0\tTarget";
    for (my $bc=0; $bc<=$#FPTmetricsFiles; $bc++)
    {
	my($barcode) = $FPTmetricsFiles[$bc] =~ m/.*_(bc\d+)_.*/g;
	print OUT "\t$barcode";
    }
    print OUT "\n";
    for (my $i=0; $i<=$#{$tilingvalues[0]}; $i++)
    {
	my $tally = $i+1;
	print OUT "$tally\t$tilingintervals[$i]";
	for (my $bc=0; $bc<=$#FPTmetricsFiles; $bc++)
	{
	    print OUT "\t$tilingvalues[$bc][$i]";
	}
	print OUT "\n";
    }
    close(OUT);
    return($FPTout);
}

#####################################
#####################################
#Compile Finger Print Tiling Coverage Metrics
sub RunFPtilingNomapqCoverageMetrics
{
    my ($bamList,$outdir) = @_;
    my @bams = @$bamList;
    my @FPTmetricsFiles = ();
    my $tilingheader;
    my @tilingvalues;
    my @tilingintervals; 
    my($FPTout) =  $bams[0] =~ /.*_bc\d+_(.*)_L\d{1,3}.*/;
    $FPTout =  $FPTout . "_ALL_tilingnomapqcoverage.txt";
    foreach my $file (@bams)
    {
	my $fptFile = $file;
	$fptFile =~ s/\.bam/\.tiling_nomapq\.covg\.sample_interval_summary/;
	#print "fpt:$fptFile\n";
	push(@FPTmetricsFiles, "$outdir/$fptFile");
    }
    open (OUT, ">$FPTout") or die $logger->fatal("Can't open output file $FPTout. Error: $!\n");
    for (my $j=0; $j<=$#FPTmetricsFiles; $j++)
    {
	@tilingintervals = ();
	open (FH, "<$FPTmetricsFiles[$j]") or die $logger->info("Can't open $FPTmetricsFiles[$j].Skip:  $!\n");
	my $topline = <FH>;
	while (my $text = <FH>)
	{
	    chomp $text;
	    my @cov = split "\t", $text;
	    push @tilingintervals, $cov[0];
	    push @{$tilingvalues[$j]}, $cov[2];
	}
	close (FH);
    }
    print OUT "0\tTarget";
    for (my $bc=0; $bc<=$#FPTmetricsFiles; $bc++)
    {
	my($barcode) = $FPTmetricsFiles[$bc] =~ m/.*_(bc\d+)_.*/g;
	print OUT "\t$barcode";
    }
    print OUT "\n";
    for (my $i=0; $i<=$#{$tilingvalues[0]}; $i++)
    {
	my $tally = $i+1;
	print OUT "$tally\t$tilingintervals[$i]";
	for (my $bc=0; $bc<=$#FPTmetricsFiles; $bc++)
	{
	    print OUT "\t$tilingvalues[$bc][$i]";
	}
	print OUT "\n";
    }
    close(OUT);
    return($FPTout);
}

sub RunTargetCoverageMetrics
{
    my ($bamList,$outdir) = @_;
    my @bams = @$bamList;
    my @HSTmetricsFiles = ();
    my $targetheader;
    my @targetvalues;
    my @targetintervals;
    tie (my %gcbiasResults, 'Tie::IxHash');
    my($HSTout) =  $bams[0] =~ /.*_bc\d+_(.*)_L\d{1,3}.*/;
    my $HST1out =  $HSTout . "_ALL_targetcoverage.txt";
    my $HST2out =  $HSTout . "_ALL_gcbias.txt";
   
    foreach my $file (@bams)
    {
	my $hstFile = $file;
	$hstFile =~ s/\.bam/\.target\.covg/;
	#print "hst:$hstFile\n";
	push(@HSTmetricsFiles, "$outdir/$hstFile");
    }
    ##Target Coverage Analysis
    open (OUT1, ">$HST1out") or die $logger->fatal("Can't open output file $HST1out. Error: $!\n");
    for (my $j=0; $j<=$#HSTmetricsFiles; $j++)
    {
	@targetintervals = ();
	open (FH, "<$HSTmetricsFiles[$j]") or die $logger->fatal("Can't open $HSTmetricsFiles[$j], $!\n");
	my $topline = <FH>;
	while (my $text = <FH>)
	{
	    chomp $text;
	    my @cov = split "\t", $text;
	    my $interval = "$cov[0]" . ":$cov[1]" . "-$cov[2]";
	    push @targetintervals, $interval;
	    push @{$targetvalues[$j]}, $cov[6];
	}
	close (FH);
    }
    print OUT1 "0\tTarget";
    for (my $bc=0; $bc<=$#HSTmetricsFiles; $bc++)
    {
	my($barcode) = $HSTmetricsFiles[$bc] =~ m/.*_(bc\d+)_.*/g;
	print OUT1 "\t$barcode";
    }
    print OUT1 "\n";
    for (my $i=0; $i<=$#{$targetvalues[0]}; $i++)
    {
	my $tally = $i+1;
	print OUT1 "$tally\t$targetintervals[$i]";
	for (my $bc=0; $bc<=$#HSTmetricsFiles; $bc++)
	{
	    print OUT1 "\t$targetvalues[$bc][$i]";
	}
	print OUT1 "\n";
    }
    close(OUT1);
    ###GC Bias Analysis
    my $GCheader = "\t";
    my @xkeys = ();
    my @gcResults =();
    open (OUT2, ">$HST2out") or die "can't open output file $HST2out, $!\n";

    for (my $j=0; $j<=$#HSTmetricsFiles; $j++)
    {
	my ($barcode) = $HSTmetricsFiles[$j] =~ m/.*_(bc\d+)_.*/g;
	$GCheader = $GCheader . "$barcode\t";
	my (@gcCalcArr1,@gcCalcArr2,@gcCalcArr3,@gcCalcArr4,@gcCalcArr5,@gcCalcArr6,@gcCalcArr7,@gcCalcArr8,@gcCalcArr9,@gcCalcArr10,@gcCalcArr11) = ();	
	my (@ncCalcArr1,@ncCalcArr2,@ncCalcArr3,@ncCalcArr4,@ncCalcArr5,@ncCalcArr6,@ncCalcArr7,@ncCalcArr8,@ncCalcArr9,@ncCalcArr10,@ncCalcArr11) = ();
	my ($gcM1,$gcM2,$gcM3,$gcM4,$gcM5,$gcM6,$gcM7,$gcM8,$gcM9,$gcM10,$gcM11) = 0;	
	my ($ncM1,$ncM2,$ncM3,$ncM4,$ncM5,$ncM6,$ncM7,$ncM8,$ncM9,$ncM10,$ncM11) = 0;
	open (FH, "<$HSTmetricsFiles[$j]") or die "can't open $HSTmetricsFiles[$j], $!\n";
	my $topline = <FH>;
	while (my $text = <FH>)
	{
	    chomp $text;
	    my @cov = split "\t", $text;
	   
	    if($cov[3] > 2)
	    {
		if($cov[5] < 0.30)
		{
		    push(@gcCalcArr1,$cov[5]);
		    push(@ncCalcArr1,$cov[7]);
		}
		if(($cov[5] >= 0.30) and ($cov[5] < 0.35))
		{
		    push(@gcCalcArr2,$cov[5]);
		    push(@ncCalcArr2,$cov[7]);
		}
		if(($cov[5] >= 0.35) and ($cov[5] < 0.40))
		{
		    push(@gcCalcArr3,$cov[5]);
		    push(@ncCalcArr3,$cov[7]);
		}
		if(($cov[5] >= 0.40) and ($cov[5] < 0.45))
		{
		    push(@gcCalcArr3,$cov[5]);
		    push(@ncCalcArr3,$cov[7]);
		}
		if(($cov[5] >= 0.45) and ($cov[5] < 0.50))
		{
		    push(@gcCalcArr4,$cov[5]);
		    push(@ncCalcArr4,$cov[7]);
		}
		if(($cov[5] >= 0.50) and ($cov[5] < 0.55))
		{
		    push(@gcCalcArr5,$cov[5]);
		    push(@ncCalcArr5,$cov[7]);
		}
		if(($cov[5] >= 0.55) and ($cov[5] < 0.60))
		{
		    push(@gcCalcArr6,$cov[5]);
		    push(@ncCalcArr6,$cov[7]);
		}
		if(($cov[5] >= 0.60) and ($cov[5] < 0.65))
		{
		    push(@gcCalcArr7,$cov[5]);
		    push(@ncCalcArr7,$cov[7]);
		}
		if(($cov[5] >= 0.65) and ($cov[5] < 0.70))
		{
		    push(@gcCalcArr8,$cov[5]);
		    push(@ncCalcArr8,$cov[7]);
		}
		if(($cov[5] >= 0.70) and ($cov[5] < 0.75))
		{
		    push(@gcCalcArr9,$cov[5]);
		    push(@ncCalcArr9,$cov[7]);
		}
		if(($cov[5] >= 0.75) and ($cov[5] < 0.80))
		{
		    push(@gcCalcArr10,$cov[5]);
		    push(@ncCalcArr10,$cov[7]);
		}
		if($cov[5] >= 0.80)
		{
		    push(@gcCalcArr11,$cov[5]);
		    push(@ncCalcArr11,$cov[7]);
		}
	    }
	    else
	    {
		next;
	    }
	}
	close (FH);
	$gcM1 = sprintf "%.4f",mean(@gcCalcArr1);
	$gcM2 = sprintf "%.4f",mean(@gcCalcArr2);
	$gcM3 = sprintf "%.4f",mean(@gcCalcArr3);
	$gcM4 = sprintf "%.4f",mean(@gcCalcArr4);
	$gcM5 = sprintf "%.4f",mean(@gcCalcArr5);
	$gcM6 = sprintf "%.4f",mean(@gcCalcArr6);
	$gcM7 = sprintf "%.4f",mean(@gcCalcArr7);
	$gcM8 = sprintf "%.4f",mean(@gcCalcArr8);
	$gcM9 = sprintf "%.4f",mean(@gcCalcArr9);
	$gcM10 = sprintf "%.4f",mean(@gcCalcArr10);
	$gcM11 = sprintf "%.4f",mean(@gcCalcArr11);
	$ncM1 = sprintf "%.4f",mean(@ncCalcArr1);
	$ncM2 = sprintf "%.4f",mean(@ncCalcArr2);
	$ncM3 = sprintf "%.4f",mean(@ncCalcArr3);
	$ncM4 = sprintf "%.4f",mean(@ncCalcArr4);
	$ncM5 = sprintf "%.4f",mean(@ncCalcArr5);
	$ncM6 = sprintf "%.4f",mean(@ncCalcArr6);
	$ncM7 = sprintf "%.4f",mean(@ncCalcArr7);
	$ncM8 = sprintf "%.4f",mean(@ncCalcArr8);
	$ncM9 = sprintf "%.4f",mean(@ncCalcArr9);
	$ncM10 = sprintf "%.4f",mean(@ncCalcArr10);
	$ncM11 = sprintf "%.4f",mean(@ncCalcArr11);
	tie(my %gcResults, "Tie::IxHash");#Hash to store all calculated data.
	%gcResults = ($gcM1=>$ncM1,$gcM2=>$ncM2,$gcM3=>$ncM3,$gcM4=>$ncM4,$gcM5=>$ncM5,$gcM6=>$ncM6,$gcM7=>$ncM7,$gcM8=>$ncM8,$gcM9=>$ncM9,$gcM10=>$ncM10,$gcM11=>$ncM11);
	#my @gcBias = ($gcM1,$gcM2,$gcM3,$gcM4,$gcM5,$gcM6,$gcM7,$gcM8,$gcM9,$gcM10,$gcM11);
	#my @normCov = ($ncM1,$ncM2,$ncM3,$ncM4,$ncM5,$ncM6,$ncM7,$ncM8,$ncM9,$ncM10,$ncM11);
	push(@gcResults,\%gcResults);
    }
    #make the final hash to print stuff
    tie(my %gcResultstable, "Tie::IxHash");
    foreach my $sample (@gcResults)
    {
	while (my ($key,$value) = each %$sample)
	{
	    if(exists $gcResultstable{$key})
	    {
		my $sampleCovs = $gcResultstable{$key};
		$sampleCovs = "$sampleCovs" . "\t$value";
		$gcResultstable{$key} = $sampleCovs;
	    }
	    else
	    {
		$gcResultstable{$key} = $value;
	    }
	}
    }
    #print the final data to the file.
    chop($GCheader);
    print OUT2 "$GCheader\n";
    while (my($pergc,$cov) = each (%gcResultstable))
    {
	print OUT2 "$pergc\t$cov\n";
    }

    close(OUT2);
    return($HST1out,$HST2out);
}

sub PlotGraphs
{
    my($FPG1out,$ECout,$outdir) = @_;
    my @notifynames = ();
    &MakeCSH($outdir);
    my ($basename) = $FPG1out =~ /(.*)_ALL_/;
    
    my $Rscript = "$RHOME/RScript";
    my $R = "$RHOME/R";

    #my $Rscript = "/ifs/e63data/bergerm1/Resources/SupportTools/R/Install/R-2.15.2/lib64/R/bin/Rscript";
    #my $R = "/ifs/e63data/bergerm1/Resources/SupportTools/R/Install/R-2.15.2/lib64/R/bin/R";
    
    open(FH,">","$outdir/PostCompileWrapper.csh") || die $logger->fatal("Cannot Open PostCompileWrapper.csh. Error: $!\n");
    print FH "#!/bin/csh\n";
    print FH "$R --slave --vanilla --args ${basename} ${outdir} ${gcbiasfile} < $loessNorm\n";
    print FH "$R --slave --vanilla --args ${basename} ${stdnormloess_nvn} ${geneIntervalAnn} ${tilingIntervalAnn} Normal,NormalGL MIN < $nvnCN\n";
    print FH "cp ${outdir}/${basename}_copynumber_segclusp.nvn.pdf ${outdir}/${basename}_normvsnorm_segclusp.pdf\n";
    print FH "$R --slave --vanilla --args ${basename} ${stdnormloess_tm} ${geneIntervalAnn} ${tilingIntervalAnn} MIN < $bestCN\n";
    close(FH); 
    eval{`chmod 755 $outdir/PostCompileWrapper.csh`;};
    if($@){
	$logger->warn("Can not change the permissions on PostCompileWrapper.sh\n");
    }
    #my $cmd = "$QSUB -q $queue -v R_LIBS_USER=$RLIBS -wd $outdir -N PCW.$$ -b y -o /dev/null -e /dev/null -l h_vmem=5G,virtual_free=5G -pe smp 1 $outdir/PostCompileWrapper.csh";
    if(defined $QSUB)
    {
    my $cmd = "$QSUB -q $queue -wd $outdir -N PCW.$$ -b y -l h_vmem=5G,virtual_free=5G -pe smp 1 $outdir/PostCompileWrapper.csh";
    $logger->debug( "COMMAND: $cmd");
    eval{
	`$QSUB -q $queue -wd $outdir -N PCW.$$ -b y -l h_vmem=5G,virtual_free=5G -pe smp 1 "$outdir/PostCompileWrapper.csh"`;
	`$QSUB -q $queue -V -wd $outdir -hold_jid PCW.$$ -N NotifyPCW.$$ -e /dev/null -o NotifyPCW.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
    };
 
        if ($@) {
        $logger->fatal(
            "PostCompileWrapper:Cannot run PostCompileWrapper Error:$@"
        );
        exit(1);
    }
    }
    else
    {
    	my $cmd = "$BSUB -q $queue -J PCW.$$ -cwd $outdir -e PCW.$$.%J.stderr -o PCW.$$.%J.stdout -We 24:00 -R \"rusage[mem=5]\" -M 10 -n 1 \"$outdir/PostCompileWrapper.csh\"";
    	$logger->debug( "COMMAND: $cmd");
    	eval{
    		`$BSUB -q $queue -J PCW.$$ -cwd $outdir -e PCW.$$.%J.stderr -o PCW.$$.%J.stdout -We 24:00 -R "rusage[mem=5]" -M 10 -n 1 "$outdir/PostCompileWrapper.csh"`;
			`$BSUB -q $queue -cwd $outdir -w "done(PCW.$$) -J NotifyPCW.$$ -e NotifyPCW.$$.%J.stderr -o NotifyPCW.$$.stat -R "rusage[mem=2]" -M 4 -n 1 "$outdir/Notify.csh"`;	
           
    	};
    	if ($@) {
        $logger->fatal(
            "PostCompileWrapper:Cannot run PostCompileWrapper Error:$@"
        );
        exit(1);
    }
    }

    push(@notifynames,"NotifyPCW.$$.stat");
    if(defined $QSUB)
    {
    my $cmd = "$QSUB -q $queue -wd $outdir -N PreparePDF.$$ -b y -o PreparePDF.$$.stdout -e PreparePDF.$$.stderr -l h_vmem=5G,virtual_free=5G -pe smp 1 $PERL $MetricsScript $basename $R";
     $logger->debug( "COMMAND: $cmd");
 eval{
    `$QSUB -q $queue -wd $outdir -N PreparePDF.$$ -b y -o PreparePDF.$$.stdout -e PreparePDF.$$.stderr -l h_vmem=5G,virtual_free=5G -pe smp 1 "perl $MetricsScript $basename $R"`;
    `$QSUB -q $queue -V -wd $outdir -hold_jid PreparePDF.$$ -N NotifyPreparePDF.$$ -e NotifyPreparePDF.$$.stderr -o NotifyPreparePDF.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
};
    if ($@) {
        $logger->fatal(
            "PreParePDF:Cannot run PreparePdf Error:$@"
        );
        exit(1);
    }
    }
    else{
    	my $cmd = "$BSUB -q $queue -J PreparePDF.$$ -cwd $outdir -e PreparePDF.$$.%J.stderr -o PreparePDF.$$.%J.stdout -We 24:00 -R \"rusage[mem=5]\" -M 10 -n 1 "$PERL $MetricsScript $basename $R\"";
    	$logger->debug( "COMMAND: $cmd");
    	eval{
    		`$BSUB -q $queue -J PreparePDF.$$ -cwd $outdir -e PreparePDF.$$.%J.stderr -o PreparePDF.$$.%J.stdout -We 24:00 -R "rusage[mem=5]" -M 10 -n 1 "$PERL $MetricsScript $basename $R"`;
			`$BSUB -q $queue -cwd $outdir -w "post_done(PreparePDF.$$)" -J NotifyPreparePDF.$$ -e NotifyPreparePDF.$$.%J.stderr -o NotifyPreparePDF.$$.stat -We 24:00 -R "rusage[mem=2]" -M 4 -n 1 "$outdir/Notify.csh"`;	
           
    	};
    	if ($@) {
        $logger->fatal(
            "PreParePDF:Cannot run PreparePdf Error:$@"
        );
        exit(1);
    }
    }
    push(@notifynames,"NotifyPreparePDF.$$.stat");
    WaitToFinish($outdir,@notifynames);
    
    return;
}
