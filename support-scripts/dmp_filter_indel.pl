#!/usr/bin/perl 
#####FilterIndel.pl#####
#Author: Ronak Shah
#Date: 03/21/2013
#LastModified: 08/06/2013
#Version:1.2
#Description: Filter result of Somatic Indel Detector using text file and vcf.
##07/29/2013
#v1.2
#RS:Added functionality for alleledepth and varaint freq to be parameter.
#RS:All filter through single filter
######################
use strict;
use Getopt::Long;
use IO::File;
use Cwd;
use Tie::IxHash; 
use File::Basename;
use Vcf;
use lib "/dmp/resources/prod/libs/mskcclib/production/";
use MSKCC_DMP_Logger;
#--This variable holds the current time 
my $now = time;
my $logger = MSKCC_DMP_Logger->get_logger('Filter_SomaticIndel_Logger');
$logger->start_local();
my($IndelTxtFile,$IndelVcfFile,$sampleName,$outdir,$totaldepth,$alleledepth,$variantfreq,$TNratio);
if (@ARGV < 3 or !GetOptions (
	    'IndelTxt|t:s'              => \$IndelTxtFile,
	    'IndelVcf|v:s'              => \$IndelVcfFile,
	    'outdir|o:s'                => \$outdir, 
	    'totaldepth|dp:i'           => \$totaldepth,,
	    'alleledepth|ad:i'          => \$alleledepth,
	    'variantfreq|vf:f'          => \$variantfreq,
            'TNratio|tnr:i'            => \$TNratio,
            'sampleName|s:s'            => \$sampleName))
	{
		Usage();
	}
if((!defined $IndelTxtFile) or (!defined $IndelVcfFile) or (!defined $sampleName))
{
    $logger->fatal("Please enter indel text and vcf files as well as sample name for the process");;
    exit(1);
}

if(!defined $outdir)
{
    $outdir = getcwd;
}
if(!defined $totaldepth)
{
    $totaldepth = 0;
}
if(!defined $alleledepth)
{
    $alleledepth = 3;
}
if(!defined $variantfreq)
{
    $variantfreq = 0.01;
}
if(!defined $TNratio)
{
    $TNratio = 5;
}

&ReadAndFilter($IndelTxtFile,$IndelVcfFile,$sampleName,$outdir);

#--Calculate total runtime (current time minus start time) 
$now = time - $now;

#--Print runtime 
printf("\n\nTotal running time: %02d:%02d:%02d\n\n", int($now / 3600), int(($now % 3600) / 60), int($now % 60));

$logger->info("SomIndelStdFilter:!!!!Done, Thanks for using the script!!!");
exit(0);

#####################################
#####################################
#How to use the script.

sub Usage
{
	print "Unknow option: @_\n" if (@_);

	print "\nUsage : FilterIndel.pl [options]
        [--IndelTxtFile|t              S tab-delimted Indel file describing details about the mutations (required)]
        [--IndelVcfFile|v              S VCF format Indel file describing details about the mutations (required)]
        [--sampleName|s                S Name of the sample (required)]
        [--totaldepth|dp               I Tumor total depth threshold for Somatic Indel Detector(default:0,optional).]
        [--alleledepth|ad              I Tumor Allele depth threshold for Somatic Indel Detector(default:3,optional).]
        [--variantfreq|vf              F Tumor variant frequency threshold for Somatic Indel Detector(default:0.01,optional).]
        [--TNratio|tnr                 I Tumor-Normal variant frequency ratio threshold for Somatic Indel Detector(default:5,optional).]
        [--outdir|o                    S Path where all the output files will be written (optional) [default:current working directory]]
	\n";

	print "For any question please contact Ronak Shah (shahr2\@mskcc.org)\n";
	print "!~!~! Comments, flames and suggestion are welcome !~!~!\n";
	exit;
}
#####################################
#####################################
#Read files and filter

sub ReadAndFilter
{
    my($txtFile,$vcfFile,$sampleName,$outdir) = @_;
    tie my %keepCoordinates, 'Tie::IxHash';
    my $basename;
    if($vcfFile =~ /\//)
    {
	$basename = basename($vcfFile);
    }
    else
    {
	$basename = $vcfFile;
    }
    my($filteredVCFfile) = $basename =~ /(.*)\.vcf/;
    $filteredVCFfile = $filteredVCFfile . "_STDfilter.vcf";
    my($filteredOutFile) = $basename =~ /(.*)\.vcf/;
    $filteredOutFile = $filteredOutFile . "_STDfilter.txt";
    open(TFH,"$txtFile") or die ($logger->fatal("Cannot open $txtFile, Error:$!"));
    while(<TFH>)
    {
	next if (($. == 1) or ($. == 2));
	chomp($_);
	my @dataCols = split("\t",$_);
	my @newDataCols = grep(s/\s*$//g, @dataCols);
	my $chromCoord = $newDataCols[0] . ":$newDataCols[1]";
	my $n_obs_cons = shift@{[split("/",pop@{[split(":",$newDataCols[4])]})]};
	my $n_obs_total = pop@{[split("/",pop@{[split(":",$newDataCols[4])]})]};
	my $n_freq = 0;
	$n_freq = sprintf("%.5f",($n_obs_cons/$n_obs_total)) if ($n_obs_total != 0);
	my $t_obs_cons = shift@{[split("/",pop@{[split(":",$newDataCols[12])]})]};
	my $t_obs_total = pop@{[split("/",pop@{[split(":",$newDataCols[12])]})]};
	my $t_freq = 0;
	$t_freq = sprintf("%.5f",($t_obs_cons/$t_obs_total))if ($t_obs_total != 0);
	my $n_freq_5times = sprintf("%.5f" , ($TNratio*($n_freq)));
        #print "$totaldepth\t$alleledepth\t$variantfreq\t$n_freq\n";
        
	if($t_freq > $n_freq_5times)
	{
	    if(($t_obs_total >= $totaldepth) and ($t_obs_cons >= $alleledepth) and ($t_freq >= $variantfreq))
	    {
		if (exists $keepCoordinates{$chromCoord})
		{
		    $logger->fatal("SomaticIndelDetector:There is a repeat.$chromCoord");
		    exit(1);
		}
		else
		{
		    $keepCoordinates{$chromCoord} = "";
		}
	    }
	    else
                {
                    #print "$chromCoord\t$n_freq\t$t_freq\n";
                    next;
                }
	}
	else
            {
                #print "$chromCoord\t$n_freq\t$t_freq\n";
                next;
            }
    }
    close(TFH);
    
    #Populate text and vcf file
    open(FVFH,">","$outdir/$filteredVCFfile") or die ($logger->fatal("FilteredVCFfile: Cannot Open $outdir/$filteredVCFfile, Error:$!"));
    open(FTFH,">","$outdir/$filteredOutFile") or die ($logger->fatal("FilteredOutFile: Cannot open $outdir/$filteredOutFile, Error:$!"));
    open(VFH,"$vcfFile") or die ($logger->fatal("SomIndelVCFfile:Cannot Open $vcfFile, Error:$!"));
    while(<VFH>)
    {
	chomp($_);
	if ($_ =~ /^#/)
	{
	    print FVFH "$_\n";
	}
	my @dataCols = split("\t",$_);	
	my @newDataCols = grep(s/\s*$//g, @dataCols);
	my $key = $newDataCols[0] . ":$newDataCols[1]"; ;
	if(exists $keepCoordinates{$key})
	{
	    print FVFH "$_\n";
	    print FTFH "$sampleName\t$newDataCols[0]\t$newDataCols[1]\t$newDataCols[3]\t$newDataCols[4]\t.\n";
	}
	else
	{
	    next;
	}
    }
    close(VFH);
    close(FVFH);
    close(FTFH);
    return;
}

