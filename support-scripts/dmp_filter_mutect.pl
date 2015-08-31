#!/usr/bin/perl
#####FilterMutect.pl#####
#Author: Ronak Shah
#Date: 03/21/2013
#LastModified: 08/06/2013
#Version:1.1
#Description: Filter result of muTect v1.14 using callstats file and vcf.
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
my $logger = MSKCC_DMP_Logger->get_logger('Filter_Mutect_Logger');
$logger->start_local();
my($MutationTxtFile,$MutationVcfFile,$totaldepth,$alleledepth,$variantfreq,$TNratio,$sampleName,$outdir);
if (@ARGV < 3 or !GetOptions (
	    'MutationTxt|t:s'              => \$MutationTxtFile,
	    'MutationVcf|v:s'              => \$MutationVcfFile,
	    'outdir|o:s'                   => \$outdir,
            'totaldepth|dp:i'              => \$totaldepth,,
	    'alleledepth|ad:i'             => \$alleledepth,
	    'variantfreq|vf:f'             => \$variantfreq,
            'TNratio|tnr:i'                => \$TNratio,
	    'sampleName|s:s'               => \$sampleName))
	{
		Usage();
	}
if((!defined $MutationTxtFile) or (!defined $MutationVcfFile) or (!defined $sampleName))
{
    $logger->fatal("Please enter indel text and vcf files as well as sample name for the process");
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
    $alleledepth = 5;
}
if(!defined $variantfreq)
{
    $variantfreq = 0.01;
}
if(!defined $TNratio)
{
    $TNratio = 5;
}
&ReadAndFilter($MutationTxtFile,$MutationVcfFile,$sampleName,$outdir);

#--Calculate total runtime (current time minus start time) 
$now = time - $now;

#--Print runtime 
printf("\n\nTotal running time: %02d:%02d:%02d\n\n", int($now / 3600), int(($now % 3600) / 60), int($now % 60));

$logger->info("MutectStdFilter:!!!!Done, Thanks for using the script!!!");
exit;

#####################################
#####################################
#How to use the script.

sub Usage
{
	print "Unknow option: @_\n" if (@_);

	print "\nUsage : FilterMutect.pl [options]
        [--MutationTxtFile|t              S tab-delimted Mutect file describing details about the mutations (required)]
        [--MutationVcfFile|v              S VCF format Mutect file describing details about the mutations (required)]
        [--sampleName|s                   S Name of the sample (required)]
        [--totaldepth|dp                  I Tumor total depth threshold for Mutect(default:0,optional).]
        [--alleledepth|ad                 I Tumor Allele depth threshold for Mutect(default:3,optional).]
        [--variantfreq|vf                 F Tumor variant frequency threshold for Mutect(default:0.01,optional).]
        [--TNratio|tnr                    I Tumor-Normal variant frequency ratio threshold for Mutect(default:5,optional).]
        [--outdir|o                       S Path where all the output files will be written (optional) [default:current working directory]]
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
    open(TFH,"$txtFile") or die($logger->fatal("MutectFilter:Cannot open $txtFile, Error:$!"));
    while(<TFH>)
    {
	next if (($. == 1) or ($. == 2));
	chomp($_);
	my @dataCols = split("\t",$_);
	my @newDataCols = grep(s/\s*$//g, @dataCols);
	my $chromCoord = $newDataCols[0] . ":$newDataCols[1]" . ":$newDataCols[3]" . ":$newDataCols[4]";
	if($newDataCols[68] =~ /KEEP/)
	{
	    #print "$chromCoord\n";
	    if (exists $keepCoordinates{$chromCoord})
	    {
		$logger->fatal("MutectStdFilter:There is a repeat.$chromCoord");
		exit(1);
	    }
	    else
	    {
		$keepCoordinates{$chromCoord} = $newDataCols[68];
	    }
	}else{
	    my $n_refAllele = $newDataCols[46];
	    my $n_altAllele = $newDataCols[47]; 
	    my $n_total = ($n_refAllele + $n_altAllele); 
	    my $n_freq = 0;
	    $n_freq = sprintf("%.5f", ($n_altAllele/$n_total)) if($n_total != 0);
	    my $n_freq_5times = $TNratio*($n_freq);
	    my $t_refAllele = $newDataCols[28] ;
	    my $t_altAllele = $newDataCols[29]; 
	    my $t_total = ($t_refAllele + $t_altAllele);
	    my $t_freq = 0;
	    my $t_freq = sprintf("%.5f", ($t_altAllele/$t_total)) if($t_total != 0);
            #print "$totaldepth\t$alleledepth\t$variantfreq\t$TNratio\n";
            
	    if($t_freq >= $n_freq_5times)
	    {
		#print "$t_freq:$n_freq\n";
		if(($t_total >= $totaldepth) and ($t_altAllele >= $alleledepth) and ($t_freq >= $variantfreq))
		{
		    my @reasonToreject = split(",",$newDataCols[67]);
		    
		    my @accepted_tags = qw(alt_allele_in_normal nearby_gap_events triallelic_site possible_contamination clustered_read_position);
		    my $count = 0;
		    foreach my $accepted_tag (@accepted_tags){
			if($newDataCols[67] =~/$accepted_tag/){
			    $count++;
			}
		    }
		    if($count != scalar (@reasonToreject)){
			next;
		    }
		    if (exists $keepCoordinates{$chromCoord}) {
			$logger->fatal("MutectStdFilter:There is a repeat.$chromCoord");
			exit(1);
		    }else{
			$keepCoordinates{$chromCoord} = $newDataCols[67];
		    }
		}   
	    }     
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
	my $key = $newDataCols[0] . ":$newDataCols[1]" . ":$newDataCols[3]" . ":$newDataCols[4]"; 
	if(exists $keepCoordinates{$key})
	{
	    my $failureReason = $keepCoordinates{$key};
	    $failureReason = "." if($failureReason eq "");
	    if(($_) =~ /REJECT/){$_ =~ s/REJECT/PASS/g};
	    print FVFH "$_\n";
	    print FTFH "$sampleName\t$newDataCols[0]\t$newDataCols[1]\t$newDataCols[3]\t$newDataCols[4]\t$failureReason\n";
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

