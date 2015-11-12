#!/usr/bin/perl 
#####GenotypeAllele.pl#####
#Author: Ronak Shah(RS),Donavan Cheng(DC)
#Date: 04/09/2013
#LastModified: 11/14/2013
#Version:1.1
#Description: Genotype each variant.
##04/04/2013
#v1.0
##04/09/2013
#v1.1
#RS:Updated to handle multiple same location different sample vcf
##04/12/2013
#v1.1
#RS:Updated to handle overlapping bed file and Mutation and Indel in same sample and same coordinates
##04/16/2013
#v1.1
#Updated to give refCount.
##04/22/2013
#v1.1
#RS:Removed the bug for multiple deletions
##08/08/2013
#v1.2
#RS:Added option of SGE queue
######################

use strict;
use Getopt::Long;
use IO::File;
use Cwd;
use Tie::IxHash; 
use File::Basename;
use MSKCC_DMP_Logger;

my $logger = MSKCC_DMP_Logger->get_logger('GENOTYPE_ALLELE');
$logger->start_local();

#--This variable holds the current time 
my $now = time;
my ($filteredMutationVcfFile,$bamFile,$deleteFiles,$QSUB,$BSUB,$refFile,$mbq,$mmq,$pathToSamtools,$outdir,$outFile,$pathToBedtools,$mpileUpOutFile,$bamId,$queue,$typeOfSample);
if (@ARGV < 3 or !GetOptions (
	    'FilteredMutationListVcf|fmv:s'     => \$filteredMutationVcfFile,
	    'BamFile|bam:s'                     => \$bamFile,
	    'RefFile|rf:s'                      => \$refFile,
	    'samtools|s:s'                      => \$pathToSamtools,
            'bedtools|b:s'                      => \$pathToBedtools,
	    'MinBaseQuality|mbq:i'              => \$mbq,
            'MinMappingQuality|mmq:i'           => \$mmq,
	    'deleteUnwantedFiles|d:i'           => \$deleteFiles,
            'outFile|of:s'                      => \$outFile,
	    'mpileUpOutFile|mof:s'              => \$mpileUpOutFile,
	    'bamId|bi:s'                        => \$bamId,
	    'qsub:s'			      	=> \$QSUB,
	    'bsub:s'			      	=> \$BSUB,
            'queue|q:s'                         => \$queue,
            'typeOfSample|tos:s'                => \$typeOfSample,
	    'outdir|o:s'                        => \$outdir))
	{
		Usage();
	}
if((!defined $filteredMutationVcfFile) or (!defined $bamFile))
{
    $logger->fatal("Mutation file or .bam file is missing. Please enter filtered Mutation file in as vcf format. Also enter the Bam file to genotype.See Usage\n");
    &Usage();
    exit;
}
if(!defined $deleteFiles)
{
    $deleteFiles = 2;
}
if(!defined $refFile)
{
    $logger->warn("Reference file was not supplied. The default (/home/shahr2/References/Homo_sapiens_assembly19.fasta) will be used.\n");
    $refFile = "/dmp/data/pubdata/hg-fasta/production/Homo_sapiens_assembly19.fasta";
}
if(!defined $pathToSamtools)
{
    $logger->warn("Path for samtools was not supplied. The default (/home/shahr2/Software/Samtools/current/samtools) will be used.\n");
    $pathToSamtools = "/home/shahr2/Software/Samtools/current/samtools";
}
if(!defined $pathToBedtools)
{
    $logger->warn("Path for bedtools was not supplied. The default (/home/shahr2/Software/BEDTools/current/bin/) will be used.\n");
    $pathToBedtools = "/home/shahr2/Software/BEDTools/current/bin/";
}
if(!defined $mbq)
{
    $logger->warn("Minimum base quality was not supplied. MBQ will be set to 5.\n");
    $mbq = 5;
}
if(!defined $mmq)
{
    $logger->warn("Minimum mapping quality was not supplied. MMQ will be set to 5.\n");
    $mmq = 5;
}
if(!defined $outdir)
{
    $logger->warn("Output directory was not supplied. The current working directory will be used.\n");
    $outdir = getcwd;
}
if(!defined $outFile)
{
    if($bamFile =~ /\//)
    {
	($outFile) = basename($bamFile);
	$outFile=~ s/\.bam/_mpileup\.alleledepth/;
    }
    else
    {
	($outFile) = $bamFile =~ /(.*)\.bam/;
	$outFile = $outFile . "_mpileup.alleledepth";
    }
}
else
{
    ($bamId) = $outFile =~ /(.*)_bc\d+/;
}
if(!defined $bamId)
{
    if($bamFile =~ /\//)
    {
	$bamId = basename($bamFile);
	($bamId) = $bamId =~ /(.*)_bc\d+/;
	if($bamId eq "")
	{
	    $bamId = basename($bamFile);
	}
    }
    else
    {
	($bamId) = $bamFile =~ /(.*)_bc\d+/;
	if($bamId eq "")
	{
	    $bamId = $bamFile;
	}
    }
}
if(!defined $mpileUpOutFile)
{
    if($bamFile =~ /\//)
    {
	($mpileUpOutFile) = basename($bamFile);
	$mpileUpOutFile=~ s/\.bam/\.mpileup/;
    }
    else
    {
	($mpileUpOutFile) = $bamFile =~ /(.*)\.bam/;
	$mpileUpOutFile = $mpileUpOutFile . ".mpileup";
    }
}
#Check  for qsub
if(!defined $QSUB)
{
    $logger->warn("Path to QSUB command not given, Assuming we will run LSF");
}
else
{
    $logger->info( "SGE QSUB:$QSUB\n");
}
#Check  for bsub
if(!defined $BSUB)
{
    $logger->warn("Path to BSUB command not given, Assuming we will run SGE");
}
else
{
    $logger->info( "LSF BSUB:$BSUB\n");
}
#Check  for queue
if(!defined $queue)
{
    $logger->warn("Name of the SGE/LSF queue not given default (all.q) will be used.\n");
    $queue = "all.q";
}
else
{
    $logger->info( "SGE/LSF Queue:$queue\n");
}
#Check  for type of Sample
if(!defined $typeOfSample)
{
    $logger->warn("Type of Sample not given default (Tumor) will be used.\n");
    $typeOfSample = "Tumor";
}
&MakeCSH($outdir);
my($bedFile,$mergedBedFile,$vcfHash,$mutationHash) = &MakeBedFile($filteredMutationVcfFile,$mpileUpOutFile,$outdir);
my($mpileupOut) = &RunMpileup($bamFile,$mergedBedFile,$mpileUpOutFile,$bamId);
my($pileupHash,$pileupstatsHash) = &AnalyzeMpileup($mpileupOut,$mutationHash);
my($mergeOutFile) = &Merge_Vcf_mpileup($vcfHash,$mutationHash,$pileupHash,$pileupstatsHash,$outdir,$outFile,$bamId);
my @deleteFiles = ();
push(@deleteFiles,$mpileupOut);
push(@deleteFiles,$bedFile);
push(@deleteFiles,$mergedBedFile);
#deletion process
if($deleteFiles == 2)
{
    &HouseKeeping($outdir,@deleteFiles);
}
if($deleteFiles == 1)
{
    $logger->info("No house keeping performed.\n");
}
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

    print "\nUsage : GenotypeAllele.pl [options]
        [--FilteredMutationVcfFile|fmv     S vcf file describing details about the mutations (required)]
        [--BamFile|bam                     S bam file to be used for genotyping (required)]
        [--RefFile|rf                      S Path to genome reference file (optional;default:/home/shahr2/References/Homo_sapiens_assembly19.fasta)]
        [--samtools|s                      S Path to samtools (optional;default:/home/shahr2/Software/Samtools/current/samtools)]
        [--bedtools|b                      S Path to bedtools (optional;default:/home/shahr2/Software/BEDTools/current/bin/)]
        [--MinBaseQuality|mbq              I Min. Base Quality Threshold (optional;default:5)]
        [--MinMappingQuality|mmq           I Min. Mapping Quality Threshold (optional;default:5)]
        [--deleteUnwantedFiles|d           I 2=>To delete files 1=> To keep files (default:2,optional)]
        [--outdir|o                        S Path where all the output files will be written (optional) [default:current working directory]]
        [--outFile|of                      S Name of the allele depth output file (optional) [default:BamFame-.bam+_mpileup.alleledepth]]
        [--bamId|bi                        S Bam Id to be used (optional) [default:bamfile name]]
        [--queue|q                         S Name of the SGE / LSF Queue where the pipeline needs to run (default:all.q,optional)]
        [--qsub                            S Path to qsub executable for SGE(default:None,optional)]
        [--bsub                            S Path to bsub executable for LSF(default:None,required)]
        [--mpileUpOutFile|mof              S Name of samtools mpileup output file (optional) [default:BamFile-.bam+.mpileup]]
        [--typeOfSample|tos                S Type of Sample (default:Tumor,optional)]
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
    if($@){
	$logger->warn("Couln't change the permissions on $filename.\n");
	return;
    }
}
###################################################
###################################################
#--Waiting for the process to finish

sub WaitToFinish
{
	my($outdir,@waitfilenames) = @_;
	$logger->info("Waiting for the Process to finish...");
	foreach my $wfile(@waitfilenames)
	{
            next if($wfile eq "NULL");
	    sleep 10 while(!(-e "$outdir/$wfile"));
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
###################################################
###################################################
#--Make Bed file and make a tracking hashes
sub MakeBedFile
{
    my($vcf,$mpileUpOutFile,$outdir) = @_;
    my $basename;
    my $intervalsFile;
    if($mpileUpOutFile =~ /\//g)
    {
	$basename = basename($mpileUpOutFile);
	if($basename =~/\.mpileup$/)
	{

	    ($intervalsFile) = $basename =~ /(.*)\.mpileup/;
	    $intervalsFile = $intervalsFile . "_intervals.bed";
	}
	else
	{
	    $intervalsFile = $basename . "_intervals.bed";
	}
    }
    else
    {
	$basename = $mpileUpOutFile;
	if($basename =~/\.mpileup$/)
	{

	    ($intervalsFile) = $basename =~ /(.*)\.mpileup/;
	    $intervalsFile = $intervalsFile . "_intervals.bed";
	}
	else
	{
	    $intervalsFile = $basename . "_intervals.bed";
	}
    }
    tie(my %vcfHash, 'Tie::IxHash');
    tie(my %altHash, 'Tie::IxHash');
    #making bed file from mutation file
    open(IN,"<",$vcf) || die "Cannot open $vcf,$!\n";
    open(OUT, ">","$outdir/$intervalsFile") ||  die $logger->fatal("Cannot open $intervalsFile,$!\n");
    while(<IN>)
    {
	next if ($_ =~ /^#/);
	chomp($_);
	my @dataCols = split("\t",$_);
	my @newDataCols = grep(s/\s*$//g, @dataCols);
	my $ref = $newDataCols[3];
	my $alt = $newDataCols[4];
	my $chrom = $newDataCols[0];
	my $start = $newDataCols[1];
	my $info = $newDataCols[7];
	my($sampleinfo,$rest) = split(";",$info);
	my($name,$sampleId) = split("=",$sampleinfo);
	my $pileupStart;
	my $pileupEnd;
	#Start for Mutation then form Indels
	$pileupStart = $start - 1;
	if((abs(length($ref))) == (abs(length($alt))))
	{
	    $pileupEnd = $start;
	    $altHash{$chrom.":".$start}{$sampleId. ":" .$alt.":".$ref} = $alt;
	}
	else
	{
	    $pileupEnd = $start + 1;
	    my $indel;
	    my $indel_len = abs(length($alt) - length($ref));
	    if(length($ref) > length($alt))#deletion
	    {
		$indel = "-".$indel_len.substr($ref,1);
	    }
	    elsif(length($alt) > length($ref))#insertion
	    {
		$indel = "+".$indel_len.substr($alt,1);
	    }
	    
	    $altHash{$chrom.":".$start}{$sampleId.":".$alt.":".$ref} = $indel; 
	    $altHash{$chrom.":".$pileupEnd}{$sampleId.":".$alt.":".$ref} = $indel;
	    
	}
	$vcfHash{$chrom.":".$start.":".$ref.":".$alt.":".$sampleId}=$_;
	my $data = "$chrom\t$pileupStart\t$pileupEnd\n";
	print OUT "$data";
    }
    close IN;
    close OUT;
    my($mergedBedFile) = $intervalsFile;
    $mergedBedFile =~ s/\.bed/_sorted_merged\.bed/;
    if((-e "$outdir/$mergedBedFile") and (-s "$outdir/$mergedBedFile" != 0)){
	$logger->info("$outdir/$mergedBedFile exists andbed file will not be created.\n");
    }
    else{
	($mergedBedFile) = &MakeProperBed($intervalsFile,$outdir);
    }
    return($intervalsFile,$mergedBedFile,\%vcfHash,\%altHash);
}
###################################################
###################################################
#--Merge Overlapping and Adjacent Entries.
sub MakeProperBed
{
    my($bedFile,$outdir) = @_;
    my $mergedBedFile = $bedFile;
    $mergedBedFile =~ s/\.bed/_sorted_merged\.bed/;
    eval{`$pathToBedtools/sortBed -i $outdir/$bedFile | $pathToBedtools/bedtools merge -d 1 > $outdir/$mergedBedFile`;};
    if($@){
	$logger->fatal("Could not merge .bed file: $outdir/$bedFile.\n");
    }
    return($mergedBedFile);
}

###################################################
###################################################
#--Run mpileup from samtools
sub RunMpileup
{
    my($bam,$bed,$mpileUpOutFile,$bamId) = @_;
    my @notifyNames = ();
    if(! $bam =~ /\//)
    {
	$bam = "$outdir/$bam";
    }
    my $command = "";
    if(defined $QSUB){
    	$command = "mpileup -f $refFile -l $bed -C 0 -A -B -Q $mbq -q $mmq -d 50000 $bam";
    }
    else
    {
    	$command = "mpileup -f $refFile -l $bed -C 0 -A -B -Q $mbq -q $mmq -d 50000 $bam > $outdir/$mpileUpOutFile";
    }
    if((-e "$outdir/$mpileUpOutFile") and (-s "$outdir/$mpileUpOutFile" != 0))
    {
	$logger->info("$outdir/$mpileUpOutFile exists and mpileup will not be ran.\n");
	return("$outdir/$mpileUpOutFile");
    }
    else
    {
	eval{
		if(defined $QSUB){
			

	my $cmd = "$QSUB -q $queue -wd $outdir -V -N RunPileup.$bamId.$$ -o $mpileUpOutFile -e RunPileup.$bamId.$$.stderr -l h_vmem=5G,virtual_free=5G -pe smp 1 -b y $pathToSamtools $command";
	$logger->debug($cmd);
	`$QSUB -q $queue -wd $outdir -V -N RunPileup.$bamId.$$ -o $mpileUpOutFile -e RunPileup.$bamId.$$.stderr -l h_vmem=5G,virtual_free=5G -pe smp 1 -b y $pathToSamtools $command`;
	`$QSUB -q $queue -V -wd $outdir -hold_jid RunPileup.$bamId.$$ -N NotifyMpileup.$bamId.$$ -e /dev/null -o NotifyMpileup.$bamId.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
           	}
           	else
           	{
           		my $cmd = "$BSUB -q $queue -J RunPileup.$bamId.$$ -cwd $outdir -e RunPileup.$bamId.$$.%J.stderr -o RunPileup.$bamId.$$.%J.stdout -We 0:59 -R \"rusage[mem=5]\" -M 8 -n 1 \"$pathToSamtools $command\"";
           		$logger->debug($cmd);
           		`$BSUB -q $queue -J RunPileup.$bamId.$$ -cwd $outdir -e RunPileup.$bamId.$$.%J.stderr -o RunPileup.$bamId.$$.%J.stdout -We 0:59 -R "rusage[mem=5]" -M 8 -n 1 "$pathToSamtools $command"`;
				`$BSUB -q $queue -cwd $outdir -w "post_done(RunPileup.$bamId.$$)" -J NotifyMpileup.$bamId.$$ -e NotifyMpileup.$bamId.$$.%J.stderr -o NotifyMpileup.$bamId.$$.stat -R "rusage[mem=2]" -M 4 -n 1 "$outdir/Notify.csh"`;	
           	}
		};
	if($@){
	    $logger->fatal("RunMpileup: Job submission failed. Error: $@\n");
	}
	push (@notifyNames,"NotifyMpileup.$bamId.$$.stat");
    }

    &WaitToFinish($outdir,@notifyNames);
    $logger->info( "Completed running mpileup on bam file\n");
    return("$outdir/$mpileUpOutFile");
}
###################################################
###################################################
#--Analyze Mpileup results
sub AnalyzeMpileup
{
    my($mpileupFile,$mutationHash) =@_;
    my(%altHash) = %$mutationHash;
    tie (my %pileupHash, 'Tie::IxHash');
    tie (my %pileupstatsHash, 'Tie::IxHash');
    open(FH,"$mpileupFile") || die $logger->fatal("Can not open $mpileupFile. Error: $!\n");
    #Working on a mpileup file
    while(<FH>)
    {
	my($chr,$pos,$ref,$dp,$str) = split("\t",$_);
	chomp($chr,$pos,$ref,$dp,$str);
	tie (my %baseHash, 'Tie::IxHash');
	tie (my %totalHash, 'Tie::IxHash');
	tie (my %RefHash, 'Tie::IxHash');
	my $index = 0;
	my $prev_index = "";
	#working on each line of pileup containg results
	while($index < length($str))
	{
	    $prev_index = $index;
	    if(substr($str,$index,1) eq ".")
	    {
		my $base = uc($ref);
		$baseHash{$base}{"+"}++;
		$totalHash{"+"}++;
		$RefHash{"+"}++;
		$index++;
	    }
	    elsif(substr($str,$index,1) eq ",")
	    {
		my $base = uc($ref);
		$baseHash{$base}{"-"}++;
		$totalHash{"-"}++;
		$RefHash{"-"}++;
		$index++;
	    }
	    elsif(substr($str,$index,1)=~/[ACGTN]/i)
	    {
		my $base = uc(substr($str,$index,1));
		if($base eq substr($str,$index,1))
		{
		    $baseHash{$base}{"+"}++;
		    $totalHash{"+"}++;
		}
		else
		{
		    $baseHash{$base}{"-"}++;
		    $totalHash{"-"}++;
		}
		$index++;
	    }
	    elsif(substr($str,$index,1) eq "+" || substr($str,$index,1) eq "-")
	    {
		my $num = "";
		my $num_i = $index+1;
		while(substr($str,$num_i,1)=~/^\d$/)
		{
		    $num .= substr($str,$num_i,1);
		    $num_i++;
		}
		
		my $indel = substr($str,$index,$num+1+length($num));
		my $base = uc($indel);
		if($base eq $indel)
		{
		    $baseHash{$base}{"+"}++;
		}
		else
		{
		    $baseHash{$base}{"-"}++;
		}
		my $cr =  $baseHash{$base}{"-"};
		my $cf =  $baseHash{$base}{"+"};
		$index += length($indel);
	    }
	    elsif(substr($str,$index,1) eq "^")
	    {
		my $map_qual = ord(substr($str,$index+1,1))-33;
		$index += 2;
	    }
	    elsif(substr($str,$index,1)=~/^\$$/)
	    {
		$index++;
	    }
	    elsif(substr($str,$index,1)=~/^\*$/)
	    {
		$baseHash{"-"}{"+"}++; 
		$totalHash{"*"}++;
		$index++;
	    }
	    if($prev_index == $index)
	    {
		last;
	    }
	}
	my $MpileupChrStart = $chr.":".$pos."";
	my $sampleId = $altHash{$MpileupChrStart};
	#Calculate Values for Samples
	while(my($id,$value) = each (%$sampleId))
	{
	    my @stats=();
	    push @stats, "DP=${dp}";
	    my $tot_f = 0;
	    my $tot_r = 0;
	    my $reftot_f = 0;
	    my $reftot_r = 0;
	    my $tot_star = 0;
	    if(exists $totalHash{"+"})
	    {
		$tot_f = $totalHash{"+"};
	    }
	    if(exists $totalHash{"-"})
	    {
		$tot_r = $totalHash{"-"};
	    }
	    if(exists $totalHash{"*"})
	    {
		$tot_star = $totalHash{"*"};
	    }
	    if(exists $RefHash{"+"})
	    {
		$reftot_f = $RefHash{"+"};
	    }
	    if(exists $RefHash{"-"})
	    {
		$reftot_r = $RefHash{"-"};
	    }
	    my $tot = $tot_f + $tot_r + $tot_star;
	    my $reftot = $reftot_f + $reftot_r;
	    push (@stats, "${reftot},${reftot_f},${reftot_r}");	
	    my $flag=0;
	    foreach my $key (sort {$a cmp $b} keys %baseHash)
	    {
		next if ($key ne $value);
		$flag=1;
		my $cnt_f=0;
		if(exists $baseHash{$key}{"+"})
		{
			$cnt_f = $baseHash{$key}{"+"};
		}
		my $cnt_r=0;
		if(exists $baseHash{$key}{"-"})
		{
		    $cnt_r = $baseHash{$key}{"-"};
		}
		my $total = $cnt_f + $cnt_r;
		push @stats, "${key}=${total}\t${cnt_f}\t${cnt_r}";
		my $vf = sprintf("%0.5f",($total/$dp));
		push @stats, "VF=${vf}";
	    }
	    if($flag==0)
	    {
		    push @stats, $value."=0\t0\t0";
		    push @stats, "VF=0";
	    }
	    $pileupstatsHash{$chr.":".$pos}{$id}=join("; ",@stats);
	}
	
    }
    close(FH);
    return(\%pileupHash,\%pileupstatsHash);
}
###################################################
###################################################
#--Merge mpileup output with vcf
sub Merge_Vcf_mpileup
{
    my($orgVcfHash,$mutationHash,$mpileupHash,$mpileupstatsHash,$outdir,$outFile,$bamId) = @_;
    my(%vcfHash) = %$orgVcfHash;
    my(%altHash) = %$mutationHash;
    my(%pileupHash) = %$mpileupHash;
    my(%pileupstatsHash) = %$mpileupstatsHash;
    open (FH,">","$outdir/$outFile") || die $logger->fatal("Cannot open $outdir/$outFile,$!\n");
    print FH "Ref_BAM\tSample\tChrom\tPOS\tRef\tAlt\tTotal_Depth\tRef_Depth\tAlt_Counts\tAlt_Freq\tRef_Forward\tRef_Reverse\tAlt_Forward\tAlt_Reverse\n";
    foreach my $key (keys %vcfHash)
    {
	my($chr,$pos,$id,$ref,$alt,$qual,$filter,$info) = split("\t",$vcfHash{$key});
	chomp($chr,$pos,$id,$ref,$alt,$qual,$filter,$info);
	#my @infoField = split(";",$info);
	my($sampleinfo,$rest) = split(";",$info);
	my($name,$sampleId) = split("=",$sampleinfo);
	my ($tf,$tff,$tfr,$dp,$ad,$aff,$afr,$vf);
	my(@g) = ();
	if((abs(length($ref))) == (abs(length($alt))))
	{
	    if(exists $pileupstatsHash{$chr.":".$pos}{$sampleId . ":" . $alt . ":" . $ref})
	    {
		($dp)=$pileupstatsHash{$chr.":".$pos}{$sampleId . ":" . $alt . ":" . $ref}=~/DP\=(\d+)/;
		@g=split('\; ',$pileupstatsHash{$chr.":".$pos}{$sampleId . ":" . $alt . ":" . $ref});
		($tf,$tff,$tfr) = split(",",$g[1]);
		my @h=split('\=',$g[2]);
		my $var = $h[0];
		($ad,$aff,$afr) = split("\t",$h[1]);
		$vf = 0;
		$vf = sprintf("%.5f",($ad/$dp)) if ($dp != 0);
	    }
	    else
	    {
		($dp,$tf,$tff,$tfr,$ad,$aff,$afr,$vf) = 0
	    }
	}
	else
	{
	    my $pos1 = $pos+1;
	    if(! exists $pileupstatsHash{$chr.":".$pos1}{$sampleId . ":" . $alt . ":" . $ref})
	    {
                if (!($typeOfSample =~ m/NTC/i)) {
                    $logger->fatal( "$bamId:Position ahead for $key did not report pileup stats");
                }
		($dp,$tf,$tff,$tfr,$ad,$aff,$afr,$vf) = 0;
	    }
	    else
	    {
		#my $pos1record = $pileupstatsHash{$chr.":".$pos1}{$sampleId};
		#my $posrecord = $pileupstatsHash{$chr.":".$pos}{$sampleId};
		#print "Pos1Indel:$pos1record\tPosIndel:$posrecord\n";
		($dp) = $pileupstatsHash{$chr.":".$pos1}{$sampleId . ":" . $alt . ":" . $ref}=~/DP\=(\d+)/;
		@g=split('\; ',$pileupstatsHash{$chr.":".$pos1}{$sampleId . ":" .$alt . ":" . $ref});
		($tf,$tff,$tfr) = split(",",$g[1]);
		@g=split('\; ',$pileupstatsHash{$chr.":".$pos}{$sampleId . ":" .$alt . ":" . $ref});
		my @h=split('\=',$g[2]);
		my $var = $h[0];
		($ad,$aff,$afr) = split("\t",$h[1]);
		
		if(length($ref) < length($alt)){ # fix added 20140108 by AZ. For insertions, RefDepth is equal to TotalDepth. However, it should be TotalDepth-AltDepth. 
			$tf = $tf-$ad;
		}
		
		$vf = 0;
		$vf = sprintf("%.5f",($ad/$dp)) if ($dp != 0);
	    } 
	}
	
	print FH "$bamId\t$sampleId\t$chr\t$pos\t$ref\t$alt\t$dp\t$tf\t$ad\t$vf\t$tff\t$tfr\t$aff\t$afr\n";
	
    }
    close(FH);
    return($outFile);
}

###################################################
###################################################
#--Delete Unwanted Files
sub HouseKeeping
{
    my($outdir,@list) = @_;
    $logger->info("Removing Unwanted files..\n");
    foreach my  $file (@list)
    {
	if($file =~ /\//)
	{
	    eval{`rm $file`;};
	    if($@){
		$logger->fatal("Can not remove $file,Error:$@");
	    }
	}
	else
	{
	    eval{`rm $outdir/$file`;};
            if($@){
                $logger->fatal("Can not remove $outdir/$file,Error:$@");
            }
	}
    }
    return;
}
