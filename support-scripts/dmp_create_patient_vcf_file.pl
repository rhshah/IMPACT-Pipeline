#!/usr/bin/perl -w
# dmp_create_patient_vcf_file.pl ---
# Author: Zehir <zehira@phoenix-h1>
# Created: 26 Nov 2013
# Version: 0.01

use warnings;
use strict;
use Getopt::Long;
use IO::File;
use Cwd;
use Tie::IxHash;
use MSKCC_DMP_Logger;

my $logger = MSKCC_DMP_Logger->get_logger('ANNOTATE_FILTER_VARIANTS');
$logger->start_local();

my ( $filteredAnnoFile, $outdir, $titleFile, $java, $IGVtools );

if (
	@ARGV < 1
	or !GetOptions(
		'annotatedFilteredVariants|v:s' => \$filteredAnnoFile,
		'outdir|o:s'                    => \$outdir,
		'TitleFile|t:s'                 => \$titleFile,
		'JAVA|j:s'                      => \$java,
		'IGVTools|it:s'                 => \$IGVtools

	)
  )
{
	Usage();
}

if ( ( !$filteredAnnoFile ) or ( !$titleFile ) ) {
	$logger->fatal("Mutation File is missing") if ( !$filteredAnnoFile );
	$logger->fatal("Title File is missing")    if ( !$titleFile );
	Usage();
	exit;
}

my ( %samples, @header );
if ( !-e "$outdir/VCFs" ) {
	`mkdir VCFs`;
}

#Get Title file information
my ( $patientIDPerSampleId, $classPerPatientIdSampleId ) = &ReadTitleFile( $titleFile, $outdir );

$logger->info("Individual VCF files are being generated......");
my %patientIDPerSampleId      = %$patientIDPerSampleId;
my %classPerPatientIdSampleId = %$classPerPatientIdSampleId;

#print "$filteredAnnoFile\n";
open IN, "<", $filteredAnnoFile || die $logger->fatal("Can not open $filteredAnnoFile: $!");
while (<IN>) {
	chomp;
	if (m/^Sample\t/) {
		@header = split("\t");
		next;
	}
	my @line   = split("\t");
	my $sample = $line[0];
	my $normal = $line[1];
	my $chr    = $line[2];
	my $pos    = $line[3];
	my $ref    = $line[4];
	my $alt    = $line[5];
	$samples{ $sample . ":" . $normal }{ $chr . ":" . $pos . ":" . $ref . ":" . $alt } = $_;
}
close IN;

foreach my $key1 ( sort keys %samples ) {
	my ( $sample, $normal ) = split( ":", $key1 );
	my $normalIndex;
	foreach my $i ( 0 .. $#header ) {
		if ( $header[$i] eq $normal ) {
			$normalIndex = $i;
		}
	}

	#determine if the sample is matched or unmatched
	my $patientID = $patientIDPerSampleId{$sample};
	my $key       = $normal . ":" . $patientID;
	my $tag;
	if ( ( exists $classPerPatientIdSampleId{$key} ) and ( $classPerPatientIdSampleId{$key} =~ /Normal/ ) ) {
		$tag = "Matched";
	}
	else {
		$tag = "Unmatched";
	}

	#print "$sample\t$normal\t$tag\t$normalIndex\t$header[$normalIndex]\n";
	my $vcfFile = $sample . "_" . $tag . "_annovar.vcf";
	open OUT, ">", "$outdir/VCFs/$vcfFile" || die $logger->fatal("Can not open $vcfFile: $!");
	my $header = <<"END";
##fileformat=VCFv4.1
##source=MuTect,SomaticIndelDetector
##FILTER=<ID=VF_T,Description="VF tumor threshold: 0.02">
##FILTER=<ID=VF_T,Description="VF_T/VF_N threshold: 5">
##FORMAT=<ID=DP,Number=3,Type=Integer,Description="Approximate read depth,followed by reference and alternative depths, genotyped">
##FORMAT=<ID=VF,Number=1,Type=Integer,Description="Variant frequency determined after genotyping">
##INFO=<ID=STRBIAS,Number=2,Type=Integer,Description="Number of reads aligned to positive and negative strand after genotyping">
##INFO=<ID=DP_ALL_N,Number=1,Type=Integer,Description="Aggregate of total coverage across all normals and alternative allele depth">
##INFO=<ID=VF_MEDIAN_N,Number=1,Type=Integer,Description="Median of variant frequencies across all normals sequenced">
##INFO=<ID=TN_VF_RATIO,Number=1,Type=Integer,Description="Ratio of variant frequency in tumor sample to median of normals">
##INFO=<ID=OCCURANCE_N,Number=1,Type=Integer,Description="Number of times the variant is observed in the normals">
##INFO=<ID=ANNOVAR,Number=1,Type=String,Description="Annotation from Annovar. GeneName, transcript ID, variant type, cDNA change, AA change, dbSNPid, cosmic ID, 1000G minor allele freq (MAF)">
##INFO=<ID=FAILURE_REASON,Number=1,Type=String,Description="Failure reason field from Mutect">
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
##contig=<ID=3,length=198022430>
##contig=<ID=4,length=191154276>
##contig=<ID=5,length=180915260>
##contig=<ID=6,length=171115067>
##contig=<ID=7,length=159138663>
##contig=<ID=8,length=146364022>
##contig=<ID=9,length=141213431>
##contig=<ID=10,length=135534747>
##contig=<ID=11,length=135006516>
##contig=<ID=12,length=133851895>
##contig=<ID=13,length=115169878>
##contig=<ID=14,length=107349540>
##contig=<ID=15,length=102531392>
##contig=<ID=16,length=90354753>
##contig=<ID=17,length=81195210>
##contig=<ID=18,length=78077248>
##contig=<ID=19,length=59128983>
##contig=<ID=20,length=63025520>
##contig=<ID=21,length=48129895>
##contig=<ID=22,length=51304566>
##contig=<ID=X,length=155270560>
##contig=<ID=Y,length=59373566>
##contig=<ID=MT,length=16569>
##contig=<ID=GL000207.1,length=4262>
##contig=<ID=GL000226.1,length=15008>
##contig=<ID=GL000229.1,length=19913>
##contig=<ID=GL000231.1,length=27386>
##contig=<ID=GL000210.1,length=27682>
##contig=<ID=GL000239.1,length=33824>
##contig=<ID=GL000235.1,length=34474>
##contig=<ID=GL000201.1,length=36148>
##contig=<ID=GL000247.1,length=36422>
##contig=<ID=GL000245.1,length=36651>
##contig=<ID=GL000197.1,length=37175>
##contig=<ID=GL000203.1,length=37498>
##contig=<ID=GL000246.1,length=38154>
##contig=<ID=GL000249.1,length=38502>
##contig=<ID=GL000196.1,length=38914>
##contig=<ID=GL000248.1,length=39786>
##contig=<ID=GL000244.1,length=39929>
##contig=<ID=GL000238.1,length=39939>
##contig=<ID=GL000202.1,length=40103>
##contig=<ID=GL000234.1,length=40531>
##contig=<ID=GL000232.1,length=40652>
##contig=<ID=GL000206.1,length=41001>
##contig=<ID=GL000240.1,length=41933>
##contig=<ID=GL000236.1,length=41934>
##contig=<ID=GL000241.1,length=42152>
##contig=<ID=GL000243.1,length=43341>
##contig=<ID=GL000242.1,length=43523>
##contig=<ID=GL000230.1,length=43691>
##contig=<ID=GL000237.1,length=45867>
##contig=<ID=GL000233.1,length=45941>
##contig=<ID=GL000204.1,length=81310>
##contig=<ID=GL000198.1,length=90085>
##contig=<ID=GL000208.1,length=92689>
##contig=<ID=GL000191.1,length=106433>
##contig=<ID=GL000227.1,length=128374>
##contig=<ID=GL000228.1,length=129120>
##contig=<ID=GL000214.1,length=137718>
##contig=<ID=GL000221.1,length=155397>
##contig=<ID=GL000209.1,length=159169>
##contig=<ID=GL000218.1,length=161147>
##contig=<ID=GL000220.1,length=161802>
##contig=<ID=GL000213.1,length=164239>
##contig=<ID=GL000211.1,length=166566>
##contig=<ID=GL000199.1,length=169874>
##contig=<ID=GL000217.1,length=172149>
##contig=<ID=GL000216.1,length=172294>
##contig=<ID=GL000215.1,length=172545>
##contig=<ID=GL000205.1,length=174588>
##contig=<ID=GL000219.1,length=179198>
##contig=<ID=GL000224.1,length=179693>
##contig=<ID=GL000223.1,length=180455>
##contig=<ID=GL000195.1,length=182896>
##contig=<ID=GL000212.1,length=186858>
##contig=<ID=GL000222.1,length=186861>
##contig=<ID=GL000200.1,length=187035>
##contig=<ID=GL000193.1,length=189789>
##contig=<ID=GL000194.1,length=191469>
##contig=<ID=GL000225.1,length=211173>
##contig=<ID=GL000192.1,length=547496>
##contig=<ID=NC_007605,length=171823>
##reference=file:///home/shahr2/References/Homo_sapiens_assembly19.fasta
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample\t$normal
END
	print OUT "$header";

	foreach my $key2 ( sort my_sort keys %{ $samples{$key1} } ) {
		my @line       = split( "\t", $samples{$key1}->{$key2} );
		#@line = sort my_sort2 (@line1);
		my $chr        = $line[2];
		my $pos        = $line[3];
		my $ref        = $line[4];
		my $alt        = $line[5];
		my $varClass   = $line[6];
		my $gene       = $line[7];
		my $exon       = $line[8];
		my $txID       = $line[11];
		my $cDNAchange = $line[12];
		my $AAchange   = $line[13];
		my $dbSNP      = $line[14];
		if ( $dbSNP eq "" ) {
			$dbSNP = ".";
		}
		my $cosmic      = $line[15];
		my $MAF         = $line[16];
		my $failReason  = $line[22];
		my $DP_N        = $line[23];
		my $DP_N_Ref    = $line[24];
		my $DP_N_Alt    = $line[25];
		my $VF_N        = $line[26];
		my $DP_T        = $line[27];
		my $DP_T_Ref    = $line[28];
		my $DP_T_Alt    = $line[29];
		my $VF_T        = $line[30];
		my $STR_T_Ref_p = $line[31];
		my $STR_T_Ref_n = $line[32];
		my $STR_T_Alt_p = $line[33];
		my $STR_T_Alt_n = $line[34];
		my $DP_N_AGG    = $line[35];
		my $VF_N_MEDIAN = $line[36];
		my $TNfreqRatio = $line[37];
		my $occurance   = $line[38];
		my $normalInfo  = $line[$normalIndex];
		print OUT
"$chr\t$pos\t$dbSNP\t$ref\t$alt\t.\t.\t.\tDP_N=$DP_N($DP_N_Ref,$DP_N_Alt);VF_N=$VF_N;DP_T=$DP_T($DP_T_Ref,$DP_T_Alt);VF_T=$VF_T;ANNOVAR=($gene|$txID|$varClass|$exon|$cDNAchange|$AAchange|$dbSNP|$MAF|$cosmic);STRBIAS_REF=($STR_T_Ref_p,$STR_T_Ref_n);STRBIAS_ALT=($STR_T_Alt_p,$STR_T_Alt_n);VF_N_MEDIAN=$VF_N_MEDIAN;TN_VF_RATIO=$TNfreqRatio;OCCURANCE_N=$occurance\t$normalInfo\n";
	}
	close OUT;
	eval { `$java -Xmx2g -jar $IGVtools sort $outdir/VCFs/$vcfFile $outdir/VCFs/$vcfFile.2` };
	if ($@) {
		$logger->fatal("Command: $java -Xmx2g -jar $IGVtools sort $outdir/VCFs/$vcfFile $outdir/VCFs/$vcfFile.2 can not be run! Error: $@");
		exit(1);
	}
	else {
		$logger->info("Command: $java -Xmx2g -jar $IGVtools sort $outdir/VCFs/$vcfFile $outdir/VCFs/$vcfFile.2 is being run!");

	}
	eval { `mv $outdir/VCFs/$vcfFile.2 $outdir/VCFs/$vcfFile` };
	if ($@) {
		$logger->fatal("Command: mv $outdir/VCFs/$vcfFile.2 $outdir/VCFs/$vcfFile can not be run! Error: $@");
		exit(1);
	}
	else {
		$logger->info("Command: mv $outdir/VCFs/$vcfFile.2 $outdir/VCFs/$vcfFile is being run!");

	}
	eval { `$java -Xmx2g -jar $IGVtools index $outdir/VCFs/$vcfFile` };
	if ($@) {
		$logger->fatal("Command: $java -Xmx2g -jar $IGVtools index $outdir/VCFs/$vcfFile can not be run! Error: $@");
		exit(1);
	}
	else {
		$logger->info("Command: $java -Xmx2g -jar $IGVtools index $outdir/VCFs/$vcfFile is being run!");

	}
}
$logger->info("VCF file generation has been completed.");

#####################################
#####################################
#How to use the script.
sub Usage {
	print "Unknow option: @_\n" if (@_);

	print "\nUsage : dmp_create_patient_vcf_file.pl [options]
        [--annotatedFilteredVariants|v S File containing annotated and filtered mutations (required and submit with full path,Ex:/SomePath/CytoOv_SomaticMutIndel.txt)]
        [--titleFile|t                 S tab-delimited title file for the samples (required and submit with full path)]
        [--outdir|o                    S Path where all the output files will be written (optional) [default:cwd]]
        [--JAVA|j						S Path where java is accessible]
		[--IGVTools|it					S Path to IGV tools]
        
        \n";

	exit;
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
	tie( my %patientIDPerSampleID,      'Tie::IxHash' );
	tie( my %classPerPatientIDSampleID, 'Tie::IxHash' );
	open( TFH, $titleFile )
	  || die $logger->fatal("Can not open $titleFile, $!");

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
	for ( my $i = 0 ; $i < scalar(@barcode) ; $i++ ) {

		#print "$[$i] => $class[$i]\n";
		$classPerPatientIDSampleID{ $sampleId[$i] . ":" . $patientId[$i] } = $class[$i];
		$patientIDPerSampleID{ $sampleId[$i] } = $patientId[$i];
	}
	return ( \%patientIDPerSampleID, \%classPerPatientIDSampleID );

}



####################################
# chromosome sort
sub my_sort{
    my ($chr_a, $t, $y, $z) = split(":", $a);
    my ($chr_b, $r, $e, $w) = split(":", $b);
    if ($chr_a eq "X") {$chr_a = 23;}
    if ($chr_a eq "Y") {$chr_a = 24;}
    if ($chr_a eq "M") {$chr_a = 25;}
    if ($chr_b eq "X") {$chr_b = 23;}
    if ($chr_b eq "Y") {$chr_b = 24;}
    if ($chr_b eq "M") {$chr_b = 25;}
    return ($chr_a <=> $chr_b || $t <=> $r);
}





__END__

=head1 NAME

dmp_create_patient_vcf_file.pl - Describe the usage of script briefly

=head1 SYNOPSIS

dmp_create_patient_vcf_file.pl [options] args

      -opt --long      Option description

=head1 DESCRIPTION

Stub documentation for dmp_create_patient_vcf_file.pl, 

=head1 AUTHOR

Zehir, E<lt>zehira@phoenix-h1E<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 by Zehir

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
