#!/usr/bin/perl -w
# dmp_annotate_filter_variants.pl --- Given a file of variants:
# 1) Annotates the variants, picking the variant annotations within the canonical exons
# 2) Filters the variants based on the variant type and MAF value
# Author: Zehir <zehira@phoenix-h1>
# Created: 19 Nov 2013
# Version: 0.01

use warnings;
use strict;
use Getopt::Long;
use IO::File;
use Cwd;
use Tie::IxHash;
use MSKCC_DMP_Logger;

my ( $mutationFile, $outdir, $deleteFiles, $titleFile, $configFile, $MAFthreshold );
my $logger = MSKCC_DMP_Logger->get_logger('ANNOTATE_FILTER_VARIANTS');
$logger->start_local();

if (
	@ARGV < 1
	or !GetOptions(
		'somaticMutIndelFile|si:s' => \$mutationFile,
		'configurationFile|c:s'    => \$configFile,
		'outdir|o:s'               => \$outdir,
		'deleteUnwantedFiles|d:i'  => \$deleteFiles,
		'TitleFile|t:s'            => \$titleFile,

	)
  )
{
	Usage();
}

if ( ( !$mutationFile ) or ( !$titleFile ) or ( !$configFile ) ) {
	$logger->fatal("Mutation File is missing")      if ( !$mutationFile );
	$logger->fatal("Title File is missing")         if ( !$titleFile );
	$logger->fatal("Configuration File is missing") if ( !$configFile );
	Usage();
	exit;
}
if ( !$deleteFiles ) {
	$logger->warn( "The option to delete intermediary files was not set, the files will be deleted by default." );
	$deleteFiles = 2;
}

if ( !$outdir ) {
	$outdir = getcwd;

}

if ( !$deleteFiles ) {
	$logger->warn( "The option to delete intermediary files was not set, the files will be deleted by default" );
	$deleteFiles = 1;
}
if ( !$outdir ) {
	$outdir = getcwd;

}
$logger->info("Variant annotation and filtering has started for $mutationFile");

#Get Title file information
my ( $patientIDPerSampleId, $classPerPatientIdSampleId ) = &ReadTitleFile( $titleFile, $outdir );

#print "MutationFile: $mutationFile\nTitle File: $titleFile\nConfigFile: $configFile\nOutdir: $outdir\n";
my ( $locationRef, $versionRef, $parameterRef ) = &GetConfiguration($configFile);
my %location  = %$locationRef;
my %version   = %$versionRef;
my %parameter = %$parameterRef;

#import file locations
my $annovar          = $location{"Annovar"};
my $annovarDB        = $location{"Annovar_db"};
my $canonicalTxFile  = $location{"Canonical_refFlat_file"};
my $translatedFolder = $location{"TranslationFolder"};
my $geneIntervalAnn = $location{"GeneIntervalAnn"};

#Run annovar annotation

my $mergedFile = $mutationFile;
my $filesToDelete;

$mergedFile =~ s/_withAlleleDepth\.txt/_withAlleleDepth_annovarAnnotated.txt/g;

if ( !-s "$mergedFile" ) {
	&RunAnnovar( $mutationFile, $outdir );

	# Merge annovar annotations with mutation file
	( $mergedFile, $filesToDelete ) = &MergeAnnotations( $mutationFile, $outdir, $canonicalTxFile, $translatedFolder );
	my @deleteFiles = @$filesToDelete;

	#deletion process
	if ( $deleteFiles == 1 ) {
		&HouseKeeping( $outdir, @deleteFiles );
	}
	if ( $deleteFiles == 2 ) {
		$logger->info("No house keeping performed");
	}
}
&FilterAnnotatedVariants( $mergedFile, $outdir, $patientIDPerSampleId, $classPerPatientIdSampleId, $geneIntervalAnn );

#####################################
#####################################
sub FilterAnnotatedVariants {
	my ( $input, $outdir, $patientIDPerSampleId, $classPerPatientIdSampleId, $geneIntervalAnn  ) = @_;
	if ( $input =~ /\// ) {
		my @f = split( '\/', $input );
		$input = $f[-1];
	}
	my $outputExonic = $input;
	$outputExonic =~ s/annovarAnnotated/annovarAnnotatedExonic/;
	my $outputNonPanelExonic = $input;
	$outputNonPanelExonic =~ s/annovarAnnotated/annovarAnnotatedNonPanelExonic/;
	my $outputDropped = $input;
	$outputDropped =~ s/annovarAnnotated/annovarAnnotatedNonExonic/;
	my $outputSilent = $input;
	$outputSilent =~ s/annovarAnnotated/annovarAnnotatedSilent/;
	my $outputNonPanelSilent = $input;
	$outputNonPanelSilent =~ s/annovarAnnotated/annovarAnnotatedNonPanelSilent/;	
	my %patientIDPerSampleId      = %$patientIDPerSampleId;
	my %classPerPatientIdSampleId = %$classPerPatientIdSampleId;

	my %geneList;	

	open EXONIC, ">", "${outdir}/${outputExonic}" || die $logger->fatal("Can not create $outputExonic file: $!");
	open NP_EXONIC, ">", "${outdir}/${outputNonPanelExonic}" || die $logger->fatal("Can not create $outputNonPanelExonic file: $!");

	#open DROPPED, ">", "$outdir/$outputDropped" || die $logger->fatal("Can not create $outputDropped file: $!;");
	open SILENT, ">", "${outdir}/${outputSilent}" || die $logger->fatal("Can not create $outputSilent file: $!;");
	open NP_SILENT, ">", "${outdir}/${outputNonPanelSilent}" || die $logger->fatal("Can not create $outputNonPanelSilent file: $!;");

	open IN, "<", "${input}" || die $logger->info("Filtering annotated variants file.");
	open INTERVAL, "<", "${geneIntervalAnn}" or die $logger->info("Filtering step annotated gene intervals file.");

	while(<INTERVAL>){
		chomp;
		my @temp = split("\t", $_);
		my @ann = split(":", $temp[2]);
		my $gene = $ann[0];

		$geneList{$gene} = 1;

	}

	while (<IN>) {
		chomp;
		our @header;
		our ( $varTypeIndex, $geneIndex, $MAFindex );
		if (/^Sample/) {
			print EXONIC "$_\n";
			print NP_EXONIC "$_\n";

			#print DROPPED "$_\n";
			print SILENT "$_\n";
			print NP_SILENT "$_\n";

			@header = split("\t");
			foreach my $i ( 0 .. scalar(@header) - 1 ) {    #find the index for these columns

				#print "$i\t$header[$i]\n";
				if ( $header[$i] eq "VariantClass" ) {
					$varTypeIndex = $i;
				}
				elsif ( $header[$i] eq "Gene" ) {
					$geneIndex = $i;
				}

			}
			next;
		}
		my @line         = split("\t");
		my $sample       = $line[0];
		my $normal       = $line[1];
		my $gene         = $line[$geneIndex];
		my $variantClass = $line[$varTypeIndex];

		if(exists $geneList{$gene}){
			if (   $variantClass eq "frameshift_deletion"
					or $variantClass eq "frameshift_insertion"
					or $variantClass eq "nonframeshift_deletion"
					or $variantClass eq "nonframeshift_insertion"
					or $variantClass eq "nonsynonymous_SNV"
					or $variantClass eq "splicing"
					or $variantClass eq "stopgain_SNV"
                                        or $variantClass eq "stoploss_SNV")
				{
					print EXONIC "$_\n";
				}
				else {
					if ( $gene eq "TERT" && $variantClass eq "upstream" ) {
						print EXONIC "$_\n";
					}
					else {
						print SILENT "$_\n";
					}
				}
		}
		else{
			if (   $variantClass eq "frameshift_deletion"
					or $variantClass eq "frameshift_insertion"
					or $variantClass eq "nonframeshift_deletion"
					or $variantClass eq "nonframeshift_insertion"
					or $variantClass eq "nonsynonymous_SNV"
					or $variantClass eq "splicing"
					or $variantClass eq "stopgain_SNV"
                                        or $variantClass eq "stoploss_SNV")
				{
					print NP_EXONIC "$_\n";
				}
			else{
				print NP_SILENT "$_\n";
			}
		}
=begin		
		if ( $gene eq "TERT" && $variantClass eq "upstream" ) {
			print EXONIC "$_\n";
			next;
		}
		if (   $variantClass eq "frameshift_deletion"
			|| $variantClass eq "frameshift_insertion"
			|| $variantClass eq "nonframeshift_deletion"
			|| $variantClass eq "nonframeshift_insertion"
			|| $variantClass eq "nonsynonymous_SNV"
			|| $variantClass eq "splicing"
			|| $variantClass eq "stopgain_SNV"
                        || $variantClass eq "stoploss_SNV")
		{
			print EXONIC "$_\n";
		}
		elsif ( $variantClass eq "synonymous_SNV" ) {
			print SILENT "$_\n";
		}
		else {
			print SILENT "$_\n";
		}
=cut

	}
	close IN;
	close INTERVAL;
	close EXONIC;
	close NP_EXONIC;
	close SILENT;
	close NP_SILENT;
}
#####################################
#####################################
#How to use the script.
sub Usage {
	print "Unknow option: @_\n" if (@_);

	print "\nUsage : AnnotateAssessFilterVarinats.pl [options]
        [--SomaticMutIndelFile|si          S File containing mutations (required and submit with full path,Ex:/SomePath/CytoOv_SomaticMutIndel.txt)]
        [--ConfigurationFile|c             S Configuration file that contains the locations for the programs and the databases (required and submit with full path)]
        [--titleFile|t                     S tab-delimited title file for the samples (required and submit with full path)]
        [--outdir|o                        S Path where all the output files will be written (optional) [default:cwd]]
        [--deleteUnwantedFiles|d           I 2=>To delete files 1=> To keep files (default:2,optional)]
        
        \n";

	exit;
}
##############################
##############################
# Run annovar on the mutation file

sub RunAnnovar {
	my ( $sample, $outDir ) = @_;
	my ($annInput) = $sample =~ /.*\/(.*).txt/;
	$logger->info("Annovar annotation has started....");
	if ( ( -e "$outDir/$annInput" ) and ( ( -s "$outDir/$annInput" ) != 0 ) ) {
		$logger->info( "Annovar input file: $annInput already exists and this process will not create a new one" );
	}
	else {
		$logger->info("Creating the annovar input file");
		open( IN, "<", "$sample" )
		  or die $logger->fatal("Can not open mutation file: $!");
		open ANN, ">", "$outDir/$annInput"
		  || die $logger->fatal("Can not creat $annInput: $!");
		while (<IN>) {
			chomp;

			#M-1608-42-T	1	9784916	C	T	KEEP
			next if ( $_ =~ /^Sample\t/ );
			my @line   = split("\t");
			my $sample = shift @line;
			my $normal = shift @line;
			my $chr    = shift @line;
			my $pos    = shift @line;
			my $ref    = shift @line;
			my $alt    = shift @line;
			my $at     = join( "|", @line );
			if ( length($ref) > length($alt) ) {    # deletion 2002 CGG C need become 2003 2004 GG -
				my $newref = substr( $ref, 1 );
				my $posEnd = $pos + length($newref);
				$pos++;
				print ANN "$chr\t$pos\t$posEnd\t$newref\t-\t$sample|$normal|$at\n";
			}
			elsif ( length($ref) < length($alt) ) {    #insertion
				my $newalt = substr( $alt, 1 );
				my $posEnd = $pos;
				print ANN "$chr\t$pos\t$posEnd\t-\t$newalt\t$sample|$normal|$at\n";
			}
			else {
				print ANN "$chr\t$pos\t$pos\t$ref\t$alt\t$sample|$normal|$at\n";
			}
		}
		close ANN;
		close IN;
	}

	if (   ( -e "$outDir/$annInput.exonic_variant_function" )
		&& ( ( -s "$outDir/$annInput.exonic_variant_function" ) != 0 ) )
	{
		$logger->warn( "Exonic variant annotation file exists. Exonic annotation will not be performed." );
	}
	else {
		eval {
			$logger->info( "Command: $annovar --geneanno -buildver hg19 $outDir/$annInput $annovarDB is being run" );
			`$annovar --geneanno -buildver hg19  $outDir/$annInput $annovarDB`;
		};
		if ($@) {
			$logger->fatal( "Command: $annovar --geneanno -buildver hg19 $outDir/$annInput $annovarDB can not be run! Error: $@" );
			exit(1);
		}

	}
	if (   ( -e "$outDir/$annInput.hg19_snp137NonFlagged_filtered" )
		&& ( ( -s "$outDir/$annInput.hg19_snp137NonFlagged_filtered" ) != 0 ) )
	{
		$logger->warn( "dnSNP variant annotation file exists. dbSNP  annotation will not be performed." );
	}
	else {
		eval {
			$logger->info( "Command: $annovar -filter -dbtype snp137NonFlagged -buildver hg19 $outDir/$annInput $annovarDB is being run!" );
			`$annovar -filter -dbtype snp137NonFlagged -buildver hg19 $outDir/$annInput $annovarDB`;
		};
		if ($@) {
			$logger->fatal( "Command: $annovar -filter -dbtype snp137NonFlagged -buildver hg19 $outDir/$annInput $annovarDB can not be run! Error: $@" );
			exit(1);
		}
	}

	if (   ( -e "$outDir/$annInput.hg19_cosmic68_filtered" )
		&& ( ( -s "$outDir/$annInput.hg19_cosmic68_filtered" ) != 0 ) )
	{
		$logger->warn( "Cosmic variant annotation file exists. Cosmic annotation will not be performed." );
	}
	else {
		eval {
			$logger->info( "Command: $annovar -filter -dbtype cosmic68 -buildver hg19 $outDir/$annInput $annovarDB is being run!" );
			`$annovar -filter -dbtype cosmic68 -buildver hg19 $outDir/$annInput $annovarDB`;
		};
		if ($@) {
			$logger->fatal( "Command: $annovar -filter -dbtype cosmic68 -buildver hg19 $outDir/$annInput $annovarDB can not be run! Error: $@" );
			exit(1);
		}
	}

	if (   ( -e "$outDir/$annInput.hg19_ALL.sites.2012_04_filtered" )
		&& ( ( -s "$outDir/$annInput.hg19_ALL.sites.2012_04_filtered" ) != 0 ) )
	{
		$logger->warn( "1000G MAF variant annotation file exists. 1000G MAF annotation will not be performed." );
	}
	else {
		eval {
			$logger->info( "Command: $annovar -filter -dbtype 1000g2012apr_all -buildver hg19 $outDir/$annInput $annovarDB is being run!" );
			`$annovar -filter -dbtype 1000g2012apr_all -buildver hg19 $outDir/$annInput $annovarDB`;
		};
		if ($@) {
			$logger->fatal( "Command: $annovar -filter -dbtype 1000g2012apr_all -buildver hg19 $outDir/$annInput $annovarDB can not be run! Error: $@" );
			exit(1);
		}
	}

	#	if (   ( -e "$outDir/$annInput.hg19_ljb2_gerp++_filtered" )
	#		&& ( ( -s "$outDir/$annInput.hg19_ljb2_gerp++_filtered" ) != 0 ) )
	#	{
	#		$logger->warn("GERP++ database variant annotation file exists. GERP++ database annotation will not be performed.");
	#	}
	#	else {
	#		eval {
	#			$logger->info("Command: $annovar -filter -dbtype ljb2_gerp++ -buildver hg19 $outDir/$annInput $annovarDB is being run!");
	#			`$annovar -filter -dbtype ljb2_gerp++ -buildver hg19 $outDir/$annInput $annovarDB`;
	#		};
	#		if ($@) {
	#			$logger->fatal("Command: $annovar -filter -dbtype ljb2_gerp++ -buildver hg19 $outDir/$annInput $annovarDB can not be run Error: $@");
	#			exit(1);
	#		}
	#	}
	#
	#	if (   ( -e "$outDir/$annInput.hg19_ljb2_lrt_filtered" )
	#		&& ( ( -s "$outDir/$annInput.hg19_ljb2_lrt_filtered" ) != 0 ) )
	#	{
	#		$logger->warn("LRT database variant annotation file exists. LRT database annotation will not be performed.");
	#	}
	#	else {
	#		eval {
	#			$logger->info("Command: $annovar -filter -dbtype ljb2_lrt -buildver hg19 $outDir/$annInput $annovarDB is being run!");
	#			`$annovar -filter -dbtype ljb2_lrt -buildver hg19 $outDir/$annInput $annovarDB`;
	#		};
	#		if ($@) {
	#			$logger->fatal("Command: $annovar -filter -dbtype ljb2_lrt -buildver hg19 $outDir/$annInput $annovarDB can not be run! Error: $@");
	#			exit(1);
	#		}
	#	}
	#
	#	if (   ( -e "$outDir/$annInput.hg19_ljb2_phylop_filtered" )
	#		&& ( ( -s "$outDir/$annInput.hg19_ljb2_phylop_filtered" ) != 0 ) )
	#	{
	#		$logger->warn("PhyloP database variant annotation file exists. PhyloP database annotation will not be performed.");
	#	}
	#	else {
	#		eval {
	#			$logger->info("Command: $annovar -filter -dbtype ljb2_phylop -buildver hg19 $outDir/$annInput $annovarDB is being run!");
	#			`$annovar -filter -dbtype ljb2_phylop -buildver hg19 $outDir/$annInput $annovarDB`;
	#		};
	#		if ($@) {
	#			$logger->fatal("Command: $annovar -filter -dbtype ljb2_phylop -buildver hg19 $outDir/$annInput $annovarDB can not be run! Error: $@");
	#			exit(1);
	#		}
	#	}
	#
	#	if (   ( -e "$outDir/$annInput.hg19_ljb2_ma_filtered" )
	#		&& ( ( -s "$outDir/$annInput.hg19_ljb2_ma_filtered" ) != 0 ) )
	#	{
	#		$logger->warn("Mutation Assesor database variant annotation file exists. Mutation Assesor database annotation will not be performed.");
	#	}
	#	else {
	#		eval {
	#			$logger->info("Command: $annovar -filter -dbtype ljb2_ma -buildver hg19 $outDir/$annInput $annovarDB is being run!");
	#			`$annovar -filter -dbtype ljb2_ma -buildver hg19 $outDir/$annInput $annovarDB`;
	#		};
	#		if ($@) {
	#			$logger->fatal("Command: $annovar -filter -dbtype ljb2_ma -buildver hg19 $outDir/$annInput $annovarDB can not be run! Error: $@");
	#			exit(1);
	#		}
	#	}
	#	if (   ( -e "$outDir/$annInput.hg19_avsift_filtered" )
	#		&& ( ( -s "$outDir/$annInput.hg19_avsift_filtered" ) != 0 ) )
	#	{
	#		$logger->warn("AVSIFT database variant annotation file exists. AVSIFT database annotation will not be performed.");
	#	}
	#	else {
	#		eval {
	#			$logger->info("Comand: $annovar -filter -dbtype avsift -buildver hg19 $outDir/$annInput $annovarDB is being run!");
	#			`$annovar -filter -dbtype avsift -buildver hg19 $outDir/$annInput $annovarDB`;
	#		};
	#		if ($@) {
	#			$logger->fatal("Command: $annovar -filter -dbtype avsift -buildver hg19 $outDir/$annInput $annovarDB can not be run! Error: $@");
	#			exit(1);
	#		}
	#	}

	return ($annInput);
}

##############################
##############################
# Merge Annotations

sub MergeAnnotations {
	my ( $sample, $annFolder, $canonical, $translatedFolder ) = @_;
	my ( %canonicalTx, %strandHash, %aa_annotHash, @deleteFiles );

	#print "$sample\n";
	my ($annInput) = $sample =~ /(.*)\.txt/;
	push @deleteFiles, $annInput;
	$logger->info("Annovar annotations are being merged to variant file......");
	open( IN, "<", $canonical )
	  or die $logger->fatal("Can not open $canonical: $!");
	while (<IN>) {
		chomp;
		next if ( $_ =~ m/\#/ );
		my @line  = split( "\t", $_ );
		my @annot = split( ":",  $line[4] );
		my $gene  = $annot[0];
		my $transcriptID = $annot[1];
		my $strand       = $line[3];
		$canonicalTx{$gene} = $transcriptID;
		$strandHash{$gene}  = $strand;
	}
	close IN;

	my ( %annExonicHash, %annVariantHash, %annDBSNPhash, %annCOSMIChash, %ann1000Ghash, %annMAhash, %annGERPhash, %annLRThash, %annPhyloPhash, %annLJBSIFThash, %annAVSIFThash, @header );

	tie( my %inputHash, 'Tie::IxHash' );
	open IN, "<", $sample || die $logger->fatal("Can not open $sample: $!");
	while (<IN>) {
		chomp;

		#M-1608-42-T	1	9784916	C	T	KEEP
		if ( $_ =~ /^Sample\t/ ) {
			@header = split("\t");
			next;
		}
		my @line     = split("\t");
		my $sample   = shift @line;
		my $normal   = shift @line;
		my $chr      = shift @line;
		my $pos      = shift @line;
		my $ref      = shift @line;
		my $alt      = shift @line;
		my $at       = join( "|", @line );
		my $key_orig = $sample . ":" . $normal . ":" . $chr . ":" . $pos . ":" . $ref . ":" . $alt;

		#print "$key_orig\n";
		$inputHash{$key_orig} = $at;
	}
	close IN;

	# Read in exonic variant annotation
	($sample) = $sample =~ /.*\/(.*)/;
	my $annExonicOutput = $sample;
	$annExonicOutput =~ s/.txt/.exonic_variant_function/;
	push @deleteFiles, $annExonicOutput;
	open IN, "<", "$annFolder/$annExonicOutput"
	  || die $logger->fatal("Can not open $annExonicOutput: $!");
	while (<IN>) {
		chomp;
		my @line     = split("\t");
		my $result   = $line[1];
		my $chr      = $line[3];
		my $startPos = $line[4];
		my $endPos   = $line[5];
		my $ref      = $line[6];
		my $alt      = $line[7];
		$result =~ s/\ /_/;
		my @annots = split( ",", $line[2] );
		my ( $gene, $transcriptID, $exon, $cDNAchange, $AAchange );

		#traverse through annovar annotations and pick the canonical one
		foreach my $annots ( sort @annots ) {
			my @e          = split( ":", $annots );
			my $gen        = $e[0];
			my $transcript = $e[1];
			next if ( $annots eq "" );
			if ( $gen eq "CDKN2A" ) {    #CDKN2A exception: depending on the tx ID either change the name to p16 or p14
				if ( $transcript eq "NM_000077" ) {
					$gene         = "CDKN2Ap16INK4A";
					$transcriptID = $e[1];
					$exon         = $e[2];
					$cDNAchange   = $e[3];
					$AAchange     = $e[4];
					( $cDNAchange, $AAchange ) = &HGVScomplianceCheck( $result, $AAchange, $cDNAchange, $gene, $ref, \%strandHash, $translatedFolder );
					my @rest   = split( "\Q|", $line[8] );
					my $sample = $rest[0];
					my $normal = $rest[1];
					my $key    = $sample . ":" . $normal . ":" . $chr . ":" . $startPos . ":" . $endPos . ":" . $ref . ":" . $alt;
					my $value  = $result . ":" . $gene . ":" . $transcriptID . ":" . $exon . ":" . $cDNAchange . ":" . $AAchange;
					push @{ $annExonicHash{$key} }, $value;
				}
				elsif ( $transcript eq "NM_058195" ) {
					$gene         = "CDKN2Ap14ARF";
					$transcriptID = $e[1];
					$exon         = $e[2];
					$cDNAchange   = $e[3];
					$AAchange     = $e[4];
					( $cDNAchange, $AAchange ) = &HGVScomplianceCheck( $result, $AAchange, $cDNAchange, $gene, $ref, \%strandHash, $translatedFolder );
					my @rest   = split( "\Q|", $line[8] );
					my $sample = $rest[0];
					my $normal = $rest[1];
					my $key    = $sample . ":" . $normal . ":" . $chr . ":" . $startPos . ":" . $endPos . ":" . $ref . ":" . $alt;
					my $value  = $result . ":" . $gene . ":" . $transcriptID . ":" . $exon . ":" . $cDNAchange . ":" . $AAchange;
					push @{ $annExonicHash{$key} }, $value;
				}
				next;
			}
			elsif ( exists( $canonicalTx{$gen} )
				&& $transcript eq $canonicalTx{$gen} )
			{
				$gene         = $e[0];
				$transcriptID = $e[1];
				$exon         = $e[2];
				$cDNAchange   = $e[3];
				$AAchange     = $e[4];

			}
		}
		if ( !defined($gene) ) {    #Some of the variants are exonic only in non-canonical transcripts. If that's the case then pick the first annotation. I'm also adding a star to the tx ID to indicate it's non-canonical so we can fix it and make it canonical or something
			                        #my @e = split( ":", $annots[0] );
			                        #$gene         = $e[0];
			                        #$transcriptID = $e[1] . "*";
			                        #$exon         = $e[2];
			                        #$cDNAchange   = $e[3];
			                        #$AAchange     = $e[4];
			                        # Not doing this anymore. if a variant is exonic only in a non-canonical tx, print nothing.
			next;
		}

		if ( $exon eq "wholegene" ) {    # Annovar labels mutations in the start codon as whole gene
			$cDNAchange = "NA";
			$AAchange   = "wholegene";
			next;
		}

		#print "$result\t$AAchange\t";
		next if ( $gene =~ /CDKN2A/ );
		if ( $gene =~ /CDKN2A/ ) {
			print "WTF?\n";
		}
		( $cDNAchange, $AAchange ) = &HGVScomplianceCheck( $result, $AAchange, $cDNAchange, $gene, $ref, \%strandHash, $translatedFolder );
		my @rest   = split( "\Q|", $line[8] );
		my $sample = $rest[0];
		my $normal = $rest[1];
		my $key    = $sample . ":" . $normal . ":" . $chr . ":" . $startPos . ":" . $endPos . ":" . $ref . ":" . $alt;
		my $value  = $result . ":" . $gene . ":" . $transcriptID . ":" . $exon . ":" . $cDNAchange . ":" . $AAchange;
		push @{ $annExonicHash{$key} }, $value;
	}
	close IN;

	# Read in all variant annotations
	my $annVariantOutput = $sample;
	$annVariantOutput =~ s/.txt/.variant_function/;
	push @deleteFiles, $annVariantOutput;
	open IN, "<", "$annFolder/$annVariantOutput"
	  || die $logger->fatal("Can not open $annVariantOutput: $!");
	while (<IN>) {
		chomp;
		my @line   = split("\t");
		my $result = $line[0];
		$result =~ s/\ /_/;
		my ( $gene, $transcriptID, $exon, $cDNAchange, $AAchange, $ann );
		if ( $result eq "splicing" ) {
			( $gene, $ann ) = $line[1] =~ /(.*)\((.*)\)/;
			if ($ann) {
				my @annots = split( ",", $ann );
				foreach my $annots ( sort @annots ) {
					my @e = split( ":", $annots );
					if ( exists( $canonicalTx{$gene} )
						&& $e[0] eq $canonicalTx{$gene} )
					{
						if ( $e[0] eq "NM_000077" ) {
							$gene       = "CDKN2Ap16INK4A";
							$cDNAchange = $e[2];
						}
						elsif ( $e[0] eq "NM_058195" ) {
							$gene       = "CDKN2Ap14ARF";
							$cDNAchange = $e[2];
						}
						$transcriptID = $e[0];
						$exon         = $e[1];
						$cDNAchange   = $e[2];

					}
				}
			}

			if ( !$cDNAchange ) {
				$result       = "splicing_noncanonical";
				$gene         = $line[1];
				$cDNAchange   = "";
				$transcriptID = "";
				$exon         = "";
			}
			if ( !$gene ) {
				$gene         = $line[1];
				$cDNAchange   = "";
				$transcriptID = "";
				$exon         = "";
			}
		}
		else {
			$cDNAchange   = "";
			$transcriptID = "";
			$exon         = "";
			$gene         = $line[1];
		}
		my $chr      = $line[2];
		my $startPos = $line[3];
		my $endPos   = $line[4];
		my $ref      = $line[5];
		my $alt      = $line[6];
		my @rest     = split( "\Q|", $line[7] );
		my $sample   = $rest[0];
		my $normal   = $rest[1];

		if ( $startPos == 55248980 ) {    # EGFR ex20 intronic insertion exception
			$transcriptID = "NM_005228";
			$exon         = "exon20";
			$cDNAchange   = "c.2284-5_2290dupTCCAGGAAGCCT";
			$result       = "splicing";
			$gene         = "EGFR";
		}
		my $key = $sample . ":" . $normal . ":" . $chr . ":" . $startPos . ":" . $endPos . ":" . $ref . ":" . $alt;

		#print "$key\n";
		$annVariantHash{$key} = $result . ":" . $gene . ":" . $transcriptID . ":" . $exon . ":" . $cDNAchange;
	}

	close IN;

	# Read in dbSNP annotation
	my $annDBSNPoutput = $sample;
	$annDBSNPoutput =~ s/txt/hg19_snp137NonFlagged_dropped/;
	push @deleteFiles, $annDBSNPoutput;
	my $a = $annDBSNPoutput;
	$a =~ s/dropped/filtered/;
	push @deleteFiles, $a;
	open IN, "<", "$annFolder/$annDBSNPoutput"
	  || die $logger->fatal("Can not open $annDBSNPoutput: $!");

	while (<IN>) {
		chomp;
		my @line     = split("\t");
		my $dbSNPid  = $line[1];
		my $chr      = $line[2];
		my $startPos = $line[3];
		my $endPos   = $line[4];
		my $ref      = $line[5];
		my $alt      = $line[6];
		my @rest     = split( "\Q|", $line[7] );
		my $sample   = $rest[0];
		my $normal   = $rest[1];
		my $key      = $sample . ":" . $normal . ":" . $chr . ":" . $startPos . ":" . $endPos . ":" . $ref . ":" . $alt;

		#print "$key\n";
		$annDBSNPhash{$key} = $dbSNPid;
	}

	close IN;

	# Read in Cosmic annotation
	my $annCOSMICoutput = $sample;
	$annCOSMICoutput =~ s/txt/hg19_cosmic68_dropped/;
	push @deleteFiles, $annCOSMICoutput;
	$a = $annCOSMICoutput;
	$a =~ s/dropped/filtered/;
	push @deleteFiles, $a;
	open IN, "<", "$annFolder/$annCOSMICoutput"
	  || die $logger->fatal("Can not open $annCOSMICoutput: $!");

	while (<IN>) {
		chomp;
		my @line     = split("\t");
		my $COSMICid = $line[1];
		my $chr      = $line[2];
		my $startPos = $line[3];
		my $endPos   = $line[4];
		my $ref      = $line[5];
		my $alt      = $line[6];
		my @rest     = split( "\Q|", $line[7] );
		my $sample   = $rest[0];
		my $normal   = $rest[1];
		my $key      = $sample . ":" . $normal . ":" . $chr . ":" . $startPos . ":" . $endPos . ":" . $ref . ":" . $alt;

		#print "$key\n";
		$annCOSMIChash{$key} = $COSMICid;
	}

	close IN;

	# Read in 1000G MAF annotation
	my $ann1000Goutput = $sample;
	$ann1000Goutput =~ s/txt/hg19_ALL.sites.2012_04_dropped/;
	push @deleteFiles, $ann1000Goutput;
	$a = $ann1000Goutput;
	$a =~ s/dropped/filtered/;
	push @deleteFiles, $a;
	open IN, "<", "$annFolder/$ann1000Goutput"
	  || die $logger->fatal("Can not open $ann1000Goutput: $!");

	while (<IN>) {
		chomp;
		my @line     = split("\t");
		my $MAF      = $line[1];
		my $chr      = $line[2];
		my $startPos = $line[3];
		my $endPos   = $line[4];
		my $ref      = $line[5];
		my $alt      = $line[6];
		my @rest     = split( "\Q|", $line[7] );
		my $sample   = $rest[0];
		my $normal   = $rest[1];
		my $key      = $sample . ":" . $normal . ":" . $chr . ":" . $startPos . ":" . $endPos . ":" . $ref . ":" . $alt;

		#print "$key\n";
		$ann1000Ghash{$key} = $MAF;
	}

	close IN;

	# Read in Mutation Assesor annotation
	#	my $annMAoutput = $sample;
	#	$annMAoutput =~ s/txt/hg19_ljb2_ma_dropped/;
	#	push @deleteFiles, $annMAoutput;
	#	$a = $annMAoutput;
	#	$a =~ s/dropped/filtered/;
	#	push @deleteFiles, $a;
	#	open IN, "<", "$annFolder/$annMAoutput"
	#	  || die $logger->fatal("Can not open $annMAoutput: $!");
	#
	#	while (<IN>) {
	#		chomp;
	#		my @line     = split("\t");
	#		my $MA       = $line[1];
	#		my $chr      = $line[2];
	#		my $startPos = $line[3];
	#		my $endPos   = $line[4];
	#		my $ref      = $line[5];
	#		my $alt      = $line[6];
	#		my @rest     = split( "\Q|", $line[7] );
	#		my $sample   = $rest[0];
	#		my $normal   = $rest[1];
	#		my $key      = $sample . ":" . $normal . ":" . $chr . ":" . $startPos . ":" . $endPos . ":" . $ref . ":" . $alt;
	#
	#		#print "$key\n";
	#		$annMAhash{$key} = $MA;
	#	}
	#	close IN;

	# Read in AV SIFT scores
	#	my $annAVSIFToutput = $sample;
	#	$annAVSIFToutput =~ s/txt/hg19_avsift_dropped/;
	#	push @deleteFiles, $annAVSIFToutput;
	#	$a = $annAVSIFToutput;
	#	$a =~ s/dropped/filtered/;
	#	push @deleteFiles, $a;
	#	open IN, "<", "$annFolder/$annAVSIFToutput"
	#	  || die $logger->fatal("Can not open $annAVSIFToutput: $!");
	#
	#	while (<IN>) {
	#		chomp;
	#		my @line     = split("\t");
	#		my $AVSIFT   = $line[1];
	#		my $chr      = $line[2];
	#		my $startPos = $line[3];
	#		my $endPos   = $line[4];
	#		my $ref      = $line[5];
	#		my $alt      = $line[6];
	#		my @rest     = split( "\Q|", $line[7] );
	#		my $sample   = $rest[0];
	#		my $normal   = $rest[1];
	#		my $key      = $sample . ":" . $normal . ":" . $chr . ":" . $startPos . ":" . $endPos . ":" . $ref . ":" . $alt;
	#
	#		#print "$key\n";
	#		$annAVSIFThash{$key} = $AVSIFT;
	#	}
	#
	#	close IN;

	# Read in Phylo P (conservation) scores
	#	my $annPhyloPoutput = $sample;
	#	$annPhyloPoutput =~ s/txt/hg19_ljb2_phylop_dropped/;
	#	push @deleteFiles, $annPhyloPoutput;
	#	$a = $annPhyloPoutput;
	#	$a =~ s/dropped/filtered/;
	#	push @deleteFiles, $a;
	#	open IN, "<", "$annFolder/$annPhyloPoutput"
	#	  || die $logger->fatal("Can not open $annPhyloPoutput: $!");
	#
	#	while (<IN>) {
	#		chomp;
	#		my @line     = split("\t");
	#		my $PhyloP   = $line[1];
	#		my $chr      = $line[2];
	#		my $startPos = $line[3];
	#		my $endPos   = $line[4];
	#		my $ref      = $line[5];
	#		my $alt      = $line[6];
	#		my @rest     = split( "\Q|", $line[7] );
	#		my $sample   = $rest[0];
	#		my $normal   = $rest[1];
	#		my $key      = $sample . ":" . $normal . ":" . $chr . ":" . $startPos . ":" . $endPos . ":" . $ref . ":" . $alt;
	#
	#		#print "$key\n";
	#		$annPhyloPhash{$key} = $PhyloP;
	#	}
	#
	#	close IN;

	# Read in LRT scores
	#	my $annLRToutput = $sample;
	#	$annLRToutput =~ s/txt/hg19_ljb2_lrt_dropped/;
	#	push @deleteFiles, $annLRToutput;
	#	$a = $annLRToutput;
	#	$a =~ s/dropped/filtered/;
	#	push @deleteFiles, $a;
	#	open IN, "<", "$annFolder/$annLRToutput"
	#	  || die $logger->fatal("Can not open $annLRToutput: $!");
	#
	#	while (<IN>) {
	#		chomp;
	#		my @line     = split("\t");
	#		my $LRT      = $line[1];
	#		my $chr      = $line[2];
	#		my $startPos = $line[3];
	#		my $endPos   = $line[4];
	#		my $ref      = $line[5];
	#		my $alt      = $line[6];
	#		my @rest     = split( "\Q|", $line[7] );
	#		my $sample   = $rest[0];
	#		my $normal   = $rest[1];
	#		my $key      = $sample . ":" . $normal . ":" . $chr . ":" . $startPos . ":" . $endPos . ":" . $ref . ":" . $alt;
	#
	#		#print "$key\n";
	#		$annLRThash{$key} = $LRT;
	#	}
	#
	#	close IN;

	# Read in GERP++ scores
	#	my $annGERPoutput = $sample;
	#	$annGERPoutput =~ s/txt/hg19_ljb2_gerp++_dropped/;
	#	push @deleteFiles, $annGERPoutput;
	#	$a = $annGERPoutput;
	#	$a =~ s/dropped/filtered/;
	#	push @deleteFiles, $a;
	#	open IN, "<", "$annFolder/$annGERPoutput"
	#	  || die $logger->fatal("Can not open $annGERPoutput: $!");
	#
	#	while (<IN>) {
	#		chomp;
	#		my @line     = split("\t");
	#		my $GERP     = $line[1];
	#		my $chr      = $line[2];
	#		my $startPos = $line[3];
	#		my $endPos   = $line[4];
	#		my $ref      = $line[5];
	#		my $alt      = $line[6];
	#		my @rest     = split( "\Q|", $line[7] );
	#		my $sample   = $rest[0];
	#		my $normal   = $rest[1];
	#		my $key      = $sample . ":" . $normal . ":" . $chr . ":" . $startPos . ":" . $endPos . ":" . $ref . ":" . $alt;
	#
	#		#print "$key\n";
	#		$annGERPhash{$key} = $GERP;
	#	}
	#
	#	close IN;
	$sample =~ s/.txt/_annovarAnnotated.txt/;

	open OUT, ">", $sample || die $logger->fatal("Can not open $sample: $!");

	my $H1         = shift @header;
	my $H2         = shift @header;
	my $H3         = shift @header;
	my $H4         = shift @header;
	my $H5         = shift @header;
	my $H6         = shift @header;
	my $H7         = $header[3];
	my $H8         = $header[4];
	my @header_new = splice( @header, 9 );
	my @Hrest      = join( "\t", @header_new );
	print OUT "$H1\t$H2\t$H3\t$H4\t$H5\t$H6\tVariantClass\tGene\tExon\tCall_Confidence\tComments\tTranscriptID\tcDNAchange\tAAchange\tdbSNP_ID\tCosmic_ID\t1000G_MAF\t$H7\t$H8\t@Hrest\n";

	# Merge Annotations...
	foreach my $key ( keys %inputHash ) {
		my ( $sample, $normal, $chr, $pos, $ref, $alt ) =
		  split( ":", $key );
		my $annotKey;
		if ( length($ref) > length($alt) ) {    # deletion
			my $newref   = substr( $ref, 1 );
			my $posEnd   = $pos + length($newref);
			my $posStart = $pos + 1;
			$annotKey = $sample . ":" . $normal . ":" . $chr . ":" . $posStart . ":" . $posEnd . ":" . $newref . ":" . "-";
		}
		elsif ( length($ref) < length($alt) ) {    #insertion
			my $newalt = substr( $alt, 1 );
			my $posEnd = $pos;
			$annotKey = $sample . ":" . $normal . ":" . $chr . ":" . $pos . ":" . $pos . ":-:" . $newalt;
		}
		else {
			$annotKey = $sample . ":" . $normal . ":" . $chr . ":" . $pos . ":" . $pos . ":" . $ref . ":" . $alt;
		}
		if ( exists( $annExonicHash{$annotKey} ) ) {
			my @a = @{ $annExonicHash{$annotKey} };

			#print "@a\n";
			if ( scalar(@a) > 1 ) {
				print scalar(@a) . "\n";

				#check to see if there is more than one annotation in the hash for a given key (this is for variants occuring in both txIDs of CDKN2A)
				foreach my $CDKN2Aannot (@a) {
					my ( $result, $gene, $txID, $exon, $cDNAchange, $AAchange ) = split( ":", $CDKN2Aannot );

					#print "$CDKN2Aannot\n";
					print OUT "$sample\t$normal\t$chr\t$pos\t$ref\t$alt\t$result\t$gene\t$exon\t\t\t$txID\t$cDNAchange\t$AAchange\t";
					if ( exists( $annDBSNPhash{$annotKey} ) ) {
						print OUT "$annDBSNPhash{$annotKey}\t";
					}
					else {
						print OUT "\t";
					}
					if ( exists( $annCOSMIChash{$annotKey} ) ) {
						print OUT "$annCOSMIChash{$annotKey}\t";
					}
					else {
						print OUT "\t";
					}
					if ( exists( $ann1000Ghash{$annotKey} ) ) {
						print OUT "$ann1000Ghash{$annotKey}\t";
					}
					else {
						print OUT "\t";
					}

					#					if ( exists( $annMAhash{$annotKey} ) ) {
					#						print OUT "$annMAhash{$annotKey}\t";
					#					}
					#					else {
					#						print OUT "\t";
					#					}
					#					if ( exists( $annAVSIFThash{$annotKey} ) ) {
					#						print OUT "$annAVSIFThash{$annotKey}\t";
					#					}
					#					else {
					#						print OUT "\t";
					#					}
					#					if ( exists( $annPhyloPhash{$annotKey} ) ) {
					#						print OUT "$annPhyloPhash{$annotKey}\t";
					#					}
					#					else {
					#						print OUT "\t";
					#					}
					#					if ( exists( $annLRThash{$annotKey} ) ) {
					#						print OUT "$annLRThash{$annotKey}\t";
					#					}
					#					else {
					#						print OUT "\t";
					#					}
					#					if ( exists( $annGERPhash{$annotKey} ) ) {
					#						print OUT "$annGERPhash{$annotKey}\t";
					#					}
					#					else {
					#						print OUT "\t";
					#					}
					my @rest           = split( "\Q|", $inputHash{$key} );
					my $failure_reason = $rest[3];
					my $call_method = $rest[4];
					my @rest_new       = splice( @rest, 9 );
					@rest_new = join( "\t", @rest_new );
					print OUT "$failure_reason\t$call_method\t@rest_new\n";
				}
				next;
			}
			else {
				my ( $result, $gene, $txID, $exon, $cDNAchange, $AAchange ) =
				  split( ":", $a[0] );
				print OUT "$sample\t$normal\t$chr\t$pos\t$ref\t$alt\t$result\t$gene\t$exon\t\t\t$txID\t$cDNAchange\t$AAchange\t";
			}
		}
		elsif ( exists( $annVariantHash{$annotKey} ) ) {
			my ( $result, $gene, $txID, $exon, $cDNAchange ) =
			  split( ":", $annVariantHash{$annotKey} );
			print OUT "$sample\t$normal\t$chr\t$pos\t$ref\t$alt\t$result\t$gene\t$exon\t\t\t$txID\t$cDNAchange\t\t";
		}
		if ( exists( $annDBSNPhash{$annotKey} ) ) {
			print OUT "$annDBSNPhash{$annotKey}\t";
		}
		else {
			print OUT "\t";
		}
		if ( exists( $annCOSMIChash{$annotKey} ) ) {
			print OUT "$annCOSMIChash{$annotKey}\t";
		}
		else {
			print OUT "\t";
		}
		if ( exists( $ann1000Ghash{$annotKey} ) ) {
			print OUT "$ann1000Ghash{$annotKey}\t";
		}
		else {
			print OUT "\t";
		}

		#		if ( exists( $annMAhash{$annotKey} ) ) {
		#			print OUT "$annMAhash{$annotKey}\t";
		#		}
		#		else {
		#			print OUT "\t";
		#		}
		#		if ( exists( $annAVSIFThash{$annotKey} ) ) {
		#			print OUT "$annAVSIFThash{$annotKey}\t";
		#		}
		#		else {
		#			print OUT "\t";
		#		}
		#		if ( exists( $annPhyloPhash{$annotKey} ) ) {
		#			print OUT "$annPhyloPhash{$annotKey}\t";
		#		}
		#		else {
		#			print OUT "\t";
		#		}
		#		if ( exists( $annLRThash{$annotKey} ) ) {
		#			print OUT "$annLRThash{$annotKey}\t";
		#		}
		#		else {
		#			print OUT "\t";
		#		}
		#		if ( exists( $annGERPhash{$annotKey} ) ) {
		#			print OUT "$annGERPhash{$annotKey}\t";
		#		}
		#		else {
		#			print OUT "\t";
		#		}
		my @rest           = split( "\Q|", $inputHash{$key} );
		my $failure_reason = $rest[3];
		my $call_method = $rest[4];
		my @rest_new       = splice( @rest, 9 );
		@rest_new = join( "\t", @rest_new );
		print OUT "$failure_reason\t$call_method\t@rest_new\n";
	}
	$logger->info("Annotations are merged...");
	close OUT;
	return ( $sample, \@deleteFiles );
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

#################################
# Get Config file and program locations - adopted from Ronak's pipeline
# Input config file

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
	return ( \%location, \%version, \%parameters );
}

###################################################
###################################################
#--Delete Unwanted Files
sub HouseKeeping {
	my ( $outdir, @list ) = @_;
	$logger->info("Removing Unwanted files...");
	foreach my $file (@list) {
		if ( $file =~ /\// ) {
			`rm $file`;
		}
		else {
			`rm $outdir/$file`;
		}
	}
	return;
}

###########################
# exon sort
sub my_sort3 {
	my ( $gene_a, $exon_a ) = split( ":", $a );
	my ( $gene_b, $exon_b ) = split( ":", $b );
	my ($ex_a) = $exon_a =~ /exon(\d+)/;
	my ($ex_b) = $exon_b =~ /exon(\d+)/;
	return ( $gene_a cmp $gene_b || $ex_a <=> $ex_b );
}

####################################
# chromosome sort
sub my_sort {
	my ( $chr_a, $t, $y, $z ) = split( ":", $a );
	my ( $chr_b, $r, $e, $w ) = split( ":", $b );
	if ( $chr_a eq "X" ) { $chr_a = 23; }
	if ( $chr_a eq "Y" ) { $chr_a = 24; }
	if ( $chr_a eq "M" ) { $chr_a = 25; }
	if ( $chr_b eq "X" ) { $chr_b = 23; }
	if ( $chr_b eq "Y" ) { $chr_b = 24; }
	if ( $chr_b eq "M" ) { $chr_b = 25; }
	return ( $chr_a <=> $chr_b );
}

####################################
# chromosome sort
sub my_sort2 {
	my ( $chr_a, $gene_a, undef, $exon1, $strand_a, $aa_range_a ) =
	  split( ":", $a );
	my ($num1) = $exon1 =~ /exon(\d+)/;
	my ( $chr_b, $gene_b, undef, $exon2, $strand_b, $aa_range_b ) =
	  split( ":", $b );
	my ($num2) = $exon2 =~ /exon(\d+)/;
	if ( $chr_a eq "X" ) { $chr_a = 23; }
	if ( $chr_a eq "Y" ) { $chr_a = 24; }
	if ( $chr_a eq "M" ) { $chr_a = 25; }
	if ( $chr_b eq "X" ) { $chr_b = 23; }
	if ( $chr_b eq "Y" ) { $chr_b = 24; }
	if ( $chr_b eq "M" ) { $chr_b = 25; }
	return ( $chr_a <=> $chr_b || $gene_a cmp $gene_b || $num1 <=> $num2 );
}

#####################################
# Reverse compliment DNA strand
sub revcomp {
	my $dna = shift @_;
	my $ret = reverse($dna);
	$ret =~ tr/ACGTacgt/TGCAtgca/;
	return ($ret);
}

# Check the HGVS compliance of the cDNA and AA annotations

sub HGVScomplianceCheck {

	my ( $result, $AAchange, $cDNAchange, $gene, $ref, $strandHash_ref, $translatedFolder ) = @_;
	my %strandHash = %$strandHash_ref;

	# Modify the output to be HGVS compliant (Modified from Donavan's script)
	my %aa_annotHash;
	if ( $result eq "nonframeshift_insertion" ) {

		# cDNA is correct
		# AA needs changing
		my ( $aa1, $num1, $new_aa ) = $AAchange =~ /p\.([A-Z])(\d+)delins([A-Z]+)/;

		if ( !$num1 || !$aa1 || !$new_aa ) {
			print "$gene\t$ref\t$AAchange\t$result\n";
		}
		my $num2 = $num1 + 1;
		open IN1, "<", "$translatedFolder/$gene.translated"
		  || die $logger->fatal("Can not open $gene.translated: $!");
		while (<IN1>) {
			chomp;
			if ( $_ =~ /^\>/ ) { last; }
			my @f = split('\t');
			$aa_annotHash{ $f[0] } = $f[2];
		}
		close IN1;
		my $aa2 = $aa_annotHash{$num2};
		if ( $strandHash{$gene} eq "+" ) {
			$new_aa = substr( $new_aa, 1 );
		}
		else {
			$new_aa = substr( $new_aa, 0, length($new_aa) - 1 );
		}
		$AAchange = "p.${aa1}${num1}_${aa2}${num2}ins${new_aa}";
		return ( $cDNAchange, $AAchange );

	}
	elsif ( $result eq "frameshift_insertion" ) {

		# cDNA: do nothing
		# Protein: do nothing
		return ( $cDNAchange, $AAchange );
	}
	elsif ( $result eq "nonframeshift_deletion" ) {

		# cDNA:
		if ( $strandHash{$gene} eq "+" ) {
			$cDNAchange .= $ref;
		}
		else {
			$cDNAchange .= revcomp($ref);
		}
		if ( $AAchange =~ /_/ ) {
			return ( $cDNAchange, $AAchange );

			# Protein: for canonical transcripts only:
			open IN2, "<", "$translatedFolder/$gene.translated"
			  || die $logger->fatal("Can not open $gene.translated: $!");
			while (<IN2>) {
				chomp;
				if ( $_ =~ /^\>/ ) { last; }
				my @f = split('\t');
				$aa_annotHash{ $f[0] } = $f[2];
			}
			close IN2;
			my ( $num1, $num2 ) = $AAchange =~ /p\.(\d+)_(\d+)del/;
			if ( !$num1 || !$num2 ) {
				print "$gene\t$ref\t$AAchange\t$result\n";
			}
			my $aa1 = $aa_annotHash{$num1};
			my $aa2 = $aa_annotHash{$num2};
			if ( $num1 == $num2 ) {
				$AAchange = "p.${aa1}${num1}del";
			}
			else {
				$AAchange = "p.${aa1}${num1}_${aa2}${num2}del";
			}
			return ( $cDNAchange, $AAchange );
		}
		else {
			return ( $cDNAchange, $AAchange );
		}
	}
	elsif ( $result eq "frameshift_deletion" ) {

		# cDNA:
		my ($b4_del) = $cDNAchange =~ /(.*)del/;
		if ( !$b4_del ) {
			print "$gene\t$ref\t$AAchange\t$result\n";
		}
		if ( $strandHash{$gene} eq "+" ) {
			$cDNAchange = $b4_del . "del" . $ref;
		}
		else {
			$cDNAchange = $b4_del . "del" . revcomp($ref);
		}

		# Protein: for canonical transcripts only:
		open IN3, "<", "$translatedFolder/$gene.translated"
		  || die $logger->fatal("Can not open $gene.translated: $!");
		while (<IN3>) {
			chomp;
			if ( $_ =~ /^\>/ ) { last; }
			my @f = split('\t');
			$aa_annotHash{ $f[0] } = $f[2];
		}
		close IN3;
		if ( !( $AAchange =~ /fs/ ) ) {
			my ($num1) = $AAchange =~ /p\.(\d+)_/;
			my $aa1 = $aa_annotHash{$num1};
			$AAchange = "p.${aa1}${num1}fs";
		}
		return ( $cDNAchange, $AAchange );
	}
	elsif ($result eq "nonsynonymous_SNV"
		|| $result eq "stopgain_SNV"
		|| $result eq "synonymous_SNV"
                || $result eq "stoploss_SNV" )
	{
		my ( $firstNTD, $position, $lastNTD ) = $cDNAchange =~ /c\.([A-Z])(\d+)([A-Z])/;
		if ( $position && $firstNTD && $lastNTD ) {
			$cDNAchange = "c." . $position . $firstNTD . ">" . $lastNTD;
		}
                $AAchange =~ s/X/*/g;
		return ( $cDNAchange, $AAchange );
	}
	else {
		return ( $cDNAchange, $AAchange );
	}
}

__END__

=head1 NAME

dmp_annotate_filter_variants.pl - Describe the usage of script briefly

=head1 SYNOPSIS

dmp_annotate_filter_variants.pl [options] args

      -opt --long      Option description

=head1 DESCRIPTION

Stub documentation for dmp_annotate_filter_variants.pl, 

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
