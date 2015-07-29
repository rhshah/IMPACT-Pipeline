#!/usr/bin/perl -w
# dmp_create_coverage_files.pl ---
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
use Data::Dumper qw(Dumper);

my $logger = MSKCC_DMP_Logger->get_logger('CREATE_COVERAGE_FILES');
$logger->start_local();

my (
	$titleFile, $outdir,$normalsInVariantCalling, $sampleSheet, $geneCoverageFile, $exonIntervals, $filteredAnnoFile, $exonCoverageFile,
	$translatedFolder, $coverage_threshold, $canonicalCoverageFile, $clinicalExons, $validationExons, $copyNumberFile, $geneIntervalAnn, $intragenicCallFile
);

if (
	@ARGV < 1
	or !GetOptions(
		'annotatedFilteredVariants|v:s'  => \$filteredAnnoFile,
		'geneLevelCopynumberFile|cn:s'   => \$copyNumberFile,
		'intragenicDeletionFile|id:s'    => \$intragenicCallFile,
		'geneIntervalAnnotated|gi:s'     => \$geneIntervalAnn,
		'outdir|o:s'                     => \$outdir,
		'TitleFile|t:s'                  => \$titleFile,
		'SampleSheet|ss:s'               => \$sampleSheet,
		'geneCoverageFile|gc:s'          => \$geneCoverageFile,
		'exonCoverageFile|ec:s'          => \$exonCoverageFile,
		'canonicalExonCoverageFile|cc:s' => \$canonicalCoverageFile,
		'clinicalExonList|ce:s'          => \$clinicalExons,
		'validationExonList|ve:s'        => \$validationExons,
		'exonIntervalFile|ei:s'          => \$exonIntervals,
		'translatedFolder|tf:s'          => \$translatedFolder,
		'coverageThreshold|ct:i'         => \$coverage_threshold,
		'normalsUsedInVariantCalling|nu:s' => \$normalsInVariantCalling,

	)
  )
{
	Usage();
}

if ( !$outdir ) {
	$outdir = getcwd;

}

# Generate per patient per gene coverage files
&GenerateGeneCoverage( $titleFile, $outdir, $geneCoverageFile, $exonIntervals );

# Generate per patient per exon coverage files
my ($validatedExonCoverage) = &GenerateExonCoverage( $filteredAnnoFile, $canonicalCoverageFile, $translatedFolder, $exonIntervals, $outdir, $titleFile, $coverage_threshold, $validationExons );

# Generate per patient clinical reports
&GenerateClinicalReports(
	$filteredAnnoFile,      $normalsInVariantCalling, $copyNumberFile, $sampleSheet,   $intragenicCallFile,    $geneCoverageFile, $geneIntervalAnn,
	$canonicalCoverageFile, $titleFile,      $clinicalExons, $validatedExonCoverage, $validationExons
);

#####################################
#####################################

sub GenerateClinicalReports {
	my (
		$variantsList,          $normalsInVariantCalling, $copyNumberFile, $sampleSheet,   $intragenicCallFile,    $geneCovFile, $geneIntervalAnn,
		$canonicalCoverageFile, $titleFile,      $clinicalExons, $validatedExonCoverage, $validationExons
	  )
	  = @_;
	$logger->info("Starting clinical report generation.");
	my %validatedExonCoverage = %$validatedExonCoverage;
	my $totalGeneNumber       = 0;
	my ( %clinicalExonsHash, %validationExonsHash, %copyNumberHash, %cytobandHash, %intragenicDeletionHash, %patientNameHash );

	open SAMPLESHEET, "<", $sampleSheet || $logger->fatal("Can not open $sampleSheet: $!");
	while (<SAMPLESHEET>) {
		next if (/FCID/);
		my @line = split(",");

		#H7HRLADXX,1,35197348-N,HUMAN,GTCCGC,Normal,N,,JS|Friedman;Cynthia|M14-1428|Female,IMPACTv3-CLIN-20140005
		my ( $tech, $name, $Mnumber, $sex ) = split( "\Q|", $line[8] );
		my $sampleName = $line[2];
		if ( $name && $Mnumber && $sex ) {
			$patientNameHash{$sampleName} = "$name $Mnumber";
		}
	}
	close SAMPLESHEET;

	open INTRAGENIC, "<", $intragenicCallFile || $logger->fatal("Can not open $intragenicCallFile : $!");
	while (<INTRAGENIC>) {
		chomp;
		next if (/sample/);
		my @line        = split("\t");
		my $sample      = $line[0];
		my $significant = $line[9];
		my $FC          = $line[7];
		my $annot       = $line[12];
		my $cytoband    = $line[11];
		my ( $gene, $txID, $exon, undef, undef, undef ) = split( ":", $annot );

		if ( $significant == 1 ) {
			$FC = sprintf( "%.2f", $FC );
			if ( $gene eq "EGFR" || $gene eq "PDGFRA" ) {
				my $a = "$gene ($txID - $cytoband) Intragenic deletion of the following exon(s):";
				push @{ $intragenicDeletionHash{$sample}{$a} }, $exon;
			}
			else {
				my $a = "$gene ($txID - $cytoband) Intragenic deletion";
				push @{ $intragenicDeletionHash{$sample}{$a} }, $exon;
			}
		}

	}
	close INTRAGENIC;

	open GENEINT, "<", $geneIntervalAnn || $logger->fata("Can not open $geneIntervalAnn: $!");
	while (<GENEINT>) {
		chomp;
		my ( $coordinates, $cytoband, $annotation ) = split(" ");
		my ( $gene, $txID, $exon, $chr, $start, $stop ) = split( ":", $annotation );
		if ( $gene && $txID ) {
			$cytobandHash{$gene} = "(" . $txID . " - " . $cytoband . ")";
		}
	}
	close GENEINT;

	# Fix for CDKN2A name change:
	$cytobandHash{"CDKN2A"} = "(NM_000077 - 9p21.3)";
	open COPYNUMBER, "<", $copyNumberFile || $logger->fatal("Can not open $copyNumberFile: $!");
	while (<COPYNUMBER>) {
		chomp;
		next if (/sample/);
		my @line        = split("\t");
		my $sample      = $line[0];
		my $gene        = $line[3];
		my $significant = $line[9];
		my $FC          = $line[7];
		my $cytoband    = $cytobandHash{$gene};
		next if ( $gene =~ /Tiling/ || $gene =~ /FP/ );

		if ( $significant == 1 ) {
			$FC = sprintf( "%.2f", $FC );
			if ( $FC > 0 ) {
				push @{ $copyNumberHash{$sample} }, "$gene $cytoband Amplification (Fold Change: $FC)";
			}
			else {
				push @{ $copyNumberHash{$sample} }, "$gene $cytoband Deletion (Fold Change: $FC)";
			}
		}
	}
	close COPYNUMBER;

	open VALEXONS, "<", $validationExons || $logger->fatal("Can not open $validationExons: $!");
	while (<VALEXONS>) {
		chomp;
		my ( $gene, $exon ) = split("\t");
		$validationExonsHash{ $gene . ":" . $exon } = 0;
	}
	close VALEXONS;

	open CLINICAL, "<", $clinicalExons || $logger->fatal("Can not open $clinicalExons: $!");
	while (<CLINICAL>) {
		chomp;
		my ( $gene, $exon, undef ) = split("\t");
		$clinicalExonsHash{ $gene . ":" . $exon } = 0;
	}
	close CLINICAL;

	open IN, "<", $geneCovFile || $logger->fatal("Can not open $geneCovFile: $! ");
	while (<IN>) {
		chomp;
		next if (/Gene/);
		$totalGeneNumber++;
	}
	close IN;

	my ( %barcodeHash, %samplesHash, %patientIDHash, %matchedHash );
	open TITLE, "<", $titleFile || $logger->fatal("Can not open $titleFile: $!");
	while (<TITLE>) {
		chomp;
		next if (/Barcode/);
		my @line = split("\t");
		$barcodeHash{ $line[0] }   = $line[2];
		$patientIDHash{ $line[2] } = $line[4];
		if ( $line[5] =~ /Tumor/ ) {
			$samplesHash{ $line[2] } = 1;
		}
	}
	close TITLE;
	
	open NORMALS, "<", $normalsInVariantCalling || $logger->fatal("Can not open $normalsInVariantCalling: $!");
		while (<NORMALS>) {
			chomp;
			next if(/NTC/);
			my @line     = split("\t");
			my $sampleID = $line[0];
			my $normalID = $line[1];
			my $patID    = $patientIDHash{$sampleID};
			my $matched;
			if ( $normalID =~ /$patID/ ) {
				$matched = "Matched";
			}
			else {
				$matched = "Unmatched";
			}
			$matchedHash{$sampleID} = $matched;
			print "$sampleID\t$normalID\t$patID\t$matched\n";
			
		}

		close NORMALS;

	my $totalExonNumber = 0;
	open COV, "<", $canonicalCoverageFile || $logger->fatal("Can not open $canonicalCoverageFile: $!");
	while (<COV>) {
		chomp;
		my @header;
		if (/Target/) {
			@header = split("\t");
			next;
		}
		$totalExonNumber++;
	}
	close COV;

	my %meanCoverageHash;
	foreach my $barcode ( sort keys %barcodeHash ) {
		my $sample   = $barcodeHash{$barcode};
		my $totalCov = 0;
		open COV, "<", $canonicalCoverageFile || $logger->fatal("Can not open $canonicalCoverageFile: $!");
		my $head = <COV>;
		chomp($head);
		my @header = split( "\t", $head );
		my $num    = scalar(@header);

		#print "@header";
		while (<COV>) {
			chomp;
			my $index;
			foreach my $i ( 0 .. $#header ) {

				#print "$barcode\t$header[$i]\n";
				if ( $barcode eq $header[$i] ) {
					$index = $i;
				}
			}
			print "$sample\t$barcode\n" if ( !$index );
			my @line = split("\t");

			#print "$line[$index]\n";
			$totalCov += $line[$index];
		}
		my $meanCov = sprintf "%d", $totalCov / $totalExonNumber;
		$meanCoverageHash{$sample} = $meanCov;
	}
	my %variantsHash;

	open VARIANTS, "<", $variantsList || $logger->fatal("Can not open $variantsList: $!");
	while (<VARIANTS>) {
		chomp;
		next if (/Sample/);
		my @line        = split("\t");
		my $sampleID    = $line[0];
		my $normalID    = $line[1];
		my $gene        = $line[7];
		my $exon        = $line[8];
		my $txID        = $line[11];
		my $variantType = $line[6];
		my $cDNAchange  = $line[12];
		my $AAchange    = $line[13];
		my $variant;

		if ( $gene eq "TERT" && $variantType eq "upstream" ) {    #check for TERT promoter mutation and construct special case
			my $pos = $line[3];
			my $ref = &complemantary( $line[4] );
			my $alt = &complemantary( $line[5] );
			$variant = $gene . " (NM_198253) promoter variant (g." . $pos . $ref . ">" . $alt . ")";
		}
		elsif ( $variantType eq "splicing" ) {
			$variant = $gene . " (" . $txID . ") " . $exon . " splicing variant (" . $cDNAchange . ")";
		}
		else {
			$variant = $gene . " (" . $txID . ") " . $exon . " " . $AAchange . " (" . $cDNAchange . ")";
		}
		if ( exists( $validationExonsHash{ $gene . ":" . $exon } ) ) {
			push @{ $variantsHash{$sampleID}{"Clinical"} }, $variant;    # create hash of arrays for clinical and investigational variants

		}
		else {
			push @{ $variantsHash{$sampleID}{"Investigational"} }, $variant;

		}
	}
	close VARIANTS;
	if ( !-e "$outdir/ClinicalReports" ) {
		`mkdir ClinicalReports`;
	}

	foreach my $sample ( sort keys %samplesHash ) {

		my $matched = $matchedHash{$sample};
		

		#print "$outdir/ClinicalReports/${sample}_clinical_report.txt\n";
		open OUT, ">", "$outdir/ClinicalReports/${sample}_clinical_report.txt" || $logger->fatal("Can not open ${sample}_clinical_report.txt for writing: $!");
		my $name = $patientNameHash{$sample};
		if (!$name){print "$sample\n";}
		print OUT "Patient ID: $sample $name \n\n";
		print OUT "----------------------------------\n\n";
		my $coverage = $meanCoverageHash{$sample};
		my ( @ClinicalVariants, @InvestigationalVariants, @copyNumberVariants );

		if ( exists( $variantsHash{$sample}{"Clinical"} ) ) {    # This check is necessary in case a sample doesn't have any variant (or any clinical variants)
			@ClinicalVariants = sort ( @{ $variantsHash{$sample}{"Clinical"} } );

			#print "$sample\t@ClinicalVariants\n";

		}
		if ( exists( $variantsHash{$sample}{"Investigational"} ) ) {
			@InvestigationalVariants = sort @{ $variantsHash{$sample}{"Investigational"} };

			#print "$sample\t@InvestigationalVariants\n";
		}
		if ( exists( $copyNumberHash{$sample} ) ) {
			@copyNumberVariants = sort my_sort @{ $copyNumberHash{$sample} };
		}
		my $totalIntragenicDeletions = 0;
		if ( exists( $intragenicDeletionHash{$sample} ) ) {
			foreach my $gene ( sort keys %{ $intragenicDeletionHash{$sample} } ) {
				$totalIntragenicDeletions++;
			}
		}
		my $numberOfVariants = scalar(@ClinicalVariants) + scalar(@InvestigationalVariants) + scalar(@copyNumberVariants) + $totalIntragenicDeletions;    # total number of variants

		#print "$sample\t$coverage\t$numberOfVariants\n";
		if ( $coverage <= 50 ) {
			print OUT <<ENDOUT;
DNA extracted from this sample did not meet the minimum quality/quantity criteria for adequate performance of this assay. Mutation analysis for selected genes will be performed by another method and will be reported as an addendum.

TEST AND METHODOLOGY:
MSK-IMPACT (Integrated Mutation Profiling of Actionable Cancer Targets) was used to identify specific mutations in $totalGeneNumber genes. The following genes contain exons that have been clinically validated:

Gene (Transcript_ID) Exon (Amino Acid Range)
AKT1 (NM_001014431) exon 3 (16-59)
ALK (NM_004304) exon 23 (1172-1215), exon 25 (1248-1279)
BRAF (NM_004333) exon 11 (439-478), exon 15 (581-620)
EGFR (NM_005228) exon 18 (688-728), exon 19 (729-761), exon 20 (762-823), exon 21 (824-875)
ERBB2 (NM_004448) exon 8 (301-341), exon 19 (737-769), exon 20 (770-831)
FGFR2 (NM_000141) exon 7 (250-313), exon 9 (362-429), exon 12 (521-558)
FGFR3 (NM_000142) exon 7 (247-310), exon 9 (359-422), exon 18 (759-807)
GNA11 (NM_002067) exon 5 (202-245)
GNAQ (NM_002072) exon 5 (202-245)
GNAS (NM_000516) exon 8 (196-220)
HRAS (NM_001130442) exon 2 (1-37), exon 3 (38-97)
IDH1 (NM_005896) exon 4 (41-138)
IDH2 (NM_002168) exon 4 (125-178)	
KIT (NM_000222) exon 9 (449-514), exon 11 (550-592), exon 13 (627-664), exon 17 (788-828)
KRAS (NM_033360) exon 2 (1-37), exon 3 (38-97), exon 4 (97-150)
NRAS (NM_002524) exon 2 (1-37), exon 3 (38-97)
PDGFRA (NM_006206) exon 12 (552-596), exon 18 (814-854)
PIK3CA (NM_006218) exon 2 (1-118), exon 5 (272-353), exon 8 (418-468), exon 10 (514-555), exon 21 (979-1069)
TP53 (NM_000546) exon 4 (33-125), exon 5 (126-187), exon 6 (187-224), exon 7 (225-261), exon 8 (261-307), exon 10 (332-367)

For the full list of genes, please refer to http://cmo/uploads/2/4/9/3/24933115/impact_clinical_genes_2.0.xlsx (MSKCC intranet website). The specific mutations are detected by hybridization capture of DNA followed by massively parallel sequencing on an Illumina HiSeq 2500 instrument.
Diagnostic sensitivity: This assay is designed to detect single nucleotide variants and small insertions and deletions (< 30bp) in protein-coding exons of the $totalGeneNumber gene panel.
Technical sensitivity: This assay may not detect mutations in validated exons if the proportion of tumor cells in the sample studied is less than 10%. This assay is at risk of false negatives when sequence coverage for an exon is below 100X.
ENDOUT

			#close OUT;
		}
		elsif ( $coverage > 50 && $numberOfVariants == 0 ) {

			#print OUT "NEGATIVE FOR VARIANTS IN THE CLINICALLY VALIDATED PANEL.\n";
			#print OUT "NEGATIVE FOR VARIANTS IN THE INVESTIGATIONAL PANEL.\n";
			print OUT "NEGATIVE FOR ALTERATIONS";
			print OUT "\n";
			if ( $matched eq "Matched" ) {
				print OUT "Mutations called against: Matched Normal\n";
				print OUT "Note: Because a patient-matched normal sample was used, this assay reports somatic variants confirmed to be absent in the matched normal.\n";
			}
			else {
				print OUT "Mutations called against: Unmatched Normal\n";
				print OUT "Note: Because a patient-matched normal sample was unavailable, this assay does not distinguish between somatic and rare germline variants\n";
			}
			print OUT "\n";
			if ( exists( $validatedExonCoverage{$sample} ) ) {
				print OUT "MEAN OVERALL COVERAGE (SEQUENCING DEPTH) IN THIS SAMPLE: ${coverage}X\n";
				print OUT
"The following clinically validated exons showed coverage of less than 100X. A false negative result cannot be excluded in these regions, especially in samples with low neoplastic cell content.\n";
				print OUT "\nGene (Transcript_ID) Exon (Amino Acid Range)\n";
				foreach my $key ( sort keys %{ $validatedExonCoverage{$sample} } ) {
					my $exons = join( ", ", @{ $validatedExonCoverage{$sample}{$key} } );
					print OUT "$key $exons\n";
				}
			}
			else {
				print OUT "MEAN OVERALL COVERAGE (SEQUENCING DEPTH) IN THIS SAMPLE: ${coverage}X\n";
				print OUT "Unless specified, all exons tested had minimum depth of coverage of 100X.\n";
			}
			print OUT <<ENDOUT;

TEST AND METHODOLOGY:
MSK-IMPACT (Integrated Mutation Profiling of Actionable Cancer Targets) was used to identify specific mutations in $totalGeneNumber genes. The following genes contain exons that have been clinically validated:

Gene (Transcript_ID) Exon (Amino Acid Range)
AKT1 (NM_001014431) exon 3 (16-59)
ALK (NM_004304) exon 23 (1172-1215), exon 25 (1248-1279)
BRAF (NM_004333) exon 11 (439-478), exon 15 (581-620)
EGFR (NM_005228) exon 18 (688-728), exon 19 (729-761), exon 20 (762-823), exon 21 (824-875)
ERBB2 (NM_004448) exon 8 (301-341), exon 19 (737-769), exon 20 (770-831)
FGFR2 (NM_000141) exon 7 (250-313), exon 9 (362-429), exon 12 (521-558)
FGFR3 (NM_000142) exon 7 (247-310), exon 9 (359-422), exon 18 (759-807)
GNA11 (NM_002067) exon 5 (202-245)
GNAQ (NM_002072) exon 5 (202-245)
GNAS (NM_000516) exon 8 (196-220)
HRAS (NM_001130442) exon 2 (1-37), exon 3 (38-97)
IDH1 (NM_005896) exon 4 (41-138)
IDH2 (NM_002168) exon 4 (125-178)	
KIT (NM_000222) exon 9 (449-514), exon 11 (550-592), exon 13 (627-664), exon 17 (788-828)
KRAS (NM_033360) exon 2 (1-37), exon 3 (38-97), exon 4 (97-150)
NRAS (NM_002524) exon 2 (1-37), exon 3 (38-97)
PDGFRA (NM_006206) exon 12 (552-596), exon 18 (814-854)
PIK3CA (NM_006218) exon 2 (1-118), exon 5 (272-353), exon 8 (418-468), exon 10 (514-555), exon 21 (979-1069)
TP53 (NM_000546) exon 4 (33-125), exon 5 (126-187), exon 6 (187-224), exon 7 (225-261), exon 8 (261-307), exon 10 (332-367)

For the full list of genes, please refer to http://cmo/uploads/2/4/9/3/24933115/impact_clinical_genes_2.0.xlsx (MSKCC intranet website). The specific mutations are detected by hybridization capture of DNA followed by massively parallel sequencing on an Illumina HiSeq 2500 instrument.
Diagnostic sensitivity: This assay is designed to detect single nucleotide variants and small insertions and deletions (< 30bp) in protein-coding exons of the $totalGeneNumber gene panel.
Technical sensitivity: This assay may not detect mutations in validated exons if the proportion of tumor cells in the sample studied is less than 10%. This assay is at risk of false negatives when sequence coverage for an exon is below 100X.
ENDOUT

			#close OUT;

		}
		else {
			my $i = 1;
			print OUT "POSITIVE FOR THE FOLLOWING ALTERATIONS:\n\n";
			if ( scalar(@ClinicalVariants) > 0 ) {

				#print OUT "POSITIVE FOR THE FOLLOWING VARIANTS IN THE CLINICALLY VALIDATED PANEL:\n";

				foreach my $index ( 0 .. scalar(@ClinicalVariants) - 1 ) {
					print OUT "$i. $ClinicalVariants[$index]\n";
					$i++;
				}
				print OUT "\n";
			}
			else {

				#print OUT "NEGATIVE FOR VARIANTS IN THE CLINICALLY VALIDATED PANEL.\n";
			}

			if ( scalar(@copyNumberVariants) > 0 ) {
				foreach my $index ( 0 .. scalar(@copyNumberVariants) - 1 ) {
					print OUT "$i. $copyNumberVariants[$index]\n";
					$i++;
				}
				print OUT "\n";
			}
			if ( $totalIntragenicDeletions > 0 ) {
				foreach my $gene ( sort keys %{ $intragenicDeletionHash{$sample} } ) {
					if ( $gene =~ "EGFR" || $gene =~ "PDGFRA" ) {
						my $a = join( ", ", @{ $intragenicDeletionHash{$sample}{$gene} } );
						print OUT "$i. $gene $a\n";
						$i++;
					}
					else {

						#my $a = join( ", ", @{ $intragenicDeletionHash{$sample}{$gene} } );
						print OUT "$i. $gene\n";
						$i++;
					}
				}
				print OUT "\n";
			}
			if ( scalar(@InvestigationalVariants) > 0 ) {

				#print OUT "POSITIVE FOR THE FOLLOWING VARIANTS IN THE INVESTIGATIONAL PANEL:\n";
				foreach my $index ( 0 .. scalar(@InvestigationalVariants) - 1 ) {
					print OUT "$i. $InvestigationalVariants[$index]\n";
					$i++;
				}
			}
			else {

				#print OUT "NEGATIVE FOR VARIANTS IN THE INVESTIGATIONAL PANEL.\n";
			}
			print OUT "\n";
			print OUT "Note: results currently for investigational use only.\n\n";
			if ( $matched eq "Matched" ) {
				print OUT "Mutations called against: Matched Normal\n";
				print OUT "Note: Because a patient-matched normal sample was used, this assay reports somatic variants confirmed to be absent in the matched normal.\n";
			}
			else {
				print OUT "Mutations called against: Unmatched Normal\n";
				print OUT "Note: Because a patient-matched normal sample was unavailable, this assay does not distinguish between somatic and rare germline variants\n";
			}
			print OUT "\n";
			if ( exists( $validatedExonCoverage{$sample} ) ) {
				print OUT "MEAN OVERALL COVERAGE (SEQUENCING DEPTH) IN THIS SAMPLE: ${coverage}X\n";
				print OUT
"The following clinically validated exons showed coverage of less than 100X. A false negative result cannot be excluded in these regions, especially in samples with low neoplastic cell content.\n";
				print OUT "\nGene (Transcript_ID) Exon (Amino Acid Range)\n";
				foreach my $key ( sort keys %{ $validatedExonCoverage{$sample} } ) {
					my $exons = join( ", ", @{ $validatedExonCoverage{$sample}{$key} } );
					print OUT "$key $exons\n";
				}
			}
			else {
				print OUT "MEAN OVERALL COVERAGE (SEQUENCING DEPTH) IN THIS SAMPLE: ${coverage}X\n";
				print OUT "Unless specified, all exons tested had minimum depth of coverage of 100X.\n";
			}
			print OUT <<ENDOUT;

TEST AND METHODOLOGY:
MSK-IMPACT (Integrated Mutation Profiling of Actionable Cancer Targets) was used to identify specific mutations in $totalGeneNumber genes. The following genes contain exons that have been clinically validated:

Gene (Transcript_ID) Exon (Amino Acid Range)
AKT1 (NM_001014431) exon 3 (16-59)
ALK (NM_004304) exon 23 (1172-1215), exon 25 (1248-1279)
BRAF (NM_004333) exon 11 (439-478), exon 15 (581-620)
EGFR (NM_005228) exon 18 (688-728), exon 19 (729-761), exon 20 (762-823), exon 21 (824-875)
ERBB2 (NM_004448) exon 8 (301-341), exon 19 (737-769), exon 20 (770-831)
FGFR2 (NM_000141) exon 7 (250-313), exon 9 (362-429), exon 12 (521-558)
FGFR3 (NM_000142) exon 7 (247-310), exon 9 (359-422), exon 18 (759-807)
GNA11 (NM_002067) exon 5 (202-245)
GNAQ (NM_002072) exon 5 (202-245)
GNAS (NM_000516) exon 8 (196-220)
HRAS (NM_001130442) exon 2 (1-37), exon 3 (38-97)
IDH1 (NM_005896) exon 4 (41-138)
IDH2 (NM_002168) exon 4 (125-178)	
KIT (NM_000222) exon 9 (449-514), exon 11 (550-592), exon 13 (627-664), exon 17 (788-828)
KRAS (NM_033360) exon 2 (1-37), exon 3 (38-97), exon 4 (97-150)
NRAS (NM_002524) exon 2 (1-37), exon 3 (38-97)
PDGFRA (NM_006206) exon 12 (552-596), exon 18 (814-854)
PIK3CA (NM_006218) exon 2 (1-118), exon 5 (272-353), exon 8 (418-468), exon 10 (514-555), exon 21 (979-1069)
TP53 (NM_000546) exon 4 (33-125), exon 5 (126-187), exon 6 (187-224), exon 7 (225-261), exon 8 (261-307), exon 10 (332-367)

For the full list of genes, please refer to http://cmo/uploads/2/4/9/3/24933115/impact_clinical_genes_2.0.xlsx (MSKCC intranet website). The specific mutations are detected by hybridization capture of DNA followed by massively parallel sequencing on an Illumina HiSeq 2500 instrument.
Diagnostic sensitivity: This assay is designed to detect single nucleotide variants and small insertions and deletions (< 30bp) in protein-coding exons of the $totalGeneNumber gene panel.
Technical sensitivity: This assay may not detect mutations in validated exons if the proportion of tumor cells in the sample studied is less than 10%. This assay is at risk of false negatives when sequence coverage for an exon is below 100X.
ENDOUT

			#close OUT;
		}
	}
	$logger->info("Clinical report generation is completed.");
}

#####################################
#####################################
sub complemantary {
	my $ntd = $_[0];
	$ntd =~ tr/ACTG/TGAC/;
	return $ntd;
}

#####################################
#####################################
# Generate per patient per gene coverage files
sub GenerateGeneCoverage {
	my ( $titleFile, $outdir, $coverageFile, $exonIntervals ) = @_;
	$logger->info("Gene based coverage files are being created.");
	my ( %barcodes, %txIDs, @header );
	if ( !-e "$outdir/Coverage" ) {
		`mkdir $outdir/Coverage`;
	}

	open IN, "<", $titleFile || die $logger->fatal("Can not open $titleFile: $!");
	while (<IN>) {
		chomp;
		next if (/Barcode/);
		my @line    = split("\t");
		my $barcode = $line[0];
		my $sample  = $line[2];
		$barcodes{$barcode} = $sample;
	}
	close IN;
	open IN2, "<", $exonIntervals || die $logger->fatal("Can not open $exonIntervals: $!");
	while (<IN2>) {
		chomp;
		next if ( /^\@SQ/ || /^\@HD/ );
		my ( $chr, $start, $stop, $strand, $annotation ) = split("\t");
		my ( $gene, $txID, $exon, $aa_range ) = split( ":", $annotation );
		$txIDs{$gene} = $txID;
	}
	close IN2;
	foreach my $barcode ( sort keys %barcodes ) {
		my $sample = $barcodes{$barcode};
		my $file   = $sample . "_gene.coverage.txt";
		open OUT,          ">", "$outdir/Coverage/$file" || die $logger->fatal("Can not open $sample.coverage.txt: $!");
		open COVERAGEFILE, "<", $coverageFile            || die $!;
		my $k = 1;
		while (<COVERAGEFILE>) {
			chomp;
			if (/Gene/) {
				@header = split("\t");
				next;
			}
			my $index;
			foreach my $i ( 0 .. $#header ) {
				if ( $barcode eq $header[$i] ) {
					$index = $i;
				}
			}
			print "$barcode\n" if ( !$index );
			my @line = split("\t");
			my $gene = $line[0];
			my $txID = $txIDs{$gene};
			if ( $gene =~ "CDKN2A" ) {
				$txID = "NM_000077";
			}
			my $cov = sprintf "%.0f", $line[$index];
			my $result;
			if ( $cov >= 100 ) { $result = "GOOD"; }
			else { $result = "POOR"; }
			print OUT "$k\t$gene\t$txID\t$cov\t$result\n";
			$k++;
		}
	}
	$logger->info("Gene based coverage files are created.");
	return;
}

#####################################
#####################################
# Generate per patient per exon coverage files

sub GenerateExonCoverage {
	my ( $filteredAnnoFile, $covFile, $translatedFolder, $exonIntervals, $outdir, $titleFile, $coverage_threshold, $validationExons ) = @_;
	my ( %barcodes, %canonicalExons, @header, $index, %coverage, %variants, %validationExonsHash, %validationExonCoverage );
	$logger->info("Exon based coverage files are being generated.");

	open VALEXONS, "<", $validationExons || die $logger->("Can not open $validationExons: $!");
	while (<VALEXONS>) {
		chomp;
		my ( $gene, $exon ) = split("\t");
		$validationExonsHash{ $gene . ":" . $exon } = 1;
	}
	close VALEXONS;

	open IN, "<", $titleFile || die $logger->fatal("Can not open $titleFile: $!");
	while (<IN>) {
		chomp;
		next if (/Barcode/);
		my @line    = split("\t");
		my $barcode = $line[0];
		my $sample  = $line[2];
		$barcodes{$barcode} = $sample;
	}
	close IN;

	open IN2, "<", $exonIntervals || die $logger->fatal("Can not open $exonIntervals: $!");
	while (<IN2>) {
		chomp;
		next if ( /^\@SQ/ || /^\@HD/ );
		my ( $chr, $start, $stop, $strand, $annotation ) = split("\t");
		my ( $gene, $txID, $exon, $aa_range ) = split( ":", $annotation );
		$canonicalExons{$chr}{ $start . ":" . $stop } = $txID . ":" . $gene . ":" . $exon . ":" . $strand . ":" . $aa_range;
	}
	close IN2;

	#print Dumper \%canonicalExons;
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
		my $gene   = $line[7];
		my $exon   = $line[8];
		push @{ $variants{$sample} }, $gene . ":" . $exon;
	}
	close IN;

	foreach my $barcode ( sort keys %barcodes ) {
		my $sample = $barcodes{$barcode};
		my $file   = $sample . "_exon.coverage.txt";

		#print "$file\n";
		#print "$covFile\n";
		open OUT, ">", "$outdir/Coverage/$file" || die $logger->fatal("Can not open $file: $!");
		$logger->info("Exon coverage file for $sample is being generated.");
		open IN, "<", $covFile || die $!;
		while (<IN>) {
			chomp;
			if (m/Target/) {
				@header = split("\t");
				next;
			}
			foreach my $i ( 0 .. $#header ) {
				if ( $barcode eq $header[$i] ) {
					$index = $i;
				}
			}

			my @line = split("\t");
			my ( $chr,   $range ) = split( ":", $line[1] );
			my ( $start, $stop )  = split( "-", $range );
			my $a = $canonicalExons{$chr}{ $start . ":" . $stop };
			my ( $txID, $gene, $exon, $strand, $aa_range ) = split( ":", $a );
			if ( $gene && $exon && $strand && $txID ) {

				# DC:
				$coverage{"$chr:$gene:$txID:$exon:$strand:$aa_range"} = $range . ":" . $line[$index];
			}
			else {
				die "Something went wrong with annotation";

				#	print OUT "$chr\t$range\t not assigned an exon\n"; #this is for when the region falls into a UTR or an intron
			}
		}

		my %genes = map { $_ => 1 } @{ $variants{$sample} };
		foreach my $key ( sort my_sort2 keys %coverage ) {

			#foreach my $key ( sort keys %coverage ) {
			my ( $chr, $gene, $txID, $exon, $strand, $aa_range ) = split( ":", $key );
			my ( $range, $cov ) = split( ":", $coverage{$key} );
			if ( exists( $validationExonsHash{ $gene . ":" . $exon } ) ) {
				if ( $cov < $coverage_threshold ) {
					push @{ $validationExonCoverage{$sample}{ $gene . " (" . $txID . ") " } }, $exon . " (" . $aa_range . ")";
				}
			}
			if ( exists( $genes{"$gene:$exon"} ) && $cov > $coverage_threshold ) {
				print OUT "$chr\t$gene\t$txID\t$exon\t$range\t$aa_range\t$strand\t$cov\tGOOD\tPOS\n";
			}
			elsif ( exists( $genes{"$gene:$exon"} ) && $cov < $coverage_threshold ) {
				print OUT "$chr\t$gene\t$txID\t$exon\t$range\t$aa_range\t$strand\t$cov\tPOOR\tPOS\n";
			}
			elsif ( !exists( $genes{"$gene:$exon"} ) && $cov > $coverage_threshold ) {
				print OUT "$chr\t$gene\t$txID\t$exon\t$range\t$aa_range\t$strand\t$cov\tGOOD\tNEG\n";
			}
			elsif ( !exists( $genes{"$gene:$exon"} ) && $cov < $coverage_threshold ) {
				print OUT "$chr\t$gene\t$txID\t$exon\t$range\t$aa_range\t$strand\t$cov\tPOOR\tINDET\n";
			}
		}
		$logger->info("Generation of exon coverage file for $sample is completed.");
	}

	close IN;
	close OUT;
	$logger->info("Exon based coverage files are generated.");
	return ( \%validationExonCoverage );
}

sub my_sort {

	#"$gene $cytoband Amplification (Fold Change: $FC)";
	my @f1 = split( " ", $a );
	my @f2 = split( " ", $b );
	my ($fc1) = $a =~ /.*\(Fold Change: (.*)\)/;
	my ($fc2) = $b =~ /.*\(Fold Change: (.*)\)/;
	if ( $fc2 > 0 ) {
		return $f1[4] cmp $f2[4] || $fc2 <=> $fc1;
	}
	else {
		return $f1[4] cmp $f2[4] || $fc1 <=> $fc2;
	}
}

sub my_sort2 {    #DC
	my @f1 = split( ':', $a );
	my @f2 = split( ':', $b );
	$f1[3] =~ s/exon//g;
	$f2[3] =~ s/exon//g;
	if ( $f1[0] eq "X" ) {
		$f1[0] = 23;
	}
	elsif ( $f1[0] eq "Y" ) {
		$f1[0] = 24;
	}
	if ( $f2[0] eq "X" ) {
		$f2[0] = 23;
	}
	elsif ( $f2[0] eq "Y" ) {
		$f2[0] = 24;
	}

	if ( $f1[0] == $f2[0] ) {
		if ( $f1[1] eq $f2[1] ) {
			return ( $f1[3] <=> $f2[3] );
		}
		else {
			return ( $f1[1] cmp $f2[1] );
		}
	}
	else {
		return ( $f1[0] <=> $f2[0] );
	}
}

#####################################
#####################################
#How to use the script.
sub Usage {
	print "Unknow option: @_\n" if (@_);

	print "\nUsage : dmp_create_patient_vcf_file.pl [options]
        [--annotatedFilteredVariants|v  S File containing annotated and filtered mutations (required and submit with full path,Ex:/SomePath/CytoOv_SomaticMutIndel.txt)]
        [--geneLevelCopynumberFile|cn  S File containing gene level copy number results (required)]
        [--intragenicDeletionFile|id    S File containing intergenic loss calls (required)]
        [--geneIntervalAnnotated|gi  S Gene interval file with annotations (txID and cytoband) (required)]
        [--titleFile|t                  S tab-delimited title file for the samples (required and submit with full path)]
        [--sampleSheet|ss                  S Sample Sheet supplied by the technologist (required)]
        [--outdir|o             	    S Path where all the output files will be written (optional) [default:cwd]]
		[--geneCoverageFile|gc			S File containing coverage information for each gene, created in step 4]
		[--exonCoverageFile|ec			S File containing coverage information for each exon, created in step 4]
		[--clinicalExonFile|ce			S File containing the list of clinical exons]
		[--validationExonFile|ve			S File containing the list of validated exons]
		[--CanonicalExonCoverageFile|cc			S File containing coverage information for each canonical exon, created in step 4]
		[--exonIntervalFile|ei			S Interval file  listing all the canonical exons and their start and end positions]
		[--translatedFolder|tf			S Path to the folder where every gene has a translation file]
		[--coverageThreshold|ct			I Threshold for an exon to be called a failure]
		[--normalsUsedInVariantCalling|nu			I The file generated after variant calling, listing all the normals used for VC for each sample]

        
        \n";

	exit;
}

__END__

=head1 NAME

dmp_create_coverage_files.pl - Describe the usage of script briefly

=head1 SYNOPSIS

dmp_create_coverage_files.pl [options] args

      -opt --long      Option description

=head1 DESCRIPTION

Stub documentation for dmp_create_coverage_files.pl, 

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
