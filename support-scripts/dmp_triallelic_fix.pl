#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use MSKCC_DMP_Logger;

my $logger = MSKCC_DMP_Logger->get_logger('IMPACT_Pipeline_Logger');
$logger->start_local();

my ($outdir, $MutationVerboseFilename, $MutationOutFilename);

if (
	@ARGV < 1
	or !GetOptions(
		'output_dir|o:s' => \$outdir,
		'mutect_raw_file|m:s'  => \$MutationVerboseFilename,
		'mutect_vcf_file|v:s'  => \$MutationOutFilename
	)
  )
{
	Usage();
}

if((!$MutationVerboseFilename) or (!$MutationOutFilename)){
	$logger->fatal("Mutect RAW file missing") if (!$MutationVerboseFilename);
	$logger->fatal("Mutect VCF file missing") if (!$MutationOutFilename);
}

	my $MutationVerboseTriAllFilename = $MutationVerboseFilename;# 
	my $MutationOutTriAllFilename = $MutationOutFilename;#
	$MutationVerboseTriAllFilename =~ s/\.txt/_TriAll\.txt/;
	$MutationOutTriAllFilename =~ s/\.vcf/_TriAll\.vcf/;
	my $line;
	my @temp;
	my $chr;
	my $coord;
	my $ref;
	my $alt;
	my $failure_reason;
	my @txt_file;
	my @vcf_file;
	my %triallelic_site;
	my %ATGC = ('A', 'A', 'T', 'T', 'G', 'G', 'C', 'C');

	open(TXT, "$outdir/$MutationVerboseFilename") or die "!ERROR!: Could not open $MutationVerboseFilename for Triallelic_SNV modification (Reading).\n";
	open(VCF, "$outdir/$MutationOutFilename") or die "!ERROR!: Could not open $MutationOutFilename for Triallelic_SNV modification (Reading).\n";

	while($line = <TXT>){
		chomp($line);
		if($line =~ /^#/){
			push(@txt_file, $line);	
			next;
		}
		if($line =~ /^contig/){
			push(@txt_file, $line);	
			next;
		}

		@temp = split("\t", $line);
		$chr = $temp[0];
		$coord = $temp[1];
		$ref = $temp[3];
		$alt = $temp[4];
		$failure_reason = $temp[67];

		if($failure_reason =~ m/triallelic_site/){
			$triallelic_site{$chr.":".$coord.":".$ref.":".$alt} = 1;
			push(@txt_file, join("\t", @temp));	
			delete $ATGC{$ref};
			delete $ATGC{$alt};	

			for (sort keys %ATGC){
				$temp[4] = $_ ; #change alternate allele to nucleotide other than ref and alt
				push(@txt_file, join("\t", @temp));	
				my $debug = join("\t", @temp);
				$logger->debug("Inserting line: $debug");
			}
			
			%ATGC = ('A', 'A', 'T', 'T', 'G', 'G', 'C', 'C');
		}
		else{
			push(@txt_file, join("\t", @temp));
		}
	}

	close(TXT);

	while($line = <VCF>){
		chomp($line);
		if($line =~ /^#/){
			push(@vcf_file, $line);	
			next;
		}

		@temp = split("\t", $line);
		$chr = $temp[0];
		$coord = $temp[1];
		$ref = $temp[3];
		$alt = $temp[4];

		if(exists $triallelic_site{$chr.":".$coord.":".$ref.":".$alt}){
			push(@vcf_file, join("\t", @temp));	
			delete $ATGC{$ref};
			delete $ATGC{$alt};	

			for (sort keys %ATGC){
				$temp[4] = $_; #change alternate allele to nucleotide other than ref and alt
				push(@vcf_file, join("\t", @temp));	
				my $debug = join("\t", @temp);
				$logger->debug("Inserting line: $debug");
			}
			
			%ATGC = ('A', 'A', 'T', 'T', 'G', 'G', 'C', 'C');
		}
		else{
			push(@vcf_file, join("\t", @temp));	
		}
	}

	close(VCF);

	open(TXT, ">$outdir/$MutationVerboseTriAllFilename") or die "!ERROR!: Could not open $MutationVerboseFilename for Triallelic_SNV modification (Writing).\n";
	open(VCF, ">$outdir/$MutationOutTriAllFilename") or die "!ERROR!: Could not open $MutationOutFilename for Triallelic_SNV modification (Writing).\n";

	foreach (@txt_file){
		print TXT $_."\n";
	}

	close(TXT);

	foreach (@vcf_file){
		print VCF $_."\n";
	}

	close(VCF);


	if((-s "$outdir/$MutationVerboseTriAllFilename") <  (-s "$outdir/$MutationVerboseFilename")){
		$logger->warn("Triallelic site fix: Something went wrong here. File size after triallelic site modification decreased.");
	}
	else{
		`rm $outdir/$MutationVerboseFilename`;
		`mv $outdir/$MutationVerboseTriAllFilename $outdir/$MutationVerboseFilename`;
	}

	if((-s "$outdir/$MutationOutTriAllFilename") <  (-s "$outdir/$MutationOutFilename")){
		$logger->warn("Triallelic site fix: Something went wrong here. File size after triallelic site modification decreased.");
	}
	else{
		`rm $outdir/$MutationOutFilename`;
		`mv $outdir/$MutationOutTriAllFilename $outdir/$MutationOutFilename`;
	}



#####################################
#####################################
#How to use the script.
sub Usage {
	print "Unknow option: @_\n" if (@_);

	print "\nUsage : $0 [options]
        [--output_dir|o 			S Path where all the output files will be written (optional) [default:cwd]]
        [--mutect_raw_file|m 		S File name of raw Mutect output (TXT)]
        [--mutect_vcf_file|v 		S File name of raw Mutect output (VCF)]
        \n";

	exit;
}