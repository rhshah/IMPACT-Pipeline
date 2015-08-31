#!/usr/bin/perl
# compareResults_try1.pl ---
# Author: Shah <shahr2@phoenix-h1>
# Created: 18 Nov 2013
# Version: 0.01
use strict;
use Getopt::Long;
use Cwd;
my (
	 $stdfile1,       $stdfile2,       $annExofile1,     $annExofile2,
	 $annStfile1,     $annStfile2,     $annExoDropfile1, $annExoDropfile2,
	 $annStDropfile1, $annStDropfile2, $outfileprefix,         $outdir,$prefix1,$prefix2
);
if (
	 @ARGV < 2
	 or !GetOptions(
					 'prefix1|p1:s' => \$prefix1,
					 'prefix2|p2:s' => \$prefix2,
					 'outfileprefix|of:s'            => \$outfileprefix,
					 'outdir|o:s'              => \$outdir
	 )
  )
{
	Usage();
}
$stdfile1 = $prefix1 . "_AllSomaticMutIndel_withAlleleDepth_annovarAnnotated.txt";
$stdfile2 = $prefix2 . "_AllSomaticMutIndel_withAlleleDepth_annovarAnnotated.txt";
$annExofile1 = $prefix1 . "_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedExonic.Filtered.txt";
$annExofile2 = $prefix2 . "_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedExonic.Filtered.txt";
$annStfile1 = $prefix1 . "_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedSilent.Filtered.txt";
$annStfile2 = $prefix2 . "_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedSilent.Filtered.txt";
$annExoDropfile1 = $prefix1 . "_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedExonic.Dropped.txt";
$annExoDropfile2 = $prefix2 . "_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedExonic.Dropped.txt";
$annStDropfile1 = $prefix1 . "_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedSilent.Dropped.txt";
$annStDropfile2 = $prefix2 . "_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedSilent.Dropped.txt";

if ( (! -e $stdfile1 ) or (! -e $stdfile2 ) )
{
	print "Input $stdfile1 & $stdfile2 file have not been defined. See Usage\n";
	Usage();
	exit;
}
if ( (! -e $annExofile1 ) or (! -e $annExofile2 ) )
{
	print "Input file $annExofile1 & $annExofile2 have not been defined. See Usage\n";
	Usage();
	exit;
}
if ( (! -e $annExoDropfile1 ) or (! -e $annExoDropfile2 ) )
{
	print "Input file $annExoDropfile1 & $annExoDropfile2 have not been defined. See Usage\n";
	Usage();
	exit;
}
if ( (! -e $annStfile1 ) or (! -e $annStfile2 ) )
{
	print "Input file $annStfile1 & $annStfile2 have not been defined. See Usage\n";
	Usage();
	exit;
}
if ( (! -e $annStDropfile1 ) or (! -e $annStDropfile2 ) )
{
	print "Input file $annStDropfile1 & $annStDropfile2 have not been defined. See Usage\n";
	Usage();
	exit;
}
if ( !defined $outfileprefix )
{
	print
"Output file is not specified, default file will be used as output file.\n";
	$outfileprefix = "VariantFileCompareReport";
}
if ( !defined $outdir )
{
	print
"Output directory is not specified, current working directory will be used as output directory.\n";
	$outdir = getcwd;
}
open OUT, ">", "$outdir/$outfileprefix\_fullStats.txt" || die "Cannot Open $outfileprefix\_fullStats.txt,$!";
open OUTCOMMON, ">", "$outdir/$outfileprefix\_CommonExonicFilterDifferences.txt"
  || die "Cannot Open $outdir/$outfileprefix\_CommonExonicFilterDifferences.txt,$!";
open OUTCOMMONMUTATIONS, ">", "$outdir/$outfileprefix\_CommonMutations.txt"
  || die "Cannot Open $outdir/$outfileprefix\_CommonMutation.txt,$!";
open OUTUNIQMUTATIONSFILE1, ">", "$outdir/$outfileprefix\_UniqueMutationsFile1.txt"
  || die "Cannot Open $outdir/$outfileprefix\_UniqueMutationsFile1.txt,$!";
open OUTUNIQMUTATIONSFILE2, ">", "$outdir/$outfileprefix\_UniqueMutationsFile2.txt"
  || die "Cannot Open $outdir/$outfileprefix\_UniqueMutationsFile2.txt,$!";
open OUTUNIQ1, ">", "$outdir/$outfileprefix\_UniqueExonicFilterFile1.txt"
  || die "Cannot Open $outdir/$outfileprefix\_UniqueExonicFilterFile1.txt,$!";
open OUTUNIQ2, ">", "$outdir/$outfileprefix\_UniqueExonicFilterFile2.txt"
  || die "Cannot Open $outdir/$outfileprefix\_UniqueExonicFilterFile2.txt,$!";
print OUT "###File 1: $stdfile1\n###File 2: $stdfile2\n";
my ( $stdhash1, $stdrealCols1, $stdrows1, $stdcols1 ) = ReadFile($stdfile1);
my ( $stdhash2, $stdrealCols2, $stdrows2, $stdcols2 ) = ReadFile($stdfile2);
my ( $annExohash1, $annExorealCols1, $annExorows1, $annExocols1 ) =
  ReadFile($annExofile1);
my ( $annExohash2, $annExorealCols2, $annExorows2, $annExocols2 ) =
  ReadFile($annExofile2);
my ( $annExoDrophash1, $annExoDroprealCols1, $annExoDroprows1,
	 $annExoDropcols1 ) = ReadFile($annExoDropfile1);
my ( $annExoDrophash2, $annExoDroprealCols2, $annExoDroprows2, $annExoDrocols2 )
  = ReadFile($annExoDropfile2);
my ( $annSthash1, $annStrealCols1, $annStrows1, $annStcols1 ) =
  ReadFile($annStfile1);
my ( $annSthash2, $annStrealCols2, $annStrows2, $annStcols2 ) =
  ReadFile($annStfile2);
my ( $annStDrophash1, $annStDroprealCols1, $annStDroprows1, $annStDropcols1 ) =
  ReadFile($annStDropfile1);
my ( $annStDrophash2, $annStDroprealCols2, $annStDroprows2, $annStDrocols2 ) =
  ReadFile($annStDropfile2);
my ( $presentInFile1, $presentInFile2, $presentInBoth ) =
  CmpStdFilter( $stdhash1, $stdrows1, $stdhash2, $stdrows2 );
CmpUniqueRecordswithFinalResults();
CmpCommonRecordswithFinalResults();
close(OUT);
close(OUTCOMMON);
close(OUTCOMMONMUTATIONS);
close(OUTUNIQ1);
close(OUTUNIQMUTATIONSFILE1);
close(OUTUNIQ2);
close(OUTUNIQMUTATIONSFILE2);
exit;
#####################################
#####################################
#How to use the script.
sub Usage
{
	print "Unknow option: @_\n" if (@_);
	print "\nUsage : compareResults.pl [options]
        [--prefix1|p1        S prefix for run 1 havaing variants (required & will be denoted as file1)]
        [--prefix2|p2        S prefix for run 2 having variants (required & will be denoted as file2)]
        [--outfile|of      S Name where of the output file (optional) [default:VariantFileCompareReport.txt]]
        [--outdir|o        S Path where all the output file will be written (optional) [default:cwd]]
	\n";
	print "For any question please contact Ronak Shah (shahr2\@mskcc.org)\n";
	print "!~!~! Comments, flames and suggestion are welcome !~!~!\n";
	exit;
}
#####################################
#####################################
#Read the file
sub ReadFile
{
	my ($file) = @_;
	my (%variantsHash);
	my $rows            = 0;
	my $cols            = 0;
	my (@line)          = ();
	my (@processedLine) = ();
	my (
		 $sample, $normal, $chr,          $start,
		 $ref,    $alt,    $varaintClass, $gene,
		 $exon,   $cdna,   $aachange,     $failureReason,
		 $ndp,    $nrc,    $nac,          $naf,
		 $tdp,    $trc,    $tac,          $taf,
		 $trefp,  $trefn,  $taltp,        $taltn,
		 $allN,   $tnfreq, $occ,          $pct
	);
	my $type;
	open IN, "<", $file || die $!;

	while (<IN>)
	{
		$rows++;
		chomp;
		if ( $. == 1 )
		{
			my @HEADER = split( "\t", $_ );
			$cols = scalar(@HEADER);
			next;
		}
		my $null        = 2;
		my $val         = 1;
		my (@inferLine) = ();
		@line = split("\t");
		foreach my $data (@line)
		{
			$data =~ s/\s//g;
			if ( !$data )
			{
				push( @inferLine, $null );
			} else
			{
				push( @inferLine, $val );
			}
		}
		my $data = join( ":", @inferLine );
		push( @processedLine, $data );
		$sample = $line[0];
		$normal = $line[1];
		$chr    = $line[2];
		$start  = $line[3];
		$ref    = $line[4];
		$ref =~ s/\s//g;
		$alt = $line[5];
		$alt =~ s/\s//g;
		$varaintClass  = $line[6];
		$gene          = $line[7];
		$exon          = $line[8];
		$cdna          = $line[12];
		$aachange      = $line[13];
		$failureReason = $line[17];
		$ndp           = $line[18];
		$nrc           = $line[19];
		$nac           = $line[20];
		$naf           = $line[21];
		$tdp           = $line[22];
		$trc           = $line[23];
		$tac           = $line[24];
		$taf           = $line[25];
		$trefp         = $line[26];
		$trefn         = $line[27];
		$taltp         = $line[28];
		$taltn         = $line[29];
		$allN          = $line[31];
		$tnfreq        = $line[32];
		( $occ, $pct ) = split( ";", $line[33] );

		if ( length($ref) == length($alt) )
		{
			$type = "SNP";
		}
		if ( length($ref) > length($alt) )
		{
			$type = "DEL";
		}
		if ( length($ref) < length($alt) )
		{
			$type = "INS";
		}
		$variantsHash{"$sample:$normal:$chr:$start:$ref:$alt:$type"} = "$_";
	}
	close(IN);
	return ( \%variantsHash, \@processedLine, $rows, $cols );
}
#####################################
#####################################
#Read the file
sub ReadStdFile
{
	my ($file) = @_;
	my (%variantsHash);
	my $rows            = 0;
	my $cols            = 0;
	my (@line)          = ();
	my (@processedLine) = ();
	my ( $sample, $normal, $chr, $start, $ref, $alt, $failureReason );
	my $type;
	open IN, "<", $file || die $!;

	while (<IN>)
	{
		$rows++;
		chomp;
		if ( $. == 1 )
		{
			my @HEADER = split( "\t", $_ );
			$cols = scalar(@HEADER);
			next;
		}
		my $null        = 2;
		my $val         = 1;
		my (@inferLine) = ();
		@line = split("\t");
		foreach my $data (@line)
		{
			$data =~ s/\s//g;
			if ( !$data )
			{
				push( @inferLine, $null );
			} else
			{
				push( @inferLine, $val );
			}
		}
		my $data = join( ":", @inferLine );
		push( @processedLine, $data );
		$sample        = $line[0];
		$normal        = $line[1];
		$chr           = $line[2];
		$start         = $line[3];
		$ref           = $line[4];
		$alt           = $line[5];
		$failureReason = $line[6];

		if ( length($ref) == length($alt) )
		{
			$type = "SNP";
		}
		if ( length($ref) > length($alt) )
		{
			$type = "DEL";
		}
		if ( length($ref) < length($alt) )
		{
			$type = "INS";
		}
		$variantsHash{"$sample:$normal:$chr:$start:$ref:$alt:$type"} =
		  "$failureReason";
	}
	close(IN);
	return ( \%variantsHash, \@processedLine, $rows, $cols );
}

sub CmpStdFilter
{
	my ( $hash1, $rows1, $hash2, $rows2 ) = @_;
	my (%vh1) = %$hash1;
	my (%vh2) = %$hash2;
	print OUT "###Variant Matching Statistics###\n";
	my %vh1only          = ();
	my %vh2only          = ();
	my %vh_1_2           = ();
	my $totalCount_1     = 0;
	my $totalCount_snp_1 = 0;
	my $totalCount_del_1 = 0;
	my $totalCount_ins_1 = 0;
	my $totalCount_2     = 0;
	my $totalCount_snp_2 = 0;
	my $totalCount_del_2 = 0;
	my $totalCount_ins_2 = 0;
	my $count_1          = 0;
	my $count_snp_1      = 0;
	my $count_del_1      = 0;
	my $count_ins_1      = 0;
	my $count_1_2        = 0;
	my $count_snp_1_2    = 0;
	my $count_del_1_2    = 0;
	my $count_ins_1_2    = 0;
	my $count_2          = 0;
	my $count_snp_2      = 0;
	my $count_del_2      = 0;
	my $count_ins_2      = 0;
	my $count_2_1        = 0;
	my $count_snp_2_1    = 0;
	my $count_del_2_1    = 0;
	my $count_ins_2_1    = 0;

	#check in file 1
	while ( my ( $keys1, $values1 ) = each(%vh1) )
	{
		$totalCount_1++;
		if ( $keys1 =~ /SNP/ )
		{
			$totalCount_snp_1++;
		}
		if ( $keys1 =~ /DEL/ )
		{
			$totalCount_del_1++;
		}
		if ( $keys1 =~ /INS/ )
		{
			$totalCount_ins_1++;
		}
		if ( exists $vh2{$keys1} )
		{
			print OUTCOMMONMUTATIONS "$values1\n";
			$count_1_2++;
			if ( $keys1 =~ /SNP/ )
			{
				$count_snp_1_2++;
			}
			if ( $keys1 =~ /DEL/ )
			{
				$count_del_1_2++;
			}
			if ( $keys1 =~ /INS/ )
			{
				$count_ins_1_2++;
			}
			$vh_1_2{$keys1} = $values1;
		} else
		{
			$vh1only{$keys1} = $values1;
			if ( $keys1 =~ /SNP/ )
			{
				$count_snp_1++;
			}
			if ( $keys1 =~ /DEL/ )
			{
				$count_del_1++;
			}
			if ( $keys1 =~ /INS/ )
			{
				$count_ins_1++;
			}
			print OUTUNIQMUTATIONSFILE1 "$values1\n";
			$count_1++;
		}
	}

	#Check in file 2
	while ( my ( $keys2, $values2 ) = each(%vh2) )
	{
		$totalCount_2++;
		if ( $keys2 =~ /SNP/ )
		{
			$totalCount_snp_2++;
		}
		if ( $keys2 =~ /DEL/ )
		{
			$totalCount_del_2++;
		}
		if ( $keys2 =~ /INS/ )
		{
			$totalCount_ins_2++;
		}
		if ( exists $vh1{$keys2} )
		{
			$count_2_1++;
			if ( $keys2 =~ /SNP/ )
			{
				$count_snp_2_1++;
			}
			if ( $keys2 =~ /DEL/ )
			{
				$count_del_2_1++;
			}
			if ( $keys2 =~ /INS/ )
			{
				$count_ins_2_1++;
			}
			$vh_1_2{$keys2} = $values2;
		} else
		{
			$vh2only{$keys2} = $values2;
			print OUTUNIQMUTATIONSFILE2 "$values2\n";
			$count_2++;
			if ( $keys2 =~ /SNP/ )
			{
				$count_snp_2++;
			}
			if ( $keys2 =~ /DEL/ )
			{
				$count_del_2++;
			}
			if ( $keys2 =~ /INS/ )
			{
				$count_ins_2++;
			}
		}
	}
	print OUT "##Number of Total Variants in File 1: $totalCount_1\n";
	print OUT "Number of Total SNP Variants in File 1: $totalCount_snp_1\n";
	print OUT "Number of Total DEL Variants in File 1: $totalCount_del_1\n";
	print OUT "Number of Total INS Variants in File 1: $totalCount_ins_1\n";
	print OUT "##Number of Total Variants in File 2: $totalCount_2\n";
	print OUT "Number of Total SNP Variants in File 2: $totalCount_snp_2\n";
	print OUT "Number of Total DEL Variants in File 2: $totalCount_del_2\n";
	print OUT "Number of Total INS Variants in File 2: $totalCount_ins_2\n";
	print OUT "##Number of Variants Present in File 1 & File 2: $count_1_2\n";
	print OUT "Number of SNP Variants in File 1 & File 2: $count_snp_1_2\n";
	print OUT "Number of DEL Variants in File 1 & File 2: $count_del_1_2\n";
	print OUT "Number of INS Variants in File 1 & File 2: $count_ins_1_2\n";
	print OUT "##Number of Variants Present in File 2 & File 1: $count_2_1\n";
	print OUT "Number of SNP Variants in File 2 & File 1: $count_snp_2_1\n";
	print OUT "Number of DEL Variants in File 2 & File 1: $count_del_2_1\n";
	print OUT "Number of INS Variants in File 2 & File 1: $count_ins_2_1\n";
	print OUT "##Number of Variants Present only in File 1: $count_1\n";
	print OUT "Number of SNP Variants only in File 1: $count_snp_1\n";
	print OUT "Number of DEL Variants only in File 1: $count_del_1\n";
	print OUT "Number of INS Variants only in File 1: $count_ins_1\n";
	print OUT "##Number of Variants Present only in File 2: $count_2\n";
	print OUT "Number of SNP Variants only in File 2: $count_snp_2\n";
	print OUT "Number of Total DEL Variants only in File 2: $count_del_2\n";
	print OUT "Number of Total INS Variants only in File 2: $count_ins_2\n";
	return ( \%vh1only, \%vh2only, \%vh_1_2 );
}

sub CmpUniqueRecordswithFinalResults
{
	my $countInExoFtfile1    = 0;
	my $countInExoFtsnpfile1 = 0;
	my $countInExoFtdelfile1 = 0;
	my $countInExoFtinsfile1 = 0;
	my $countInExoFtfile2    = 0;
	my $countInExoFtsnpfile2 = 0;
	my $countInExoFtdelfile2 = 0;
	my $countInExoFtinsfile2 = 0;
	my $countInExoDpfile1    = 0;
	my $countInExoDpsnpfile1 = 0;
	my $countInExoDpdelfile1 = 0;
	my $countInExoDpinsfile1 = 0;
	my $countInExoDpfile2    = 0;
	my $countInExoDpsnpfile2 = 0;
	my $countInExoDpdelfile2 = 0;
	my $countInExoDpinsfile2 = 0;
	my $countInStFtfile1     = 0;
	my $countInStFtsnpfile1  = 0;
	my $countInStFtdelfile1  = 0;
	my $countInStFtinsfile1  = 0;
	my $countInStFtfile2     = 0;
	my $countInStFtsnpfile2  = 0;
	my $countInStFtdelfile2  = 0;
	my $countInStFtinsfile2  = 0;
	my $countInStDpfile1     = 0;
	my $countInStDpsnpfile1  = 0;
	my $countInStDpdelfile1  = 0;
	my $countInStDpinsfile1  = 0;
	my $countInStDpfile2     = 0;
	my $countInStDpsnpfile2  = 0;
	my $countInStDpdelfile2  = 0;
	my $countInStDpinsfile2  = 0;
	my %cmpExoFtfile1        = ();
	my %cmpExoFtfile2        = ();
	print OUT "###Unique Variant Statistics###\n";

	#Compare File1 with Exonic Filtered for File 1
	#print OUT "File 1 Records Going to Exonic Filtered of File 1:\n";
	while ( my ( $keys1, $values1 ) = each(%$presentInFile1) )
	{
		if ( exists $$annExohash1{$keys1} )
		{
			my $newval = $$annExohash1{$keys1};
			#my $record = "$keys1:$newval";
			#$record =~ s/:/\t/g;
			print OUTUNIQ1 "$newval\n";
			$countInExoFtfile1++;
			if ( $keys1 =~ /SNP/ )
			{
				$countInExoFtsnpfile1++;
			}
			if ( $keys1 =~ /DEL/ )
			{
				$countInExoFtdelfile1++;
			}
			if ( $keys1 =~ /INS/ )
			{
				$countInExoFtinsfile1++;
			}
		}
	}
	print OUT
"##Number of Records of File 1 Going to Exonic Filtered of File 1: $countInExoFtfile1\n";
	print OUT
"Number of SNP Records of File 1 Going to Exonic Filtered of File 1: $countInExoFtsnpfile1\n";
	print OUT
"Number of DEL Records of File 1 Going to Exonic Filtered of File 1: $countInExoFtdelfile1\n";
	print OUT
"Number of INS Records of File 1 Going to Exonic Filtered of File 1: $countInExoFtinsfile1\n";

	#Compare File 2 with Exonic Filtered for File 2
	#print OUT "File 2 Records Going to Exonic Filtered of File 2:\n";
	while ( my ( $keys2, $values2 ) = each(%$presentInFile2) )
	{
		if ( exists $$annExohash2{$keys2} )
		{
			my $newval = $$annExohash2{$keys2};
			#my $record = "$keys2:$newval";
			#$record =~ s/:/\t/g;
			print OUTUNIQ2 "$newval\n";
			$countInExoFtfile2++;
			if ( $keys2 =~ /SNP/ )
			{
				$countInExoFtsnpfile2++;
			}
			if ( $keys2 =~ /DEL/ )
			{
				$countInExoFtdelfile2++;
			}
			if ( $keys2 =~ /INS/ )
			{
				$countInExoFtinsfile2++;
			}
		}
	}
	print OUT
"##Number of Records of File 2 Going to Exonic Filtered of File 2: $countInExoFtfile2\n";
	print OUT
"Number of SNP Records of File 2 Going to Exonic Filtered of File 2: $countInExoFtsnpfile2\n";
	print OUT
"Number of DEL Records of File 2 Going to Exonic Filtered of File 2: $countInExoFtdelfile2\n";
	print OUT
"Number of INS Records of File 2 Going to Exonic Filtered of File 2: $countInExoFtinsfile2\n";

	#Compare File 1 with Exonic Dropped for File 1
	#print OUT "File 1 Records Going to Exonic Dropped of File 1:\n";
	while ( my ( $keys1, $values1 ) = each(%$presentInFile1) )
	{
		if ( exists $$annExoDrophash1{$keys1} )
		{
			my $newval = $$annExoDrophash1{$keys1};

			#print OUT "$keys1:$newval\n";
			$countInExoDpfile1++;
			if ( $keys1 =~ /SNP/ )
			{
				$countInExoDpsnpfile1++;
			}
			if ( $keys1 =~ /DEL/ )
			{
				$countInExoDpdelfile1++;
			}
			if ( $keys1 =~ /INS/ )
			{
				$countInExoDpinsfile1++;
			}
		}
	}
	print OUT
"##Number of Records of File 1 Going to Exonic Dropped of File 1: $countInExoDpfile1\n";
	print OUT
"Number of SNP Records of File 1 Going to Exonic Dropped of File 1: $countInExoDpsnpfile1\n";
	print OUT
"Number of DEL Records of File 1 Going to Exonic Dropped of File 1: $countInExoDpdelfile1\n";
	print OUT
"Number of INS Records of File 1 Going to Exonic Dropped of File 1: $countInExoDpinsfile1\n";

	#Compare File 2 with Exonic Filtered for File 2
	#print OUT "File 2 Records Going to Exonic Dropped of File 2:\n";
	while ( my ( $keys2, $values2 ) = each(%$presentInFile2) )
	{
		if ( exists $$annExoDrophash2{$keys2} )
		{
			my $newval = $$annExoDrophash2{$keys2};

			#print OUT "$keys2:$newval\n";
			$countInExoDpfile2++;
			if ( $keys2 =~ /SNP/ )
			{
				$countInExoDpsnpfile2++;
			}
			if ( $keys2 =~ /DEL/ )
			{
				$countInExoDpdelfile2++;
			}
			if ( $keys2 =~ /INS/ )
			{
				$countInExoDpinsfile2++;
			}
		}
	}
	print OUT
"##Number of Records of File 2 Going to Exonic Dropped of File 2: $countInExoDpfile2\n";
	print OUT
"Number of SNP Records of File 2 Going to Exonic Dropped of File 2: $countInExoDpsnpfile2\n";
	print OUT
"Number of DEL Records of File 2 Going to Exonic Dropped of File 2: $countInExoDpdelfile2\n";
	print OUT
"Number of INS Records of File 2 Going to Exonic Dropped of File 2: $countInExoDpinsfile2\n";

	#Compare File1 with Silent Filtered for File 1
	#print OUT "File 1 Records Going to Silent Filtered of File 1:\n";
	while ( my ( $keys1, $values1 ) = each(%$presentInFile1) )
	{
		if ( exists $$annSthash1{$keys1} )
		{
			my $newval = $$annSthash1{$keys1};

			#print OUT "$keys1:$newval\n";
			$countInStFtfile1++;
			if ( $keys1 =~ /SNP/ )
			{
				$countInStFtsnpfile1++;
			}
			if ( $keys1 =~ /DEL/ )
			{
				$countInStFtdelfile1++;
			}
			if ( $keys1 =~ /INS/ )
			{
				$countInStFtinsfile1++;
			}
		}
	}
	print OUT
"##Number of Records of File 1 Going to Silent Filtered of File 1: $countInStFtfile1\n";
	print OUT
"Number of SNP Records of File 1 Going to Silent Filtered of File 1: $countInStFtsnpfile1\n";
	print OUT
"Number of DEL Records of File 1 Going to Silent Filtered of File 1: $countInStFtdelfile1\n";
	print OUT
"Number of INS Records of File 1 Going to Silent Filtered of File 1: $countInStFtinsfile1\n";

	#Compare File1 with Silent Filtered for File 2
	#print OUT "File 1 Records Going to Silent Filtered of File 2:\n";
	while ( my ( $keys2, $values2 ) = each(%$presentInFile2) )
	{
		if ( exists $$annSthash2{$keys2} )
		{
			my $newval = $$annSthash2{$keys2};

			#print OUT "$keys2:$newval\n";
			$countInStFtfile2++;
			if ( $keys2 =~ /SNP/ )
			{
				$countInStFtsnpfile2++;
			}
			if ( $keys2 =~ /DEL/ )
			{
				$countInStFtdelfile2++;
			}
			if ( $keys2 =~ /INS/ )
			{
				$countInStFtinsfile2++;
			}
		}
	}
	print OUT
"Number of Records of File 2 Going to Silent Filtered of File 2: $countInStFtfile2\n";
	print OUT
"Number of SNP Records of File 2 Going to Silent Filtered of File 2: $countInStFtsnpfile2\n";
	print OUT
"Number of DEL Records of File 2 Going to Silent Filtered of File 2: $countInStFtdelfile2\n";
	print OUT
"Number of INS Records of File 2 Going to Silent Filtered of File 2: $countInStFtinsfile2\n";

	#Compare File 1 with Silent Dropped for File 1
	#print OUT "File 1 Records Going to Silent Dropped of File 1:\n";
	while ( my ( $keys1, $values1 ) = each(%$presentInFile1) )
	{
		if ( exists $$annStDrophash1{$keys1} )
		{
			my $newval = $$annStDrophash1{$keys1};

			#print OUT "$keys1:$newval\n";
			$countInStDpfile1++;
			if ( $keys1 =~ /SNP/ )
			{
				$countInStDpsnpfile1++;
			}
			if ( $keys1 =~ /DEL/ )
			{
				$countInStDpdelfile1++;
			}
			if ( $keys1 =~ /INS/ )
			{
				$countInStDpinsfile1++;
			}
		}
	}
	print OUT
"##Number of Records of File 1 Going to Silent Dropped of File 1: $countInStDpfile1\n";
	print OUT
"Number of SNP Records of File 1 Going to Silent Dropped of File 1: $countInStDpsnpfile1\n";
	print OUT
"Number of DEL Records of File 1 Going to Silent Dropped of File 1: $countInStDpdelfile1\n";
	print OUT
"Number of INS Records of File 1 Going to Silent Dropped of File 1: $countInStDpinsfile1\n";

	#Compare File 2 with Silent Dropped for File 2
	#print OUT "File 2 Records Going to Silent Dropped of File 2:\n";
	while ( my ( $keys2, $values2 ) = each(%$presentInFile2) )
	{
		if ( exists $$annStDrophash2{$keys2} )
		{
			my $newval = $$annStDrophash2{$keys2};

			#print OUT "$keys2:$newval\n";
			$countInStDpfile2++;
			if ( $keys2 =~ /SNP/ )
			{
				$countInStDpsnpfile2++;
			}
			if ( $keys2 =~ /DEL/ )
			{
				$countInStDpdelfile2++;
			}
			if ( $keys2 =~ /INS/ )
			{
				$countInStDpinsfile2++;
			}
		}
	}
	print OUT
"##Number of Records of File 2 Going to Silent Dropped of File 2: $countInStDpfile2\n";
	print OUT
"Number of SNP Records of File 2 Going to Silent Dropped of File 2: $countInStDpsnpfile2\n";
	print OUT
"Number of DEL Records of File 2 Going to Silent Dropped of File 2: $countInStDpdelfile2\n";
	print OUT
"Number of INS Records of File 2 Going to Silent Dropped of File 2: $countInStDpinsfile2\n";
}

sub CmpCommonRecordswithFinalResults
{
	my $countInExoFtfile1    = 0;
	my $countInExoFtsnpfile1 = 0;
	my $countInExoFtdelfile1 = 0;
	my $countInExoFtinsfile1 = 0;
	my $countInExoFtfile2    = 0;
	my $countInExoFtsnpfile2 = 0;
	my $countInExoFtdelfile2 = 0;
	my $countInExoFtinsfile2 = 0;
	my $countInExoDpfile1    = 0;
	my $countInExoDpsnpfile1 = 0;
	my $countInExoDpdelfile1 = 0;
	my $countInExoDpinsfile1 = 0;
	my $countInExoDpfile2    = 0;
	my $countInExoDpsnpfile2 = 0;
	my $countInExoDpdelfile2 = 0;
	my $countInExoDpinsfile2 = 0;
	my $countInStFtfile1     = 0;
	my $countInStFtsnpfile1  = 0;
	my $countInStFtdelfile1  = 0;
	my $countInStFtinsfile1  = 0;
	my $countInStFtfile2     = 0;
	my $countInStFtsnpfile2  = 0;
	my $countInStFtdelfile2  = 0;
	my $countInStFtinsfile2  = 0;
	my $countInStDpfile1     = 0;
	my $countInStDpsnpfile1  = 0;
	my $countInStDpdelfile1  = 0;
	my $countInStDpinsfile1  = 0;
	my $countInStDpfile2     = 0;
	my $countInStDpsnpfile2  = 0;
	my $countInStDpdelfile2  = 0;
	my $countInStDpinsfile2  = 0;
	my %cmpExoFtfile1        = ();
	my %cmpExoFtfile2        = ();
	print OUT "###Common Variant Statistics###\n";

	#Compare File1 with Exonic Filtered for File 1
	#print OUT "File 1 Records Going to Exonic Filtered of File 1:\n";
	while ( my ( $keys1, $values1 ) = each(%$presentInBoth) )
	{
		if ( exists $$annExohash1{$keys1} )
		{
			my $newval = $$annExohash1{$keys1};
			$cmpExoFtfile1{$keys1} = $newval;

			#print OUT "$keys1:$newval\n";
			$countInExoFtfile1++;
			if ( $keys1 =~ /SNP/ )
			{
				$countInExoFtsnpfile1++;
			}
			if ( $keys1 =~ /DEL/ )
			{
				$countInExoFtdelfile1++;
			}
			if ( $keys1 =~ /INS/ )
			{
				$countInExoFtinsfile1++;
			}
		}
	}
	print OUT
"##Number of Records of File 1 Going to Exonic Filtered of File 1: $countInExoFtfile1\n";
	print OUT
"Number of SNP Records of File 1 Going to Exonic Filtered of File 1: $countInExoFtsnpfile1\n";
	print OUT
"Number of DEL Records of File 1 Going to Exonic Filtered of File 1: $countInExoFtdelfile1\n";
	print OUT
"Number of INS Records of File 1 Going to Exonic Filtered of File 1: $countInExoFtinsfile1\n";

	#Compare File 2 with Exonic Filtered for File 2
	#print OUT "File 2 Records Going to Exonic Filtered of File 2:\n";
	while ( my ( $keys2, $values2 ) = each(%$presentInBoth) )
	{
		if ( exists $$annExohash2{$keys2} )
		{
			my $newval = $$annExohash2{$keys2};
			$cmpExoFtfile2{$keys2} = $newval;

			#print OUT "$keys2:$newval\n";
			$countInExoFtfile2++;
			if ( $keys2 =~ /SNP/ )
			{
				$countInExoFtsnpfile2++;
			}
			if ( $keys2 =~ /DEL/ )
			{
				$countInExoFtdelfile2++;
			}
			if ( $keys2 =~ /INS/ )
			{
				$countInExoFtinsfile2++;
			}
		}
	}
	print OUT
"##Number of Records of File 2 Going to Exonic Filtered of File 2: $countInExoFtfile2\n";
	print OUT
"Number of SNP Records of File 2 Going to Exonic Filtered of File 2: $countInExoFtsnpfile2\n";
	print OUT
"Number of DEL Records of File 2 Going to Exonic Filtered of File 2: $countInExoFtdelfile2\n";
	print OUT
"Number of INS Records of File 2 Going to Exonic Filtered of File 2: $countInExoFtinsfile2\n";
	if ( $countInExoFtfile1 != $countInExoFtfile2 )
	{
		
		my @common = ();
		my @uncommon = ();
		foreach ( keys %cmpExoFtfile1 )
		{
			if (exists $cmpExoFtfile2{$_}){
				push( @common, $_ );
			}
			else{
				push( @uncommon, $_ );
			}
		}
		foreach ( keys %cmpExoFtfile2 )
		{
			if (exists $cmpExoFtfile1{$_}){
				push( @common, $_ );
			}
			else{
				push( @uncommon, $_ );
			}
		}
		foreach my $ckey (@uncommon)
		{
			if ( not exists $cmpExoFtfile1{$ckey} )
			{
				my $record = $$annExohash2{$ckey};
				#$record = $ckey . ":" . $record;
				#$record =~ s/:/\t/g;
				print OUTCOMMON "File2\t$record\n";
			}
			if ( not exists $cmpExoFtfile2{$ckey} )
			{
				my $record = $$annExohash1{$ckey};
				#$record = $ckey . ":" . $record;
				#$record =~ s/:/\t/g;
				print OUTCOMMON "File1\t$record\n";
			}
		}
	}

	#Compare File 1 with Exonic Dropped for File 1
	#print OUT "File 1 Records Going to Exonic Dropped of File 1:\n";
	while ( my ( $keys1, $values1 ) = each(%$presentInBoth) )
	{
		if ( exists $$annExoDrophash1{$keys1} )
		{
			my $newval = $$annExoDrophash1{$keys1};

			#print OUT "$keys1:$newval\n";
			$countInExoDpfile1++;
			if ( $keys1 =~ /SNP/ )
			{
				$countInExoDpsnpfile1++;
			}
			if ( $keys1 =~ /DEL/ )
			{
				$countInExoDpdelfile1++;
			}
			if ( $keys1 =~ /INS/ )
			{
				$countInExoDpinsfile1++;
			}
		}
	}
	print OUT
"##Number of Records of File 1 Going to Exonic Dropped of File 1: $countInExoDpfile1\n";
	print OUT
"Number of SNP Records of File 1 Going to Exonic Dropped of File 1: $countInExoDpsnpfile1\n";
	print OUT
"Number of DEL Records of File 1 Going to Exonic Dropped of File 1: $countInExoDpdelfile1\n";
	print OUT
"Number of INS Records of File 1 Going to Exonic Dropped of File 1: $countInExoDpinsfile1\n";

	#Compare File 2 with Exonic Filtered for File 2
	#print OUT "File 2 Records Going to Exonic Dropped of File 2:\n";
	while ( my ( $keys2, $values2 ) = each(%$presentInBoth) )
	{
		if ( exists $$annExoDrophash2{$keys2} )
		{
			my $newval = $$annExoDrophash2{$keys2};

			#print OUT "$keys2:$newval\n";
			$countInExoDpfile2++;
			if ( $keys2 =~ /SNP/ )
			{
				$countInExoDpsnpfile2++;
			}
			if ( $keys2 =~ /DEL/ )
			{
				$countInExoDpdelfile2++;
			}
			if ( $keys2 =~ /INS/ )
			{
				$countInExoDpinsfile2++;
			}
		}
	}
	print OUT
"##Number of Records of File 2 Going to Exonic Dropped of File 2: $countInExoDpfile2\n";
	print OUT
"Number of SNP Records of File 2 Going to Exonic Dropped of File 2: $countInExoDpsnpfile2\n";
	print OUT
"Number of DEL Records of File 2 Going to Exonic Dropped of File 2: $countInExoDpdelfile2\n";
	print OUT
"Number of INS Records of File 2 Going to Exonic Dropped of File 2: $countInExoDpinsfile2\n";

	#Compare File1 with Silent Filtered for File 1
	#print OUT "File 1 Records Going to Silent Filtered of File 1:\n";
	while ( my ( $keys1, $values1 ) = each(%$presentInBoth) )
	{
		if ( exists $$annSthash1{$keys1} )
		{
			my $newval = $$annSthash1{$keys1};

			#print OUT "$keys1:$newval\n";
			$countInStFtfile1++;
			if ( $keys1 =~ /SNP/ )
			{
				$countInStFtsnpfile1++;
			}
			if ( $keys1 =~ /DEL/ )
			{
				$countInStFtdelfile1++;
			}
			if ( $keys1 =~ /INS/ )
			{
				$countInStFtinsfile1++;
			}
		}
	}
	print OUT
"##Number of Records of File 1 Going to Silent Filtered of File 1: $countInStFtfile1\n";
	print OUT
"Number of SNP Records of File 1 Going to Silent Filtered of File 1: $countInStFtsnpfile1\n";
	print OUT
"Number of DEL Records of File 1 Going to Silent Filtered of File 1: $countInStFtdelfile1\n";
	print OUT
"Number of INS Records of File 1 Going to Silent Filtered of File 1: $countInStFtinsfile1\n";

	#Compare File1 with Silent Filtered for File 2
	#print OUT "File 1 Records Going to Silent Filtered of File 2:\n";
	while ( my ( $keys2, $values2 ) = each(%$presentInBoth) )
	{
		if ( exists $$annSthash2{$keys2} )
		{
			my $newval = $$annSthash2{$keys2};

			#print OUT "$keys2:$newval\n";
			$countInStFtfile2++;
			if ( $keys2 =~ /SNP/ )
			{
				$countInStFtsnpfile2++;
			}
			if ( $keys2 =~ /DEL/ )
			{
				$countInStFtdelfile2++;
			}
			if ( $keys2 =~ /INS/ )
			{
				$countInStFtinsfile2++;
			}
		}
	}
	print OUT
"Number of Records of File 2 Going to Silent Filtered of File 2: $countInStFtfile2\n";
	print OUT
"Number of SNP Records of File 2 Going to Silent Filtered of File 2: $countInStFtsnpfile2\n";
	print OUT
"Number of DEL Records of File 2 Going to Silent Filtered of File 2: $countInStFtdelfile2\n";
	print OUT
"Number of INS Records of File 2 Going to Silent Filtered of File 2: $countInStFtinsfile2\n";

	#Compare File 1 with Silent Dropped for File 1
	#print OUT "File 1 Records Going to Silent Dropped of File 1:\n";
	while ( my ( $keys1, $values1 ) = each(%$presentInBoth) )
	{
		if ( exists $$annStDrophash1{$keys1} )
		{
			my $newval = $$annStDrophash1{$keys1};

			#print OUT "$keys1:$newval\n";
			$countInStDpfile1++;
			if ( $keys1 =~ /SNP/ )
			{
				$countInStDpsnpfile1++;
			}
			if ( $keys1 =~ /DEL/ )
			{
				$countInStDpdelfile1++;
			}
			if ( $keys1 =~ /INS/ )
			{
				$countInStDpinsfile1++;
			}
		}
	}
	print OUT
"##Number of Records of File 1 Going to Silent Dropped of File 1: $countInStDpfile1\n";
	print OUT
"Number of SNP Records of File 1 Going to Silent Dropped of File 1: $countInStDpsnpfile1\n";
	print OUT
"Number of DEL Records of File 1 Going to Silent Dropped of File 1: $countInStDpdelfile1\n";
	print OUT
"Number of INS Records of File 1 Going to Silent Dropped of File 1: $countInStDpinsfile1\n";

	#Compare File 2 with Silent Dropped for File 2
	#print OUT "File 2 Records Going to Silent Dropped of File 2:\n";
	while ( my ( $keys2, $values2 ) = each(%$presentInBoth) )
	{
		if ( exists $$annStDrophash2{$keys2} )
		{
			my $newval = $$annStDrophash2{$keys2};

			#print OUT "$keys2:$newval\n";
			$countInStDpfile2++;
			if ( $keys2 =~ /SNP/ )
			{
				$countInStDpsnpfile2++;
			}
			if ( $keys2 =~ /DEL/ )
			{
				$countInStDpdelfile2++;
			}
			if ( $keys2 =~ /INS/ )
			{
				$countInStDpinsfile2++;
			}
		}
	}
	print OUT
"##Number of Records of File 2 Going to Silent Dropped of File 2: $countInStDpfile2\n";
	print OUT
"Number of SNP Records of File 2 Going to Silent Dropped of File 2: $countInStDpsnpfile2\n";
	print OUT
"Number of DEL Records of File 2 Going to Silent Dropped of File 2: $countInStDpdelfile2\n";
	print OUT
"Number of INS Records of File 2 Going to Silent Dropped of File 2: $countInStDpinsfile2\n";
}
