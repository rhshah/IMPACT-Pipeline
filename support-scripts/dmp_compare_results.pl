#!/usr/bin/perl 
# compareResults_try1.pl --- 
# Author: Shah <shahr2@phoenix-h1>
# Created: 18 Nov 2013
# Version: 0.01


use strict;
use Getopt::Long;
use Cwd;

my ($file1,$file2,$outfile,$outdir);

if (@ARGV < 2 or !GetOptions (
    'file1|f1:s'   	=> \$file1,
    'file2|f2:s'   	=> \$file2,
    'outfile|of:s'   	=> \$outfile,
    'outdir|o:s'        => \$outdir))
    {
        Usage();
    }

if ((!defined $file1) or (!defined $file2)){
    print "Input file have not been defined. See Usage\n";
    Usage();
    exit;
    
}
if (!defined $outfile){
    print "Output file is not specified, default file will be used as output file.\n";
    $outfile = "VariantFileCompareReport.txt";
}
if (!defined $outdir){
    print "Output directory is not specified, current working directory will be used as output directory.\n";
    $outdir = getcwd;
}
open OUT, ">" , "$outdir/$outfile" || die "Cannot Open $outfile,$!";

my ($hash1,$realCols1,$rows1,$cols1) = ReadFile($file1);
my ($hash2,$realCols2,$rows2,$cols2) = ReadFile($file2);
my @values1 = @$realCols1;
my @values2 = @$realCols2;
my $flag = 1;
my $null = 2;

print OUT "###Dimensional Statisticcs###\n";
my $moreRows ;
if ($rows1 == $rows2) {
    print OUT "NumberOfRows=PASS(file1:$rows1,file2:$rows2)\n";
    $moreRows = $rows1;
    
}
else {
    print OUT "NumberOfRows=FAIL(file1:$rows1,file2:$rows2)\n";
    $flag = 2;
    if ($rows1 > $rows2) {
        $moreRows = $rows1;
    }
    else {
        $moreRows = $rows2;
    }
    
}
if ($cols1 == $cols2) {
    
    print OUT "NumberOfCols=PASS(file1:$cols1,file2:$cols2)\n";
    
}
else {
    print OUT "NumberOfCols=FAIL(file1:$cols1,file2:$cols2)\n";
    $flag = 2;
    
}

for(my $cnt=0;$cnt<$moreRows;$cnt++) {
    
        my @numbers1 = split(":",$values1[$cnt]);
        my @numbers2 = split(":",$values2[$cnt]);
        my $subcols1 = scalar(@numbers1);
        my $subcols2 = scalar(@numbers2);
        my $moreCols;
        my $line = $cnt + 1;
        if ($subcols1 == $subcols2) {
            $moreCols = $subcols1;
            
        }
        else {
           
            if ($subcols1 > $subcols2) {
                $moreCols = $subcols1;
            }
            else {
                $moreCols = $subcols2;
            }
            
        }
        for(my $i=0;$i<$moreCols;$i++) {
            if ($numbers1[$i] != $numbers2[$i]) {
                if ($numbers1[$i] == $null) {
                    print OUT "ColumnsMissingValues=FAIL(LineNumber=$line;Column=$i;file=file1;NosColumnFile1:$subcols1,NosColumnsFile2:$subcols2)\n";
                    $flag = 2;
                    
                }
                if ($numbers2[$i] == $null) {
                    print OUT "ColumnsMissingValues=FAIL(LineNumber=$line;Column=$i;file=file2;ColumnFile1:$subcols1,ColumnsFile2:$subcols2)\n";
                    $flag = 2
                }
            }
        }
    }

if ($flag == 2) {
    print OUT "NumberofColsHavingValues=FAIL\n";
    
}
else {
    print OUT "NumberofColsHavingValues=PASS\n";
    
}     

if ($flag == 1) {
    print OUT "DimensionsMatchStatus=PASS.\n";
    
}
else {
    print OUT "DimensionsMatchStatus=FAIL.\n";
        
    
}
AnalyzeVariants($hash1,$rows1,$hash2,$rows2);
close(OUT);

exit;

#####################################
#####################################
#How to use the script.
sub Usage{
    print "Unknow option: @_\n" if (@_);
    
    print "\nUsage : compareResults.pl [options]
        [--file1|f1        S File 1 havaing variants (required & will be noted as file1)]
        [--file2|f2        S File 2 having variants (required & will be denoted as file2)]
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
sub ReadFile{
    my ($file) = @_;
    my (%variantsHash);
    my $rows = 0;
    my $cols = 0;
    
    my(@line) = ();
    my (@processedLine) = ();
    
    my ($sample,$normal,$chr,$start,$ref,$alt,$varaintClass,$gene,$exon,$cdna,$aachange,$ndp,$nrc,$nac,$naf,$tdp,$trc,$tac,$taf,$trefp,$trefn,$taltp,$taltn,$allN,$tnfreq,$occ,$pct);
    
    open IN, "<", $file || die $!;
    
    while (<IN>) {
        $rows++;
        chomp;
        if($. == 1){
            my @HEADER = split("\t",$_);
            $cols = scalar(@HEADER);
            next;              
        }
        my $null = 2;
        my $val = 1;
        my(@inferLine) = ();
        @line = split("\t");
        foreach my $data (@line) {
            $data =~ s/\s//g;
            if (! $data) {
                push(@inferLine,$null);
               
            }
            else {
                  push(@inferLine,$val);  
            }
        }
        my $data = join(":",@inferLine);
        push (@processedLine,$data);
        $sample = $line[0];
        $normal = $line[1];
        $chr = $line[2];
        $start = $line[3];
        $ref = $line[4];
        $alt = $line[5];
        $varaintClass = $line[6];
        $gene = $line[7];
        $exon = $line[8];
        $cdna = $line[12];
        $aachange = $line[13];
        $ndp = $line[23];
        $nrc = $line[24];
        $nac = $line[25];
        $naf = $line[26];
        $tdp = $line[27];
        $trc = $line[28];
        $tac = $line[29];
        $taf = $line[30];
        $trefp = $line[31];
        $trefn = $line[32];
        $taltp = $line[33];
        $taltn = $line[34];
        $allN = $line[35];
        $tnfreq = $line[36];
        ($occ,$pct) = split(";",$line[38]);
        
        
        $variantsHash{"$sample:$normal:$chr:$start:$ref:$alt"} = "$varaintClass:$gene:$exon:$cdna:$aachange:$ndp:$nrc:$nac:$naf:$tdp:$trc:$tac:$taf:$trefp:$trefn:$taltp:$taltn:$allN:$tnfreq:$occ:$pct";
        
    }
    close(IN);
    
    return(\%variantsHash,\@processedLine,$rows,$cols)
}
#####################################
#####################################
#Analyse the files
sub AnalyzeVariants {
    my($hash1,$rows1,$hash2,$rows2) = @_;
    my(%vh1) = %$hash1;
    my(%vh2) = %$hash2;
    my $flag = 1;
    print OUT "###Variant Matching Statistics###\n";
    my %vh3 = ();
    
    if ($rows1 >= $rows2) {
        
        while (my($keys1,$values1) = each (%vh1)) {
            my($variantClass1,$gene1,$exon1,$cdna1,$aachange1,$ndp1,$nrc1,$nac1,$naf1,$tdp1,$trc1,$tac1,$taf1,$trefp1,$trefn1,$taltp1,$taltn1,$allN1,$tnfreq1,$occ1,$pct1) = split(":",$values1);
            
            if (exists $vh2{$keys1}) {
                my($values2) = $vh2{$keys1};
                my($variantClass2,$gene2,$exon2,$cdna2,$aachange2,$ndp2,$nrc2,$nac2,$naf2,$tdp2,$trc2,$tac2,$taf2,$trefp2,$trefn2,$taltp2,$taltn2,$allN2,$tnfreq2,$occ2,$pct2) = split(":",$values2);
                if (!(($ndp1 == $ndp2) and ($nrc1 == $nrc2) and ($nac1 == $nac2) and ($naf1 == $naf2))){
                    print OUT "NormalVariantValues=FAIL($keys1)\n";
                    my $flag = 2;
                }
                
                if (!(($tdp1 == $tdp2) and ($trc1 == $trc2) and ($tac1 == $tac2) and ($taf1 == $taf2))) {
                    print OUT "TumorVariantValues=FAIL($keys1)\n";
                    my $flag = 2;
                }
                
                
                if (!(($trefp1 == $trefp2) and ($trefn1 == $trefn2) and ($taltp1 == $taltp2) and ($taltn1 == $taltn2))) {
                    print OUT "TumorVariant_Forward_Reverese_Values=FAIL($keys1)\n";
                    my $flag = 2;
                }
                    
                if (!(($allN1 eq $allN2) and ($tnfreq1 == $tnfreq2) and ($occ1 == $occ2) and ($pct1 == $pct2))) {
                        print OUT "AllNormalVariantValues=FAIL($keys1)\n";               
                        my $flag = 2;
                    }
                    
                if (!(($variantClass1 eq $variantClass2) and ($gene1 eq $gene2) and ($exon1 eq $exon2) and ($cdna1 eq $cdna2) and ($aachange1 eq $aachange2))) {
                    print OUT "VariantAnnotationValues=FAIL($keys1)\n";               
                    my $flag = 2;
                    
                }
            }
            
            else {
                print OUT "PresentOnlyInFile1=FAIL($keys1)\n";
                my $flag = 2;
                
            }
        }       
    }
    if ($rows2 > $rows1) {
        
        while (my($keys1,$values1) = each (%vh2)) {
            my($variantClass1,$gene1,$exon1,$cdna1,$aachange1,$ndp1,$nrc1,$nac1,$naf1,$tdp1,$trc1,$tac1,$taf1,$trefp1,$trefn1,$taltp1,$taltn1,$allN1,$tnfreq1,$occ1,$pct1) = split(":",$values1);
            
            if (exists $vh1{$keys1}) {
                my($values2) = $vh1{$keys1};
                my($variantClass2,$gene2,$exon2,$cdna2,$aachange2,$ndp2,$nrc2,$nac2,$naf2,$tdp2,$trc2,$tac2,$taf2,$trefp2,$trefn2,$taltp2,$taltn2,$allN2,$tnfreq2,$occ2,$pct2) = split(":",$values2);
                if (!(($ndp1 == $ndp2) and ($nrc1 == $nrc2) and ($nac1 == $nac2) and ($naf1 == $naf2))){
                    print OUT "NormalVariantValues=FAIL($keys1)\n";
                    my $flag = 2;
                }
                
                if (!(($tdp1 == $tdp2) and ($trc1 == $trc2) and ($tac1 == $tac2) and ($taf1 == $taf2))) {
                    print OUT "TumorVariantValues=FAIL($keys1)\n";
                    my $flag = 2;
                }
                
                
                if (!(($trefp1 == $trefp2) and ($trefn1 == $trefn2) and ($taltp1 == $taltp2) and ($taltn1 == $taltn2))) {
                    print OUT "TumorVariant_Forward_Reverese_Values=FAIL($keys1)\n";
                    my $flag = 2;
                }
                    
                if (!(($allN1 eq $allN2) and ($tnfreq1 == $tnfreq2) and ($occ1 == $occ2) and ($pct1 == $pct2))) {
                        print OUT "AllNormalVariantValues=FAIL($keys1)\n";               
                        my $flag = 2;
                    }
                    
                if (!(($variantClass1 eq $variantClass2) and ($gene1 eq $gene2) and ($exon1 eq $exon2) and ($cdna1 eq $cdna2) and ($aachange1 eq $aachange2))) {
                    print OUT "VariantAnnotationValues=FAIL($keys1)\n";               
                    my $flag = 2;
                    
                }
            }
            
            else {
                print OUT "PresentOnlyInFile2=FAIL($keys1)\n";
                my $flag = 2;
                
            }
        }       
    }
        
    if ($flag == 2) {
        print OUT "VariantsMatchStatus=FAIL\n";
        
    }
    else {
        print OUT "VariantsMatchStatus=PASS\n";
        
    }    
    
    
    return;
    
}

__END__

=head1 NAME

compareResults_try1.pl - Describe the usage of script briefly

=head1 SYNOPSIS

compareResults_try1.pl [options] args

      -opt --long      Option description

=head1 DESCRIPTION

Stub documentation for compareResults_try1.pl, 

=head1 AUTHOR

Shah, E<lt>shahr2@phoenix-h1E<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 by Shah

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
