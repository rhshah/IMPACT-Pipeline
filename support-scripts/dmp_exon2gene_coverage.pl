#!/usr/bin/perl
use strict;
use warnings;

##########
## Input = parsed DepthOfCoverage (column 1 = interval; column 2 = total coverage)
## Input = v1_hg19_gene_coords.txt
##
## Output = mean coverage per gene
##########

if ($#ARGV!=1) {die "usage is 'perl exon2genecovg.pl [DepthOfCoverage.exons.out] v1_hg19_gene_coords.txt'\n";}
my $filename = shift;
my $loci_file = shift;

my @target_covg;
open (FH, "<$filename") or die "can't open $filename\n";
my $topline = <FH>;
while (my $text = <FH>) {
    chomp $text;
    my @line = split "\t", $text;
    push @target_covg, $line[1];
}
close (FH);

my %gene_covg;
my %gene_length;
my %gene_start;
my $target_counter = 0;
open (FH2, "<$loci_file") or die "can't open loci\n";
while (my $text = <FH2>) {
    chomp $text;
    my @line = split "\t", $text;
    my $length = $line[4]-$line[3]+1;
    my $start_coord = "chr$line[2]:$line[3]";
    my $bases = $target_covg[$target_counter];
    my $gene = $line[1];
    push @{$gene_covg{$gene}}, $bases;
    push @{$gene_length{$gene}}, $length;
    push @{$gene_start{$gene}}, $start_coord;
    $target_counter++;
}
close (FH2);

my $key;
foreach $key (sort keys %gene_covg) {
    my $totalbases = 0;
    my $totallength = 0;
    for (my $i=0; $i<=$#{$gene_covg{$key}}; $i++) {
	if (defined $gene_covg{$key}[$i] ) {
            $totalbases += $gene_covg{$key}[$i];
        }
        else {
            $totalbases += 0;
        }
        if (defined $gene_length{$key}[$i] ){
            $totallength += $gene_length{$key}[$i];
        }
        else {
            $totallength += 0;
        }
    }
    my $totalcovg = $totalbases/$totallength;
    my $coord = $gene_start{$key}[0];
    print "$key\t$totalcovg\t$coord\n";
}
