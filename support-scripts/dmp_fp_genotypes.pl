# Modified by: Donavan Cheng
# 2013 Apr 22nd

#!/usr/bin/perl
use strict;
use warnings;
use MSKCC_DMP_Logger;

if($#ARGV!=1) {die "usage is 'perl FP_genotypes_v1.1.pl XXX.FP.counts.txt XXX_FP_tiling_genotypes.txt'\n";}

my $logger = MSKCC_DMP_Logger->get_logger('FP_GENOTYPES');
$logger->start_local();
$logger->info("Fingerprint genotyping is being analyzed\n");

my $file = $ARGV[0];
my $genotypes_file = $ARGV[1];

#my @snp;
my %snpHash=();
my %coord;
open (FH0, "<${genotypes_file}") or die $logger->fatal("Can't open $genotypes_file. Error: $!\n");
while (my $text = <FH0>) {
    chomp $text;
    my @line = split "\t", $text;
    $coord{"$line[0]"} = 0;    
    @{$snpHash{$line[0]}} = split (/\//, $line[1]);
}
close (FH0);

open (FH, "<$file") or die $logger->info("Can not open $file. Error: $!\n");
my $topline = <FH>;
while (my $text = <FH>) {
    my $gen1=0; my $gen2=0;
    my $gen1freq=0; my $gen2freq=0;
    chomp $text;
    my @line = split "\t", $text;
    if(! exists $coord{$line[0]}){ next;}

    print "$line[0]\t";
    my @counts = split " ", $line[4];
    if (@{$snpHash{$line[0]}}[0] eq "A") {print "$counts[0] "; my @array = split ":", $counts[0]; $gen1="A"; $gen1freq=$array[1];}
    elsif (@{$snpHash{$line[0]}}[0] eq "C") {print "$counts[1] "; my @array = split ":", $counts[1]; $gen1="C"; $gen1freq=$array[1];}
    elsif (@{$snpHash{$line[0]}}[0] eq "G") {print "$counts[2] "; my @array = split ":", $counts[2]; $gen1="G"; $gen1freq=$array[1];}
    elsif (@{$snpHash{$line[0]}}[0] eq "T") {print "$counts[3] "; my @array = split ":", $counts[3]; $gen1="T"; $gen1freq=$array[1];}

    if (@{$snpHash{$line[0]}}[1] eq "A") {print "$counts[0]\t"; my @array = split ":", $counts[0]; $gen2="A"; $gen2freq=$array[1];}
    elsif (@{$snpHash{$line[0]}}[1] eq "C") {print "$counts[1]\t"; my @array = split ":", $counts[1]; $gen2="C"; $gen2freq=$array[1];}
    elsif (@{$snpHash{$line[0]}}[1] eq "G") {print "$counts[2]\t"; my @array = split ":", $counts[2]; $gen2="G"; $gen2freq=$array[1];}
    elsif (@{$snpHash{$line[0]}}[1] eq "T") {print "$counts[3]\t"; my @array = split ":", $counts[3]; $gen2="T"; $gen2freq=$array[1];}
	
    if ($gen1freq == 0 && $gen2freq == 0) {
	print "--\t0\n";
    }
	
    elsif ($gen1freq < $gen2freq) {
	my $fraction = $gen1freq/($gen1freq+$gen2freq);
	if ($fraction < 0.10) {print "$gen2$gen2\t$fraction\n";}
	else {print "$gen1$gen2\t$fraction\n";}
    }
    else {
	my $fraction = $gen2freq/($gen1freq+$gen2freq);
	if ($fraction < 0.10) {print "$gen1$gen1\t$fraction\n";}
	else {print "$gen1$gen2\t$fraction\n";}
    }

}

close (FH);

$logger->info("Finger print analysis is completed.\n");
