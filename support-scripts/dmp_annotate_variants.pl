#!/usr/bin/perl -w
##########AnnotateAssessFilterVariants.pl########
#Author: Ronak Shah
#Date: 31/01/2013
#LastModified: 05/01/2013
#Version:1.2
#Description: Annotating and Assesing Varinats
#             using Oncotator & Mutation Assessor
##04/12/2013
#Updatd to work with v1.5 of the pipeline
##04/30/2013
#Added functionality to filter variants.
#
##05/01/2013
#v1.2
#Made it comaptible with the data sets.
#
##05/14/2013
#v1.3
#Made exception for TERT.
#
##07/23/2013
# Ahmet Zehir
#v1.4
# Changed annotation from Oncotator to Annovar
# Added coverage file generation (both gene and exon based)
# Added two-tiered filtering scheme based on hotspots
# Added logger functionality
#
##08/22/2013
# Added functionality to generate clinical reports
#############################
use strict;
use Getopt::Long;
use IO::File;
use Cwd;
use Tie::IxHash;
use MSKCC_DMP_Logger;

#--This variable holds the current time
my $now = time;
my ($mutationFile, $outdir, $deleteFiles, $titleFile, $configFile, $exonCoverageFile, $geneCoverageFile);

if (@ARGV < 1 or !GetOptions (
	    'somaticMutIndelFile|si:s'       => \$mutationFile,
	    'configurationFile|c:s'          => \$configFile,
	    'outdir|o:s'                     => \$outdir,
            'exonCoverageFile|ec:s'          => \$exonCoverageFile,
	    'geneCoverageFile|gc:s'          => \$geneCoverageFile,
            'deleteUnwantedFiles|d:i'        => \$deleteFiles,
	    'TitleFile|t:s'                  => \$titleFile))
{
    Usage();
}

my $logger = MSKCC_DMP_Logger->get_logger('ANNOTATE_ASSESS_FILTER_VARIANTS');

#print "MutationFile: $mutationFile\nTitle File: $titleFile\nConfigFile: $configFile\nOutdir: $outdir\n";
my ($locationRef, $versionRef, $parameterRef) = &GetConfiguration($configFile);
my %location = %$locationRef;
my %version = %$versionRef;
my %parameter = %$parameterRef; 

#import file locations
my $annovar = $location{"Annovar"};
my $annovarDB = $location{"Annovar_db"};
my $canonicalTxFile = $location{"Canonical_refFlat_file"};
my $java = $location{"JAVA_1_6"};
my $IGVtools = $location{"IGVtools"};
my $translatedFolder = $location{"TranslationFolder"};
my $hotspots = $location{"HotSpot_mutations"};
my $exonIntervals = $location{"Exon_Interval_table_cv1"};
my $tumor_suppressor = $location{"Tumor_supressor_list"};



#import parameters
my $tumorFreq = $parameter{"Tfreq_AnnotationFilter"};
my $VF_threshold_hotspot = $parameter{"VF_threshold_hotspot"};
my $tumornormalFreqRatio = $parameter{"TNfreqRatio_AnnotationFilter"};
my $MAFthreshold = $parameter{"MAFthreshold_AnnotationFilter"};
my $AD_Tthreshold_indel_high = $parameter{"AD_Tthreshold_indel_high_AnnotationFilter"};
my $AD_Tthreshold_indel_low = $parameter{"AD_Tthreshold_indel_low_AnnotationFilter"};
my $VF_Tthreshold_indel_high = $parameter{"VF_Tthreshold_indel_high_AnnotationFilter"};
my $VF_Tthreshold_indel_low = $parameter{"VF_Tthreshold_indel_low_AnnotationFilter"};
my $VF_Tthreshold_indel_hotspot_high = $parameter{"VF_Tthreshold_indel_hotspot_high_AnnotationFilter"};
my $VF_Tthreshold_indel_hotspot_low = $parameter{"VF_Tthreshold_indel_hotspot_low_AnnotationFilter"};

my $coverage_threshold = $parameter{"Coverage_threshold_darwin_report"};

if((! $mutationFile) or (! $titleFile) or (! $configFile) or (! $translatedFolder) or (! $hotspots))
   {
       $logger->fatal("Mutation File is missing") if(!$mutationFile);
       $logger->fatal("Title File is missing") if(!$titleFile);
       $logger->fatal("Configuration File is missing") if(!$configFile);
       $logger->fatal("Translated Folder is missing") if(!$translatedFolder);
       $logger->fatal("Hotspots is missing") if(!$hotspots);
       Usage();
       exit;
   }
if(! $deleteFiles)
{
    $logger->warn("The option to delete intermediary files was not set, the files will be deleted by default");
    $deleteFiles = 2;
}

if (! $outdir)
{
    $outdir = getcwd;

}


$logger->info("Variant annotation and filtering has started for $mutationFile");


#Get Title file information
my($patientIDPerSampleId,$classPerPatientIdSampleId) = &ReadTitleFile($titleFile,$outdir);

#Run annovar annotation
&RunAnnovar($mutationFile, $outdir);

# Merge annovar annotations with mutation file
my($mergedFile, $filesToDelete) = &MergeAnnotations($mutationFile, $outdir, $canonicalTxFile, $translatedFolder);

# Filter the merged annotation file 
my($filteredAnnoFile, $bedFile) = &FilterAnnotations($mergedFile, $tumor_suppressor, $outdir, $VF_threshold_hotspot, $tumorFreq, $tumornormalFreqRatio, $MAFthreshold, $AD_Tthreshold_indel_high, $AD_Tthreshold_indel_low, $VF_Tthreshold_indel_high, $VF_Tthreshold_indel_low, $VF_Tthreshold_indel_hotspot_high,$VF_Tthreshold_indel_hotspot_low, $hotspots, $patientIDPerSampleId,$classPerPatientIdSampleId);

# Generate per patient VCF Files
&GenerateVCFfiles($filteredAnnoFile, $outdir, $tumorFreq, $tumornormalFreqRatio, $patientIDPerSampleId,$classPerPatientIdSampleId);

# Generate per patient per gene coverage files
&GenerateGeneCoverage($titleFile, $outdir, $geneCoverageFile, $exonIntervals);

# Generate per patient per exon coverage files
&GenerateExonCoverage($filteredAnnoFile, $exonCoverageFile, $translatedFolder, $exonIntervals, $outdir, $titleFile, $coverage_threshold);

# Generate clinical reports
&GenerateClinicalReports($outdir, $hotspots);

#$logger->info("The Annotated Results of Variants are written in: $outdir/$filteredAnnoFile");
my @deleteFiles = @$filesToDelete;


#deletion process
if($deleteFiles == 1)
{

    &HouseKeeping($outdir,@deleteFiles);
}
if($deleteFiles == 2)
{
    $logger->info("No house keeping performed");
}


#--Calculate total runtime (current time minus start time) 
$now = time - $now;

#--Print runtime 
$logger->info("Total running time: %02d:%02d:%02d\n\n", int($now / 3600), int(($now % 3600) / 60), int($now % 60));
$logger->info("Done, Thanks for choosing us to annotate your variants...\n");
exit;
#####################################
#####################################
#How to use the script.
sub Usage
{
    print "Unknow option: @_\n" if (@_);

	print "\nUsage : AnnotateAssessFilterVarinats.pl [options]
        [--SomaticMutIndelFile|si          S File containing mutations (required and submit with full path,Ex:/SomePath/CytoOv_SomaticMutIndel.txt)]
        [--ConfigurationFile|c             S Configuration file that contains the locations for the programs and the databases (required and submit with full path)]
        [--titleFile|t                     S tab-delimited title file for the samples (required and submit with full path)]
        [--outdir|o                        S Path where all the output files will be written (optional) [default:cwd]]
        [--exonCoverageFile|ec             S Path where the all exon coverage file is located (full path]
        [--geneCoverageFile|gc             S Path where the gene coverage file is located (full path)]
        [--deleteUnwantedFiles|d           I 2=>To delete files 1=> To keep files (default:2,optional)]
        \n";

	print "For any question please contact Ronak Shah (shahr2\@mskcc.org)\n";
	print "!~!~! Comments, flames and suggestion are welcome !~!~!\n";
	exit;
}

##############################
##############################
# Run annovar on the mutation file

sub RunAnnovar{
    my ($sample, $outDir) = @_;
    my ($annInput) = $sample =~ /.*\/(.*).txt/;
    $logger->info("Annovar annotation has started....");
    if((-e "$outDir/$annInput") and ((-s "$outDir/$annInput") != 0)) {
	$logger->info( "Annovar input file: $annInput already exists and this process will not create a new one");
    }else{
	$logger->info("Creating the annovar input file");
	open(IN, "<", "$sample") or die $logger->fatal("Can not open mutation file: $!");
	open ANN, ">", "$outDir/$annInput" ||  die $logger->fatal("Can not creat $annInput: $!");
	while (<IN>) {
	    chomp;
	    #M-1608-42-T	1	9784916	C	T	KEEP
	    next if($_=~/^Sample\t/);
	    my @line = split("\t");
	    my $sample = shift @line;
	    my $normal = shift @line;
	    my $chr = shift @line;
	    my $pos = shift @line;
	    my $ref = shift @line;
	    my $alt = shift @line;
	    my $at = join("|", @line);
	    if(length($ref) > length ($alt)) {# deletion 2002 CGG C need become 2003 2004 GG -
		my $newref = substr($ref, 1);
		my $posEnd = $pos+length($newref);
		$pos++;
		print ANN "$chr\t$pos\t$posEnd\t$newref\t-\t$sample|$normal|$at\n";
	    }elsif(length($ref) < length($alt)){ #insertion
		my $newalt = substr($alt, 1);
		my $posEnd = $pos;
		print ANN "$chr\t$pos\t$posEnd\t-\t$newalt\t$sample|$normal|$at\n";
	    }else {
		print ANN "$chr\t$pos\t$pos\t$ref\t$alt\t$sample|$normal|$at\n";
	    }
	}
	close ANN;
	close IN;
    }
    eval{`$annovar --geneanno -buildver hg19 $outDir/$annInput $annovarDB`};
    if($@){
	$logger->fatal("Command: $annovar --geneanno -buildver hg19 $outDir/$annInput $annovarDB can not be run! Error: $@");
	exit(1);
    }else{
	if((-e "$outDir/$annInput.exonic_variant_function") && ((-s "$outDir/$annInput.exonic_variant_function") != 0)){
	    $logger->warn("Exonic variant annotation file exists. Exonic annotation will not be performed.");
	}else{
	    $logger->info("Command: $annovar --geneanno -buildver hg19 $outDir/$annInput $annovarDB is being run");
	    `$annovar --geneanno -buildver hg19 $outDir/$annInput $annovarDB`;
	}
    }

    eval{`$annovar -filter -dbtype snp137NonFlagged -buildver hg19 $outDir/$annInput $annovarDB`};
    if($@){
	$logger->fatal("Command: $annovar -filter -dbtype snp137NonFlagged -buildver hg19 $outDir/$annInput $annovarDB can not be run! Error: $@");
	exit(1);
    }else{
	if((-e "$outDir/$annInput.hg19_snp137NonFlagged_filtered") && ((-s "$outDir/$annInput.hg19_snp137NonFlagged_filtered") != 0)){
	    $logger->warn("dnSNP variant annotation file exists. dbSNP  annotation will not be performed.");
	}else{
	    $logger->info("Command: $annovar -filter -dbtype snp137NonFlagged -buildver hg19 $outDir/$annInput $annovarDB is being run!");
	    `$annovar -filter -dbtype snp137NonFlagged -buildver hg19 $outDir/$annInput $annovarDB`;
	}
    }

    eval{`$annovar -filter -dbtype cosmic64 -buildver hg19 $outDir/$annInput $annovarDB`};
    if($@){
	$logger->fatal("Command: $annovar -filter -dbtype cosmic64 -buildver hg19 $outDir/$annInput $annovarDB can not be run! Error: $@");
	exit(1);
    }else{
	if((-e "$outDir/$annInput.hg19_cosmic64_filtered") && ((-s "$outDir/$annInput.hg19_cosmic64_filtered") != 0)){
	    $logger->warn("Cosmic variant annotation file exists. Cosmic annotation will not be performed.");
	}else{
	    $logger->info("Command: $annovar -filter -dbtype cosmic64 -buildver hg19 $outDir/$annInput $annovarDB is being run!");
	    `$annovar -filter -dbtype cosmic64 -buildver hg19 $outDir/$annInput $annovarDB`;
	}
    }

    eval{`$annovar -filter -dbtype 1000g2012apr_all -buildver hg19 $outDir/$annInput $annovarDB`};
    if($@){
	$logger->fatal("Command: $annovar -filter -dbtype 1000g2012apr_all -buildver hg19 $outDir/$annInput $annovarDB can not be run! Error: $@");
	exit(1);
    }else{
	if((-e "$outDir/$annInput.hg19_ALL.sites.2012_04_filtered") && ((-s "$outDir/$annInput.hg19_ALL.sites.2012_04_filtered") != 0)){
	    $logger->warn("1000G MAF variant annotation file exists. 1000G MAF annotation will not be performed.");
	}else{
	    $logger->info("Command: $annovar -filter -dbtype 1000g2012apr_all -buildver hg19 $outDir/$annInput $annovarDB is being run!");
	    `$annovar -filter -dbtype 1000g2012apr_all -buildver hg19 $outDir/$annInput $annovarDB`;
	}
    }

    eval{`$annovar -filter -dbtype ljb2_gerp++ -buildver hg19 $outDir/$annInput $annovarDB`};
    if($@){
	$logger->fatal("Command: $annovar -filter -dbtype ljb2_gerp++ -buildver hg19 $outDir/$annInput $annovarDB can not be run Error: $@");
	exit(1);
    }else{
	if((-e "$outDir/$annInput.hg19_ljb2_gerp++_filtered") && ((-s "$outDir/$annInput.hg19_ljb2_gerp++_filtered") != 0)){
	    $logger->warn("GERP++ database variant annotation file exists. GERP++ database annotation will not be performed.");
	}else{
	    $logger->info("Command: $annovar -filter -dbtype ljb2_gerp++ -buildver hg19 $outDir/$annInput $annovarDB is being run!");
	    `$annovar -filter -dbtype ljb2_gerp++ -buildver hg19 $outDir/$annInput $annovarDB`;
	}
    }

    eval{`$annovar -filter -dbtype ljb2_lrt -buildver hg19 $outDir/$annInput $annovarDB`};
    if($@){
	$logger->fatal("Command: $annovar -filter -dbtype ljb2_lrt -buildver hg19 $outDir/$annInput $annovarDB can not be run! Error: $@");
	exit(1);
    }else{
	if((-e "$outDir/$annInput.hg19_ljb2_lrt_filtered") && ((-s "$outDir/$annInput.hg19_ljb2_lrt_filtered") != 0)){
	    $logger->warn("LRT database variant annotation file exists. LRT database annotation will not be performed.");
	}else{
	    $logger->info("Command: $annovar -filter -dbtype ljb2_lrt -buildver hg19 $outDir/$annInput $annovarDB is being run!");
	    `$annovar -filter -dbtype ljb2_lrt -buildver hg19 $outDir/$annInput $annovarDB`;
	}
    }

    eval{`$annovar -filter -dbtype ljb2_phylop -buildver hg19 $outDir/$annInput $annovarDB`};
    if($@){
	$logger->fatal("Command: $annovar -filter -dbtype ljb2_phylop -buildver hg19 $outDir/$annInput $annovarDB can not be run! Error: $@");
	exit(1);
    }else{
	if((-e "$outDir/$annInput.hg19_ljb2_phylop_filtered") && ((-s "$outDir/$annInput.hg19_ljb2_phylop_filtered") != 0)){
	    $logger->warn("PhyloP database variant annotation file exists. PhyloP database annotation will not be performed.");
	}else{
	    $logger->info("Command: $annovar -filter -dbtype ljb2_phylop -buildver hg19 $outDir/$annInput $annovarDB is being run!");
	    `$annovar -filter -dbtype ljb2_phylop -buildver hg19 $outDir/$annInput $annovarDB`;
	}
    }

    eval{`$annovar -filter -dbtype ljb2_ma -buildver hg19 $outDir/$annInput $annovarDB`};
    if($@){
	$logger->fatal("Command: $annovar -filter -dbtype ljb2_ma -buildver hg19 $outDir/$annInput $annovarDB can not be run! Error: $@");
	exit(1);
    }else{
	if((-e "$outDir/$annInput.hg19_ljb2_ma_filtered") && ((-s "$outDir/$annInput.hg19_ljb2_ma_filtered") != 0)){
	    $logger->warn("Mutation Assesor database variant annotation file exists. Mutation Assesor database annotation will not be performed.");
	}else{
	    $logger->info("Command: $annovar -filter -dbtype ljb2_ma -buildver hg19 $outDir/$annInput $annovarDB is being run!");
	    `$annovar -filter -dbtype ljb2_ma -buildver hg19 $outDir/$annInput $annovarDB`;
	}
    }

    eval{`$annovar -filter -dbtype avsift -buildver hg19 $outDir/$annInput $annovarDB`};
    if($@){
	$logger->fatal("Command: $annovar -filter -dbtype avsift -buildver hg19 $outDir/$annInput $annovarDB can not be run! Error: $@");
	exit(1);
    }else{
	if((-e "$outDir/$annInput.hg19_avsift_filtered") && ((-s "$outDir/$annInput.hg19_avsift_filtered") != 0)){
	    $logger->warn("AVSIFT database variant annotation file exists. AVSIFT database annotation will not be performed.");
	}else{
	    $logger->info("Comand: $annovar -filter -dbtype avsift -buildver hg19 $outDir/$annInput $annovarDB is being run!");
	    `$annovar -filter -dbtype avsift -buildver hg19 $outDir/$annInput $annovarDB`;
	}
    }

    return($annInput);
}

##############################
##############################
# Merge Annotations

sub MergeAnnotations {
    my ($sample, $annFolder, $canonical, $translatedFolder) = @_;
    my (%canonicalTx, %strandHash, %aa_annotHash, @deleteFiles);
    print "$sample\n";
    my ($annInput) = $sample =~ /(.*)\.txt/;
    push @deleteFiles, $annInput;
    $logger->info("Annovar annotations are being merged to variant file......");
    open(IN, "<", $canonical) or die $logger->fatal("Can not open $canonical: $!");
    while (<IN>) {
	chomp;
	next if($_ =~ m/\#/);
	my @line = split("\t", $_);
	my $gene = $line[0];
	my $transcriptID = $line[1];
	my $strand = $line[3];
	$canonicalTx{$gene} = $transcriptID;
	$strandHash{$gene}= $strand;
    }
    close IN;
    
    my (%annExonicHash, %annVariantHash, %annDBSNPhash, %annCOSMIChash, %ann1000Ghash, %annMAhash, %annGERPhash, %annLRThash, %annPhyloPhash, %annLJBSIFThash, %annAVSIFThash, @header);

    tie (my %inputHash, 'Tie::IxHash');
    open IN, "<", $sample || die $logger->fatal("Can not open $sample: $!");
    while (<IN>) {
	chomp;
	#M-1608-42-T	1	9784916	C	T	KEEP
	if($_=~/^Sample\t/){
	    @header = split("\t");
	    next;
	}
	my @line = split("\t");
	my $sample = shift @line;
	my $normal = shift @line;
	my $chr = shift @line;
	my $pos = shift @line;
	my $ref = shift @line;
	my $alt = shift @line;
	my $at = join("|", @line);
      	my $key_orig = $sample.":".$normal.":".$chr.":".$pos.":".$ref.":".$alt;
	#print "$key_orig\n";
	$inputHash{$key_orig} = $at;
    }
    close IN;
    
    # Read in exonic variant annotation
    ($sample) = $sample =~ /.*\/(.*)/;
    my $annExonicOutput = $sample;
    $annExonicOutput =~ s/.txt/.exonic_variant_function/;
    push @deleteFiles, $annExonicOutput;
    open IN, "<", "$annFolder/$annExonicOutput" || die $logger->fatal("Can not open $annExonicOutput: $!");
    while(<IN>){
	chomp;
	my @line = split("\t");
	my $result = $line[1];
	my $chr = $line[3];
	my $startPos = $line[4];
	my $endPos = $line[5];
	my $ref = $line[6];
	my $alt = $line[7];
	$result =~ s/\ /_/;
	my @annots = split(",", $line[2]);
	my ($gene, $transcriptID, $exon, $cDNAchange, $AAchange);
	foreach my $annots (sort @annots){
	    my @e = split(":", $annots);
	    my $gen = $e[0];
	    my $transcript = $e[1];
	    if (exists($canonicalTx{$gen}) && $transcript eq $canonicalTx{$gen}) {
		$gene = $e[0];
		$transcriptID = $e[1];
		$exon = $e[2];
		$cDNAchange = $e[3];
		$AAchange = $e[4];
	    }
	}
	if(!defined($gene)){#Some of the variants are exonic only in non-canonical transcripts. If that's the case then pick the first annotation. I'm also adding a star to the tx ID to indicate it's non-canonical so we can fix it and make it canonical or something
	    my @e = split(":", $annots[0]);
	    $gene = $e[0];
	    $transcriptID = $e[1]."*";
	    $exon = $e[2];
	    $cDNAchange = $e[3];
	    $AAchange = $e[4];
	}

	if($exon eq "wholegene"){
	    $cDNAchange = "NA";
	    $AAchange = "wholegene";
	}

	# Modify the output to be HGVS compliant (Modified from Donavan's script)
	if($result eq "nonframeshift_insertion"){
	    # cDNA is correct
	    # AA needs changing
	    my ($aa1,$num1,$new_aa)=$AAchange=~/p\.([A-Z])(\d+)delins([A-Z]+)/;
	    my $num2 = $num1+1;
	    open IN1, "<", "$translatedFolder/$gene.translated" || die $logger->fatal("Can not open $gene.translated: $!");
	    while(<IN1>){
		chomp;
		if($_=~/^\>/){last;}
		my @f = split('\t');
		$aa_annotHash{$f[0]}=$f[2];
	    }
	    close IN1;
	    my $aa2 = $aa_annotHash{$num2};
	    if($strandHash{$gene} eq "+"){
		$new_aa=substr($new_aa,1);
	    }else{
		$new_aa=substr($new_aa,0,length($new_aa)-1);
	    }
	    $AAchange="p.${aa1}${num1}_${aa2}${num2}ins${new_aa}";

	}elsif($result eq "frameshift_insertion"){
	    # cDNA: do nothing
	    # Protein: do nothing

	}elsif($result eq "nonframeshift_deletion"){
	    # cDNA:
	    if($strandHash{$gene} eq "+"){
		$cDNAchange.=$ref;
	    }else{
		$cDNAchange.=revcomp($ref);
	    }
	    # Protein: for canonical transcripts only:
	    open IN2, "<", "$translatedFolder/$gene.translated" || die $logger->fatal("Can not open $gene.translated: $!");
	    while(<IN2>){
		chomp;
		if($_=~/^\>/){last;}
		my @f = split('\t');
		$aa_annotHash{$f[0]}=$f[2];
	    }
	    close IN2;
	    my ($num1,$num2)=$AAchange=~/p\.(\d+)_(\d+)del/;
	    my $aa1 = $aa_annotHash{$num1};
	    my $aa2 = $aa_annotHash{$num2};
	    if($num1 == $num2){
		$AAchange = "p.${aa1}${num1}del";
	    }else{
		$AAchange = "p.${aa1}${num1}_${aa2}${num2}del";
	    }

	}elsif($result eq "frameshift_deletion"){
	    # cDNA: 
	    my ($b4_del) = $cDNAchange=~/(.*)del/;
	    if($strandHash{$gene} eq "+"){
		$cDNAchange=$b4_del."del".$ref;
	    }else{
		$cDNAchange=$b4_del."del".revcomp($ref);
	    }
	    # Protein: for canonical transcripts only:
	    open IN3, "<", "$translatedFolder/$gene.translated" || die $logger->fatal("Can not open $gene.translated: $!");
	    while(<IN3>){
		chomp;
		if($_=~/^\>/){last;}
		my @f = split('\t');
		$aa_annotHash{$f[0]}=$f[2];
	    }
	    close IN3;
	    if(!($AAchange=~/fs/)){
		my ($num1)=$AAchange=~/p\.(\d+)_/;
		my $aa1 = $aa_annotHash{$num1};
		$AAchange = "p.${aa1}${num1}fs";
	    }
	}
	my @rest = split("\Q|", $line[8]);
	my $sample = $rest[0];
	my $normal = $rest[1];
	my $key = $sample.":".$normal.":".$chr.":".$startPos.":".$endPos.":".$ref.":".$alt;
	#print "$key\n";
	$annExonicHash{$key} = $result.":".$gene.":".$transcriptID.":".$exon.":".$cDNAchange.":".$AAchange;
    }
    close IN;

    # Read in all variant annotations
    my $annVariantOutput = $sample;
    $annVariantOutput =~ s/.txt/.variant_function/;
    push @deleteFiles, $annVariantOutput;
    open IN, "<", "$annFolder/$annVariantOutput" || die $logger->fatal("Can not open $annVariantOutput: $!");
    while(<IN>){
	chomp;
	my @line = split("\t");
	my $result = $line[0];
	$result =~ s/\ /_/;
	my ($gene, $transcriptID, $exon, $cDNAchange, $AAchange, $ann);
	if($result eq "splicing"){
	    ($gene, $ann) = $line[1] =~ /(.*)\((.*)\)/;
	    if($ann){
		my @annots = split(",", $ann);
		foreach my $annots (sort @annots){
		    my @e = split(":", $annots);
		    $transcriptID = $e[0];
		    $exon = $e[1];
		    $cDNAchange = $e[2];
		}
	    }else{
		$cDNAchange = "";
		$transcriptID = "";
		$exon = "";
		$gene = $line[1];
	    }
	}else{
	    $cDNAchange = "";
	    $transcriptID = "";
	    $exon = "";
	    $gene = $line[1];
	}
	my $chr = $line[2];
	my $startPos = $line[3];
	my $endPos = $line[4];
	my $ref = $line[5];
	my $alt = $line[6];
	my @rest = split("\Q|", $line[7]);
	my $sample = $rest[0];
	my $normal = $rest[1];
	my $key = $sample.":".$normal.":".$chr.":".$startPos.":".$endPos.":".$ref.":".$alt;
	#print "$key\n";
	$annVariantHash{$key} = $result.":".$gene.":".$transcriptID.":".$exon.":".$cDNAchange;
    }
    
    close IN;
    

    # Read in dbSNP annotation
    my $annDBSNPoutput = $sample;
    $annDBSNPoutput =~ s/txt/hg19_snp137NonFlagged_dropped/;
    push @deleteFiles, $annDBSNPoutput;
    my $a = $annDBSNPoutput;
    $a =~ s/dropped/filtered/;
    push @deleteFiles, $a;
    open IN, "<", "$annFolder/$annDBSNPoutput" || die $logger->fatal("Can not open $annDBSNPoutput: $!");
    while(<IN>){
	chomp;
	my @line = split("\t");
	my $dbSNPid = $line[1];
	my $chr = $line[2];
	my $startPos = $line[3];
	my $endPos = $line[4];
	my $ref = $line[5];
	my $alt = $line[6];
	my @rest = split("\Q|", $line[7]);
	my $sample = $rest[0];
	my $normal = $rest[1];
	my $key = $sample.":".$normal.":".$chr.":".$startPos.":".$endPos.":".$ref.":".$alt;
	#print "$key\n";
	$annDBSNPhash{$key} = $dbSNPid;
    } 
    
    close IN;
    
    
    # Read in Cosmic annotation
    my $annCOSMICoutput = $sample;
    $annCOSMICoutput =~ s/txt/hg19_cosmic64_dropped/;
    push @deleteFiles, $annCOSMICoutput;
    $a = $annCOSMICoutput;
    $a =~ s/dropped/filtered/;
    push @deleteFiles, $a;
    open IN, "<", "$annFolder/$annCOSMICoutput" || die $logger->fatal("Can not open $annCOSMICoutput: $!");
    while(<IN>){
	chomp;
	my @line = split("\t");
	my $COSMICid = $line[1];
	my $chr = $line[2];
	my $startPos = $line[3];
	my $endPos = $line[4];
	my $ref = $line[5];
	my $alt = $line[6];
	my @rest = split("\Q|", $line[7]);
	my $sample = $rest[0];
	my $normal = $rest[1];
	my $key = $sample.":".$normal.":".$chr.":".$startPos.":".$endPos.":".$ref.":".$alt;
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
    open IN, "<", "$annFolder/$ann1000Goutput" || die $logger->fatal("Can not open $ann1000Goutput: $!");
    while(<IN>){
	chomp;
	my @line = split("\t");
	my $MAF = $line[1];
	my $chr = $line[2];
	my $startPos = $line[3];
	my $endPos = $line[4];
	my $ref = $line[5];
	my $alt = $line[6];
	my @rest = split("\Q|", $line[7]);
	my $sample = $rest[0];
	my $normal = $rest[1];
	my $key = $sample.":".$normal.":".$chr.":".$startPos.":".$endPos.":".$ref.":".$alt;
	#print "$key\n";
	$ann1000Ghash{$key} = $MAF;
    } 
    
    close IN;

    # Read in Mutation Assesor annotation
    my $annMAoutput = $sample;
    $annMAoutput =~ s/txt/hg19_ljb2_ma_dropped/;
    push @deleteFiles, $annMAoutput;
    $a = $annMAoutput;
    $a =~ s/dropped/filtered/;
    push @deleteFiles, $a;
    open IN, "<", "$annFolder/$annMAoutput" || die $logger->fatal("Can not open $annMAoutput: $!");
    while(<IN>){
	chomp;
	my @line = split("\t");
	my $MA = $line[1];
	my $chr = $line[2];
	my $startPos = $line[3];
	my $endPos = $line[4];
	my $ref = $line[5];
	my $alt = $line[6];
	my @rest = split("\Q|", $line[7]);
	my $sample = $rest[0];
	my $normal = $rest[1];
	my $key = $sample.":".$normal.":".$chr.":".$startPos.":".$endPos.":".$ref.":".$alt;
	#print "$key\n";
	$annMAhash{$key} = $MA;
    } 
    close IN;
    
    # Read in AV SIFT scores
    my $annAVSIFToutput = $sample;
    $annAVSIFToutput =~ s/txt/hg19_avsift_dropped/;
    push @deleteFiles, $annAVSIFToutput;
    $a = $annAVSIFToutput;
    $a =~ s/dropped/filtered/;
    push @deleteFiles, $a;
    open IN, "<", "$annFolder/$annAVSIFToutput" || die $logger->fatal("Can not open $annAVSIFToutput: $!");
    while(<IN>){
	chomp;
	my @line = split("\t");
	my $AVSIFT = $line[1];
	my $chr = $line[2];
	my $startPos = $line[3];
	my $endPos = $line[4];
	my $ref = $line[5];
	my $alt = $line[6];
	my @rest = split("\Q|", $line[7]);
	my $sample = $rest[0];
	my $normal = $rest[1];
	my $key = $sample.":".$normal.":".$chr.":".$startPos.":".$endPos.":".$ref.":".$alt;
	#print "$key\n";
	$annAVSIFThash{$key} = $AVSIFT;
    } 
    
    close IN;
    
    # Read in Phylo P (conservation) scores
    my $annPhyloPoutput = $sample;
    $annPhyloPoutput =~ s/txt/hg19_ljb2_phylop_dropped/;
    push @deleteFiles, $annPhyloPoutput;
    $a = $annPhyloPoutput;
    $a =~ s/dropped/filtered/;
    push @deleteFiles, $a;
    open IN, "<", "$annFolder/$annPhyloPoutput" || die $logger->fatal("Can not open $annPhyloPoutput: $!");
    while(<IN>){
	chomp;
	my @line = split("\t");
	my $PhyloP = $line[1];
	my $chr = $line[2];
	my $startPos = $line[3];
	my $endPos = $line[4];
	my $ref = $line[5];
	my $alt = $line[6];
	my @rest = split("\Q|", $line[7]);
	my $sample = $rest[0];
	my $normal = $rest[1];
	my $key = $sample.":".$normal.":".$chr.":".$startPos.":".$endPos.":".$ref.":".$alt;
	#print "$key\n";
	$annPhyloPhash{$key} = $PhyloP;
    } 
    
    close IN;
    
    # Read in LRT scores
    my $annLRToutput = $sample;
    $annLRToutput =~ s/txt/hg19_ljb2_lrt_dropped/;
    push @deleteFiles, $annLRToutput;
    $a = $annLRToutput;
    $a =~ s/dropped/filtered/;
    push @deleteFiles, $a;
    open IN, "<", "$annFolder/$annLRToutput" || die $logger->fatal("Can not open $annLRToutput: $!");
    while(<IN>){
	chomp;
	my @line = split("\t");
	my $LRT = $line[1];
	my $chr = $line[2];
	my $startPos = $line[3];
	my $endPos = $line[4];
	my $ref = $line[5];
	my $alt = $line[6];
	my @rest = split("\Q|", $line[7]);
	my $sample = $rest[0];
	my $normal = $rest[1];
	my $key = $sample.":".$normal.":".$chr.":".$startPos.":".$endPos.":".$ref.":".$alt;
	#print "$key\n";
	$annLRThash{$key} = $LRT;
    } 
    
    close IN;
    
    
    # Read in GERP++ scores
    my $annGERPoutput = $sample;
    $annGERPoutput =~ s/txt/hg19_ljb2_gerp++_dropped/;
    push @deleteFiles, $annGERPoutput;
    $a = $annGERPoutput;
    $a =~ s/dropped/filtered/;
    push @deleteFiles, $a;
    open IN, "<", "$annFolder/$annGERPoutput" || die $logger->fatal("Can not open $annGERPoutput: $!");
    while(<IN>){
	chomp;
	my @line = split("\t");
	my $GERP = $line[1];
	my $chr = $line[2];
	my $startPos = $line[3];
	my $endPos = $line[4];
	my $ref = $line[5];
	my $alt = $line[6];
	my @rest = split("\Q|", $line[7]);
	my $sample = $rest[0];
	my $normal = $rest[1];
	my $key = $sample.":".$normal.":".$chr.":".$startPos.":".$endPos.":".$ref.":".$alt;
	#print "$key\n";
	$annGERPhash{$key} = $GERP;
    } 
    
    close IN;
    $sample =~s/.txt/_annovarAnnotated.txt/;

    open OUT, ">", $sample || die $logger->fatal("Can not open $sample: $!");
    
    my $H1 = shift @header;
    my $H2 = shift @header;
    my $H3 = shift @header;
    my $H4 = shift @header;
    my $H5 = shift @header;
    my $H6 = shift @header;
    my $H7 = $header[3];;
    my @header_new = splice(@header, 9);
    my @Hrest = join("\t", @header_new);
    print OUT "$H1\t$H2\t$H3\t$H4\t$H5\t$H6\tVariantClass\tGene\tExon\tCall_Confidence\tComments\tTranscriptID\tcDNAchange\tAAchange\tdbSNP_ID\tCosmic_ID\t1000G_MAF\tMutation_Assesor\tAV_SIFT_score\tPhyloP_score\tLRT_score\tGERP++_score\t$H7\t@Hrest\n";
    
    
    # Merge Annotations...
    foreach my $key (keys %inputHash){
	my ($sample, $normal, $chr, $pos, $ref, $alt) = split(":", $key);
	my $annotKey;
	if(length($ref) > length ($alt)) {# deletion
	    my $newref = substr($ref, 1);
	    my $posEnd = $pos+length($newref);
	    my $posStart = $pos+1;
	    $annotKey = $sample.":".$normal.":".$chr.":".$posStart.":".$posEnd.":".$newref.":"."-";
	    #print "$key\t$annotKey\n";	
	}elsif(length($ref) < length($alt)){ #insertion
	    my $newalt = substr($alt, 1);
	    my $posEnd = $pos;
	    $annotKey = $sample.":".$normal.":".$chr.":".$pos.":".$pos.":-:".$newalt;
	    
	}else {
	    $annotKey = $sample.":".$normal.":".$chr.":".$pos.":".$pos.":".$ref.":".$alt;
	}
	if(exists($annExonicHash{$annotKey})){
	    my ($result, $gene, $txID, $exon, $cDNAchange, $AAchange) = split(":", $annExonicHash{$annotKey});
	    print OUT "$sample\t$normal\t$chr\t$pos\t$ref\t$alt\t$result\t$gene\t$exon\t\t\t$txID\t$cDNAchange\t$AAchange\t";
	}elsif(exists($annVariantHash{$annotKey})){
	    my ($result, $gene, $txID, $exon, $cDNAchange) = split(":", $annVariantHash{$annotKey});
	    print OUT "$sample\t$normal\t$chr\t$pos\t$ref\t$alt\t$result\t$gene\t$exon\t\t\t$txID\t$cDNAchange\t\t";
	}
	if(exists($annDBSNPhash{$annotKey})){
	    print OUT "$annDBSNPhash{$annotKey}\t";
	}else{
	    print OUT "\t";
	}
	if(exists($annCOSMIChash{$annotKey})) {
	    print OUT "$annCOSMIChash{$annotKey}\t";
	}else{
	    print OUT "\t";
	}
	if(exists($ann1000Ghash{$annotKey})){
	    print OUT "$ann1000Ghash{$annotKey}\t";
	}else{
	    print OUT "\t";
	}
	if(exists($annMAhash{$annotKey})){
	    print OUT "$annMAhash{$annotKey}\t";
	}else{
	    print OUT "\t";
	}
	if(exists($annAVSIFThash{$annotKey})) {
	    print OUT "$annAVSIFThash{$annotKey}\t";
	}else{
	    print OUT "\t";
	}
	if(exists($annPhyloPhash{$annotKey})){
	    print OUT "$annPhyloPhash{$annotKey}\t";
	}else{
	    print OUT "\t";
	}
	if(exists($annLRThash{$annotKey})){
	    print OUT "$annLRThash{$annotKey}\t";
	}else{
	    print OUT "\t";
	}
	if(exists($annGERPhash{$annotKey})){
	    print OUT "$annGERPhash{$annotKey}\t";
	}else{
	    print OUT "\t";
	}
	my @rest = split("\Q|", $inputHash{$key});
	my $failure_reason =  $rest[3];
	my @rest_new = splice(@rest, 9);
	@rest_new = join("\t", @rest_new);
	print OUT "$failure_reason\t@rest_new\n";
    }
    $logger->info("Annotations are merged...");
    return($sample, \@deleteFiles);
}

#####################################
#####################################
# Filter annotations based on several filters


sub FilterAnnotations  {
    my ($mergedFile, $tumorSupressors, $outdir, $VF_threshold_hotspot, $VF_Tthreshold, $TNfreqRatioThreshold, $MAFthreshold, $AD_Tthreshold_indel_high, $AD_Tthreshold_indel_low, $VF_Tthreshold_indel_high, $VF_Tthreshold_indel_low, $VF_Tthreshold_indel_hotspot_high,$VF_Tthreshold_indel_hotspot_low, $hotspots, $patientIDPerSampleId,$classPerPatientIdSampleId) = @_;
    $logger->info("Variants are being filtered...");
    my ($sample) = $mergedFile =~ /(.*)_annovarAnnotated.txt/;
    my $filteredFile = $sample."_AnnotatedFiltered.txt";
    my $bedFile = $sample."_AnnotatedFiltered.bed";
    my $droppedFile = $sample."_AnnotatedDropped.txt";
    my %patientIDPerSampleId = %$patientIDPerSampleId;
    my %classPerPatientIdSampleId = %$classPerPatientIdSampleId;
    my (%hotspotHash, %hotspotExonHash, %hotspotIndelHash, %tumorSupressorHash);

    open HOTSPOT, "<", $hotspots || die  $logger->fatal("Can not open $hotspots: $!");
    while(<HOTSPOT>){
	chomp;
	my @line = split("\t");
	next if($line[0] eq "Index");
	my ($mut, undef) = split(" ", $line[2]);
	my $gene = $line[1];
	next if($mut=~/Validation/);
	if ($mut=~/_indel/){
	    $mut=~s/_indel//;
	    $hotspotExonHash{$gene.":".$mut} = $mut;
	}elsif($mut=~/fs/ || $mut =~ /_/ || $mut =~ /del/ || $mut =~/>/){
	    $hotspotIndelHash{$gene.":".$mut} = $mut;
	}else{
	    my ($mut2) = $mut =~/([A-Z]\d+)\>*[A-Za-z]*.*/; # storing only the original amino acid and the position
	    $hotspotHash{$gene.":".$mut2} = $mut2 if($mut2);
	    #print "$gene\t$mut\t$mut2\n";
	}
    }


    open IN, "<", "$tumorSupressors" || die $logger->fatal("Can not open $tumorSupressors: $!");
    while(<IN>){
	chomp;
	$tumorSupressorHash{$_} = 1;
    }
    close IN;

    open FILTERED, ">", "$outdir/$filteredFile" || die $logger->fatal("Can not open $filteredFile: $!");
    open BEDFILE, ">", "$outdir/$bedFile" || die $logger->fatal("Can not open $bedFile: $!");
    open DROPPED, ">", "$outdir/$droppedFile" || die $logger->fatal("Can not open $droppedFile: $!");
    open VARIANTS,"<", "$outdir/$mergedFile" || die $logger->fatal("Can not open $mergedFile: $!");
    while(<VARIANTS>){
	chomp;
	my $header;
	if(m/^Sample\t/){
	    print FILTERED "$_\n";
	    print DROPPED "$_\n";
	    next;
	}
	my @line = split("\t");
	my $sample = $line[0];
	my $normal = $line[1];
	my $chr = $line[2];
	my $pos = $line[3];
	my $ref = $line[4];
	my $alt = $line[5];
	# Filter SNVs
	if (length($ref) == length($alt) && length($ref) == 1){
	    #determine if the sample is matched or unmatched
	    my $patientID = $patientIDPerSampleId{$sample};
	    my $key= $normal . ":" .  $patientID;
	    my $tag;
	    if((exists $classPerPatientIdSampleId{$key}) and ($classPerPatientIdSampleId{$key} =~ /Normal/)){
		$tag = "Matched";
	    }else{
		$tag = "Unmatched";
	    }
	    my $gene = $line[7];
	    my $variantClass = $line[6];
	    if ($gene eq "TERT" && $variantClass eq "upstream"){
		print FILTERED "$_\n";
		print BEDFILE "$chr\t$pos\t$pos\t$sample:$normal\t\t+\n";
		next;
	    }
	    # Filter for variant class
	    if ($variantClass eq "frameshift_deletion" || $variantClass eq "frameshift_insertion" || $variantClass eq "nonframeshift_deletion" || $variantClass eq "nonframeshift_insertion" || $variantClass eq "nonsynonymous_SNV" || $variantClass eq "splicing" || $variantClass eq "stopgain_SNV") {
		my ($AAchange) = $line[13] =~ /p\.([A-Z]\d+).*/;  # Grabbing only G12 instead of G12D (includes fs or * notations)
		$AAchange = "NA" if ($variantClass eq "splicing"); # This is to fix the error that comes up when splicing doesn't have a protein annotations

		my $AD_T = $line[20];
		my $VF_T = $line[30];
		my $TNfreqRatio = $line[37];
		my $MAF = $line[16];
		if(!$MAF) {$MAF =0;}
		if ($TNfreqRatio == 0){
		    $TNfreqRatio = $TNfreqRatioThreshold;
		}

		# Filter for VF_T and TNfreqRatio: two-tiered filtering. 
		if(exists($hotspotHash{$gene.":".$AAchange})){ # isHotSpot???
		    if($tag eq "Unmatched"){
			if ($VF_T >= $VF_threshold_hotspot && $TNfreqRatio >= $TNfreqRatioThreshold && $MAF <  $MAFthreshold) {
			    print FILTERED "$_\n";
			    print BEDFILE "$chr\t$pos\t$pos\t$sample:$normal\t\t+\n";
			}else{
			    print DROPPED "$_\n";
			}
		    }elsif($tag eq "Matched"){ # There are some matched variants with MAF > 0.01 
			if ($VF_T >= $VF_threshold_hotspot && $TNfreqRatio >= $TNfreqRatioThreshold){
			    print FILTERED "$_\n";
			    print BEDFILE "$chr\t$pos\t$pos\t$sample:$normal\t\t+\n";
			}else{
			    print DROPPED "$_\n";
			}
		    }
		}elsif(exists($tumorSupressorHash{$gene}) && ($variantClass eq "frameshift_deletion" || $variantClass eq "frameshift_insertion" || $variantClass eq "stopgain_SNV")){ #IsTumorSupressor and is frameshift????
		    if($tag eq "Unmatched"){
			if ($VF_T >= $VF_threshold_hotspot && $TNfreqRatio >= $TNfreqRatioThreshold && $MAF <  $MAFthreshold) {
			    print FILTERED "$_\n";
			    print BEDFILE "$chr\t$pos\t$pos\t$sample:$normal\t\t+\n";
			}else{
			    print DROPPED "$_\n";
			}
		    }elsif($tag eq "Matched"){ # There are some matched variants with MAF > 0.01 
			if ($VF_T >= $VF_threshold_hotspot && $TNfreqRatio >= $TNfreqRatioThreshold){
			    print FILTERED "$_\n";
			    print BEDFILE "$chr\t$pos\t$pos\t$sample:$normal\t\t+\n";
			}else{
			    print DROPPED "$_\n";
			}
		    }
		}else{ # For the rest
		    # Filter for MAF if the sample is unmatched
		    if($tag eq "Unmatched") {
			if ($VF_T >= $VF_Tthreshold && $TNfreqRatio >= $TNfreqRatioThreshold && $MAF <  $MAFthreshold) {
			    print FILTERED "$_\n";
			    print BEDFILE "$chr\t$pos\t$pos\t$sample:$normal\t\t+\n";
			}else{
			    print DROPPED "$_\n";
			}
		    }elsif($tag eq "Matched"){ # There are some matched variants with MAF > 0.01 
			if ($VF_T >= $VF_Tthreshold && $TNfreqRatio >= $TNfreqRatioThreshold){
			    print FILTERED "$_\n";
			    print BEDFILE "$chr\t$pos\t$pos\t$sample:$normal\t\t+\n";
			}else{
			    print DROPPED "$_\n";
			}
		    }
		}
		}else{
		    print DROPPED "$_\n";
		}
	    # Filter INDELs
	}elsif(length($ref) > length($alt) || length($ref) < length($alt)){
	    #determine if the sample is matched or unmatched
	    my $patientID = $patientIDPerSampleId{$sample};
	    my $key= $normal . ":" .  $patientID;
	    my $tag;
	    if((exists $classPerPatientIdSampleId{$key}) and ($classPerPatientIdSampleId{$key} =~ /Normal/)){
		$tag = "Matched";
	    }else{
		$tag = "Unmatched";
	    }
	    my $exon = $line[8];
	    my $gene = $line[7];
	    my $variantClass = $line[6];
	    if ($gene eq "TERT" && $variantClass eq "upstream"){
		print FILTERED "$_\n";
		print BEDFILE "$chr\t$pos\t$pos\t$sample:$normal\t\t+\n";
		next;
	    }
	    # Filter for variant class
	    if ($variantClass eq "frameshift_deletion" || $variantClass eq "frameshift_insertion" || $variantClass eq "nonframeshift_deletion" || $variantClass eq "nonframeshift_insertion" || $variantClass eq "nonsynonymous_SNV" || $variantClass eq "splicing" || $variantClass eq "stopgain_SNV") {
		my ($AAchange) = $line[13] =~ /p\.(.*)/; #INDEL annotation can be varied. Need to make sure the formats are the same between Donavan's file and HGVS compliant Annovar output
		$AAchange = "NA" if($variantClass eq "splicing");
		my $AD_T = $line[29];
		my $VF_T = $line[30];
		my $TNfreqRatio = $line[37];
		my $MAF = $line[16];
		#print "\t$line[13]\t$AAchange\t$AD_T\t$VF_T\t$MAF\t$TNfreqRatio\n";
		if(!$MAF) {$MAF =0;}
		if ($MAF > 0.01){
			print DROPPED "$_\n";
			next;
		}
		if ($TNfreqRatio == 0){
		    $TNfreqRatio = $TNfreqRatioThreshold;
		}
		# Print hotspotINDELs directly
		if(exists($hotspotIndelHash{$gene.":".$AAchange})){
		    if($tag eq "Unmatched"){
			if(($AD_T >= $AD_Tthreshold_indel_high && $VF_T >= $VF_Tthreshold_indel_hotspot_low && $MAF < $MAFthreshold) || ($AD_T >= $AD_Tthreshold_indel_low && $VF_T >= $VF_Tthreshold_indel_hotspot_high && $MAF < $MAFthreshold)){
			    print FILTERED "$_\n";
			    print BEDFILE "$chr\t$pos\t$pos\t$sample:$normal\t\t+\n";
			}else{
			    print DROPPED "$_\n";
			}
		    }elsif($tag eq "Matched"){
			if(($AD_T >= $AD_Tthreshold_indel_high && $VF_T >= $VF_Tthreshold_indel_hotspot_low) || ($AD_T >= $AD_Tthreshold_indel_low && $VF_T >= $VF_Tthreshold_indel_hotspot_high)){
			    print FILTERED "$_\n";
			    print BEDFILE "$chr\t$pos\t$pos\t$sample:$normal\t\t+\n";
			}else{
			    print DROPPED "$_\n";
			}
		    }
		}elsif(exists($tumorSupressorHash{$gene}) && ($variantClass eq "frameshift_deletion" || $variantClass eq "frameshift_insertion" || $variantClass eq "stopgain_SNV")){ #IsTumorSupressor and is frameshift????
		    if($tag eq "Unmatched"){
			if(($AD_T >= $AD_Tthreshold_indel_high && $VF_T >= $VF_Tthreshold_indel_hotspot_low && $MAF < $MAFthreshold) || ($AD_T >= $AD_Tthreshold_indel_low && $VF_T >= $VF_Tthreshold_indel_hotspot_high && $MAF < $MAFthreshold)){
			    print FILTERED "$_\n";
			    print BEDFILE "$chr\t$pos\t$pos\t$sample:$normal\t\t+\n";
			}else{
			    print DROPPED "$_\n";
			}
		    }elsif($tag eq "Matched"){
			if(($AD_T >= $AD_Tthreshold_indel_high && $VF_T >= $VF_Tthreshold_indel_hotspot_low) || ($AD_T >= $AD_Tthreshold_indel_low && $VF_T >= $VF_Tthreshold_indel_hotspot_high)){
			    print FILTERED "$_\n";
			    print BEDFILE "$chr\t$pos\t$pos\t$sample:$normal\t\t+\n";
			}else{
			    print DROPPED "$_\n";
			}
		    }
		}elsif(exists($hotspotExonHash{$gene.":".$exon})){ # Should we add a length filter here? I have it in my notes as indels around 3-12bp
		    print FILTERED "$_\n";
		    print BEDFILE "$chr\t$pos\t$pos\t$sample:$normal\t\t+\n";
		    next;
		}elsif(($AD_T >= $AD_Tthreshold_indel_high && $VF_T >= $VF_Tthreshold_indel_low) || ($AD_T >= $AD_Tthreshold_indel_low && $VF_T >= $VF_Tthreshold_indel_high)){ # if the INDEL is not hotspot then higher stringency
		    print FILTERED "$_\n";
		    print BEDFILE "$chr\t$pos\t$pos\t$sample:$normal\t\t+\n";
		}else{
		    print DROPPED "$_\n";
		    #print "$sample\t$gene\t$AD_T\t$AD_Tthreshold_indel_high\t$AD_Tthreshold_indel_low\t$VF_T\t$VF_Tthreshold_indel_low\t$VF_Tthreshold_indel_high\n";
		}
	    }else{
		print DROPPED "$_\n";
	    }
	}else{
	    #die "it shouldn't come here!\n";
	}
    }
    close IN;
    close FILTERED;
    close BEDFILE;
    close DROPPED;
    $logger->info("Variants are filtered based on the two-tired filtering scheme. Filtered variants can be found in $filteredFile and variants that did not meet the requried criteria can be found in $droppedFile");
    return($filteredFile, $bedFile);
}


####################################
####################################
# Generate per patient VCF files

sub GenerateVCFfiles {
    my ($filteredAnnoFile, $outdir, $VF_Tthreshold, $TNfreqRatioThreshold, $patientIDPerSampleId,$classPerPatientIdSampleId) = @_;
    my (%samples, @header);
    if(!-e "$outdir/VCFs"){
	`mkdir VCFs`;
    }
    $logger->info("Individual VCF files are being generated......");
    my %patientIDPerSampleId = %$patientIDPerSampleId;
    my %classPerPatientIdSampleId = %$classPerPatientIdSampleId;
    print "$filteredAnnoFile\n";
    open IN, "<", $filteredAnnoFile || die $logger->fatal("Can not open $filteredAnnoFile: $!");
    while(<IN>){
	chomp;
	if(m/^Sample\t/){
	    @header = split("\t");
	    next;
	}
	my @line = split("\t");
	my $sample = $line[0];
	my $normal = $line[1];
	my $chr = $line[2];
	my $pos = $line[3];
	my $ref = $line[4];
	my $alt = $line[5];
	$samples{$sample.":".$normal}{$chr.":".$pos.":".$ref.":".$alt} = $_;
    }
    close IN;
    foreach my $key1 (sort keys %samples){
	my ($sample, $normal) = split(":", $key1);
	my $normalIndex;
	foreach my $i (0..$#header){
	    if($header[$i] eq $normal){
		$normalIndex = $i;
	    }
	}
	#determine if the sample is matched or unmatched
	my $patientID = $patientIDPerSampleId{$sample};
	my $key= $normal . ":" .  $patientID;
	my $tag;
	if((exists $classPerPatientIdSampleId{$key}) and ($classPerPatientIdSampleId{$key} =~ /Normal/)){
	    $tag = "Matched";
	}else{
	    $tag = "Unmatched";
	}
	#print "$sample\t$normal\t$tag\t$normalIndex\t$header[$normalIndex]\n";
	my $vcfFile = $sample."_".$tag."_annovar.vcf";
	open OUT, ">", "$outdir/VCFs/$vcfFile" || die $logger->fatal("Can not open $vcfFile: $!");
	my $header=  <<"END";
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

	foreach my $key2 (sort my_sort  keys %{$samples{$key1}}){
	    my @line = split("\t", $samples{$key1}->{$key2});
	    my $chr = $line[2];
	    my $pos = $line[3];
	    my $ref = $line[4];
	    my $alt = $line[5];
	    my $varClass = $line[6];
	    my $gene = $line[7];
	    my $exon = $line[8];
	    my $txID = $line[11];
	    my $cDNAchange = $line[12];
	    my $AAchange = $line[13];
	    my $dbSNP = $line[14];
	    if($dbSNP eq ""){
		$dbSNP = ".";
	    }
	    my $cosmic = $line[15];
	    my $MAF = $line[16];
	    my $failReason = $line[22];
	    my $DP_N = $line[23];
	    my $DP_N_Ref = $line[24];
	    my $DP_N_Alt = $line[25];
	    my $VF_N = $line[26];
	    my $DP_T = $line[27];
	    my $DP_T_Ref = $line[28];
	    my $DP_T_Alt = $line[29];
	    my $VF_T = $line[30];
	    my $STR_T_Ref_p = $line[31];
	    my $STR_T_Ref_n = $line[32];
	    my $STR_T_Alt_p = $line[33];
	    my $STR_T_Alt_n = $line[34];
	    my $DP_N_AGG = $line[35];
	    my $VF_N_MEDIAN = $line[36];
	    my $TNfreqRatio = $line[37];
	    my $occurance = $line[38];
	    my $normalInfo = $line[$normalIndex];
	    print OUT "$chr\t$pos\t$dbSNP\t$ref\t$alt\t.\t.\t.\tDP_N=$DP_N($DP_N_Ref,$DP_N_Alt);VF_N=$VF_N;DP_T=$DP_T($DP_T_Ref,$DP_T_Alt);VF_T=$VF_T;ANNOVAR=($gene|$txID|$varClass|$exon|$cDNAchange|$AAchange|$dbSNP|$MAF|$cosmic);STRBIAS_REF=($STR_T_Ref_p,$STR_T_Ref_n);STRBIAS_ALT=($STR_T_Alt_p,$STR_T_Alt_n);VF_N_MEDIAN=$VF_N_MEDIAN;TN_VF_RATIO=$TNfreqRatio;OCCURANCE_N=$occurance\t$normalInfo\n";
	}
	close OUT;
	eval{`$java -Xmx2g -jar $IGVtools sort $outdir/VCFs/$vcfFile $outdir/VCFs/$vcfFile.2`};
	if($@){
	    $logger->fatal("Command: $java -Xmx2g -jar $IGVtools sort $outdir/VCFs/$vcfFile $outdir/VCFs/$vcfFile.2 can not be run! Error: $@");
	    exit(1);
	}else{
	    $logger->info("Command: $java -Xmx2g -jar $IGVtools sort $outdir/VCFs/$vcfFile $outdir/VCFs/$vcfFile.2 is being run!");
	    
	}
	eval{`mv $outdir/VCFs/$vcfFile.2 $outdir/VCFs/$vcfFile`};
	if($@){
	    $logger->fatal("Command: mv $outdir/VCFs/$vcfFile.2 $outdir/VCFs/$vcfFile can not be run! Error: $@");
	    exit(1);
	}else{
	    $logger->info("Command: mv $outdir/VCFs/$vcfFile.2 $outdir/VCFs/$vcfFile is being run!");
	   
	}
	eval{`$java -Xmx2g -jar $IGVtools index $outdir/VCFs/$vcfFile`};
	if($@){
	    $logger->fatal("Command: $java -Xmx2g -jar $IGVtools index $outdir/VCFs/$vcfFile can not be run! Error: $@");
	    exit(1);
	}else{
	    $logger->info("Command: $java -Xmx2g -jar $IGVtools index $outdir/VCFs/$vcfFile is being run!");
	   
	}
    }
    $logger->info("VCF file generation has been completed.");
}


###################################
###################################
# Generate per patient gene coverage files


sub GenerateGeneCoverage{
    my ($titleFile, $outdir, $coverageFile, $exonIntervals) = @_;
    $logger->info("Gene based coverage files are being created.");
    my (%barcodes, %txIDs, @header);
    if(!-e "$outdir/Coverage"){
	`mkdir $outdir/Coverage`;
    }
    open IN, "<", $titleFile || die $logger->fatal("Can not open $titleFile: $!");
    while(<IN>){
	chomp;
	next if(/Barcode/);
	my @line = split("\t");
	my $barcode = $line[0];
	my $sample = $line[2];
	$barcodes{$barcode} = $sample;
    }
    close IN;
    open IN2, "<", $exonIntervals || die $logger->fatal("Can not open $exonIntervals: $!");
    while(<IN2>){
	chomp;
	next if(/^\@SQ/ || /^\@HD/);
	my ($chr, $start, $stop, $strand, $annotation) = split("\t");
	my ($gene, $txID, $exon, $aa_range) = split(":", $annotation);
	$txIDs{$gene} = $txID;
    }
    close IN2;
    foreach my $barcode (sort keys %barcodes){
	my $sample = $barcodes{$barcode};
	my $file = $sample."_gene.coverage.txt";
	open OUT, ">", "$outdir/Coverage/$file" || die $logger->fatal("Can not open $sample.coverage.txt: $!");
	open IN, "<", $coverageFile || die $!;
	my $k = 1;
	while(<IN>){
	    chomp;
	    if(/Gene/) {
		@header = split("\t"); 
		next;
	    }
	    my $index;
	    foreach my $i (0..$#header){
		if($barcode eq $header[$i]){
		    $index = $i;
		}
	    }
	    my @line = split("\t");
	    my $gene = $line[0];
	    my $txID = $txIDs{$gene};
	    my $cov = sprintf "%.0f", $line[$index];
	    my $result;
	    if($cov >= 100){$result = "GOOD";}else{$result = "POOR";}
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

sub GenerateExonCoverage{
    my($filteredAnnoFile, $covFile, $translatedFolder, $exonIntervals, $outdir, $titleFile, $coverage_threshold) = @_;
    my(%barcodes, %canonicalExons, @header, $index, %coverage, %variants);
    $logger->info("Exon based coverage files are being generated.");
    open IN, "<", $titleFile || die $logger->fatal("Can not open $titleFile: $!");
    while(<IN>){
	chomp;
	next if(/Barcode/);
	my @line = split("\t");
	my $barcode = $line[0];
	my $sample = $line[2];
	$barcodes{$barcode} = $sample;
    }
    close IN;

    open IN2, "<", $exonIntervals || die $logger->fatal("Can not open $exonIntervals: $!");
    while(<IN2>){
	chomp;
	next if(/^\@SQ/ || /^\@HD/);
	my ($chr, $start, $stop, $strand, $annotation) = split("\t");
	my ($gene, $txID, $exon, $aa_range) = split(":", $annotation);
	$canonicalExons{$chr}{$start.":".$stop.":".$strand} = $txID.":".$gene.":".$exon.":".$strand.":".$aa_range;
    }
    close IN2;
    
    open IN, "<", $filteredAnnoFile || die $logger->fatal("Can not open $filteredAnnoFile: $!");
    while(<IN>){
	chomp;
	if(m/^Sample\t/){
	    @header = split("\t");
	    next;
	}
	my @line = split("\t");
	my $sample = $line[0];
	my $normal = $line[1];
	my $chr = $line[2];
	my $pos = $line[3];
	my $ref = $line[4];
	my $alt = $line[5];
	my $gene = $line[7];
	my $exon = $line[8];
	push @{$variants{$sample}}, $gene.":".$exon;
    }
    close IN;

    foreach my $barcode (sort keys %barcodes){
	my $sample = $barcodes{$barcode};
	my $file = $sample."_exon.coverage.txt";
	print "$file\n";
	print "$covFile\n";
	open OUT, ">", "$outdir/Coverage/$file" || die $logger->fatal("Can not open $file: $!");
	open IN, "<", $covFile || die $!;
	while(<IN>){
	    chomp;
	    if(m/Target/){
		@header = split("\t");
		next;
	    }
	    foreach my $i (0..$#header){
		if($barcode eq $header[$i]){
		    $index = $i;
		}
	    }
	    my @line = split("\t");
	    my ($chr, $range) = split(":", $line[1]);
	    my ($start, $stop) = split("-", $range);	
	    my ($txID, $gene, $exon, $strand, $exonStart, $exonStop, $aa_range);
	    foreach my $exonalRange (sort keys %{$canonicalExons{$chr}}){
		($exonStart, $exonStop, $strand) = split(":", $exonalRange);
		if ($strand eq "+"){
		    #----ES*---->----->----ES-----
		    if($exonStart <= $start && $stop <= $exonStop){
			($txID, $gene,$exon, $strand, $aa_range) = split(":", $canonicalExons{$chr}->{$exonalRange});
			
			#---->---ES*----ES--->-----
		    }elsif($start <= $exonStart && $exonStop <= $stop){
			($txID, $gene, $exon, $strand, $aa_range) = split(":", $canonicalExons{$chr}->{$exonalRange});
			
			#---->---ES*--->---ES------
		    }elsif($start <= $exonStart && $stop <= $exonStop && $exonStart < $stop){
			($txID, $gene, $exon, $strand, $aa_range) = split(":", $canonicalExons{$chr}->{$exonalRange});
			
			#---ES*--->----ES--->------
		    }elsif($exonStart <= $start && $exonStop <= $stop && $start < $exonStop){
			($txID, $gene, $exon, $strand, $aa_range) = split(":", $canonicalExons{$chr}->{$exonalRange});
		    }
		}elsif($strand eq "-"){
		    #-----ES---<------<----ES*-----
		    if ($exonStop <= $stop && $start <= $exonStart) {
			($txID, $gene,$exon, $strand, $aa_range) = split(":", $canonicalExons{$chr}->{$exonalRange});
			#-----<----ES-----<------ES*-----
		    }elsif($stop <= $exonStop && $exonStop < $start  && $start <= $exonStart){
			($txID, $gene,$exon, $strand, $aa_range) = split(":", $canonicalExons{$chr}->{$exonalRange});

			#------ES----<------ES*------<-------
		    }elsif($exonStop <= $stop && $stop < $exonStart && $exonStart <= $start){
			($txID, $gene,$exon, $strand, $aa_range) = split(":", $canonicalExons{$chr}->{$exonalRange});

			#-----<----ES------ES*-----<------
		    }elsif($stop <= $exonStop && $exonStart <= $start){
			($txID, $gene,$exon, $strand, $aa_range) = split(":", $canonicalExons{$chr}->{$exonalRange});
		    }
		}else{
		    die "Strand info is dead wrong!\n";
		}
	    }
	    if ($gene && $exon && $strand && $txID) {
		$coverage{"$chr:$gene:$txID:$exon:$strand:$aa_range"} = $range.":".$line[$index];
	    } else {
	#	print OUT "$chr\t$range\t not assigned an exon\n"; #this is for when the region falls into a UTR or an intron
	    }
	}
	my %genes = map { $_ => 1 } @{$variants{$sample}};
	foreach my $key (sort my_sort2 keys %coverage){
	    my ($chr, $gene, $txID, $exon, $strand, $aa_range) = split(":", $key);
	    my ($range, $cov) = split(":", $coverage{$key});
	    if(exists($genes{"$gene:$exon"}) && $cov > $coverage_threshold){
		print OUT "$chr\t$gene\t$txID\t$exon\t$range\t$aa_range\t$strand\t$cov\tGOOD\tPOS\n";
	    }elsif(exists($genes{"$gene:$exon"}) && $cov < $coverage_threshold){
		print OUT "$chr\t$gene\t$txID\t$exon\t$range\t$aa_range\t$strand\t$cov\tPOOR\tPOS\n";
	    }elsif(!exists($genes{"$gene:$exon"}) && $cov > $coverage_threshold){
		print OUT "$chr\t$gene\t$txID\t$exon\t$range\t$aa_range\t$strand\t$cov\tGOOD\tNEG\n";
	    }elsif(!exists($genes{"$gene:$exon"}) && $cov < $coverage_threshold){
		print OUT "$chr\t$gene\t$txID\t$exon\t$range\t$aa_range\t$strand\t$cov\tPOOR\tINDET\n";
	    }
	}
    }
    
    close IN;
    close OUT;
    $logger->info("Exon based coverage files are generated.");
}

####################################
## Generate clinical reports
sub GenerateClinicalReports {
    my ($outdir, $hotspots) = @_;

    my(%hotspotExonHash, %hotspotIndelHash, %hotspotHash, %clinicalPanelHash, %samples, @header, @vcfFiles, $clinicalReport);
    open HOTSPOT, "<", $hotspots || die  $logger->fatal("Can not open $hotspots: $!");
    while(<HOTSPOT>){
	chomp;
	my @line = split("\t");
	next if($line[0] eq "Index");
	my ($mut, undef) = split(" ", $line[2]);
	my $gene = $line[1];
	if(/Validation/){
	    $mut =~ s/_indel//;
	    $clinicalPanelHash{$gene.":".$mut} = 1;
	    #print "$gene\t$mut\n";
	}
	if ($mut=~/exon_indel/){
	    $mut=~s/_indel//;
	    $hotspotExonHash{$gene.":".$mut} = $mut;
	}elsif($mut=~/fs/ || $mut =~ /_/ || $mut =~ /del/ || $mut =~/>/){
	    
	    $hotspotIndelHash{$gene.":".$mut} = $mut;
	}else{
	    my ($mut2) = $mut =~/([A-Z]\d+)\>*[A-Za-z]*.*/; # storing only the original amino acid and the position
	    $hotspotHash{$gene.":".$mut2} = $mut2 if($mut2);
	    #print "$gene\t$mut\t$mut2\n";
	}
    }
    close HOTSPOT;

    if(!-e "$outdir/ClinicalReports"){
	eval{`mkdir $outdir/Clinical_Reports`};
	if($@){
	    logger->fatal("Can not create the folder: Clinical_Reports. Reason: $@\n");
	}else{
	    `mkdir $outdir/Clinical_Reports`;
	}
    }

    @vcfFiles = glob("$outdir/VCFs/*.vcf");
    foreach my $vcfFile (@vcfFiles){
	my( %clinicalCoverageHash, %investigationalCoverageHash, %geneTXid, @investigationalMutString, @hotspotMutString);
	my ($sample) = $vcfFile =~ /.*\/(.*)_(Matched|Unmatched)_annovar\.vcf/;
	#print "$sample\n";
	my $covFile = $sample."_exon.coverage.txt";
	open IN, "<", "$outdir/Coverage/$covFile" || die $logger->fatal("Can not open $covFile due to $!");
	while(<IN>){
	    chomp;
	    my ($chr, $gene, $txID, $exon, $ntdRange, $aaRange, $strand, $cov, $threshold, $detected) = split("\t");
	    $geneTXid{$gene} = $txID;
	    if(exists($clinicalPanelHash{$gene.":".$exon})){
		$clinicalCoverageHash{$gene.":".$exon} = $aaRange.":".$cov.":".$threshold.":".$detected;
	    }else{
		$investigationalCoverageHash{$gene.":".$exon} = $aaRange.":".$cov.":".$threshold.":".$detected;
	    }
	}
	close IN;
	
	open VCF, "<", $vcfFile || die $logger->fatal("Can not open $vcfFile: $!");
	while(<VCF>){
	    chomp;
	    next if(/^#/);
	    my @line = split("\t");
	    my $ref = $line[3];
	    my $alt = $line[4];
	    my ($annovar) = $line[8] =~ /.*ANNOVAR\=\((.*)\)/;
	    my ($gene, $txID, $mutType, $exon, $cDNAchange, $aaChange, undef, undef, undef) = split("\Q|", $annovar);
	    if($gene eq "TERT" && $mutType eq "upstream"){
		my $text = "$gene promoter";
		push @investigationalMutString, $text;
		next;
	    }
	    if($mutType eq "splicing"){
		my $text = "$gene $exon splicing ($cDNAchange)";
		push @investigationalMutString, $text;
		next;
	    }
	    if(length($ref) == length($alt) && length($ref) == 1){
		my ($mut) = $aaChange =~ /p\.([A-Z]\d+)[A-Z]+/;
		if(exists($hotspotHash{$gene.":".$mut}) && exists($clinicalPanelHash{$gene.":".$exon})){
		    push @hotspotMutString, "$gene $exon $aaChange ($cDNAchange)";
		}else{
		    push @investigationalMutString, "$gene $exon $aaChange ($cDNAchange)";
		}
	    }else{
		if(exists($hotspotIndelHash{$gene.":".$aaChange}) && exists($hotspotExonHash{$gene.":".$exon})) {
		    push  @hotspotMutString, "$gene $exon $aaChange ($cDNAchange)";
		}else{
		    push @investigationalMutString, "$gene $exon $aaChange ($cDNAchange)";
		}
	    }
	}
	close VCF;

	
	my $clinicalReport = $sample."_clinical_report.txt";
	my ($pat_id) = $sample =~ /(.*)-*T*.*/;
	open OUT, ">", "$outdir/Clinical_Reports/$clinicalReport" || die $logger->fatal("Can't open $clinicalReport beacuse $!");
	print OUT "Patient ID : $pat_id\n";
	
	print OUT "\n";
	print OUT "Test performed: IMPACT\n";
	print OUT "---------------\n";
	print OUT "Next generation sequencing for specific mutations in 340 cancer related genes. Regions covered by the test are available upon request.\n";
	print OUT "\n";
	print OUT "Methodology:\n";
	print OUT "The specific mutations are detected by preparing a library of DNA sequences corresponging to each exon of the 340 genes. The prepared library is then sequenced on a HiSeq 2500 instrument\n";
	if ($vcfFile =~ /Unmatched/){
	    print OUT "The sample did not have a matching normal sample, therefore another normal sample was used in variant calling.\n";
	}else{
	    print OUT "The sample had a matching normal which was consequently used in the variant calling.\n";
	}
	print OUT "\n";
	print OUT "Clinical results\n";
	print OUT "\n";
	print OUT "Coding variants of known significance:\n";
	print OUT "\n";
	if(scalar(@hotspotMutString) == 0){
	    print OUT "No mutations detected.\n";
	}else{
	    foreach my $i (0..$#hotspotMutString){
		my $index = $i + 1;
		print OUT "$index. $hotspotMutString[$i] mutation is DETECTED.\n";
	    }
	}
	
	print OUT "\n";
	print OUT "Investigational results\n";
	print OUT "\n";
	print OUT "Coding variants of uncertain significance:\n";
	print OUT "\n";
	if(scalar (@investigationalMutString) == 0){
	    print OUT "No variants detected.\n";
	}else{
	    foreach my $i (0 .. $#investigationalMutString){
		my $index = $i + 1;
		print OUT "${index}. $investigationalMutString[$i] mutation is DETECTED\n";
	    }
	}
	print OUT "\n";
	
	print OUT "Note: Non-coding variants are available upon request.\n";
	print OUT "\n";
	print OUT "Diagnostic sensitivity: This assay is designed to detect single nucleotide variants only within defined regions. Nucleotide insertions, deletions and gene rearrangements will not be detected\n";
	print OUT "\n";
	print OUT "Technical sensitivity: This assay may not detect certain mutations if the proportion of tumor cells in the sample studied is less than 20%\n";
	print OUT "\n";
	print OUT "Limitations\n";
	print OUT "Mutations outside of the defined regions are not detected. Non-coding variants are available upon request.Variants may not be detected in regions with low coverage. The following regions showed coverage of less than 100X and are therefore interpreted as indeterminate:\n";
	print OUT "\n";
	print OUT "Depth of coverage: <100X=low, >100X=good\n";
	print OUT "\n";
	print OUT "Clinical Panel\n";
	print OUT "Gene\tTranscript\tExon\tAA\tCoverage\tResult\n";
	foreach my $key (sort keys %clinicalCoverageHash){
	    my ($gene, $exon) = split(":", $key);
	    my ($aaRange, $cov, $threshold, $detected) = split(":", $clinicalCoverageHash{$key});
	    print OUT "$gene\t$geneTXid{$gene}\t$exon\t(aa $aaRange)\t$cov\t$detected\n";
	}
	print OUT "\nInvestigational Panel\n";
	my $count = 0;
	foreach my $key (sort my_sort3 keys %investigationalCoverageHash){
	    my ($aaRange, $cov, $threshold, $detected) = split(":", $investigationalCoverageHash{$key});
	    if($cov < 100){
		$count++;
	    }
	}
	if ($count == 0){
	    print OUT "NONE\n";
	}else{
	    print OUT "Gene\tTranscript\tExon\tAA\tCoverage\tResult\n";
	    foreach my $key (sort my_sort3 keys %investigationalCoverageHash){
		my ($gene, $exon) = split(":", $key);
		my ($aaRange, $cov, $threshold, $detected) = split(":", $investigationalCoverageHash{$key});
		if($cov < 100){
		    print OUT "$gene\t$geneTXid{$gene}\t$exon\t(aa $aaRange)\t$cov\t$detected\n";
		}
	    }
	}
    }
}

###########################
# exon sort
sub my_sort3 {
    my($gene_a, $exon_a) = split(":", $a);
    my($gene_b, $exon_b) = split(":", $b);
    my ($ex_a) = $exon_a =~ /exon(\d+)/;
    my ($ex_b) = $exon_b =~ /exon(\d+)/;
    return($gene_a cmp $gene_b || $ex_a <=>$ex_b);
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
    return ($chr_a <=> $chr_b);
}


####################################
# chromosome sort
sub my_sort2{
    my ($chr_a, $gene_a,  undef, $exon1, $strand_a, $aa_range_a) = split(":", $a);
    my ($num1) = $exon1 =~ /exon(\d+)/;
    my ($chr_b, $gene_b,  undef, $exon2, $strand_b, $aa_range_b) = split(":", $b);
    my ($num2) = $exon2 =~ /exon(\d+)/;
    if ($chr_a eq "X") {$chr_a = 23;}
    if ($chr_a eq "Y") {$chr_a = 24;}
    if ($chr_a eq "M") {$chr_a = 25;}
    if ($chr_b eq "X") {$chr_b = 23;}
    if ($chr_b eq "Y") {$chr_b = 24;}
    if ($chr_b eq "M") {$chr_b = 25;}
    return($chr_a <=> $chr_b || $gene_a cmp $gene_b || $num1 <=> $num2);
}


#####################################
# Reverse compliment DNA strand
sub revcomp{
    my $dna = shift @_;
    my $ret = reverse($dna);
    $ret =~ tr/ACGTacgt/TGCAtgca/;
    return($ret);
}

#####################################
#####################################
#Read data related to samples as well as barcodes from title file.

sub ReadTitleFile
{
    my($titleFile,$outdir) = @_;
    my @barcode = ();
    my @pool = ();
    my @sampleId = ();
    my @collabId = ();
    my @patientId = ();
    my @class = ();
    my @sampleType = ();
    my @inputNg = ();
    my @libraryYeild = ();
    my @poolInput = ();
    my @baitVersion = ();
    my @fof = ();
    my @newfof = ();
    tie (my %patientIDPerSampleID, 'Tie::IxHash');
    tie (my %classPerPatientIDSampleID, 'Tie::IxHash');
    open(TFH,$titleFile)||die $logger->fatal("Can not open $titleFile, $!");
    while(<TFH>)
    {
	next if($. == 1);
	my @dataCols = split("\t",$_);
	my @newDatacols = grep(s/\s*$//g, @dataCols);#remove whitespace if any
	push(@barcode,$newDatacols[0]);
	push(@pool,$newDatacols[1]);
	push(@sampleId,$newDatacols[2]);
	push(@collabId,$newDatacols[3]);
	push(@patientId,$newDatacols[4]);
	push(@class,$newDatacols[5]);
	push(@sampleType,$newDatacols[6]);
	push(@inputNg,$newDatacols[7]);
	push(@libraryYeild,$newDatacols[8]);
	push(@poolInput,$newDatacols[9]);
	push(@baitVersion,$newDatacols[10]);
    }
    close(TFH);
    for(my $i = 0; $i < scalar(@barcode); $i++)
    {
	#print "$[$i] => $class[$i]\n";
	$classPerPatientIDSampleID{$sampleId[$i] . ":" . $patientId[$i]} = $class[$i];
	$patientIDPerSampleID{$sampleId[$i]} = $patientId[$i];
    }
    return(\%patientIDPerSampleID,\%classPerPatientIDSampleID);

}
###################################################
###################################################
#--Delete Unwanted Files
sub HouseKeeping
{
    my($outdir,@list) = @_;
    $logger->info("Removing Unwanted files...");
    foreach my  $file (@list)
    {
	if($file =~ /\//)
	{
	    `rm $file`;
	}
	else
	{
	    `rm $outdir/$file`;
	}
    }
	return;
}




#################################
# Get Config file and program locations - adopted from Ronak's pipeline
# Input config file


sub GetConfiguration
{
    my($config_file) = @_;
    my @data = ();
    tie (my %config,'Tie::IxHash');
    tie (my %location,'Tie::IxHash');
    tie (my %version,'Tie::IxHash');
    tie (my %parameters,'Tie::IxHash');
    $logger->info("Reading the configuration file"); 
    # change the default Input Record Separator.
    $/= ">";
    open(CONFIG,"$config_file") or die($logger->fatal("Cannot open $config_file. Error: $!"));
    my $junk = <CONFIG>;
    while(<CONFIG>)
    {
	next if ($_ =~ /^#/);
	chomp($_);
	my ($defLine, @configLines) = split /\n/, $_;
	if($defLine =~ /Locations/)
	{
	    foreach my $config (@configLines)
	    {
		next if($config =~ /^#/);
		@data = split("=",$config,2);
		$data[0] =~ s/\s//g;
		$data[1] =~ s/\s//g;
		$location{$data[0]} = $data[1];
	    }
	}
	if($defLine =~ /Parameters/)
	{
	    foreach my $config (@configLines)
	    {
		next if($config =~ /^#/);
		@data = split("=",$config,2);
		$data[0] =~ s/\s//g;
		$data[1] =~ s/\s//g;
		$parameters{$data[0]} = $data[1];
	    }
	}
	if($defLine =~ /Versions/)
	{
	    foreach my $config (@configLines)
	    {
		next if($config =~ /^#/);
		@data = split("=",$config,2);
		$data[0] =~ s/\s//g;
		$data[1] =~ s/\s//g;
		$version{$data[0]} = $data[1];
	    }
	}
    }
    close(CONFIG);
    $logger->info("Completed reading the configuration file");
    # Change back the Input Record Separator to default.
    $/= "\n";
    return(\%location, \%version, \%parameters);
}
