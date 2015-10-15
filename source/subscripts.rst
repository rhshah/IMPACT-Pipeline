==========================
Description of sub-scripts
==========================

--------------------------
Script in the bin folder
--------------------------

Compiling QC metrics, Generating QC-Report and running copynumber analysis
==========================================================================

**dmp_compile_qc_metrics.pl [options]**

        --bamList | -i       S File of files having list of all bam files (required)
		
        --titleFile | -t     S tab-delimited title file for the samples (required and submit with full path)
		
        --AllMetrics | -am   S Path to AllMetrics script (required)
		
        --LoessNorm | -ln    S Path to Loess Normalization Script  (required)
		
        --BestCN | -cn       S Path to Best Copy Number Script (required)
		
        --GCBias | -gcb      S Path to GC bias file (required)
		
        --HistNorm | -his    S Path to Directory with all historical normal files (required)
		
        --queue | -q         S Name of the Sun Grd Engine Queue where the pipeline needs to run (required)
		
        --qsub            S Path to qsub executable for SGE(default:None,optional)
		
        --bsub            S Path to bsub executable for LSF(default:None,required)
		
        --metricsScript | -ms         S Name of the script used to generate .html and .pdf files
		
        --outdir | -o        S Path where all the output files will be written (optional) [default:cwd]
		
		
Genotyping Variants across multiple samples
===========================================

**dmp_genotype_allele.pl [options]**

        --FilteredMutationVcfFile | -fmv     S vcf file describing details about the mutations (required)
		
        --BamFile | -bam                     S bam file to be used for genotyping (required)
		
        --RefFile | -rf                      S Path to genome reference file (required)
		
        --samtools | -s                      S Path to samtools (required)
		
        --bedtools | -b                      S Path to bedtools (required)
		
        --MinBaseQualit | -mbq              I Min. Base Quality Threshold (optional;default:5)
		
        --MinMappingQuality | -mmq           I Min. Mapping Quality Threshold (optional;default:5)
		
        --deleteUnwantedFiles | -d           I 2=>To delete files 1=> To keep files (default:2,optional)
		
        --outdir | -o                        S Path where all the output files will be written (optional;default:current working directory)
		
        --outFile | -of                      S Name of the allele depth output file (optional;default:BamFame-.bam+_mpileup.alleledepth)
		
        --bamId | -bi                        S Bam Id to be used (optional;default:bamfile name)
		
        --queue | -q                         S Name of the SGE / LSF Queue where the pipeline needs to run (required)
		
        --qsub                            S Path to qsub executable for SGE(default:None,optional)
		
        --bsub                            S Path to bsub executable for LSF(default:None,required)
		
        --mpileUpOutFile | -mof              S Name of samtools mpileup output file (optional;default:BamFile-.bam+.mpileup)
		
        --typeOfSample | -tos                S Type of Sample (optional;default:Tumor;canbe Tumor or Normal)


Running Indel Realignment using ABRA
====================================

**usage: Run_AbraRealignment.py [options]**

Run ABRA Indel Realignment

arguments:
  -h, --help            show this help message and exit
  -i BamFile.list, --bamList BamFile.list
                        Full path to the tumor bam files as a fof.
  -p PatientID, --patientId PatientID
                        Id of the Patient for which the bam files are to be
                        realigned
  -v, --verbose         make lots of noise [default]
  -t 5, --threads 5     Number of Threads to be used to run ABRA
  -d, --mdp             Threshold for downsampling depth to run ABRA
  -k [43 [43 ...]], --kmers [43 [43 ...]]
                        Number of k-mers to be used to run ABRA; Multiple
                        k-mers are separated by space
  -temp /somepath/tmpdir, --temporaryDirectory /somepath/tmpdir
                        Full Path to temporary directory
  -r /somepath/Homo_Sapeins_hg19.fasta, --referenceFile /somepath/Homo_Sapeins_hg19.fasta
                        Full Path to the reference file with the bwa index.
  -a /somepath/ABRA.jar, --abraJar /somepath/ABRA.jar
                        Full Path to the ABRA jar file.
  -tr /somepath/targetRegion.bed, --targetRegion /somepath/targetRegion.bed
                        Full Path to the target region bed file
  -j /somepath/java, --javaPATH /somepath/java
                        Path to java executable.
  -b /somepath/bin, --bwaPATH /somepath/bin
                        Path to the bin of bwa executable.
  -q all.q or clin.q, --queue all.q or clin.q
                        Name of the SGE queue
  -o /somepath/output, --outDir /somepath/output
                        Full Path to the output dir.
  -qsub /somepath/qsub, --qsubPath /somepath/qsub
                        Full Path to the qsub executables of SGE.
  -bsub /somepath/bsub, --bsubPath /somepath/bsub
                        Full Path to the bsub executables of LSF.


Call Indels > 25 bp using Pindel
================================

**usage: Run_Pindel.py [options]**

Run Pindel for Long Indels & MNPS (32bp-350bp)

optional arguments:
  -h, --help            show this help message and exit
  -i pindel.conf, --pindelConfig pindel.conf
                        Full path to the pindel configuration
  -pId PatientID, --patientId PatientID
                        Id of the Patient for which the bam files are to be
                        realigned
  -v, --verbose         make lots of noise [default]
  -t 5, --threads 5     Number of Threads to be used to run Pindel
  -r /somepath/Homo_Sapeins_hg19.fasta, --referenceFile /somepath/Homo_Sapeins_hg19.fasta
                        Full Path to the reference file with the bwa index.
  -p /somepath/pindel/bin, --pindelDir /somepath/pindel/bin
                        Full Path to the Pindel executables.
  -chr ALL, --chromosomes ALL
                        Which chr/fragment. Pindel will process reads for one
                        chromosome each time. ChrName must be the same as in
                        reference sequence and in read file.
  -q all.q or clin.q, --queue all.q or clin.q
                        Name of the SGE queue
  -o /somepath/output, --outDir /somepath/output
                        Full Path to the output dir.
  -op TumorID, --outPrefix TumorID
                        Id of the Tumor bam file which will be used as the
                        prefix for Pindel output files
  -qsub /somepath/qsub, --qsubPath /somepath/qsub
                        Full Path to the qsub executables of SGE.
  -bsub /somepath/bsub, --bsubPath /somepath/bsub
                        Full Path to the bsub executables of LSF.


------------------------------------
Script in the support-scripts folder
------------------------------------

Calculate Intervals from BAM file that have some minimum coverage
=================================================================

**usage: Run_FindCoveredInterval.py [options]**

This will run find covered interval program from GATK.

optional arguments:
  -h, --help            show this help message and exit
  -i BamFile.list, --bamList BamFile.list
                        Full path to the tumor bam files as a fof.
  -of OutFilePrefix, --outFilePrefix OutFilePrefix
                        Output Covered Interval File Prefix for the bam files.
  -v, --verbose         make lots of noise [default]
  -t 5, --threads 5     Number of Threads to be used to run
                        FindCoveredIntervals
  -dp 20, --totaldepth 20
                        Total depth threshold
  -mbq 20, --minbasequality 20
                        Threshold for minimum base quality for Running Find
                        Covered Interval
  -mmq 20, --minmappingquality 20
                        Threshold for minimum mapping quality for Running Find
                        Covered Interval
  -r /somepath/Homo_Sapeins_hg19.fasta, --referenceFile /somepath/Homo_Sapeins_hg19.fasta
                        Full Path to the reference file with the bwa index.
  -g /somepath/GenomeAnalysisTK.jar, --gatkJar /somepath/GenomeAnalysisTK.jar
                        Full Path to the GATK jar file.
  -j /somepath/java, --javaPATH /somepath/java
                        Path to java executable.
  -q all.q or clin.q, --queue all.q or clin.q
                        Name of the SGE queue
  -o /somepath/output, --outDir /somepath/output
                        Full Path to the output dir.
  -qsub /somepath/qsub, --qsubPath /somepath/qsub
                        Full Path to the qsub executables of SGE.
  -bsub /somepath/bsub, --bsubPath /somepath/bsub
                        Full Path to the bsub executables of LSF.

Annotating Variants All Merged Variant File
============================================

**dmp_annotate_variants.pl [options]**

        --SomaticMutIndelFile | -si          S File containing mutations (required and submit with full path,Ex:/SomePath/Some_SomaticMutIndel.txt)
		
        --ConfigurationFile | -c             S Configuration file that contains the locations for the programs and the databases (required and submit with full path)
		
        --titleFile | -t                     S tab-delimited title file for the samples (required and submit with full path)
		
        --outdir | -o                        S Path where all the output files will be written (optional;default:cwd)
		
		--exonCoverageFile | -ec             S Path where the all exon coverage file is located (full path]
		
		--geneCoverageFile | -gc             S Path where the gene coverage file is located (full path)]
		
        --deleteUnwantedFiles | -d           I 2=>To delete files 1=> To keep files (default:2,optional)


Filter Variants after annotation
=================================

**dmp_filter_genotyped_variants.pl [options]**

    --input | -i                                            S  File containing mutations with genotype information (required)
	
    --hotspots | -h										    S  File containing the list of hotspots (required)
	
    --clinicalExons | -ce								    S  File containing the list of clinical exons (required)
	
    --titleFile | -t										S  Title file (required)
	
    --minimumDPforSNVs | -dp_snv                            I  Minumum accepted DP for novel SNVs (default: 20)
	
    --minimumADforSNVs | -ad_snv                            I  Minimum accepted AD for novel SNVs (default: 10)
	
    --minimumVFforSNVs | -vf_snv                            F  Minimum accepted VF for novel SNVs (default: 0.05)
	
    --minimumDPforSNVhotspot | -dp_snvHS                    I  Minumum accepted DP for Hotspot SNVs (default: 20)
	
    --minimumADforSNVhotspot | -ad_snvHS                    I  Minimum accepted AD for Hotspot SNVs (default: 8)
	
    --minimumVFforSNVhotspot | -vf_snvHS                    F  Minimum accepted VF for Hotspot SNVs (default: 0.02)
	
    --minimumDPforINDELs | -dp_indel                        I  Minumum accepted DP for novel INDELs (default: 20)
	
    --minimumADforINDELs | -ad_indel                        I  Minimum accepted AD for novel INDELs (default: 10)
    
	--minimumVFforINDELs | -vf_indel                        I  Minimum accepted VF for novel INDELs (default: 0.05)
    
	--minimumDPforINDELhotspot | -dp_indelHS                F  Minumum accepted DP for Hotspot INDELs (default: 20)
    
	--minimumADforINDELhotspot | -ad_indelHS                I  Minimum accepted AD for Hotspot INDELs (default: 8)
    
	--minimumVFforINDELhotspot | -vf_indelHS                F  Minimum accepted VF for Hotspot INDELs (default: 0.02)
    
	--minimumOccurrencePercent | -occurrence                S  Minimum accepted value of occurrence in other normals, in percent (default: 20)
    
	--TNfreqRatioThreshold | -tn_ratio					    S  Minimum value for VFt/VFn value (default: 5)
    
	--MAFthreshold | -mt                 					F Minimum accepted MAF values for unmatched variant calls (default : 0.01)
    
