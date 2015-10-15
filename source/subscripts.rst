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


**Running Indel Realignment using ABRA**
========================================

usage: Run_AbraRealignment.py [options]

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
