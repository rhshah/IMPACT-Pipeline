=====
Usage
=====

Quick Usage
===========
**RunIlluminaProcess.pl -c configuration.txt -sc configuration_sc.txt -d /path/to/fastq/files -o /path/to/output/directory [options]**
	
	--config | -c                        S Path to configration file(required)
	
	--svConfig | -sc                     S Path to structural variant configration file(optional)
	
	--symLinkFlag | -sf           	   I Flag for Keeping or removing the symolic links(1:Remove;2:Keep)(default:2)
	
	--dataDirectory | -d                 S Path where all the files to be processed are located (required)
	
	--outputDirectory | -o               S Path where all the output files will be written (required)
	
Assuming you have setup the configuration file properly and you have SampleSheet.csv and title_file.txt in the **dataDirectory**you can run:
``perl RunIlluminaProcess.pl -c configuration.txt -sc configuration_sc.txt -d /path/to/fastq/files -o /path/to/output/directory``

Detailed Usage
==============

The behaviour of the prgram depends on the inputs in the configuration file:

In the configuration file the **Process** variable in section **>Parameters** tells pipeline following:

+---------+-----------------------------------------------------------------------------------------------+
| Process | Things Pipeline will do                                                                       |
+=========+===============================================================================================+
| 1       | Merge Fastq 																				  |
+---------+-----------------------------------------------------------------------------------------------+
| 2       | Trimming, Mapping & sorting of SAM file giving you a BAM file								  |
+---------+-----------------------------------------------------------------------------------------------+
| 3       | Mark Duplicates, Indel Reaglinment, Base Quality Recalibration 								  |
+---------+-----------------------------------------------------------------------------------------------+
| 4       | Metrics Calculation, QC Report Genaration and launching IMPACT-SV if given -sc flag specified |
+---------+-----------------------------------------------------------------------------------------------+
| 5       | Variant Calling 																		      |
+---------+-----------------------------------------------------------------------------------------------+
| 6       | Variant Filtering and Genotyping 															  |
+---------+-----------------------------------------------------------------------------------------------+
| 7       | Variant Annotation and Variant Filtering 													  |
+---------+-----------------------------------------------------------------------------------------------+


1. To run the complete pipeline. Set the following in the configuration file:
	
	:Process: 1,2,3,4,5,6,7

2. To run from **Process 1**. Set the following in the configuration file:
	
	:ListOfFile: fastq.list #where fastq.list contains all the fastq files to be proces, this needs to be an even number as it automatically pairs them.
	:Process: 2,3,4,5,6,7
	
3. To run from the **Process 3 to 7**. Set the following in the configuration file:
	
	:ListOfFile: SortedBam.list #where SortedBam.list contains all the sorted bam files from Process 2 to be processed
	:Process: 3,4,5,6,7
	
4. To run from the **Process 4 to 7**. Set the following in the configuration file:
	:ListOfFile: RecalibratedBam.list #where Recalibrated.list contains all the recalibrated bam files from Process 3 to be processed
	:Process: 4,5,6,7

5. To run from the **Process 5 to 7**. Set the following in the configuration file:
	:ListOfFile: RecalibratedBam.list #where Recalibrated.list contains all the recalibrated bam files from Process 3 to be processed
	:Process: 5,6,7
	
	**Note:** For this to be sucessfull you should hve the files from step 4 in the **outputDirectory** 
	
6. To run from the **Process 6 to 7**. Set the following in the configuration file:
	:ListOfFile: RecalibratedBam.list #where Recalibrated.list contains all the recalibrated bam files from Process 3 to be processed
	:Process: 6,7
	
	**Note:** For this to be sucessfull you should hve the files from step 5 in the **outputDirectory**
	
7.  To run from the **Process 7**. Set the following in the configuration file:
	:ListOfFile: RecalibratedBam.list #where Recalibrated.list contains all the recalibrated bam files from Process 3 to be processed
	:Process: 7
	
	**Note:** For this to be sucessfull you should hve the files from step 6 in the **outputDirectory**
	
If you want to run each Process separetly that is also possible but you need to make sure that files from previous procss are present in the **outputDirectory**