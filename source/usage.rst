=====
Usage
=====

Quick Usage
===========
**RunIlluminaProcess.pl [options]**
	
	--config | -c                        S Path to configration file(required)
	
	--svConfig | -sc                     S Path to structural variant configration file(optional)
	
	--symLinkFlag | -sf           	   I Flag for Keeping or removing the symolic links(1:Remove;2:Keep)(default:2)
	
	--dataDirectory | -d                 S Path where all the files to be processed are located (required)
	
	--outputDirectory | -o               S Path where all the output files will be written (required)
	
Assuming you have setup the configuration file properly and you have SampleSheet.csv and title_file.txt in the **dataDirectory** you can run:

``nohup perl RunIlluminaProcess.pl -c configuration.txt -sc configuration_sc.txt -d /path/to/fastq/files -o /path/to/output/directory``

Detailed Usage
==============

The behaviour of the prgram depends on the inputs in the configuration file:

In the configuration file the **Process** variable in section **>Parameters** tells pipeline following:

What does each number represent
-------------------------------

+---------+-----------------------------------------------------------------------------------------------+
| Process | Things Pipeline will do                                                                       |
+=========+===============================================================================================+
| 1       | Merge Fastq 										  |										  
+---------+-----------------------------------------------------------------------------------------------+
| 2       | Trimming, Mapping & sorting of SAM file giving you a BAM file				  |				  
+---------+-----------------------------------------------------------------------------------------------+
| 3       | Mark Duplicates, Indel Realignment, Base Quality Recalibration 				  |
+---------+-----------------------------------------------------------------------------------------------+
| 4       | Metrics Calculation, QC Report Genaration and launching IMPACT-SV if given -sc flag specified |
+---------+-----------------------------------------------------------------------------------------------+
| 5       | Variant Calling 										  |
+---------+-----------------------------------------------------------------------------------------------+
| 6       | Variant Filtering and Genotyping 						                  |
+---------+-----------------------------------------------------------------------------------------------+
| 7       | Variant Annotation and Variant Filtering 							  |
+---------+-----------------------------------------------------------------------------------------------+


Using different Process to run Pipeline
---------------------------------------

1. To run the complete pipeline. Set the following in the configuration file:
	
	:Process: 1,2,3,4,5,6,7

2. To run from **Process 1**. Set the following in the configuration file:
	
	:ListOfFilesListOfFiles: fastq.list #where fastq.list contains all the fastq files to be proces, this needs to be an even number as it automatically pairs them.
	:Process: 2,3,4,5,6,7
	
3. To run from the **Process 3 to 7**. Set the following in the configuration file:
	
	:ListOfFiles: SortedBam.list (where SortedBam.list contains all the sorted bam files from Process 2 to be processed)
	:Process: 3,4,5,6,7
	
4. To run from the **Process 4 to 7**. Set the following in the configuration file:
	:ListOfFiles: RecalibratedBam.list (where Recalibrated.list contains all the recalibrated bam files from Process 3 to be processed)
	:Process: 4,5,6,7

5. To run from the **Process 5 to 7**. Set the following in the configuration file:
	:ListOfFiles: RecalibratedBam.list (where Recalibrated.list contains all the recalibrated bam files from Process 3 to be processed)
	:Process: 5,6,7
	
	**Note:** For this to be sucessfull you should hve the files from step 4 in the **outputDirectory**
	
6. To run from the **Process 6 to 7**. Set the following in the configuration file:
	:ListOfFiles: RecalibratedBam.list (where Recalibrated.list contains all the recalibrated bam files from Process 3 to be processed)
	:Process: 6,7
	
	**Note:** For this to be sucessfull you should hve the files from step 5 in the **outputDirectory**
	
7.  To run from the **Process 7**. Set the following in the configuration file:
	:ListOfFiles: RecalibratedBam.list (where Recalibrated.list contains all the recalibrated bam files from Process 3 to be processed)
	:Process: 7
	
	**Note:** For this to be sucessfull you should hve the files from step 6 in the **outputDirectory**
	
If you want to run each Process separetly that is also possible but you need to make sure that files from previous procss are present in the **outputDirectory**

Shell Script to run pipeline
----------------------------
There is also a **helper shell script (Run_Pipeline.sh)** in the bin directory which will help to run the framework.
Which looks like this:

.. code-block:: sh

	##Run_Pipeline.sh
	#author:Ronak H Shah
	#v1.0.1
	##Path where the fastq are stored
	export DATADIR=<Path for Data Directory>
	##Path where the output shoud be written
	export OUTDIR=<Path To Output Directory>
	##Path to the IMPACT-Pipeline script
	export PipelineScript=<Path to IMPACT-Pipeline Script>
	##Path to Perl installation
	export Perl=<Path to Perl>
	##Project associated with the Run
	export ProjectName=<ProjectName>
	##Path to working directory where you will write the LSF/SGE outputs
	export WorkingDir=<Path to write sge/lsf files>
	##Path to configfile for running main IMPACT pipeline
	export CONFIGFILE=<Path To Pipeline Configuration File>
	##Path to structural variants pipeline configuration file
	export SV_ConfigFile=<Path to SV detection configuration file>

	##Run both IMPACT & SV Process on Luna LSF
	echo bsub -q sol -cwd ${WorkingDir} -J ${ProjectName} -e${ProjectName}.stderr -o ${ProjectName}.stdout -We 24:00 -R "rusage[mem=2]" -M 4 \"${Perl} ${PipelineScript} -c ${CONFIGFILE} -sc {$SV_ConfigFile} -d ${DATADIR} -o ${OUTDIR}\"

	bsub -q sol -cwd ${WorkingDir} -J ${PoolName} -e${ProjectName}.stderr -o ${ProjectName}.stdout -We 24:00 -R "rusage[mem=2]" -M 4 \"${Perl} ${PipelineScript} -c ${CONFIGFILE} -sc {$SV_ConfigFile} -d ${DATADIR} -o ${OUTDIR}\"

	##Run both IMPACT & SV Process on SGE
	echo qsub -q test.q -wd ${WorkingDir} -N ${ProjectName} -l hvmem=2G,virtual_free=2G -pe smp 1 -b y \"${Perl} ${PipelineScript} -c ${CONFIGFILE} -sc ${SV_ConfigFile} -d ${DATADIR} -o ${OUTDIR}\"
	qsub -q test.q -wd ${WorkingDir} -N ${ProjectName} -l hvmem=2G,virtual_free=2G -pe smp 1 -b y \"${Perl} ${PipelineScript} -c ${CONFIGFILE} -sc ${SV_ConfigFile} -d ${DATADIR} -o ${OUTDIR}\"
	
**Note** Comment out the commands according to the cluster type.
