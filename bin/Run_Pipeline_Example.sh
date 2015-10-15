##Run_Pipeline_Example.sh
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

##Run both IMPACT-Pipeline & SV Process on LSF
echo bsub -q sol -cwd ${WorkingDir} -J ${ProjectName} -e${ProjectName}.stderr -o ${ProjectName}.stdout -We 24:00 -R "rusage[mem=2]" -M 4 \"${Perl} ${PipelineScript} -c ${CONFIGFILE} -sc {$SV_ConfigFile} -d ${DATADIR} -o ${OUTDIR}\"
bsub -q sol -cwd ${WorkingDir} -J ${PoolName} -e${ProjectName}.stderr -o ${ProjectName}.stdout -We 24:00 -R "rusage[mem=2]" -M 4 \"${Perl} ${PipelineScript} -c ${CONFIGFILE} -sc {$SV_ConfigFile} -d ${DATADIR} -o ${OUTDIR}\"
##Run IMPACT-Pipeline on LSF
echo bsub -q sol -cwd ${WorkingDir} -J ${ProjectName} -e${ProjectName}.stderr -o ${ProjectName}.stdout -We 24:00 -R "rusage[mem=2]" -M 4 \"${Perl} ${PipelineScript} -c ${CONFIGFILE} -d ${DATADIR} -o ${OUTDIR}\"
bsub -q sol -cwd ${WorkingDir} -J ${PoolName} -e${ProjectName}.stderr -o ${ProjectName}.stdout -We 24:00 -R "rusage[mem=2]" -M 4 \"${Perl} ${PipelineScript} -c ${CONFIGFILE} -d ${DATADIR} -o ${OUTDIR}\"

##Run both IMPACT-Pipeline & SV Process on SGE
echo qsub -q test.q -wd ${WorkingDir} -N ${ProjectName} -l hvmem=2G,virtual_free=2G -pe smp 1 -b y \"${Perl} ${PipelineScript} -c ${CONFIGFILE} -sc ${SV_ConfigFile} -d ${DATADIR} -o ${OUTDIR}\"
qsub -q test.q -wd ${WorkingDir} -N ${ProjectName} -l hvmem=2G,virtual_free=2G -pe smp 1 -b y \"${Perl} ${PipelineScript} -c ${CONFIGFILE} -sc ${SV_ConfigFile} -d ${DATADIR} -o ${OUTDIR}\"
##Run both IMPACT-Pipeline on SGE
echo qsub -q test.q -wd ${WorkingDir} -N ${ProjectName} -l hvmem=2G,virtual_free=2G -pe smp 1 -b y \"${Perl} ${PipelineScript} -c ${CONFIGFILE} -d ${DATADIR} -o ${OUTDIR}\"
qsub -q test.q -wd ${WorkingDir} -N ${ProjectName} -l hvmem=2G,virtual_free=2G -pe smp 1 -b y \"${Perl} ${PipelineScript} -c ${CONFIGFILE} -d ${DATADIR} -o ${OUTDIR}\"
