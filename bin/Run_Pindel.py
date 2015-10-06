'''
Created on 07/31/2014
@Ronak Shah

'''

import argparse
import sys
import time
import os.path
import stat
from subprocess import Popen
import shlex
import shutil
from datetime import date

def main():
   parser = argparse.ArgumentParser(prog='Run_Pindel.py', description='Run Pindel for Long Indels & MNPS (32bp-350bp)', usage='%(prog)s [options]')
   parser.add_argument("-i", "--pindelConfig", action="store", dest="config", required=True, metavar='pindel.conf', help="Full path to the pindel configuration") 
   parser.add_argument("-pId", "--patientId", action="store", dest="patientId", required=True, metavar='PatientID', help="Id of the Patient for which the bam files are to be realigned")
   parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=True, help="make lots of noise [default]")
   parser.add_argument("-t", "--threads", action="store", dest="threads", required=True, metavar='5', help="Number of Threads to be used to run Pindel")
   parser.add_argument("-r", "--referenceFile", action="store", dest="ref", required=True, metavar='/somepath/Homo_Sapeins_hg19.fasta', help="Full Path to the reference file with the bwa index.")
   parser.add_argument("-p", "--pindelDir", action="store", dest="PINDEL", required=True, metavar='/somepath/pindel/bin', help="Full Path to the Pindel executables.")
   parser.add_argument("-chr", "--chromosomes", action="store", dest="chr", required=True, metavar='ALL', help="Which chr/fragment. Pindel will process reads for one chromosome each time. ChrName must be the same as in reference sequence and in read file.")
   parser.add_argument("-q", "--queue", action="store", dest="queue", required=False, metavar='all.q or clin.q', help="Name of the SGE queue")
   parser.add_argument("-o", "--outDir", action="store", dest="outdir", required=True, metavar='/somepath/output', help="Full Path to the output dir.")
   parser.add_argument("-op", "--outPrefix", action="store", dest="outprefix", required=True, metavar='TumorID', help="Id of the Tumor bam file which will be used as the prefix for Pindel output files")
   parser.add_argument("-qsub", "--qsubPath", action="store", dest="qsub", required=False, metavar='/somepath/qsub', help="Full Path to the qsub executables of SGE.")
   parser.add_argument("-bsub", "--bsubPath", action="store", dest="bsub", required=False, metavar='/somepath/bsub', help="Full Path to the bsub executables of LSF.")
   # parser.add_argument("-tr", "--targetRegion", action="store", dest="targetRegion", required=True, metavar='/somepath/targetRegion.bed', help="Full Path to bedfile for target region.")
   
   args = parser.parse_args()
   if(args.verbose):
       print "I have Started the run for Pindel."
   (wd, outvcf, tag) = ProcessArgs(args)
   RunPindel(args, wd, outvcf, tag)
   if(args.verbose):
       print "I have finished the run for Pindel."   
   
def ProcessArgs(args):
    if(args.verbose):
        print "I am currently processing the arguments.\n"
    #Check qsub and bsub
    if(args.qsub and args.bsub):
       print "Please give either qsub or bsub arguments. Script does not allow usage of both\n"
       sys.exit(1)           
    if((not args.qsub) and (not args.bsub)):
       print "Please give either qsub or bsub arguments. None are provided\n"
       sys.exit(1)      
    tumorBam = '' 
    with open(args.config, 'r') as filecontent:
        for line in filecontent:
            if(args.patientId not in line):
                continue
            else:
                data = line.rstrip('\n').split('\t')
                tumorBam = os.path.basename(data[0]).rstrip('.bam')
                break
    outvcf = tumorBam + ".pindel.detailed.TN.matched.vcf"
    SampleDirName = args.patientId
    staticDir = "PindelAnalysis"
    AnalysisDir = os.path.join(args.outdir, staticDir)
    SampleAnalysisDir = os.path.join(AnalysisDir, SampleDirName)
    try:
        os.mkdir(AnalysisDir)
    except OSError:
        if(args.verbose):
            print "Dir:", AnalysisDir, " exists thus we wont be making it\n"
    if os.path.isdir(SampleAnalysisDir):
            if(args.verbose):
                print "Dir:", SampleAnalysisDir, " exists and we wont run the analysis\n"
            pindeltag = 0
    else:
        os.mkdir(SampleAnalysisDir)
        pindeltag = 1
    if(args.verbose):
        print "I am done processing the arguments.\n"   
    
    return(SampleAnalysisDir, outvcf, pindeltag)

def RunPindel(args, wd, vcfoutName, tag):
    myPid = os.getpid()
    day = date.today()
    today = day.isoformat()
    today = today.replace("-", "")
    pindel = os.path.join(args.PINDEL, "pindel")
    pindel2vcf = os.path.join(args.PINDEL, "pindel2vcf")
    vcfOutPath = os.path.join(args.outdir, vcfoutName)
    cl_cmd = ''
    mem = int(args.threads) * 6
    maxmem = int(mem)+6
    # myPid = str(myPid)
    if(args.verbose):
        print "I am running Pindel for ", args.patientId, " using SGE/LSF"
    
    # Setting Job for SGE   
    if(tag):
        cmd = pindel + " -i " + args.config + " -f " + args.ref + " -c " + args.chr + " -o " + args.outprefix + " -r false -t false -I false"
        # print "CMD==>",cmd,"\n"
        if(args.qsub):
            cl_cmd = args.qsub + " -q " + args.queue + " -N " + "Pindel_" + args.patientId + "_" + str(myPid) + " -o " + "Pindel_" + args.patientId + "_" + str(myPid) + ".stdout" + " -e " + "Pindel_" + args.patientId + "_" + str(myPid) + ".stderr" + " -V -l h_vmem=6G,virtual_free=6G -pe smp " + args.threads + " -wd " + wd + " -sync y " + " -b y " + cmd 
        else:
            cl_cmd = args.bsub + " -q " + args.queue + " -J " + "Pindel_" + args.patientId + "_" + str(myPid) + " -o " + "Pindel_" + args.patientId + "_" + str(myPid) + ".stdout" + " -e " + "Pindel_" + args.patientId + "_" + str(myPid) + ".stderr" + " -We 24:00 -R \"rusage[mem=" + str(mem) + "]\" -M " + str(maxmem) + " -n " + args.threads + " -cwd " + wd + " -K " + cmd 
        print "Cluster_CMD==>", cl_cmd , "\n"
        cl_args = shlex.split(cl_cmd)
        proc = Popen(cl_args)
        proc.wait()
        retcode = proc.returncode
    else:
        retcode = 1
    if(retcode >= 0):
        if(args.verbose):
            print "I have finished running Pindel for ", args.patientId, " using SGE/LSF"
        if(os.path.isfile(vcfOutPath)):
            retcode = 1
        else:
            p2v_cmd = pindel2vcf + " --pindel_output_root " + args.outprefix + " --reference " + args.ref + " --reference_name hg19 --reference_date " + today + " --vcf " + vcfOutPath + " -b true --gatk_compatible"         
            if(args.qsub):
                p2v_cl = args.qsub + " -q " + args.queue + " -N " + "Pindel2Vcf_" + args.patientId + "_" + str(myPid) + " -o " + "Pindel2Vcf_" + args.patientId + "_" + str(myPid) + ".stdout" + " -e " + "Pindel2Vcf_" + args.patientId + "_" + str(myPid) + ".stderr" + " -V -l h_vmem=8G,virtual_free=8G -pe smp 1 " + " -wd " + wd + " -sync y " + " -b y " + p2v_cmd 
            else:
                p2v_cl = args.bsub + " -q " + args.queue + " -J " + "Pindel2Vcf_" + args.patientId + "_" + str(myPid) + " -o " + "Pindel2Vcf_" + args.patientId + "_" + str(myPid) + ".stdout" + " -e " + "Pindel2Vcf_" + args.patientId + "_" + str(myPid) + ".stderr" +  " -We 24:00 -R \"rusage[mem=8]\" -M 16 -n 1" + " -cwd " + wd + " -K " + p2v_cmd 
           
            print "P2V cmd: ", p2v_cl
            p2v_args = shlex.split(p2v_cl)
            proc = Popen(p2v_args)
            proc.wait()
            retcode = proc.returncode
        if(retcode >= 0):
            if(args.verbose):
                print "Converted Pindel Output to VCF\n"
        else:
            if(args.verbose):
                    print "Pindel2Vcf is either still running or it errored out with return code", retcode, "\n"        
    else:
        if(args.verbose):
            print "Pindel is either still running or it errored out with return code", retcode, "\n"    
               
    
if __name__ == "__main__":
    start_time = time.time()  
    main()
    end_time = time.time()
    print("Elapsed time was %g seconds" % (end_time - start_time))      
    
