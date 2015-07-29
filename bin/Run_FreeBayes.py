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
   parser = argparse.ArgumentParser(prog='Run_FreeBayes.py', description='Run FreeBayes for Long Indels & MNPS (32bp-350bp)', usage='%(prog)s [options]')
   parser.add_argument("-pId", "--patientId", action="store", dest="patientId", required=True, metavar='PatientID', help="Id of the Patient for which the bam files are to be realigned")
   parser.add_argument("-tbam", "--tumorBAM", action="store", dest="tbam", required=True, metavar='tbam', help="Full Path to tumor bam file")
   parser.add_argument("-nbam", "--normalBAM", action="store", dest="nbam", required=True, metavar='nbam', help="Full Path to normal bam file")
   parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=True, help="make lots of noise [default]")
   parser.add_argument("-t", "--threads", action="store", dest="threads", required=True, metavar='1', help="Number of Threads to be used to run freebayes")
   parser.add_argument("-r", "--referenceFile", action="store", dest="ref", required=True, metavar='/somepath/Homo_Sapeins_hg19.fasta', help="Full Path to the reference file with the bwa index.")
   parser.add_argument("-floc", "--freebayes", action="store", dest="FREEBAYES", required=True, metavar='/somepath/bin/freebayes', help="Full Path to the freebayes executables.")
   parser.add_argument("-q", "--queue", action="store", dest="queue", required=False, metavar='all.q or clin.q', help="Name of the SGE queue")
   parser.add_argument("-o", "--outDir", action="store", dest="outdir", required=True, metavar='/somepath/output', help="Full Path to the output dir.")
   parser.add_argument("-of", "--outfile", action="store", dest="outfile", required=True, metavar='outputfilename', help="Outputfile name")
   parser.add_argument("-qsub", "--qsubPath", action="store", dest="qsub", required=True, metavar='/somepath/qsub', help="Full Path to the qsub executables of SGE.")
   parser.add_argument("-mapQ", "--mappingquality", action="store", dest="MAPQ", required=True, metavar='20', help="Mapping Quality Threshold")
   parser.add_argument("-baseQ", "--basequality", action="store", dest="BASEQ", required=True, metavar='20', help="BASE Quality Threshold")
   parser.add_argument("-mac", "--minimumAlternateCount", action="store", dest="MAC", required=True, metavar='2', help="Minimum Alternate Allele COunt")
   parser.add_argument("-maf", "--minimumAlternateFrequnecy", action="store", dest="MAF", required=True, metavar='0.01', help="Minimum Alternate Allele Frequency")
   parser.add_argument("-lai", "--leftAlignIndels", action="store_true", dest="lai", default=False, help="Pass if you wish to left align indels")
   
   args = parser.parse_args()
   if(args.verbose):
       print "I have Started the run for Freebayes."
   #(wd)= ProcessArgs(args)
   #if(wd.__ne__("NULL")):
   RunFreebayes(args)
   if(args.verbose):
       print "I have finished the run for Freebayes."   
   
def ProcessArgs(args):
    if(args.verbose):
        print "I am currently processing the arguments.\n"
    SampleDirName = args.patientId
    staticDir = "FreeBayesAnalysis"
    AnalysisDir = os.path.join(args.outdir,staticDir)
    SampleAnalysisDir = os.path.join(AnalysisDir,SampleDirName)
    if os.path.isdir(AnalysisDir):
        if(args.verbose):
            print "Dir:", AnalysisDir, " exists thus we wont be making it\n"
    else:
        os.mkdir(AnalysisDir)
        
    if os.path.isdir(SampleAnalysisDir):
            if(args.verbose):
                print "Dir:", SampleAnalysisDir," exists and we wont run the analysis\n"
            #return("NULL")
    else:
        os.mkdir(SampleAnalysisDir)
    if(args.verbose):
        print "I am done processing the arguments.\n"    
    return(SampleAnalysisDir)

def RunFreebayes(args):
    myPid = os.getpid()
    day = date.today()
    today = day.isoformat()
    today = today.replace("-","")
    #myPid = str(myPid)
    if(args.verbose):
        print "I am running freebayes for ", args.patientId, " using SGE"
    #Setting Job for SGE   
    cmd = args.FREEBAYES + " -b " + args.tbam + " -b " + args.nbam + " -f " + args.ref + " -v " + args.outfile + " -I " + " -X " + " -O " + " -m " + args.MAPQ + " -q " + args.BASEQ + " -F " + args.MAF + " -C " + args.MAC + " -J --genotype-qualities"
    #print "CMD==>",cmd,"\n"
    qsub_cmd = args.qsub + " -q " + args.queue + " -N " + "Freebayes_"+args.patientId+"_"+str(myPid) + " -o " + "Freebayes_"+ args.patientId+"_"+str(myPid)+".stdout" + " -e " + "Freebayes_"+args.patientId+"_"+str(myPid)+".stderr" + " -V -l h_vmem=6G,virtual_free=6G -pe smp " + args.threads + " -wd " + args.outdir + " -sync y " + " -b y " + cmd 
    print "QSUB_CMD==>", qsub_cmd , "\n"
    qsub_args = shlex.split(qsub_cmd)
    proc = Popen(qsub_args)
    proc.wait()
    retcode = proc.returncode
    if(retcode >= 0):
        if(args.verbose):
            print "I have finished running Freebayes for ", args.patientId, " using SGE"        
    else:
        if(args.verbose):
            print "Freebayes is either still running or it errored out with return code", retcode,"\n"    
               
    
if __name__ == "__main__":
    start_time = time.time()  
    main()
    end_time = time.time()
    print("Elapsed time was %g seconds" % (end_time - start_time))      
    