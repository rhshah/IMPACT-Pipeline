'''
Created on 03/30/2015
Description: This will run find covered interval program from GATK.
@Ronak H Shah

'''

import argparse
import os
import sys
import time
import stat
from subprocess import Popen
import shlex
import shutil
import pybedtools

def main():
   parser = argparse.ArgumentParser(prog='Run_FindCoveredInterval.py', description='This will run find covered interval program from GATK.', usage='%(prog)s [options]')
   parser.add_argument("-i", "--bamList", action="store", dest="bamList", required=True, metavar='BamFile.list', help="Full path to the tumor bam files as a fof.") 
   parser.add_argument("-of", "--outFilePrefix", action="store", dest="outFilePrefix", required=True, metavar='OutFilePrefix', help="Output Covered Interval File Prefix for the bam files.")
   parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=True, help="make lots of noise [default]")
   parser.add_argument("-t", "--threads", action="store", dest="threads", required=True, metavar='5', help="Number of Threads to be used to run FindCoveredIntervals")
   parser.add_argument("-dp", "--totaldepth", action="store", dest="dp", default=20, metavar='20', help="Total depth threshold")
   parser.add_argument("-mbq", "--minbasequality", action="store", dest="mbq", default=20, metavar='20', help="Threshold for minimum base quality for Running Find Covered Interval")
   parser.add_argument("-mmq", "--minmappingquality", action="store", dest="mmq", default=20, metavar='20', help="Threshold for minimum mapping quality for Running Find Covered Interval")
   parser.add_argument("-r", "--referenceFile", action="store", dest="ref", required=True, metavar='/somepath/Homo_Sapeins_hg19.fasta', help="Full Path to the reference file with the bwa index.")
   parser.add_argument("-g", "--gatkJar", action="store", dest="GATK", required=True, metavar='/somepath/GenomeAnalysisTK.jar', help="Full Path to the GATK jar file.")
   parser.add_argument("-j", "--javaPATH", action="store", dest="JAVA", required=True, metavar='/somepath/java', help="Path to java executable.")
   parser.add_argument("-q", "--queue", action="store", dest="queue", required=True, metavar='all.q or clin.q', help="Name of the SGE queue")
   parser.add_argument("-o", "--outDir", action="store", dest="outdir", required=True, metavar='/somepath/output', help="Full Path to the output dir.")
   parser.add_argument("-qsub", "--qsubPath", action="store", dest="qsub", required=True, metavar='/somepath/qsub', help="Full Path to the qsub executables of SGE.")
   args = parser.parse_args()
   if(args.verbose):
       print "Running the FindCoveredIntervals for the re-alignment process."
   (listFile) = RunFindCoveredIntervals(args)
   (bedFile) = ListToBed(listFile,args)
   (listFile) = Bed2List(bedFile, args)
   if(args.verbose):
       print "Done running FindCoveredIntervals for the re-alignment process."  
   
def RunFindCoveredIntervals(args):
    outFile = args.outdir + "/" + args.outFilePrefix + "_covered.list"
    outFileSrt = args.outdir + "/" + args.outFilePrefix + "_covered_srt.list"
    if(os.path.isfile(outFileSrt)):
       if(args.verbose):
           print outFileSrt, " already exists\n"
    else:
        cmd = args.JAVA + " -Xmx20g -jar " + args.GATK + " -T FindCoveredIntervals -R " + args.ref + " -I " + args.bamList + " -minBQ " + str(args.mbq) + " -minMQ " + str(args.mmq) + " -cov " + str(args.dp) + " -o " + outFile + " -rf FailsVendorQualityCheck -rf BadMate -rf UnmappedRead -rf BadCigar"
        myPid = os.getpid()
        qsub_cmd = args.qsub + " -q " + args.queue + " -N " + "FindCoveredIntervals_" + str(myPid) + " -o " + "FindCoveredIntervals_" + str(myPid) + ".stdout" + " -e " + "FindCoveredIntervals_"+ str(myPid) +".stderr" + " -V -l h_vmem=8G,virtual_free=8G -pe smp " + str(args.threads) + " -wd " + args.outdir + " -sync y " + " -b y " + cmd 
        if(args.verbose):
            print "QSUB_CMD==>", qsub_cmd , "\n"
            qsub_args = shlex.split(qsub_cmd)
            proc = Popen(qsub_args)
        proc.wait()
        retcode = proc.returncode
        if(retcode >= 0):
            if(args.verbose):
                print "Finished Running FindCoveredIntervals for ", myPid, " using SGE"
        else:
            if(args.verbose):
                print "FindCoveredIntervals is either still running or it errored out with return code", retcode,"\n" 
    return(outFile)

def ListToBed(file,args):
    outFile = args.outdir + "/" + args.outFilePrefix + "_covered.bed"
    outFileSrt = args.outdir + "/" + args.outFilePrefix + "_covered_srt.bed"
    if(os.path.isfile(outFileSrt)):
        if(args.verbose):
            print outFileSrt, " already exists\n"
    else:
        outHandle = open(outFile,"w")
        if(os.stat(file).st_size == 0):
            outHandle.write("1\t963754\t963902\n")
        else:
            with open(file,'r') as filecontent:
                for line in filecontent:
                    data = line.rstrip('\n').split(":")
                    chr = data[0]
                    if "-" in data[1]:
                        (st,en) = data[1].split("-")
                    else:
            			st = data[1]
            			en = int(data[1]) + 1
                    outHandle.write(str(chr) + "\t" + str(st) + "\t" + str(en) + "\n")
        outHandle.close()
        bedtool = pybedtools.BedTool(outFile)
        stbedtool = bedtool.sort()
        mbedtool = stbedtool.merge(d=50)
        c = mbedtool.saveas(outFileSrt)
        if(os.path.isfile(outFileSrt)):
            os.remove(outFile)
    return(outFileSrt)
    
def Bed2List(file,args):
    outFile = args.outdir + "/" + args.outFilePrefix + "_covered_srt.list"
    if(os.path.isfile(outFile)):
        if(args.verbose):
            print outFile, " already exists\n"
    else:
        outHandle = open(outFile,"w")
        if(os.stat(file).st_size == 0):
            outHandle.write("1:963754-963902\n")
        else:
            with open(file,'r') as filecontent:
                 for line in filecontent:
                     data = line.rstrip('\n').split("\t")
                     outHandle.write(str(data[0]) + ":" + str(data[1]) + "-" + str(data[2]) + "\n")
        outHandle.close()
    return(outFile)

if __name__ == "__main__":
    start_time = time.time()  
    main()
    end_time = time.time()
    print("Elapsed time was %g seconds" % (end_time - start_time))
    
