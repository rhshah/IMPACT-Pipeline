'''
Created on 07/24/2014
Modified on 03/31/2015
Note: Now it requires the libIntelDeflater.so file as -Dsamjdk.intel_deflater_so_path in command line
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

def main():
   parser = argparse.ArgumentParser(prog='Run_AbraRealignment.py', description='Run ABRA Indel Realignment', usage='%(prog)s [options]')
   parser.add_argument("-i", "--bamList", action="store", dest="bamList", required=True, metavar='BamFile.list', help="Full path to the tumor bam files as a fof.") 
   parser.add_argument("-p", "--patientId", action="store", dest="patientId", required=True, metavar='PatientID', help="Id of the Patient for which the bam files are to be realigned")
   parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=True, help="make lots of noise [default]")
   parser.add_argument("-t", "--threads", action="store", dest="threads", required=True, metavar='5', help="Number of Threads to be used to run ABRA")
   parser.add_argument("-d", "--mdp", action="store_true", dest="dp", default=1000, help="Threshold for downsampling depth to run ABRA")
   parser.add_argument("-k", "--kmers", action="store", nargs='*', dest="kmers", required=True, metavar='43', help="Number of k-mers to be used to run ABRA; Multiple k-mers are separated by space")
   parser.add_argument("-temp", "--temporaryDirectory", action="store", dest="tmpDir", required=True, metavar='/somepath/tmpdir', help="Full Path to temporary directory")
   parser.add_argument("-r", "--referenceFile", action="store", dest="ref", required=True, metavar='/somepath/Homo_Sapeins_hg19.fasta', help="Full Path to the reference file with the bwa index.")
   parser.add_argument("-a", "--abraJar", action="store", dest="ABRA", required=True, metavar='/somepath/ABRA.jar', help="Full Path to the ABRA jar file.")
   parser.add_argument("-tr", "--targetRegion", action="store", dest="targetRegion", required=True, metavar='/somepath/targetRegion.bed', help="Full Path to the target region bed file")
   parser.add_argument("-j", "--javaPATH", action="store", dest="JAVA", required=True, metavar='/somepath/java', help="Path to java executable.")
   parser.add_argument("-b", "--bwaPATH", action="store", dest="BWA", required=True, metavar='/somepath/bin', help="Path to the bin of bwa executable.")
   parser.add_argument("-q", "--queue", action="store", dest="queue", required=False, metavar='all.q or clin.q', help="Name of the SGE queue")
   parser.add_argument("-o", "--outDir", action="store", dest="outdir", required=True, metavar='/somepath/output', help="Full Path to the output dir.")
   parser.add_argument("-qsub", "--qsubPath", action="store", dest="qsub", required=False, metavar='/somepath/qsub', help="Full Path to the qsub executables of SGE.")
   parser.add_argument("-bsub", "--bsubPath", action="store", dest="bsub", required=False, metavar='/somepath/bsub', help="Full Path to the bsub executables of LSF.")
   
   args = parser.parse_args()
   print "Running the ABRA re-alignment process."
   (inBamList,outBamList,kmers,tmpdir)= ProcessArgs(args)
   RunABRA(args,inBamList,outBamList,kmers,tmpdir)
   RunHouseKeeping(args,tmpdir)
   print "Done running ABRA realignment."        

def ProcessArgs(args):
    if(args.verbose):
        print "Lets make some noise while doing re-alignment"
    inBamFiles = []
    outBamFiles = []
    inBamList = ""
    outBamList = ""
    kmers = ""
    if(args.verbose):
        print "Going to see how many bam files are there for us to run the analysis."
    #Check qsub and bsub
    if(args.qsub and args.bsub):
       print "Please give either qsub or bsub arguments. Script does not allow usage of both\n"
       sys.exit(1)           
   if((not args.qsub) and (not args.bsub)): 
       print "Please give either qsub or bsub arguments. None are provided\n"
       sys.exit(1)       
    
    #Open fof of bam and make a list
    with open(args.bamList, 'r') as filecontent:
        for line in filecontent:
            data = line.rstrip('\n')
            infile = data
            outfile = data.replace('.bam','_IR.bam')
            inBamFiles.append(infile)
            outBamFiles.append(outfile)
    #Check if there are more then one files to run 
    if (inBamFiles.__len__() > 1):
        inBamList = ",".join(inBamFiles)
    else:
        inBamList = "".join(inBamFiles)
    #Check if there are more then one files to run    
    if (outBamFiles.__len__() > 1):
        outBamList = ",".join(outBamFiles)
    else:
        outBamList = "".join(outBamFiles)
    if(args.verbose):
        print "Going to run analysis for", inBamFiles.__len__(), "bam files."
    
    #Check how many kemrs the analysis should be ran for
    if(args.verbose):
        print "Checking the number of k-mers given for analysis."
    if (args.kmers.__len__() > 1):
        kmers = ",".join(args.kmers)
    else:
        kmers = "".join(args.kmers)
    
    #Make temporary directory specific to this samples
    myPid = os.getpid()
    if(args.verbose):
        print "Checking to see if there are any existing temporary directory for same Patient ID. We will try to remove it,","\n","If we are unsuccessful we would exit."
    tmpDir_name = "abra_temp_" + str(myPid) + "_"
    tmpDir_name = tmpDir_name + args.patientId
    tmpdir = os.path.join(args.tmpDir,tmpDir_name)
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir,onerror = remShut)
    
    return(inBamList,outBamList,kmers,tmpdir)
   
def RunABRA(args,inBamList,outBamList,kmers,tmpdir):
    myPid = os.getpid()
    #myPid = str(myPid)
    if(args.verbose):
        print "Running ABRA realignment for ", args.patientId, " using SGE"
        
    abra_dir, abra_filename = os.path.split(args.ABRA)
    #libdir = abra_dir.replace("target","lib")
    #libFile = libdir + "/libIntelDeflater.so"
    #Setting Job for SGE   
    #cmd = args.JAVA + " -Dsamjdk.intel_deflater_so_path=" + libFile + " -Xmx40g -jar " + args.ABRA + " --in " + inBamList + " --out " + outBamList + " --ref " + args.ref + " --targets " + args.targetRegion + " --threads "  + args.threads + " --mad " + str(args.dp) + " --kmer " + kmers + " --working " + tmpdir
    cmd = args.JAVA + " -Xmx40g -jar " + args.ABRA + " --in " + inBamList + " --out " + outBamList + " --ref " + args.ref + " --targets " + args.targetRegion + " --threads "  + args.threads + " --mad " + str(args.dp) + " --kmer " + kmers + " --working " + tmpdir
    #print "CMD==>",cmd,"\n"
    os.environ["PATH"] = os.environ["PATH"] + ":" + args.BWA
    cl_cmd = ''
    mem = int(args.threads) * 5
    maxmem = int(mem)+5
    if(args.qsub):
        cl_cmd = args.qsub + " -q " + args.queue + " -N " + "ABRA_"+args.patientId+"_"+str(myPid) + " -o " + "ABRA_"+ args.patientId+"_"+str(myPid)+".stdout" + " -e " + "ABRA_"+args.patientId+"_"+str(myPid)+".stderr" + " -V -l h_vmem=5G,virtual_free=5G -pe smp " + args.threads + " -wd " + args.outdir + " -sync y " + " -b y " + cmd 
    else:
        cl_cmd = args.bsub + " -q " + args.queue + " -J " + "ABRA_"+args.patientId+"_"+str(myPid) + " -o " + "ABRA_"+ args.patientId+"_"+str(myPid)+".stdout" + " -e " + "ABRA_"+args.patientId+"_"+str(myPid)+".stderr" + " -We 24:00 -R \"rusage[mem=" + str(mem) + "]\" -M " + str(maxmem) + " -n " + args.threads + " -cwd " + args.outdir + " -K " + cmd
    print "CLUSTER_CMD==>", cl_cmd , "\n"
    cl_args = shlex.split(cl_cmd)
    proc = Popen(cl_args)
    proc.wait()
    retcode = proc.returncode
    if(retcode >= 0):
        if(args.verbose):
            print "Finished Running ABRA realignment for ", args.patientId, " using SGE/LSF"
    else:
        if(args.verbose):
            print "ABRA is either still running or it errored out with return code", retcode,"\n" 
def RunHouseKeeping(args,tmpdir):
    if(args.verbose):
        print "Removing Temporary Directory:", tmpdir,"\n"
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir,onerror = remShut)
def remShut(*args):
    #onerror returns a tuple containing function, path and exception info
    func, path, _ = args 
    os.chmod(path, stat.S_IWRITE)
    os.remove(path)
    
if __name__ == "__main__":
    start_time = time.time()  
    main()
    end_time = time.time()
    print("Elapsed time was %g seconds" % (end_time - start_time))