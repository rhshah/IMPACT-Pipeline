'''
Created on 08/20/2014
@Ronak Shah

'''
import argparse
import glob
import os
import sys
import time
import stat
from subprocess import Popen
from subprocess import call
import shlex
import shutil
from pprint import pprint
import re

def main():
   parser = argparse.ArgumentParser(prog='dmp_run_pipeline_in_batch.py', description='Run DMP IMPACT Pipeline In Batch', usage='%(prog)s [options]')
   parser.add_argument("-i", "--dmsRunInfo", action="store", dest="dmsRunInfo", required=True, metavar='RunIlluminaProcess.pl', help="Full path to RunIlluminaProcess.pl.") 
   parser.add_argument("-z", "--pipeline", action="store", dest="pipeline", required=True, metavar='RunInformation.txt', help="Full path to the file containing dms database dump for each Run.") 
   parser.add_argument("-c", "--dmpConf", action="store", dest="dmpConf", required=True, metavar='dmp_impact.conf', help="Full Path to dmp_impact.conf file")
   parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=True, help="make lots of noise [default]")
   parser.add_argument("-sc", "--svconf", action="store", dest="svconfig", required=False, metavar='configuration_cv3.txt', help="Full Path to configuration_cv3.txt file")
   parser.add_argument("-p", "--process", action="store", nargs='*', dest="process", required=True, metavar='1', help="Number of process to be used to run pipeline; Multiple process are separated by space")
   parser.add_argument("-o", "--outDir", action="store", dest="outdir", required=True, metavar='/somepath/output', help="Full Path to the output dir.")
   args = parser.parse_args()
   #Check how many process the analysis should be ran for
   if(args.verbose):
       print "Checking the number of process given for analysis."

   if (args.process.__len__() > 1):
       process = ",".join(args.process)
   else:
       process = "".join(args.process)
    
   if(process.startswith('3')):
       RunForThree(args)

def RunForThree(args):
    if(args.verbose):
        print "Going to Run IMPACT pipeline for all Runs in ", args.dmsRunInfo,"\n"
    #Open fof of bam and make a list
    with open(args.dmsRunInfo, 'r') as filecontent:
        for line in filecontent:
            data = line.rstrip('\n').split(",")
            fastqLocation = data[0]
            poolName = data[1]
            BamLocation = data[2]
            poolOutput = args.outdir + "/" + poolName
            if(args.verbose):
                print "Will try to run the process on " , poolName,".\n"
            #make the output dir
            if os.path.isdir(poolOutput):
                if(args.verbose):
                    print "The output directory ", poolOutput, " exists & thus we will skip the run ", poolName
                continue
            else:    
                #print "test\n"#
                if(args.verbose):
                    print"Making directory & changing directory to", poolOutput," to run the Pipeline.\n"
                os.mkdir(poolOutput)
                call(['chmod', '755', poolOutput])
                os.chdir(poolOutput)
            
            grepPattern = BamLocation + "/*.bam"
            filelist = glob.glob(grepPattern)
            #Srt bam list
            if(args.verbose):
                print "Making the InputBam.list file\n"
            InputBamList = poolOutput + "/InputBam.list" 
            InputBam = open(InputBamList,'w')
            call(['chmod', '755', InputBamList])
            for bamFile in filelist:
                baiFile = bamFile.replace('.bam','.bai')
                fileName = os.path.basename(bamFile)
                baseName = re.search('(.*)_MD.*', fileName).group(1)
                srtBam = baseName + ".bam"
                srtBai = baseName + ".bai"
                mdBam = baseName + "_MD.bam"
                mdBai = baseName + "_MD.bai"
                destsrtBam = poolOutput + "/" + srtBam
                destsrtBai = poolOutput + "/" + srtBai
                InputBam.write("%s\n" % destsrtBam)
                destmdBam = poolOutput + "/" + mdBam
                destmdBai = poolOutput + "/" + mdBai
                print bamFile,"\n"
                print baiFile,"\n"
                print destsrtBam,"\n"
                print destmdBam,"\n"
                print destsrtBai,"\n"
                print destmdBai,"\n"
                #Softlink files
                if os.path.isfile(bamFile):
                    #print "test1\n"
                    if(args.verbose):
                        print "Making symbolic soft links for bam files.\n"
                    os.symlink(bamFile, destsrtBam)
                    os.symlink(bamFile, destmdBam)
                    
                if os.path.isfile(baiFile):
                    #print "test2\n"
                    if(args.verbose):
                        print "Making symbolic soft links for index bai files.\n"
                    os.symlink(baiFile, destsrtBai)
                    os.symlink(baiFile, destmdBai)
            
            
            InputBam.close()
            #sys.exit()
            dstDmpConf = ""
            dstSVconf = ""     
            if(args.dmpConf):    
                #assert not os.path.isabs(args.dmpConf)
                if(args.verbose):
                    print "Copying the pipeline configuration file.\n"
                dstDmpConf =  os.path.join(poolOutput, os.path.basename(args.dmpConf))    
                shutil.copy(args.dmpConf, dstDmpConf)
            if(args.svconfig):    
                #assert not os.path.isabs(args.svconfig)
                if(args.verbose):
                    print "Copying the SV configuration file.\n"
                dstSVconf =  os.path.join(poolOutput, os.path.basename(args.svconfig))
                shutil.copy(args.svconfig, dstSVconf)
            if(args.verbose):
                print "Making cmd to run pipeline.\n"
            if(args.svconfig):  
                cmd = args.pipeline + " -c " + dstDmpConf + " -sc " + dstSVconf + " -d " + fastqLocation + " -o " + poolOutput
            else:
                cmd = args.pipeline + " -c " + dstDmpConf + " -d " + fastqLocation + " -o " + poolOutput
            print "cmd:",cmd,"\n"
            print "Running Pipeline for pool " ,poolName,".\n"
            cmd_args = shlex.split(cmd)
            proc = Popen(cmd_args)
            proc.wait()
            print "ReturnCode:", proc.returncode
            print "Finished Running Pipleine for pool ", poolName, ".\n"
            #sys.exit()

if __name__ == "__main__":
    start_time = time.time()  
    main()
    end_time = time.time()
    print("Elapsed time was %g seconds" % (end_time - start_time))