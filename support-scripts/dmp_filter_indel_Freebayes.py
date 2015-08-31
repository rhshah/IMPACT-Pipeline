'''
Created on 07/31/2014
@Ronak Shah

'''
from __future__ import division
import argparse
import sys
import time
import os.path
import stat
from subprocess import Popen
import shlex
import shutil
from datetime import date
import vcf
import copy
def main():
   parser = argparse.ArgumentParser(prog='dmp_filter_indel.py', description='Filter Indels from the output of freebayes', usage='%(prog)s [options]')
   parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=True, help="make lots of noise [default]")
   parser.add_argument("-i", "-inputVcf", action="store", dest="inputVcf", required=True, metavar='SomeID.vcf', help="Input vcf freebayes file which needs to be filtered")
   parser.add_argument("-tsn", "--tsampleName", action="store", dest="tsampleName", required=True, metavar='SomeName', help="Name of the tumor Sample")
   parser.add_argument("-dp", "--totaldepth", action="store", dest="dp", required=True, metavar='0', help="Tumor total depth threshold")
   parser.add_argument("-ad", "--alleledepth", action="store", dest="ad", required=True, metavar='3', help="Tumor allele depth threshold")
   parser.add_argument("-tnr", "--tnRatio", action="store", dest="tnr", required=True, metavar='5', help="Tumor-Normal variant frequency ratio threshold ")
   parser.add_argument("-vf", "--variantfrequency", action="store", dest="vf", required=True, metavar='0.01', help="Tumor variant frequency threshold ")
   parser.add_argument("-o", "--outDir", action="store", dest="outdir", required=True, metavar='/somepath/output', help="Full Path to the output dir.")
   parser.add_argument("-vlb", "--vcflibBin", action="store", dest="vlb", required=True, metavar='/dmp/resources/dev2/tools/bio/vcflibs/current/bin/', help="Full Path to the vcf lib bin")
   parser.add_argument("-qsub", "--qsubPath", action="store", dest="qsub", required=True, metavar='/somepath/qsub', help="Full Path to the qsub executables of SGE.")
   parser.add_argument("-q", "--queue", action="store", dest="queue", required=False, metavar='all.q or clin.q', help="Name of the SGE queue")
   
   args = parser.parse_args()
   if(args.verbose):
       print "I have Started the run for doing standard filter."
   #(wd)= ProcessArgs(args)
   #if(wd.__ne__("NULL")):
   #(decomposedVCF) = DecomposeVCF(args)
   (bmsVCF) = BreakVCF(args)
   (stdfilterVCF)= RunStdFilter(args,bmsVCF)
   if(args.verbose):
       print "I have finished the run for doing standard filter."   
def BreakVCF(args):
    myPid = os.getpid()
    vcf_out = os.path.basename(args.inputVcf)
    vcf_out = os.path.splitext(vcf_out)[0]
    vcf_out = vcf_out + "_BMS.vcf"
    vcfbm = args.vlb + "/vcfbreakmulti"
    if os.path.isfile(vcf_out) and os.access(vcf_out, os.R_OK):
        if(args.verbose):
            print "File:",vcf_out ," exists and is readable"
        return(vcf_out)
    else:
        cmd = vcfbm + " " + args.inputVcf
        qsub_cmd = args.qsub + " -q " + args.queue + " -N " + "FreebayesBMS_"+args.tsampleName+"_"+str(myPid) + " -o " + vcf_out + " -e " + "Freebayes_"+args.tsampleName+"_"+str(myPid)+".stderr" + " -V -l h_vmem=6G,virtual_free=6G -pe smp 1 " + " -wd " + args.outdir + " -sync y " + " -b y " + cmd
        if(args.verbose):
            print "BMS cmd: ", qsub_cmd
            qsub_args = shlex.split(qsub_cmd)
            proc = Popen(qsub_args)
            proc.wait()
    return(vcf_out)
def DecomposeVCF(args):
    myPid = os.getpid()
    vcf_out = os.path.basename(args.inputVcf)
    vcf_out = os.path.splitext(vcf_out)[0]
    vcf_out = vcf_out + "_decompose.vcf"
    vcfdc = args.vlb + "/vcfallelicprimitives"
    if os.path.isfile(vcf_out) and os.access(vcf_out, os.R_OK):
        if(args.verbose):
            print "File:",vcf_out ," exists and is readable"
        return(vcf_out)
    else:
        cmd = vcfdc + " " + args.inputVcf + " -t SplitComplex -g -k"
        qsub_cmd = args.qsub + " -q " + args.queue + " -N " + "FreebayesDecompose_"+args.tsampleName+"_"+str(myPid) + " -o " + vcf_out + " -e " + "Decompose_"+args.tsampleName+"_"+str(myPid)+".stderr" + " -V -l h_vmem=6G,virtual_free=6G -pe smp 1 " + " -wd " + args.outdir + " -sync y " + " -b y " + cmd
        if(args.verbose):
            print "Decompose cmd: ", qsub_cmd
            qsub_args = shlex.split(qsub_cmd)
            proc = Popen(qsub_args)
            proc.wait()
    return(vcf_out)
def RunStdFilter(args,bmsVCF):
    vcf1_out = os.path.basename(args.inputVcf)
    vcf1_out = os.path.splitext(vcf1_out)[0]
    txt1_out = vcf1_out
    vcf_out = vcf1_out + "_STDfilter.vcf"
    txt_out = txt1_out + "_STDfilter.txt"
    vcfc_out = vcf1_out + "_STDfilterComplex.vcf"
    txtc_out = txt1_out + "_STDfilterComplex.txt"
    vcf_reader = vcf.Reader(open(bmsVCF, 'r'))    
    allsamples = vcf_reader.samples
    sample1 = allsamples[0]
    sample2 = allsamples[1]
    if(sample1 == args.tsampleName):
        nsampleName = sample2
    else:
        nsampleName = sample1
    vcf_writer = vcf.Writer(open(vcf_out, 'w'), vcf_reader)
    vcfc_writer = vcf.Writer(open(vcfc_out, 'w'), vcf_reader)
    txt_fh = open(txt_out, "wb") 
    txtc_fh = open(txtc_out, "wb") 
    for record in vcf_reader:
        tcall = record.genotype(args.tsampleName)
        fbAlt = list(str(record.ALT[0]))
        fbRef = list(str(record.REF))
        recordType = record.INFO['TYPE'][0]
        print "recordType: ",recordType
        if(recordType != "complex"):
            print "ORG=>" , record.REF, record.ALT[0], record.CHROM, record.POS, "\n"
            (newRef, newAlt) = ConvertRefAlt(fbRef, fbAlt, record.CHROM, record.POS,recordType)
            record.REF = newRef
            record.ALT[0] = newAlt
            print "NEW=>" , record.REF, record.ALT[0], record.CHROM, record.POS, "\n"
            
        if(tcall['DP'] != None):
            tdp = int(tcall['DP'])
        else:
            tdp = 0
        if(tcall['AO'] != None):
            tad = tcall['AO']
            try:
                tad = int(tcall['AO'])
            except TypeError:
                tad = int(tad[0])           
        else:
            tad = 0
         
        if(tdp != 0):
            tvf = int(tad) / float(tdp)
        else:
            tvf = 0
            # print "I am here " , float(tvf) 
        # if(tdp > tad):    
            # print tcall['DP'], "==", tdp,"\t" ,tcall['AO'], "==", tad, "\tTVF==", tvf
        ncall = record.genotype(nsampleName)
        if(ncall):
            if(ncall['DP'] != None):
                ndp = int(ncall['DP'])
            else:
                ndp = 0
            if(ncall['AO'] != None):
                nad = ncall['AO']
                try:
                    nad = int(ncall['AO'])
                except TypeError:
                    nad = int(nad[0])
            else:
                nad = 0
            if(ndp != 0):
                nvf = nad / ndp
            else:
                nvf = 0
        else:
            ndp = 0
            nad = 0
            nvf = 0.0
             
        nvfRF = int(args.tnr) * nvf      
        
        if(tvf > nvfRF):
            if(recordType != "complex"):
                if((tdp >= int(args.dp)) & (tad >= int(args.ad)) & (tvf >= float(args.vf))):
                        vcf_writer.write_record(record)
                        txt_fh.write(args.tsampleName + "\t" + record.CHROM + "\t" + str(record.POS) + "\t" + newRef + "\t" + newAlt + "\t" + "." + "\n")
            else:
                if((tdp >= int(20)) & (tad >= int(8)) & (tvf >= float(0.05))):
                    vcfc_writer.write_record(record)
                    txtc_fh.write(args.tsampleName + "\t" + record.CHROM + "\t" + str(record.POS) + "\t" + str(record.REF) + "\t" + str(record.ALT[0]) + "\t" + "." + "\n")
    
    txt_fh.close()
    txtc_fh.close()
    return(vcf_out)

def ConvertRefAlt(fbRef,fbAlt,chr,pos,recordType):
    orgRef = copy.copy(fbRef)
    orgAlt = copy.copy(fbAlt)
    lenRef = len(orgRef)
    lenAlt = len(orgAlt)
    count = 0
    #print "ORG1",orgRef,orgAlt
    if(lenRef > lenAlt):
        while(fbRef):
            #print "count=",count
            if(fbAlt):
                itemAlt = fbAlt.pop()
            else:
                break
            if(fbRef):
                itemRef = fbRef.pop()
            else:
                break
            if((itemAlt != '') & (itemRef != '')):
                if(itemAlt == itemRef):
                    count = count+1
                    continue
                else:
                    break
        newRef = copy.copy(orgRef)
        newAlt = copy.copy(orgAlt)
        if(lenAlt == 1):
            newRef = newRef
        elif(lenAlt == 2):
            del newRef[-(count)]
        else:
            del newRef[lenRef-count:]
        del newAlt[1:]
        #print "Del",lenAlt,lenRef,count ,chr,pos, orgRef,orgAlt,"".join(newRef),"".join(newAlt)
        return("".join(newRef),"".join(newAlt))
    if(lenAlt > lenRef):
        while (fbRef):
            print "count=",count
            itemAlt = fbAlt.pop()
            itemRef = fbRef.pop()
            if(itemAlt == itemRef):
                count = count+1
                continue
            else:
                break
        newRef = copy.copy(orgRef)
        newAlt = copy.copy(orgAlt)
        del newRef[1:]
        del newAlt[lenAlt-count:]
        #print "Ins",count,chr,pos , orgRef,orgAlt,"".join(newRef),"".join(newAlt)
        return("".join(newRef),"".join(newAlt))
def list_duplicates_of(seq,item):
    start_at = -1
    locs = []
    while True:
        try:
            loc = seq.index(item,start_at+1)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return locs
if __name__ == "__main__":
    start_time = time.time()  
    main()
    end_time = time.time()
    print("Elapsed time was %g seconds" % (end_time - start_time))   
