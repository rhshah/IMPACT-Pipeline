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
   parser = argparse.ArgumentParser(prog='dmp_filter_pindel.py', description='Filter Indels from the output of pindel', usage='%(prog)s [options]')
   parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=True, help="make lots of noise [default]")
   parser.add_argument("-i", "-inputVcf", action="store", dest="inputVcf", required=True, metavar='SomeID.vcf', help="Input vcf freebayes file which needs to be filtered")
   parser.add_argument("-tsn", "--tsampleName", action="store", dest="tsampleName", required=True, metavar='SomeName', help="Name of the tumor Sample")
   parser.add_argument("-dp", "--totaldepth", action="store", dest="dp", required=True, metavar='0', help="Tumor total depth threshold")
   parser.add_argument("-ad", "--alleledepth", action="store", dest="ad", required=True, metavar='3', help="Tumor allele depth threshold")
   parser.add_argument("-tnr", "--tnRatio", action="store", dest="tnr", required=True, metavar='5', help="Tumor-Normal variant frequency ratio threshold ")
   parser.add_argument("-vf", "--variantfrequency", action="store", dest="vf", required=True, metavar='0.01', help="Tumor variant frequency threshold ")
   parser.add_argument("-o", "--outDir", action="store", dest="outdir", required=True, metavar='/somepath/output', help="Full Path to the output dir.")
   parser.add_argument("-min", "--min_var_len", action="store", dest="min", required=False, metavar='25', help="Minimum length of the Indels")
   parser.add_argument("-max", "--max_var_len", action="store", dest="max", required=False, metavar='500', help="Max length of the Indels")
   #parser.add_argument("-qsub", "--qsubPath", action="store", dest="qsub", required=True, metavar='/somepath/qsub', help="Full Path to the qsub executables of SGE.")
   #parser.add_argument("-q", "--queue", action="store", dest="queue", required=False, metavar='all.q or clin.q', help="Name of the SGE queue")
   
   args = parser.parse_args()
   if(args.verbose):
       print "I have Started the run for doing standard filter."
   (stdfilterVCF)= RunStdFilter(args)
   if(args.verbose):
       print "I have finished the run for doing standard filter."   
def RunStdFilter(args):
    vcf_out = os.path.basename(args.inputVcf)
    vcf_out = os.path.splitext(vcf_out)[0]
    txt_out = vcf_out
    vcf_out = vcf_out + "_STDfilter.vcf"
    txt_out = txt_out + "_STDfilter.txt"
    vcf_reader = vcf.Reader(open(args.inputVcf, 'r'))   
    vcf_writer = vcf.Writer(open(vcf_out, 'w'), vcf_reader)
    txt_fh = open(txt_out, "wb") 
    allsamples = vcf_reader.samples
    sample1 = allsamples[0]
    sample2 = allsamples[1]
    if(sample1 == args.tsampleName):
        nsampleName = sample2
    else:
        nsampleName = sample1 
    for record in vcf_reader:
        tcall = record.genotype(args.tsampleName)
        recordType = record.INFO['SVTYPE']
        recordLen = abs(int(record.INFO['SVLEN']))
        #print tcall, "recordType: ",recordType, "recordLen: ", recordLen
            
        if(tcall['AD'] != None):
            trd = int(tcall['AD'][0])
            tad = int(tcall['AD'][1])
        else:
            trd = 0
            tad = 0
        tdp = trd + tad
        if(tdp != 0):
            tvf = int(tad) / float(tdp)
        else:
            tvf = 0
        #print "Tdata: ",trd,tad
        ncall = record.genotype(nsampleName)
        if(ncall):
            if(ncall['AD'] != None):
               nrd = int(ncall['AD'][0])
               nad = int(ncall['AD'][1])
            else:
                nrd = 0
                nad = 0
            ndp = nrd + nad
            if(ndp != 0):
                nvf = nad / ndp
            else:
                nvf = 0
            #print "Ndata: ",trd,tad
             
            nvfRF = int(args.tnr) * nvf      
            #print recordLen, args.min, args.max
            if((recordLen >= int(args.min)) and (recordLen <= int(args.max))):
                if(tvf > nvfRF):
                    if((tdp >= int(args.dp)) & (tad >= int(args.ad)) & (tvf >= float(args.vf))):
                        if(recordType != "RPL"): 
                        #print args.tsampleName , "\t" , record.CHROM , "\t" , str(record.POS) , "\t" , str(record.REF) , "\t" + str(record.ALT[0]) , "\t" , "." + "\n" 
                            vcf_writer.write_record(record)
                            txt_fh.write(args.tsampleName + "\t" + record.CHROM + "\t" + str(record.POS) + "\t" + str(record.REF) + "\t" + str(record.ALT[0]) + "\t" + "." + "\n")
        
    txt_fh.close()
    return(vcf_out)

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
