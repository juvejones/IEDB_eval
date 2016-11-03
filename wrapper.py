#!/usr/bin/python

####################################################################################
### THIS MODULE IS MODIFIED FOR MHCII PREDICTION SUBMISSION AND DATA INTEGRATION ###
####################################################################################
import subprocess,re,os,sys,getopt

import package
from package import PepData, Dir
from package import predict, processing, fileroutine, container

__location__ = os.path.realpath(
    os.path.join(os.getcwd(),os.path.dirname(__file__)))

global methods
methods =[
          #"consensus",
          #"NetMHCIIpan",
          #"nn_align",
          #"smm_align",
          "comblib",
          "tepitope"
          ]
         
global dirname 
dirname = Dir()       

def main(argv):
    try:
        opts, args = getopt.getopt(argv,"hl:i:t:o:d:",["log=","idir=","tdir=","odir=","ddir="])
    except getopt.GetoptError:
        print('wrapper.py -l <logfile> -i <input_dir> -t <tmp_dir> -o <output_dir> -d <data_dir>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('wrapper.py -l <logfile> -i <input_dir> -t <tmp_dir> -o <output_dir> -d <data_dir>')
            sys.exit()
        elif opt in ("-l", "--log"):
            dirname.logfile = arg
        elif opt in ("-i", "--idir"):
            dirname.input_dir = arg
        elif opt in ("-t", "--tdir"):
            dirname.tmp_dir = arg
        elif opt in ("-o", "--odir"):
            dirname.output_dir = arg
        elif opt in ("-d", "--ddir"):
            dirname.data_dir = arg
    print('Log:%s   Input:%s   Tmp:%s   Output:%s   Data:%s' 
        % (dirname.logfile, dirname.input_dir, dirname.tmp_dir, dirname.output_dir, dirname.data_dir))

    if len(sys.argv) < 5:
        print('need specify names of directories')
        sys.exit(2)
     
if __name__ == "__main__":
    main(sys.argv[1:])
    indexF = open(os.path.join(__location__,dirname.logfile),'r')
    alleles = fileroutine.findAllele(indexF)
    print(alleles)
    indexF.close()
    fileroutine.InputSeq(alleles, dirname)
    predict.runPredict(methods, alleles, dirname)
