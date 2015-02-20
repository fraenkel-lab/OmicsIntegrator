#!/usr/local/bin/python
'''
This script will build a TF-binding matrix from human dnase clusters, enabling you to analyze
any type of gene expression data
'''

__author__='Sara JG Gosline'
__email__='sgosline@mit.edu'

import os
from optparse import OptionParser

if __name__=='__main__':
    parser=OptionParser()
    opts,args=parser.parse_args()
    
    #garnet requires a configuration file that has all the data
    garnet_out='mcf7_output'
    garnet_conf='mcf7_garnet.cfg' #provided config file
    gcmd='python ../../scripts/garnet.py --outdir=%s %s'%(garnet_out,garnet_conf) #command
    
    
    ##now ready to run commands
        
    print gcmd
    os.system(gcmd)

    
