#!/usr/local/bin/python
'''
A test script to run through the sample dataset provided with garnet-forest
'''

__author__='Sara JG Gosline'
__email__='sgosline@mit.edu'

import os
from optparse import OptionParser

if __name__=='__main__':
    parser=OptionParser()
    parser.add_option('--forest-only',dest='forest_only',action='store_true',default=False,help='Set this flag to run forest on phospho-proteomic data only.  DEFAULT:%default')
    parser.add_option('--msgsteiner',dest='msgsteiner',type='string',help='Path to msgsteiner9 code, be sure to include!')
    opts,args=parser.parse_args()
    
    #garnet requires a configuration file that has all the data
    garnet_conf='tgfb_garnet.cfg' #provided config file
    gcmd='python ../scripts/garnet.py '+garnet_conf #command
    
    #forest requires more inputs
    forest_conf='tgfb_forest.cfg' #provided config file
    dummy_nodes='Tgfb_interactors.txt' #use proteins that interact with TgfB as 'dummyNodes'
    phos_weights='Tgfb_phos.txt' #proteins with changes in phosphorylation levels
    edge_file='../../data/iref_mitab_miscore_2013_08_12_interactome.txt' #interactome

    msgsteinerpath=opts.msgsteiner ##WE NEED MSGSTEINER9 INSTALLED!!!
    forest_out='tgfb_forest_output'
    
    ##now ready to run commands
    if not opts.forest_only:
        print gcmd
        os.system(gcmd)
        garnet_output=''
        garnet_beta='2'
        output=forest_out+'_with_garnet'
        if not os.path.exists(output): ##BE SURE TO MAKE DIRECTORY RESULT
            os.system('mkdir '+output)
        fcmd='python ../../scripts/forest.py --prize=%s --edge=%s --conf=%s --garnet=%s --garnetBeta=%s --outpath=%s --msgpath=%s --dummyMode=%s'%(phos_weights,edge_file,forest_conf,garnet_output,garnet_beta,output,msgsteinerpath,dummy_nodes)
        print fcmd
        os.system(fcmd)

    else:
        if not os.path.exists(forest_out):##BE SURE TO MAKE DIRECTORY FOR RESULT
            os.system('mkdir '+forest_out)
        fcmd='python ../../scripts/forest.py --prize=%s --edge=%s --conf=%s  --outpath=%s --msgpath=%s --dummyMode=%s'%(phos_weights,edge_file,forest_conf,forest_out,msgsteinerpath,dummy_nodes)
        print fcmd
        os.system(fcmd)

    
