#!/usr/local/bin/python
'''
Large integration test script.
'''

import os, subprocess, filecmp, shutil, shlex
from optparse import OptionParser

def test_integration(msgsteiner):

    parser = OptionParser()

    parser.add_option('--msgpath',dest='msgsteiner',type='string',help='Path to msgsteiner9 code, be sure to include!')

    phos_weights = os.path.join(os.path.dirname(__file__), '..', 'example', 'a549', 'Tgfb_phos.txt')

    #provided config file
    garnet_conf = os.path.join(os.path.dirname(__file__), '..', 'example', 'a549', 'tgfb_garnet.cfg')

    #forest requires more inputs
    forest_conf = os.path.join(os.path.dirname(__file__), '..', 'example', 'a549', 'tgfb_forest.cfg') #provided config file

    edge_file = os.path.join(os.path.dirname(__file__), '..', 'data', 'iref_mitab_miscore_2013_08_12_interactome.txt') #interactom

    #WE NEED MSGSTEINER9 INSTALLED!!!
    msgsteinerpath = msgsteiner
    if msgsteinerpath == None:
        print 'Please provide path to msgsteiner using --msgpath option'
        assert 0
	
    forest_out = 'temp'

    #Arbitrary value
    seed = 2

    # Create a tmp directory for output
    script_dir = os.path.dirname(__file__)
    abs_file_path = os.path.join(script_dir, forest_out)
    os.makedirs(abs_file_path)

    forest_path = os.path.join(script_dir, '..', 'scripts', 'forest.py')

    fcmd='python %s --prize=%s --edge=%s --conf=%s  --outpath=%s --msgpath=%s --seed=%s'%(
        forest_path, phos_weights, edge_file, forest_conf, abs_file_path, msgsteinerpath, seed)
    subprocess.call(shlex.split(fcmd), shell=False)	

    curr_dir = os.path.dirname(__file__)
    results = filecmp.cmpfiles(abs_file_path, os.path.join(curr_dir, 'integration_test_standard'), 
		['result_augmentedForest.sif', 'result_dummyForest.sif', 'result_edgeattributes.tsv',
		'result_info.txt', 'result_nodeattributes.tsv', 'result_optimalForest.sif'], 
		shallow=False)

	
    shutil.rmtree(os.path.join(curr_dir, 'temp'))

    if len(results[0]) != 6:
        assert 0
    else:
        assert 1
