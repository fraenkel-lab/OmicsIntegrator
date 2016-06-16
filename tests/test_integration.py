#!/usr/local/bin/python
'''
Large integration test script.
'''

import os, subprocess, filecmp, shutil, shlex, tempfile

# msgsteiner is a fixture object
def test_integration(msgsteiner):
    # msgsteiner is parsed using conftest.py
    assert msgsteiner is not None, 'Please provide path to msgsteiner using --msgpath option'
    
    # Forest inputs
    forest_conf = os.path.join(os.path.dirname(__file__), '..', 'example', 'a549', 'tgfb_forest.cfg') #provided config file
    phos_weights = os.path.join(os.path.dirname(__file__), '..', 'example', 'a549', 'Tgfb_phos.txt')

    edge_file = os.path.join(os.path.dirname(__file__), '..', 'data', 'iref_mitab_miscore_2013_08_12_interactome.txt') #interactom

    # Create a tmp directory for output	
    forest_out = tempfile.mkdtemp()

    #Arbitrary value
    seed = 2
    
    try:
        
        script_dir = os.path.dirname(__file__)
        forest_path = os.path.join(script_dir, '..', 'scripts', 'forest.py')
    
        fcmd='python %s --prize=%s --edge=%s --conf=%s  --outpath=%s --msgpath=%s --seed=%s'%(
            forest_path, phos_weights, edge_file, forest_conf, forest_out, msgsteiner, seed)
        subprocess.call(shlex.split(fcmd), shell=False)	

        curr_dir = os.path.dirname(__file__)
        match, mismatch, errors = filecmp.cmpfiles(forest_out, os.path.join(curr_dir, 'integration_test_standard'), 
                                   ['result_augmentedForest.sif', 'result_dummyForest.sif', 'result_edgeattributes.tsv',
                                    'result_info.txt', 'result_nodeattributes.tsv', 'result_optimalForest.sif'], 
                                   shallow=False)

        if len(match) != 6:
            print 'Mismatching files: ', mismatch
            print 'Errors: ', errors
            assert 0, 'Not all Forest output files match'

    except IOError:
        print 'IO error'
    finally:
        shutil.rmtree(forest_out)
