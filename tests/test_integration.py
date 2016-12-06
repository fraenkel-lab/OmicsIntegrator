#!/usr/local/bin/python
'''
End-to-end integration test.
'''

import os, subprocess, filecmp, shutil, shlex, tempfile

# msgsteiner is a fixture object
def test_forest_integration(msgsteiner):
    '''
    Forest integration test.  Runs Forest on the example data in the a549
    directory and compares the output files with stored reference files.
    '''
    # msgsteiner is parsed using conftest.py
    assert msgsteiner is not None, 'Please provide path to msgsteiner using --msgpath option'
    
    # Forest inputs
    curr_dir = os.path.dirname(__file__)
    forest_conf = os.path.join(curr_dir, '..', 'example', 'a549', 'tgfb_forest.cfg') #provided config file
    phos_weights = os.path.join(curr_dir, '..', 'example', 'a549', 'Tgfb_phos.txt')

    edge_file = os.path.join(curr_dir, '..', 'data', 'iref_mitab_miscore_2013_08_12_interactome.txt') #interactom

    # Create a tmp directory for output
    forest_out = tempfile.mkdtemp()

    #Arbitrary value
    seed = 2
    
    try:
        forest_path = os.path.join(curr_dir, '..', 'scripts', 'forest.py')
    
        fcmd='python %s --prize=%s --edge=%s --conf=%s  --outpath=%s --msgpath=%s --seed=%s'%(
            forest_path, phos_weights, edge_file, forest_conf, forest_out, msgsteiner, seed)
        subprocess.call(shlex.split(fcmd), shell=False)	

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

def test_garnet_integration():
    '''
    Garnet integration test.  Runs Garnet on the example FOXA1 data in the a549
    directory and compares the output files with stored reference files.
    '''
    # Garnet config file
    curr_dir = os.path.dirname(__file__)
    garnet_conf = os.path.join(curr_dir, '..', 'example', 'a549', 'tgfb_foxa1_garnet.cfg')

    # Create a tmp directory for output
    garnet_out = tempfile.mkdtemp()

    try:
        garnet_path = os.path.join(curr_dir, '..', 'scripts', 'garnet.py')

        garnet_cmd='python %s --outdir=%s %s' % (garnet_path, garnet_out, garnet_conf)
        subprocess.call(shlex.split(garnet_cmd), shell=False)

        # Test all of the Garnet output files
        output_files = ['events_to_genes.fsa',
                        'events_to_genes.xls',
                        'events_to_genes_with_motifs.pkl',
                        'events_to_genes_with_motifs.tgm',
                        'events_to_genes_with_motifs.txt',
                        'events_to_genes_with_motifs_geneids.txt',
                        'events_to_genes_with_motifs_tfids.txt',
                        'events_to_genes_with_motifsregression_results.html',
                        'events_to_genes_with_motifsregression_results.tsv',
                        'events_to_genes_with_motifsregression_results_FOREST_INPUT.tsv']
        match, mismatch, errors = filecmp.cmpfiles(garnet_out, os.path.join(curr_dir, 'integration_test_standard'),
                                   output_files, shallow=False)

        if len(match) != len(output_files):
            print 'Mismatching files: ', mismatch
            print 'Errors: ', errors
            assert 0, 'Not all Forest output files match'

    except IOError:
        print 'IO error'
    finally:
        shutil.rmtree(garnet_out)
