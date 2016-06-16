'''
Test the Forest w parameter.

The test generates a configuration file for a particular value of w,
runs Forest using that configuration file on a small network and prize file,
and finally tests whether the optimal subnetwork matches the expected
subnetwork.  Different values of w are tested to check
whether the parameter has the desired effect.
'''

import os, sys, subprocess, shlex, tempfile, shutil, pytest

# Create the path to forest relative to the test_load_graph.py path
cur_dir = os.path.dirname(__file__)
path = os.path.abspath(os.path.join(cur_dir, '..', 'scripts'))
if not path in sys.path:
    sys.path.insert(1, path)
del path

from forest_util import loadGraph

def write_conf(tmpf, w):
    '''Write a configuration file
    
    INPUT:
    tmpf - a temporary file
    w - the value of w
    '''
    if w is not None:
        tmpf.write('w = %f\n' % w)
    tmpf.write('b = 1\n')
    tmpf.write('D = 5\n')
    tmpf.write('mu = 0\n')

class TestW:      

    def shared_w_test(self, msgsteiner, w=None):
        '''Run Forest on the test network for a particular w
        
        INPUT:
        w - the value of w (default None)
        msgsteiner - the path to msgsteiner
        OUTPUT:
        graph - the DiGraph object for the optimal Steiner forest
        '''
        assert msgsteiner is not None, 'Please provide path to msgsteiner using --msgpath option'

        # For reproducibility, though it is likely not needed for this
        # small test case        
        seed = 2016
        
        conf_filename= tempfile.NamedTemporaryFile(suffix='.txt', delete=False)
        try:
            # Write the configuration file with the specified w
            write_conf(conf_filename, w)
        finally:
            conf_filename.close()
        
        # Location of the prize and network files
        network_filename = os.path.join(cur_dir, 'small_forest_tests', 'w_test_network.txt')
        prize_filename = os.path.join(cur_dir, 'small_forest_tests', 'w_test_prizes.txt')

        # Create a tmp directory for output	
        forest_out = tempfile.mkdtemp()
        
        try:
            forest_path = os.path.join(cur_dir, '..', 'scripts', 'forest.py')
            forest_cmd = 'python %s --prize=%s --edge=%s --conf=%s  --outpath=%s --msgpath=%s --seed=%s' % \
                (forest_path, prize_filename, network_filename, conf_filename.name, forest_out, msgsteiner, seed)
            subprocess.call(shlex.split(forest_cmd), shell=False)	
    
            # Only test the optimal Forest to see if w has the intended effect,
            # ignore the other output files
            opt_forest = os.path.join(forest_out, 'result_optimalForest.sif')
            assert os.path.isfile(opt_forest), 'Forest did not generate output files'
            graph = loadGraph(opt_forest)
            
        finally:
            # Remove here because delete=False above
            os.remove(conf_filename.name)
            # Remove the Forest output directory and all files
            shutil.rmtree(forest_out)
            
        return graph

    def test_w_025(self, msgsteiner):
        ''' Run Forest with w=0.25 and check optimal subnetwork

        INPUT:
        msgsteiner - fixture object with the value of --msgpath parsed by conftest.py
        '''
        graph = self.shared_w_test(msgsteiner, 0.25)
        
        # Check that the DiGraph has the expected properties
        # Undirected edges are loaded as a pair of directed edges
        assert graph.order() == 4, "Unexpected number of nodes"
        assert graph.size() == 4, "Unexpected number of edges"
        
    def test_w_1(self, msgsteiner):
        ''' Run Forest with w=1 and check optimal subnetwork

        INPUT:
        msgsteiner - fixture object with the value of --msgpath parsed by conftest.py
        '''
        graph = self.shared_w_test(msgsteiner, 1)
        
        # Check that the DiGraph has the expected properties
        # Undirected edges are loaded as a pair of directed edges
        assert graph.order() == 3, "Unexpected number of nodes"
        assert graph.size() == 3, "Unexpected number of edges"
        
    def test_w_2(self, msgsteiner):
        ''' Run Forest with w=2 and check optimal subnetwork

        INPUT:
        msgsteiner - fixture object with the value of --msgpath parsed by conftest.py
        '''
        graph = self.shared_w_test(msgsteiner, 2)
        
        # Check that the DiGraph has the expected properties
        assert graph.order() == 0, "Unexpected number of nodes"
        assert graph.size() == 0, "Unexpected number of edges"

    def test_w_missing(self, msgsteiner):
        ''' Run Forest with w missing from the configuration file and check
        that Forest does not produce an output subnetwork

        INPUT:
        msgsteiner - fixture object with the value of --msgpath parsed by conftest.py
        '''
        with pytest.raises(Exception) as excinfo:
            self.shared_w_test(msgsteiner)
        # Forest will not raise an Exception, it simply exits if w is missing
        # Have to check for missing output files to determine whether it exited
        assert 'Forest did not generate output files' in excinfo.value.message
    
    # Could update Forest to require positive w and test for that
