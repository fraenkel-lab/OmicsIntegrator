'''
Test the Forest w parameter.

The test generates a configuration file for a particular value of w,
runs Forest using that configuration file on a small network and prize file,
and finally tests whether the optimal subnetwork matches the expected
subnetwork.  Different values of w are tested to check
whether the parameter has the desired effect.
'''
import os, sys, pytest, copy

# import repo's tests utilities
cur_dir = os.path.dirname(__file__)
path = os.path.abspath(os.path.join(cur_dir, '..', 'tests'))
if not path in sys.path:
    sys.path.insert(1, path)
del path
import test_util

# Set arguments used in all forest tests:
# Define all but the w parameter here; vary beta for the tests
conf_params = {
    'b': 5,
    'D': 5,
    'mu': 0,
    'g': 0
}
# Location of the prize, network, and root node files
forest_opts = {
  'prize': os.path.join(cur_dir, 'small_forest_tests', 'w_test_prizes.txt'),
  'edge': os.path.join(cur_dir, 'small_forest_tests', 'w_test_network.txt'),
  'dummyMode': os.path.join(cur_dir, 'small_forest_tests', 'w_test_roots.txt')
}

# All tests should in theory pass with beta = 1 but because msgsteiner
# does not converge to the global optimum, when w = 1 it returns an
# empty network instead of the optimal three node network
class TestW:      
    def test_w_025(self, msgsteiner):
        ''' Run Forest with w=0.25 and check optimal subnetwork

        INPUT:
        msgsteiner - fixture object with the value of --msgpath parsed by conftest.py
        '''
        params = copy.deepcopy(conf_params)
        params['w'] = 0.25
        graph = test_util.run_forest(msgsteiner, params, forest_opts)
        
        # Check that the DiGraph has the expected properties
        # Undirected edges are loaded as a pair of directed edges
        assert graph.order() == 4, "Unexpected number of nodes"
        assert graph.size() == 4, "Unexpected number of edges"

        # Check that the DiGraph has the expected edges
        assert graph.has_edge('A','C')
        assert graph.has_edge('C','A')
        assert graph.has_edge('B','D')
        assert graph.has_edge('D','B')


    def test_w_1(self, msgsteiner):
        ''' Run Forest with w=1 and check optimal subnetwork

        INPUT:
        msgsteiner - fixture object with the value of --msgpath parsed by conftest.py
        '''
        params = copy.deepcopy(conf_params)
        params['w'] = 1
        graph = test_util.run_forest(msgsteiner, params, forest_opts)
        
        # Check that the DiGraph has the expected properties
        # Undirected edges are loaded as a pair of directed edges
        assert graph.order() == 3, "Unexpected number of nodes"
        assert graph.size() == 3, "Unexpected number of edges"

        # Check that the DiGraph has the expected edges
        assert graph.has_edge('A','C')
        assert graph.has_edge('C','A')
        assert graph.has_edge('C','D')


    def test_w_10(self, msgsteiner):
        '''tw =  Run Forest with w=10 and check optimal subnetwork

        INPUT:
        msgsteiner - fixture object with the value of --msgpath parsed by conftest.py
        '''
        params = copy.deepcopy(conf_params)
        params['w'] = 10
        graph = test_util.run_forest(msgsteiner, params, forest_opts)

        # Check that the DiGraph has the expected properties
        assert graph.order() == 0, "Unexpected number of nodes"
        assert graph.size() == 0, "Unexpected number of edges"


    def test_w_missing(self, msgsteiner):
        ''' Run Forest with w missing from the configuration file and check
        that Forest does not produce an output subnetwork

        INPUT:
        msgsteiner - fixture object with the value of --msgpath parsed by conftest.py
        '''
        params = copy.deepcopy(conf_params)
        with pytest.raises(Exception) as excinfo:
            graph = test_util.run_forest(msgsteiner, params, forest_opts)
        # Forest will not raise an Exception, it simply exits if w is missing
        # Have to check for missing output files to determine whether it exited
        assert 'Forest did not generate output files' in excinfo.value.message
    
    # Could update Forest to require positive w and test for that
