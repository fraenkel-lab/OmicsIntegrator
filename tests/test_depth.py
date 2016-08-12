import os, sys, pytest, copy

# import repo's tests utilities
cur_dir = os.path.dirname(__file__)
path = os.path.abspath(os.path.join(cur_dir, '..', 'tests'))
if not path in sys.path:
    sys.path.insert(1, path)
del path
import test_util

# Set arguments used in all forest tests:
# Define all but the D (depth) parameter here; vary D for the tests
conf_params = {
    'w': 0,
    'mu': 0,
    'g': 0,
    'b': 5
}
# Location of the prize, network, and root node files
forest_opts = {
  'prize': os.path.join(cur_dir, 'small_forest_tests', 'w_test_prizes.txt'),
  'edge': os.path.join(cur_dir, 'small_forest_tests', 'w_test_network.txt'),
  'dummyMode': os.path.join(cur_dir, 'small_forest_tests', 'w_test_roots.txt')
}

class TestDepth:
    '''
    Test various values of the depth parameter for the following network:
      
      A    B
      |    |
      C -> D
    '''
    def test_depth_0(self, msgsteiner):
        ''' 
        Run Forest with depth = 0; msgsteiner fails and test_util.run_forest raises an error

        INPUT:
        msgsteiner - fixture object with the value of --msgpath parsed by conftest.py
        '''
        params = copy.deepcopy(conf_params)
        params['D'] = 0
        with pytest.raises(Exception) as excinfo:
            graph = test_util.run_forest(msgsteiner, params, forest_opts)
        assert 'Forest did not generate output files' in excinfo.value.message

    def test_depth_2(self, msgsteiner):
        ''' 
        Run Forest with depth = 2; while this depth is valid, it results in no edges
          because the depth parameter is used as a strict inequality and there is a dummy node
        '''
        params = copy.deepcopy(conf_params)
        params['D'] = 2
        graph = test_util.run_forest(msgsteiner, params, forest_opts)

        assert graph.order() == 0, "Unexpected number of nodes"
        assert graph.size() == 0, "Unexpected number of edges"

    def test_depth_3(self, msgsteiner):
        '''
        Run forest with depth = 3; this is the minimum value of D for which any edges are 
          in the result network
        '''
        params = copy.deepcopy(conf_params)
        params['D'] = 3
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

    def test_depth_10(self, msgsteiner):
        '''
        Run forest with depth = 10; the results are the same as in test_depth_3
        '''
        params = copy.deepcopy(conf_params)
        params['D'] = 3
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
