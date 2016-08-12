import os, sys, pytest, copy

# import repo's tests utilities
cur_dir = os.path.dirname(__file__)
path = os.path.abspath(os.path.join(cur_dir, '..', 'tests'))
if not path in sys.path:
    sys.path.insert(1, path)
del path
import test_util

# Set arguments used in all forest tests:
# Define all but the D (depth) and w parameter here; vary D,w for the tests
conf_params = {
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

class TestDepthW:
    '''
    Test the interaction of the depth parameter with the w parameter in the following network:
      
      A    B
      |    |
      C -> D

    '''
    def test_depth_3_w_1(self, msgsteiner):
        '''
        When the solver is constrained by depth, it still picks up net positive nodes.
          This results in a different network than if the solver is not constrained by depth.

        Note behavior of the following tests for comparison (which we do not repeat in this class):
          TestW#test_w_1, which has w = 1 and D = 5
            we are not constrained by depth and we pick A - C -> D
          TestDepth#test_depth_10, which has w = 0 and D = 10
            there is no forest penalty, so the optimal network is the same as this one
        '''
        params = copy.deepcopy(conf_params)
        params['D'] = 3
        params['w'] = 1
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
