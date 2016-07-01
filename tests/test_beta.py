import os, sys, pytest, copy

# import repo's tests utilities
cur_dir = os.path.dirname(__file__)
path = os.path.abspath(os.path.join(cur_dir, '..', 'tests'))
if not path in sys.path:
    sys.path.insert(1, path)
del path
import test_util

# Set arguments used in all forest tests:
# Define all but the b (beta) parameter here; vary beta for the tests
conf_params = {
    'w': 0,
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

class TestBeta:
    '''
    Test various values of the beta parameter:
    p'(v) = beta * p(v) - mu * deg(v)
    min f'(F) = sum_{v not in V_F} p'(v) + sum_{e in E_F} c(e) + w * K

    Higher values of beta increase the importance of prizes relative to their degree penalty, and
    also increase the importance of prizes relative to the cost of edges

    Use the following test network:

    A    B
     \    \
      C -> D

    p(A) = 0
    p(B) = 0
    p(C) = 1
    p(D) = 1
    c(AC) = 0.25
    c(BD) = 0.25
    c(CD) = 0.75
    '''
    def test_beta_0(self, msgsteiner):
        ''' 
        Run Forest with beta = 0 and check optimal subnetwork

        In p'(v) = beta * p(v) - mu * deg(v), beta = 0 implies that an empty graph is the optimal network

        INPUT:
        msgsteiner - fixture object with the value of --msgpath parsed by conftest.py
        '''
        params = copy.deepcopy(conf_params)
        params['b'] = 0
        graph = test_util.run_forest(msgsteiner, params, forest_opts)

        assert graph.order() == 0, "Unexpected number of nodes"
        assert graph.size() == 0, "Unexpected number of edges"
        
    def test_beta_024(self, msgsteiner):
        '''
        Run Forest with beta = 0.24

        This value of beta is just below the 0.25 value of beta which sets the prize benefit 
        of C and D equal to the cost of obtaining those prizes through edges AC and BD: we expect
        an empty network again
        '''
        params = copy.deepcopy(conf_params)
        params['b'] = 0.24
        graph = test_util.run_forest(msgsteiner, params, forest_opts)

        assert graph.order() == 0, "Unexpected number of nodes"
        assert graph.size() == 0, "Unexpected number of edges"
  
    def test_beta_026(self, msgsteiner):
        '''
        Run Forest with beta = 0.26

        See test_beta_024; this value of beta makes it worth while to obtain prizes for C and D;
        we expect the following network:

        A   B
         \   \
          C   D

        note an undirected edge ab in a digraph is represented by the network containing both edges ab and ba
        '''
        params = copy.deepcopy(conf_params)
        params['b'] = 0.26
        graph = test_util.run_forest(msgsteiner, params, forest_opts)

        assert graph.order() == 4
        assert graph.size() == 4

        assert graph.has_edge('A', 'C')
        assert graph.has_edge('C', 'A')
        assert graph.has_edge('B', 'D')
        assert graph.has_edge('D', 'B')
