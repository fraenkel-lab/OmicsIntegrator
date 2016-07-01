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
    'mu': 2,
    'g': 0
}
# Location of the prize, network, and root node files
forest_opts = {
  'prize': os.path.join(cur_dir, 'small_forest_tests', 'beta_mu_test_prizes.txt'),
  'edge': os.path.join(cur_dir, 'small_forest_tests', 'beta_mu_test_network.txt'),
  'dummyMode': os.path.join(cur_dir, 'small_forest_tests', 'beta_mu_test_roots.txt')
}

class TestBetaMu:
    '''
    Test various values of the beta parameter with respect to its interaction with a non-zero mu
    p'(v) = beta * p(v) - mu * deg(v)
    min f'(F) = sum_{v not in V_F} p'(v) + sum_{e in E_F} c(e) + w * K

    The mu parameter is typically used in practice to exclude well-studied gene regulatory elements
    that appear as high-confidence hubs in bioinformatics databases.

    Use the following test network:

      B - A - E
       \ / \ /
        C   D

    p(A) = 6
    p(B) = 5
    p(C) = 5
    p(D) = 5
    p(E) = 5
    c(AB) = 0.9
    c(AC) = 0.9
    c(AD) = 0.9
    c(AE) = 0.9
    c(BC) = 0.1
    c(DE) = 0.1
    '''
    def test_beta_1(self, msgsteiner):
        ''' 
        In p'(v) = beta * p(v) - mu * deg(v), beta = 1 is too small to overcome the hub penalty with mu = 2:
          p'(A) = 6 - 8 = -2

        We expect forest to use the more costly edges BC and CD in its network instead
        '''
        params = copy.deepcopy(conf_params)
        params['b'] = 1
        graph = test_util.run_forest(msgsteiner, params, forest_opts)

        try:
            assert graph.order() == 4, "Unexpected number of nodes"
            assert graph.size() == 4, "Unexpected number of edges"
    
            assert graph.has_edge('B', 'C')
            assert graph.has_edge('C', 'B')
            assert graph.has_edge('D', 'E')
            assert graph.has_edge('E', 'D')
        except AssertionError as e:
            test_util.print_graph(graph)
            raise e
        
    def test_beta_2(self, msgsteiner):
        '''
        See test_beta_1; beta = 2 is enough to overcome the hub penalty so that the hub at A is chosen:
          p'(A) = 12 - 8 = 4
        '''
        params = copy.deepcopy(conf_params)
        params['b'] = 2
        graph = test_util.run_forest(msgsteiner, params, forest_opts)

        try:
            assert graph.order() == 3, "Unexpected number of nodes"
            assert graph.size() == 4, "Unexpected number of edges"
    
            assert graph.has_edge('A', 'C')
            assert graph.has_edge('C', 'A')
            assert graph.has_edge('A', 'E')
            assert graph.has_edge('E', 'A')
        except AssertionError as e:
            test_util.print_graph(graph)
            raise e
