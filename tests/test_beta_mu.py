import os, sys, pytest, copy
from numpy import isclose

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
    'w': 1,
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
    c(e) = 1 - c'(e) [c' is "confidence" while c is "cost"]

    The mu parameter is typically used in practice to exclude well-studied proteins
    that appear as high-confidence hubs in bioinformatics databases.

    Use the following test network (with directed edges as mentioned in c' below):

      B - A - D
       \ / \ /
        C   E

    p(A) = 5
    p(B) = 6
    p(C) = 6
    p(D) = 6
    p(E) = 6
    c'(AB) = 0.9
    c'(AC) = 0.9
    c'(AD) = 0.9
    c'(AE) = 0.9
    c'(BC) = 0.1
    c'(DE) = 0.1
    '''
    def test_beta_1(self, msgsteiner):
        ''' 
        In p'(v) = beta * p(v) - mu * deg(v), beta = 1 is too small to overcome the hub penalty with mu = 2:
          p'(A) = 1*5 - 2*4 = -3

        We expect forest to use the more costly edges BC and CD in its network instead
        '''
        params = copy.deepcopy(conf_params)
        params['b'] = 1
        graph, objective = test_util.run_forest(msgsteiner, params, forest_opts)

        try:
            assert graph.order() == 4, "Unexpected number of nodes"
            assert graph.size() == 2, "Unexpected number of edges"
    
            assert graph.has_edge('B', 'C')
            assert graph.has_edge('D', 'E')
            
            # Check that the optimal forest has the correct objective function
            # value, using isclose to allow for minor floating point variation
            # Objective function: 0.8
            # Excluded prizes: -3.0
            # Edge costs: 1.8
            # Number of trees * w: 2 * 1.0 = 2.0
            assert isclose(0.8, objective, rtol=0, atol=1e-5), 'Incorrect objective function value'

        except AssertionError as e:
            test_util.print_graph(graph)
            raise e
        
    def test_beta_2(self, msgsteiner):
        '''
        See test_beta_1; beta = 2 is enough to overcome the hub penalty so that the hub at A is chosen:
          p'(A) = 2*5 - 2*4 = 2
        '''
        params = copy.deepcopy(conf_params)
        params['b'] = 2
        graph, objective = test_util.run_forest(msgsteiner, params, forest_opts)

        try:
            assert graph.order() == 5, "Unexpected number of nodes"
            assert graph.size() == 4, "Unexpected number of edges"
    
            assert graph.has_edge('A', 'B')
            assert graph.has_edge('A', 'C')
            assert graph.has_edge('A', 'D')
            assert graph.has_edge('A', 'E')
            
            # Check that the optimal forest has the correct objective function
            # value, using isclose to allow for minor floating point variation
            # Objective function: 1.4
            # Excluded prizes: 0
            # Edge costs: 0.4
            # Number of trees * w: 1 * 1.0 = 1.0
            assert isclose(1.4, objective, rtol=0, atol=1e-5), 'Incorrect objective function value'
            
        except AssertionError as e:
            test_util.print_graph(graph)
            raise e
