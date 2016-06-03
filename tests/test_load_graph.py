'''
Test the forest_util loadGraph function
'''

import os, sys, tempfile, pytest

# Create the path to forest relative to the test_load_graph.py path
path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'scripts'))
if not path in sys.path:
    sys.path.insert(1, path)
del path

from forest_util import loadGraph

def write_valid_graph(tmpf):
    '''Write a valid .sif file
    
    INPUT: tmpf - a temporary file
    '''
    tmpf.write('A\tpp\tB\n')
    tmpf.write('B\tpd\tC\n')
    tmpf.write('C\tpd\tD\n')
    tmpf.write('D\tpd\tC\n')
    tmpf.write('A\tpp\tE\n')

class TestLoadGraph:
    
    def test_valid_graph(self):
        sifFilename= tempfile.NamedTemporaryFile(suffix='.sif', delete=False)
        try:
            # Write the network file
            write_valid_graph(sifFilename)
        finally:
            sifFilename.close()
            
        # Load the DiGraph
        graph = loadGraph(sifFilename.name)
        
        # Check that the DiGraph has the expected properties
        assert graph.order() == 5, "Unexpected number of nodes"
        assert graph.size() == 7, "Unexpected number of edges"
        
        # Check that the DiGraph has the expected edges
        assert graph.has_edge('A','B')
        assert graph.has_edge('B','A')
        assert graph.has_edge('B','C')
        assert graph.has_edge('C','D')
        assert graph.has_edge('D','C')
        assert graph.has_edge('A','E')
        assert graph.has_edge('E','A')
        
        # Remove because delete=False above
        os.remove(sifFilename.name)

    def test_invalid_edge_type(self):
        sifFilename= tempfile.NamedTemporaryFile(suffix='.sif', delete=False)
        try:
            # Write the network file
            write_valid_graph(sifFilename)
            # Add another line with an invalid edge
            sifFilename.write('A\tdd\tF\n')
        finally:
            sifFilename.close()
            
        # Load the DiGraph and expect an Exception for the invalid edge
        with pytest.raises(Exception):
            loadGraph(sifFilename.name)
                
        # Remove because delete=False above
        os.remove(sifFilename.name)

    def test_invalid_columns(self):
        sifFilename= tempfile.NamedTemporaryFile(suffix='.sif', delete=False)
        try:
            # Write the network file
            write_valid_graph(sifFilename)
            # Add another line with an extra column
            sifFilename.write('A\tdd\tF\t0.85\n')
        finally:
            sifFilename.close()
            
        # Load the DiGraph and expect an Exception for the extra column
        with pytest.raises(Exception):
            loadGraph(sifFilename.name)
                
        # Remove because delete=False above
        os.remove(sifFilename.name)
