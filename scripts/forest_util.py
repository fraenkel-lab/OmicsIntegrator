'''
Utility functions for Forest and network analysis
'''
import networkx as nx

def loadGraph(sifFile):
    '''Parses a .sif file from Forest.
    
    pp is undirected, pd directed.
    networkx currently does not support mixed graphs, so for all pp edges,
    a directed edge both to and from each node in the entry is created.
    
    INPUT: sifFile - the filename of the Forest output .sif file
    
    OUTPUT: a networkx DiGraph representing the network from Forest
    '''
    with open(sifFile, 'r') as sf:
        G = nx.DiGraph()
        for line in sf:
            edge = line.split()
            assert len(edge) == 3, '.sif file must contain two nodes and ' \
                'edge type on each line'
            assert edge[1] == 'pp' or edge[1] == 'pd', 'Edge type must be pp or pd'
            
            # Undirected edge is pair of directed edges
            if edge[1] == "pp":
                G.add_edge(edge[0], edge[2])
                G.add_edge(edge[2], edge[0])
            # Must be a directed edge if not undirected
            else:
                G.add_edge(edge[0], edge[2])
    return G
