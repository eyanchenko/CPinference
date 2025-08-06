import cpnet
import numpy as np
import networkx as nx

def cp_cur(A):
    
    algorithm = cpnet.LapCore()
    G = nx.from_numpy_array(A)
    algorithm.detect(G)
    x = algorithm.get_coreness()
    
    return list(x.values())
