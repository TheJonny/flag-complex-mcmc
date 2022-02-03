from data_importer import simplex, clique, bbp, c_elegans, densifier, random_with_p
from scipy.sparse import coo_matrix

import numpy as np

def join_graphs(a,b):
    '''join 2 graphs a and b given by dense numpy adjacency matrices'''
    return np.block([
        [a, np.zeros(shape=(a.shape[0], b.shape[1]))],
        [np.zeros(shape=(b.shape[0], a.shape[0])), b]
        ])

def save_unweighted_flag(fname, graph):
    assert(graph.shape[0] == graph.shape[1])
    N = graph.shape[0]
    content  = "dim 0:\n"
    content += ("1 "*N)[:-1] + "\n"
    content += "dim 1:\n"
    coo_graph = coo_matrix(graph)
    for i,j in zip(coo_graph.row, coo_graph.col):
        content += f"{i} {j} 1\n"

    with open(fname, 'w') as f:
        f.write(content)

if __name__ == '__main__':
    ex_00 = simplex(3)

    ex_01 = simplex(3)
    ex_01[0,3] = 1

    ex_02 = simplex(3)
    ex_02[2,3] = 1

    ex_03 = clique(3)

    ex_04 = densifier(
                [0, 0, 1, 3, 3],
                [1, 2, 2, 1, 2])

    ex_05 = densifier(
                [0, 0, 1, 1, 3],
                [1, 2, 2, 3, 2])
    
    ex_06 = join_graphs(ex_04, ex_05)

    ex_07 = densifier(
                [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 3, 4, 5, 6, 7, 8, 9],
                [1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 2, 2, 2, 2, 2, 2, 2]
            )

    ex_08 = c_elegans()

    ex_09 = bbp(0)

    ex_20 = random_with_p(100,0.05)     #A
    ex_21 = random_with_p(1000,0.05)    #B
    ex_22 = random_with_p(10000,0.05)   #C

    ex_23 = random_with_p(10000,0.0005) #same amount of edges as B
    ex_24 = random_with_p(10000,0.005)  #10 times the edges as B


    examples = [ex_00, ex_01, ex_02, ex_03,
            ex_04, ex_05, ex_06,
            ex_07,
            ex_08, ex_09,
            ex_20, ex_21, ex_22, ex_23, ex_24, 
            ]

    for i, ex in enumerate(examples):
        print(i, ex)
        save_unweighted_flag(f'/tmp/{i:02}.flag', ex)

