from data_importer import simplex, clique, bbp, c_elegans, densifier, random_with_p
from scipy.sparse import coo_matrix

import numpy as np

def join_graphs(a,b):
    '''join 2 graphs a and b given by dense numpy adjacency matrices'''
    return np.block([
        [a, np.zeros(shape=(a.shape[0], b.shape[1]))],
        [np.zeros(shape=(b.shape[0], a.shape[0])), b]
        ])

def seoify(g):
    '''Perform SearchEngineOptimization on the graph g by turning all double edges to single edges'''
    (N,N) = g.shape
    for i in range(N):
        for j in range(N):
            if g[i,j] == 1 and g[j,i] == 1:
                if np.random.uniform() < 0.5:
                    g[i,j] = 0
                else:
                    g[j,i] = 0
    return g


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

    # smallest example with H_2(G) = 1.
    ex_10 = densifier(
                [0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4],
                [1, 2, 3, 5, 2, 3, 4, 4, 5, 4, 5, 5]
            )

    # like ex_10, but only cycles and thus s_2 = 0.
    # TODO
    ex_11 = densifier(
                [0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4],
                [1, 2, 3, 5, 2, 3, 4, 4, 5, 4, 5, 5]
            )

    ex_12 = join_graphs(ex_10, simplex(2))

    ex_20 = random_with_p(100,0.05)     #A
    ex_21 = random_with_p(1000,0.05)    #B
    ex_22 = random_with_p(10000,0.05)   #C

    ex_23 = random_with_p(10000,0.0005) #same amount of edges as B
    ex_24 = random_with_p(10000,0.005)  #10 times the edges as B

    ### single edge only graphs
    ex_30 = seoify(c_elegans())
    ex_31 = seoify(random_with_p(300,0.05))


    examples = [ex_00, ex_01, ex_02, ex_03,
            ex_04, ex_05, ex_06,
            ex_07,
            ex_08, ex_09,
            ex_10, ex_11, ex_12,
            ex_20, ex_21, ex_22, ex_23, ex_24, 
            ex_30, ex_31,
            ]

    for i, ex in enumerate(examples):
        print(i, ex)
        save_unweighted_flag(f'/tmp/{i:02}.flag', ex)

