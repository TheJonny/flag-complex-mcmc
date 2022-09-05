from data_importer import simplex, clique, bbp, c_elegans, densifier, random_with_p
from scipy.sparse import coo_matrix

import numpy as np

OUTDIR = '/tmp/'

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
    save_unweighted_flag(f'{OUTDIR}/00.flag', ex_00)

    ex_01 = simplex(3)
    ex_01[0,3] = 1
    save_unweighted_flag(f'{OUTDIR}/01.flag', ex_01)

    ex_02 = simplex(3)
    ex_02[2,3] = 1
    save_unweighted_flag(f'{OUTDIR}/02.flag', ex_02)

    ex_03 = clique(3)
    save_unweighted_flag(f'{OUTDIR}/03.flag', ex_03)

    ex_04 = densifier(
                [0, 0, 1, 3, 3],
                [1, 2, 2, 1, 2])
    save_unweighted_flag(f'{OUTDIR}/04.flag', ex_04)

    ex_05 = densifier(
                [0, 0, 1, 1, 3],
                [1, 2, 2, 3, 2])
    save_unweighted_flag(f'{OUTDIR}/05.flag', ex_05)
    
    ex_06 = join_graphs(ex_04, ex_05)
    save_unweighted_flag(f'{OUTDIR}/06.flag', ex_06)

    ex_07 = densifier(
                [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 3, 4, 5, 6, 7, 8, 9],
                [1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 2, 2, 2, 2, 2, 2, 2]
            )
    save_unweighted_flag(f'{OUTDIR}/07.flag', ex_07)


    save_unweighted_flag(f'{OUTDIR}/c_elegans.flag', c_elegans())
    save_unweighted_flag(f'{OUTDIR}/bbp0.flag', bbp(0))
    save_unweighted_flag(f'{OUTDIR}/bbp0_l13.flag', bbp(0, allowed_neuron_types=['L1_DAC', 'L1_DLAC', 'L1_HAC', 'L1_NGC-DA', 'L1_NGC-SA', 'L1_SLAC', 'L23_BP', 'L23_BTC', 'L23_ChC', 'L23_DBC',
        'L23_LBC', 'L23_MC', 'L23_NBC', 'L23_NGC', 'L23_PC', 'L23_SBC']))        # only layers 1-3
    save_unweighted_flag(f'{OUTDIR}/bbp0_l14.flag', bbp(0, allowed_neuron_types=['L1_DAC', 'L1_DLAC', 'L1_HAC', 'L1_NGC-DA', 'L1_NGC-SA', 'L1_SLAC', 'L23_BP', 'L23_BTC', 'L23_ChC', 'L23_DBC',
        'L23_LBC', 'L23_MC', 'L23_NBC', 'L23_NGC', 'L23_PC', 'L23_SBC', 'L4_BP', 'L4_BTC', 'L4_ChC', 'L4_DBC', 'L4_LBC',
        'L4_MC', 'L4_NBC', 'L4_NGC', 'L4_PC', 'L4_SBC', 'L4_SP', 'L4_SS']))      #only layers 1-4



    ex_20 = random_with_p(100,0.05)     #A
    save_unweighted_flag(f'{OUTDIR}/20.flag', ex_20)
    ex_21 = random_with_p(1000,0.05)    #B
    save_unweighted_flag(f'{OUTDIR}/21.flag', ex_21)
    ex_22 = random_with_p(10000,0.05)   #C
    save_unweighted_flag(f'{OUTDIR}/22.flag', ex_22)

    ex_23 = random_with_p(10000,0.0005) #same amount of edges as B
    save_unweighted_flag(f'{OUTDIR}/23.flag', ex_23)
    ex_24 = random_with_p(10000,0.005)  #10 times the edges as B
    save_unweighted_flag(f'{OUTDIR}/24.flag', ex_24)


