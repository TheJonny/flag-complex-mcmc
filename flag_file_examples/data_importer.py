import numpy as np


#
# biological data - statistical reconstructions
#
def bbp(i, allowed_neuron_types='all'):
    import os
    if not os.path.exists(f'data/bbp/average/cons_locs_pathways_mc{i}_Column.h5'):
        raise RuntimeError(
        """
        The connectome for BBP was not found. Please go to
        https://bbp.epfl.ch/nmc-portal/downloads.html
        to sign their form and place that data in 
        data/bbp/average/cons_locs_pathways_mc0_Column.h5
        """
       )
    import h5py
    h = h5py.File(f'data/bbp/average/cons_locs_pathways_mc{i}_Column.h5')
    neuron_types = list(h['connectivity'].keys())
    if type(allowed_neuron_types) is list:
        neuron_types = allowed_neuron_types
    if allowed_neuron_types in ['exc', 'inh']:
        exc_neuron_types =  [x for x  in neuron_types if 'PC' in x] + ['L4_SS', 'L4_SP']
        if allowed_neuron_types == 'exc':
            neuron_types = exc_neuron_types
        else:
            neuron_types = set(neuron_types) - set(exc_neuron_types)
    N_per_nt = [h[f'populations/{nt}/locations'].shape[0] for nt in neuron_types]
    N = sum(N_per_nt)

    c = np.zeros((N,N), dtype='bool')
    inds = np.cumsum([0]+N_per_nt)
    for i, nti in enumerate(neuron_types):
        for j, ntj in enumerate(neuron_types):
            c[inds[i]:inds[i+1], inds[j]:inds[j+1]] = h[f'connectivity/{nti}/{ntj}/cMat'][:]
    return c

    
#
# biological data - dense reconstructions
#
def c_elegans():
    # 279 neurons with 2194 directed synapses taken from https://github.com/lrvarshney/elegans which represent
    # the chemical network in C.elegans according to
    # https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1001066
    import os
    if not os.path.exists('data/c.elegans/A_sendjoint.mat'):
        from urllib.request import urlretrieve
        os.makedirs('data/c.elegans/', exist_ok=True)
        urlretrieve('https://github.com/lrvarshney/elegans/raw/master/A_sendjoint.mat', 'data/c.elegans/A_sendjoint.mat')
    from scipy.io import loadmat
    A = loadmat('data/c.elegans/A_sendjoint.mat')['Ac']
    return (A != 0).toarray()

#
# artifical examples
#
def simplex(d):
    '''returns a d+1 weight matrix corresponding to a simplex of dimension d'''
    return np.tril(np.ones((d+1,d+1)), k=-1)


def clique(d):
    '''returns a d+1 weight matrix corresponding to a clique of size d+1: A graph with d+1 nodes and all posible edges
    except reflexive ones'''
    c = np.ones((d+1,d+1))
    np.fill_diagonal(c, 0)
    return c



#
# 0 models for comparison
#
def random_like(c, exact=False):
    '''return a connectome matrix of shape like A, with a global connection density like A and empty main diagonal'''
    assert c.ndim == 2 and c.shape[0] == c.shape[1]

    D = c.shape[0]
    num_nonzero = (c != 0).sum()

    if exact:
        ones = np.zeros(D * (D - 1))
        ones[:num_nonzero] = 1
        np.random.shuffle(ones)

        c0 = np.zeros((D, D))
        c0[np.where(~np.eye(D, dtype=bool))] = ones

        assert c0.sum() == num_nonzero
        assert (np.diagonal(c0) == 0).all()

    else:
        p = num_nonzero / D**2
        c0 = np.random.binomial(1, p=p, size=c.shape)
        np.fill_diagonal(c0, 0)

    return c0


def random_with_p(N, p):
    '''returns a binary (N, N) matrix with connection probability of p on every edge except the diagonal'''
    c = np.random.uniform(size=(N, N)) < p * (N**2) / (N**2 - N)
    np.fill_diagonal(c, 0)
    return c


def random_spatial(i=0, N=1000, p=0.02):
    import pickle
    with open(f'data/random_spatial/random_spatial_N{N}_p{p}_{i:02}.pkl', 'rb') as f:
        x = pickle.load(f)
    return x



###
# helper

def densifier(li, lj):
    #assume min vertex is 0
    N = max(li+lj) + 1
    r = np.zeros((N,N))
    for i,j in zip(li,lj):
        r[i,j] = 1
    return r
