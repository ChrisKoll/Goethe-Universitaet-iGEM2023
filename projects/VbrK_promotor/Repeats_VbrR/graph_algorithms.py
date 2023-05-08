import numpy as np


#  G represents a Graph - here an edge matrix (size: vertices*vertices)

def bottelneck(G, w):
    new_G = G
    new_G[new_G > w] = 0
    return new_G


def MIS(G):
    I = set()
    V = set(np.arange(0, G.shape[0], 1))
    deg = np.sum(G > 0, axis=0)
    idx = np.arange(0, G.shape[0], 1)
    while len(V) > 0:
        v = int(np.argmin(deg))
        if v in V:
            deg[v] = 10000
        else:
            break
        I.add(v)
        V.remove(v)
        N = idx[G[v, :] > 0]
        for n in N:
            V.remove(n)

    return I


def min_max_k_clustering(G, w, k):
    i = 0
    cont = True
    I = set()

    while cont:
        i = i + 1
        Gb = bottelneck(G, w[i])
        I = MIS(Gb)
        if I <= k:
            cont = False

    return I


def minimum_clique_partitioning(G):
    i = 0
    Q = set()  # set of clustered vertices
    V = set(np.arange(0, G.shape[0], 1))
    C = list()

    while len(V.difference(Q)) > 0:
        C.append(set())
        Vp = V.difference(Q)
        deg = np.sum(G > 0, axis=0)
        idx = np.arange(0, G.shape[0], 1)

        while len(Vp) > 0:
            v = int(np.argmin(deg))
            deg[v] = 10000
            if v not in Vp:
                continue

            C[i].add(v)
            Vp.remove(v)
            N = set(idx[G[v, :] > 0])
            Vp = Vp.intersection(N)

        Q = Q.union(C[i])
        i += 1

    return C
