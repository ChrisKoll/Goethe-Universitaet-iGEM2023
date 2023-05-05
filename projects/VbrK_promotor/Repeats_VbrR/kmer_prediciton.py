import numpy as np
import matplotlib.pyplot as plt

WINDOWSIZE = 6  # as described in paper


def hamming_dist(x, y):
    if len(x) != len(y):
        return np.NAN
    else:
        res = 0
        for i in range(len(x)):
            if ord(x[i]) != ord(y[i]):
                res += 1
        return res


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

# counting kmers using Hash table
kmers = dict()
with open("PromoterOnly.fasta", "r") as file:
    file.readline()
    sequence = file.readline()

    for i in range(len(sequence) - WINDOWSIZE):
        try:
            kmers[sequence[i:i + WINDOWSIZE]] += 1
        except KeyError:
            kmers[sequence[i:i + WINDOWSIZE]] = 1

# Calculate distance matrix / network
dist = np.zeros((len(kmers.keys()), len(kmers.keys())))
kmers_list = list()
for i, key_i in enumerate(kmers.keys()):
    kmers_list.append(key_i)
    for j, key_j in enumerate(kmers.keys()):
        dist[i, j] = hamming_dist(key_i, key_j)
        

# Plot results
dist = (dist + 1) * (1 - np.eye(dist.shape[0]))
kmers_list = np.array(kmers_list)
maxindset = np.array(list(MIS(bottelneck(dist, 2))))
print(kmers_list[maxindset])

kmers_sorted = dict(sorted(kmers.items(), key=lambda item: item[1], reverse=True))
print(kmers_sorted)
keys = list(kmers_sorted.keys())
values = list(kmers_sorted.values())
plt.bar(keys[0:14], values[0:14], width=0.75)
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()
