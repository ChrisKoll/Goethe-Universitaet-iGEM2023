import numpy as np
import logomaker
import seaborn as sns
import matplotlib.pyplot as plt
from graph_algorithms import *

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

def count_occ(kmers_dict, kmers):
    occs = list()
    for kmer in kmers:
        occs.extend(kmers_dict[kmer])
    return check_occs(occs)

def check_occs(occs:list):
    duplicates = set()
    for i in range(len(occs)-1):
        for j in range(1, len(occs)):
            if occs[i][1] > occs[j][0] > occs[i][0] or occs[j][1] > occs[i][0] > occs[j][0]:
                duplicates.add(i)
                duplicates.add(j)

    dups_sorted = sorted(duplicates, reverse=True)
    for i in dups_sorted:
        occs.pop(i)

    return occs

kmers = dict()
with open("PromoterOnly.fasta", "r") as file:
    file.readline()
    sequence = file.readline()

    for i in range(len(sequence) - WINDOWSIZE):
        try:
            kmers[sequence[i:i + WINDOWSIZE]].append((i, i + WINDOWSIZE))
        except KeyError:
            kmers[sequence[i:i + WINDOWSIZE]] = [(i, i + WINDOWSIZE)]

######################
# Cliquen Clustering #
######################

dist = np.zeros((len(kmers.keys()), len(kmers.keys())))
kmers_list = list()
for i, key_i in enumerate(kmers.keys()):
    kmers_list.append(key_i)
    for j, key_j in enumerate(kmers.keys()):
        dist[i, j] = hamming_dist(key_i, key_j)

dist = (dist + 1) * (1 - np.eye(dist.shape[0]))
sns.clustermap(dist)
plt.show()

max_3_graph = bottelneck(dist, 2)
sns.clustermap(max_3_graph)
plt.show()

clusters = minimum_clique_partitioning(max_3_graph)
kmers_list = np.array(kmers_list)

cluster_counts = dict()
for i in range(len(clusters)):
    c_kmers = kmers_list[np.array(list(clusters[i])).astype(int)]
    counts = count_occ(kmers, c_kmers)
    cluster_counts["+".join(c_kmers)] = len(counts)

#################
# Visualization #
#################

# Plot kmers (unique)
kmers_sorted = dict(sorted(kmers.items(), key=lambda item: len(item), reverse=True))
keys = list(kmers_sorted.keys())
values = [len(occ) for occ in list(kmers_sorted.values())]
plt.bar(keys[0:14], values[0:14], width=0.75)
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()

# Plot clusters
clusters_counts = dict(sorted(cluster_counts.items(), key=lambda item: item[1], reverse=True))
keys = list(clusters_counts.keys())
values = list(clusters_counts.values())
plt.bar(keys[0:14], values[0:14], width=0.75)
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()
