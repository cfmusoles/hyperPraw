# Hypergraph generator based on the algorithm proposed in https://github.com/HariniA/EDRW/blob/master/synthetic_dataset_generation.pdf
# and # FENNEL: two models for constructing graphs: 
#       n vertices, k clusters, p_interconnectivity, p_intraconnectivity
#       n vertices, d powerlaw distribution: then use Chung-Lu method to create an instance of the corresponding power law graph
# Paper reference: Extended Discriminative Random Walk: A Hypergraph Approach to Multi-View Multi-Relational Transductive Learning, IJCAI 2015
# We do not care about classes since we are not labelling the vertices
# Based on a power-law distribution (scale free network), generate a specified number of hyperedges from a list of vertices

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

# hypergraph parameters
hypergraph_name = "random_power_law"
num_vertices = 10000
num_hyperedges = 10000
num_clusters = 10
cluster_density = [1.0/num_clusters for _ in range(num_clusters)]       # probability that a vertex belongs to each cluster
p_intraconnectivity = 0.95
gamma = 1.8   # determines the skewness of the distribution (higher values more skewed to the left). Must be >= 0 (gamma == 0 is the uniform distribution)
max_hyperedge_degree = 30
min_hyperedge_degree = 2
show_distribution = True


hypergraph_name = hypergraph_name + ".hgr"

#power law distributions
def truncated_power_law(a, m):
    x = np.arange(1, m+1, dtype='float')
    pmf = 1/x**a
    pmf /= pmf.sum()
    return stats.rv_discrete(values=(range(1, m+1), pmf))

# will sample from 1 to max_hyperedge_degree with probability power law degree
distribution = truncated_power_law(a=gamma, m=max_hyperedge_degree)
#sample = distribution.rvs(size=num_hyperedges)
#plt.hist(sample, bins=np.arange(max_hyperedge_degree)+0.5)
#plt.show()

# assign vertices to clusters (randomly)
clusters = [[] for _ in range(num_clusters)]
vertices = list(range(num_vertices))
for k in range(num_clusters):
        # how many vertices will belong to this cluster
        nodes = min(np.random.binomial(num_vertices,cluster_density[k]),len(vertices))
        selected = np.random.choice(vertices,nodes,replace=False)
        vertices = [v for v in vertices if v not in selected]
        clusters[k].extend(selected.tolist())
# assign remaining vertices randomly
for v in vertices:
        k = np.random.randint(0,num_clusters)
        clusters[k].append(v)


# build hypergraph one hyperedge at a time
hypergraph = []

with open(hypergraph_name,"w+") as f:
        # write header (NUM_HYPEREDGES NUM_VERTICES)
        f.write("{} {}\n".format(num_hyperedges,num_vertices))

        for he_id in range(num_hyperedges):
                # draw hedge degree from power law distribution (with specific min value)
                hyperedge_degree = (distribution.rvs(size=1))[0] + min_hyperedge_degree-1
                # draw vertices ids randomly
                # needs to account for premade clusters
                # choose a cluster based on cluster density
                k = np.random.choice(num_clusters,p=cluster_density)
                # whether each vertex should be local or not
                nodes = [p_intraconnectivity > np.random.random_sample() for _ in range(hyperedge_degree)]
                # select vertices from cluster k if inter_node[x] is True, from any other cluster otherwise
                vertices = []
                for p in nodes:
                        c = k
                        if not p:
                                # different cluster to k
                                while k == c:
                                        c = np.random.randint(0,num_clusters)
                        vertex = np.random.choice(clusters[c])
                        vertices.append(vertex)                             
                
                #vertices = np.random.choice(num_vertices,hyperedge_degree,replace=False)

                #write hyperedge to file
                vertices = [str(v+1) for v in sorted(vertices)]
                f.write(" ".join(vertices))
                f.write("\n")
                
                # for tracking distribution
                if show_distribution:
                        hypergraph.append(len(vertices))

if show_distribution:
        plt.hist(hypergraph, bins=np.arange(max_hyperedge_degree)+0.5)
        plt.show()





