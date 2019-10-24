# Hypergraph generator based on the algorithm proposed in https://github.com/HariniA/EDRW/blob/master/synthetic_dataset_generation.pdf
# and # FENNEL: two models for constructing graphs: 
#       n vertices, k clusters, p_interconnectivity, p_intraconnectivity
#       n vertices, d powerlaw distribution: then use Chung-Lu method to create an instance of the corresponding power law graph
# Paper reference: Extended Discriminative Random Walk: A Hypergraph Approach to Multi-View Multi-Relational Transductive Learning, IJCAI 2015
# We do not care about classes since we are not labelling the vertices
# Based on a power-law distribution (scale free network), generate a specified number of hyperedges from a list of vertices

# evaluate clusterings with Normalized Mutual Information http://dmml.asu.edu/users/xufei/Papers/ICDM2010.pdf or Sillouette 

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

# hypergraph parameters
hypergraph_name = "random_power_law"
export_folder = '../resources/'
num_vertices = 100000
num_hyperedges = 100000
num_clusters = 12
cluster_density = [1.0/num_clusters for _ in range(num_clusters)]       # probability that a vertex belongs to each cluster
p_intraconnectivity = 1.0
hyperedge_gamma = 1.8                                # determines the skewness of the distribution (higher values more skewed to the left). Must be >= 0 (gamma == 0 is the uniform distribution)
max_hyperedge_degree = 100
min_hyperedge_degree = 10
vertex_degree_power_law = True          # whether drawing vertex ids is done using power law distribution (much slower)
vertex_gamma = 1.1
show_distribution = True
store_clustering = True


hypergraph_name = hypergraph_name + ".hgr"

#power law distributions
def truncated_power_law(a, m):
    x = np.arange(1, m+1, dtype='float')
    pmf = 1/x**a
    pmf /= pmf.sum()
    return stats.rv_discrete(values=(range(1, m+1), pmf))

# will sample from 1 to max_hyperedge_degree with probability power law degree
distribution = truncated_power_law(a=hyperedge_gamma, m=max_hyperedge_degree)
#sample = distribution.rvs(size=num_hyperedges)
#plt.hist(sample, bins=np.arange(max_hyperedge_degree)+0.5)
#plt.show()

# assign vertices to clusters (randomly)
print("Generating {} clusters".format(num_clusters))
vertices = set(range(num_vertices))
nodes = np.random.binomial(num_vertices,p=cluster_density,size=num_clusters)
clusters = [[] for _ in range(num_clusters)]
if store_clustering:
        partitioning = np.array([0 for _ in range(num_vertices)])
for k in range(num_clusters):
        # how many vertices will belong to this cluster
        n = min(nodes[k],len(vertices))
        selected = np.random.choice(list(vertices),size=n,replace=False)
        vertices = vertices.difference(set(selected)) #vertices = [v for v in vertices if v not in selected]
        clusters[k].extend(selected.tolist())
        if store_clustering:
                partitioning[selected] = k
# assign remaining vertices randomly
for v in vertices:
        k = np.random.randint(0,num_clusters)
        clusters[k].append(v)

# store clustering
if store_clustering:
        print("Store clustering in file...")
        with open(hypergraph_name + '_clustering',"w") as f:
                for v in partitioning:
                        f.write(str(v) + '\n')
        partitioning = []

# precalculate power law distributions per cluster
if vertex_degree_power_law:
        print("Calculating power law distributions")
        v_dist_per_cluster = []
        for k in range(len(clusters)):
                v_dist_per_cluster.append(truncated_power_law(a=vertex_gamma, m=len(clusters[k])))                 

# build hypergraph one hyperedge at a time
hypergraph = []
print("Generating {} hyperedges".format(num_hyperedges))
with open(export_folder + hypergraph_name,"w+") as f:
        # write header (NUM_HYPEREDGES NUM_VERTICES)
        f.write("{} {}\n".format(num_hyperedges,num_vertices))
        # draw hedge degree from power law distribution (with specific min value)
        hyperedge_degrees = (distribution.rvs(size=num_hyperedges)) + min_hyperedge_degree-1
                
        for he_id in range(num_hyperedges):
                # draw vertices ids randomly
                # needs to account for premade clusters
                # choose a cluster based on cluster density
                k = np.random.choice(num_clusters,p=cluster_density)
                # whether each vertex should be local or not
                nodes = [p_intraconnectivity > np.random.random_sample() for _ in range(hyperedge_degrees[he_id])]
                # select vertices from cluster k if inter_node[x] is True, from any other cluster otherwise
                vertices = set()
                for p in nodes:
                        c = k
                        if not p:
                                # different cluster to k
                                while k == c:
                                        c = np.random.randint(0,num_clusters)
                        while True:
                                if not vertex_degree_power_law:
                                        rand_index = np.random.randint(0,len(clusters[c]))
                                        vertex = clusters[c][rand_index] 
                                else:
                                        # select vertices based on power law distribution too
                                        rand_index = v_dist_per_cluster[c].rvs() - 1
                                        vertex = clusters[c][rand_index]
                                # ensure the same vertex is not added twice to a hyperedge
                                if vertex not in vertices:
                                        break
                        vertices.add(vertex) 
                
                #write hyperedge to file
                vertices = list(vertices)
                vertices = [str(v+1) for v in sorted(vertices)]
                f.write(" ".join(vertices))
                f.write("\n")
                
                # for tracking distribution
                if show_distribution:
                        hypergraph.append(len(vertices))

if show_distribution:
        plt.hist(hypergraph, bins=np.arange(max_hyperedge_degree)+0.5)
        plt.show()





