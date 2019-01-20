import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

hgraph_file = "../resources/sat14_11pipe_q0_k.cnf.hgr" #2D_54019_highK.mtx.hgr #sparsine.mtx.hgr sat14_11pipe_q0_k.cnf.hgr

# data structure to hold hedge sizes: 1d list corresponding to each hedge size
# data structure to hold vertex degree: 1d list corresponding to each vertex degree
with open(hgraph_file) as f:
    print(f.readline()) #skip header info
    hedges = [x.split(" ") for x in f]
    hedge_sizes = [len(x) for x in hedges]
    hedges = [int(vid.replace("\n","")) for hedge in hedges for vid in hedge] #np.array([int(x.replace("\n","")) for x in hedges]).flatten()
    unique, vertex_counts = np.unique(hedges,return_counts=True)

# plotting size of hyperedges histogram
mybins=np.logspace(0,np.log10(max(hedge_sizes)),num=50)
g = sns.distplot(hedge_sizes,kde=False,bins=mybins,hist_kws={'edgecolor':'black'})
g.set_xscale('log')
g.set_yscale('log')
g.set_title("Hyperedges sizes")
g.set_xlabel("Size of hyperedge")
g.set_ylabel("Count")
plt.show()


# plotting vertex degrees histogram
mybins=np.logspace(0,np.log10(len(unique)),num=100)
g = sns.distplot(vertex_counts,kde=False,bins=mybins,hist_kws={'edgecolor':'black'})
g.set_xscale('log')
g.set_yscale('log')
g.set_title("Vertex degrees")
g.set_xlabel("Vertex degree")
g.set_ylabel("Count")
plt.show()

