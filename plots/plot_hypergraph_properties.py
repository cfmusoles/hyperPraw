import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math

folder = '../resources/'
hgraphs = ["small_dense_powerlaw.hgr","small_dense_uniform.hgr","large_sparse_powerlaw.hgr","large_sparse_uniform.hgr"]
        #['sat14_itox_vc1130.cnf.dual.hgr','2cubes_sphere.mtx.hgr','ABACUS_shell_hd.mtx.hgr','sparsine.mtx.hgr',
            #'pdb1HYS.mtx.hgr','sat14_10pipe_q0_k.cnf.primal.hgr','sat14_E02F22.cnf.hgr','webbase-1M.mtx.hgr',
            #'ship_001.mtx.hgr','sat14_atco_enc1_opt1_05_21.cnf.dual.hgr']

image_format = 'pdf'

plot_graphs = True

for hgraph_file in hgraphs:
    # data structure to hold hedge sizes: 1d list corresponding to each hedge size
    # data structure to hold vertex degree: 1d list corresponding to each vertex degree
    print('\nHypergraph: ' + hgraph_file)
    with open(folder + hgraph_file) as f:
        print('Hyperedges - Vertices')
        print(f.readline()) #skip header info
        hedges = [x.rstrip().split(" ") for x in f]
        hedge_sizes = [len(x) for x in hedges]
        print('Total NNZ')
        print(np.sum(hedge_sizes))
        print('Average cardinality')
        print(np.mean(hedge_sizes))
        hedges = [int(vid.replace("\n","")) for hedge in hedges for vid in hedge] #np.array([int(x.replace("\n","")) for x in hedges]).flatten()
        unique, vertex_counts = np.unique(hedges,return_counts=True)
        print('Number of non-zero-neighbour vertices')
        print(len(unique))

    if plot_graphs:
        # plotting size of hyperedges histogram
        mybins=np.linspace(min(hedge_sizes)-1,max(hedge_sizes)+1,num=max(hedge_sizes)-min(hedge_sizes)+2)
        g = sns.distplot(hedge_sizes,kde=False,bins=mybins,hist_kws={'edgecolor':'black'})
        g.set_xscale('linear')
        g.set_yscale('linear')
        g.set_title("Hyperedges sizes")
        g.set_xlabel("Size of hyperedge")
        g.set_ylabel("Count")
        plt.savefig(hgraph_file + "_hedge_size." + image_format,format=image_format,dpi=1000)
        plt.show()


        # plotting vertex degrees histogram
        mybins=np.linspace(min(vertex_counts)-1,max(vertex_counts)+1,num=max(vertex_counts)-min(vertex_counts)+2)
        #mybins=np.logspace(1,math.log10(max(vertex_counts)),num=max(vertex_counts)-min(vertex_counts))
        g = sns.distplot(vertex_counts,kde=False,bins=mybins,hist_kws={'edgecolor':'black'})
        g.set_xscale('linear')
        g.set_yscale('linear')
        g.set_title("Vertex degrees")
        g.set_xlabel("Vertex degree")
        g.set_ylabel("Count")
        plt.savefig(hgraph_file + "_vtx_degree." + image_format,format=image_format,dpi=1000)
        plt.show()

