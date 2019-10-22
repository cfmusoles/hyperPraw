import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import collections

folder = "../results/stream/"
hgraph_files = ["2cubes_sphere.mtx.hgr","ABACUS_shell_hd.mtx.hgr","sparsine.mtx.hgr","sat14_10pipe_q0_k.cnf.primal.hgr","webbase-1M.mtx.hgr"]# ["atmosmodj.mtx.hgr","kkt_power.mtx.hgr"]
experiment_prefix = 'stream_rHDRF_default_4'
partition = 'rHDRF'
num_processes = 72

storePlot = False
image_format = 'pdf'
image_names = [partition + "_" + str(i) for i,v in enumerate(hgraph_files)]

for i,hgraph in enumerate(hgraph_files):
    #load partitioning scheme
    partitioning = np.genfromtxt(folder + experiment_prefix + '_' + hgraph + '_' + partition + '_partitioning__' + str(num_processes),skip_header=0,delimiter=",")
    partitioning = [int(i)  for i in partitioning]

    counts = collections.Counter(partitioning)
    df = pd.DataFrame(list(counts.items()),columns=['Partition','Elements'],index=counts.keys())

    #for key in counts.keys():
    #    print("Partition {} has {} elements".format(key,counts[key]))
    
    fig, ax = plt.subplots()
    df.boxplot(column='Elements',ax=ax)
    
    plt.show()

    