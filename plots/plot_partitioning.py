import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import collections

folder = "../results/baseline/"
hgraph_files = ["atmosmodj.mtx.hgr","kkt_power.mtx.hgr"]
experiment_prefix = 'baseline_rHDRF_bandwidth_1_'
partition = 'rHDRF'
num_processes = 144

storePlot = False
image_format = 'pdf'
image_names = ["rHDRF_" + str(i) for i,v in enumerate(hgraph_files)]

for i,hgraph in enumerate(hgraph_files):
    #load partitioning scheme
    partitioning = np.genfromtxt(folder + experiment_prefix + hgraph + '_' + partition + '_partitioning__' + str(num_processes),skip_header=0,delimiter=",")
    partitioning = [int(i)  for i in partitioning]

    counts = collections.Counter(partitioning)
    df = pd.DataFrame(list(counts.items()),columns=['Partition','Elements'],index=counts.keys())

    #for key in counts.keys():
    #    print("Partition {} has {} elements".format(key,counts[key]))
    
    fig, ax = plt.subplots()
    df.boxplot(column='Elements',ax=ax)
    
    plt.show()

    