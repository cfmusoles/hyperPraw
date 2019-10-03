## Plot element degree (vertex degree for edge partitioning, hedge degree for vertex partitioning) vs whether or not the element has been replicated (i.e. partitioned)
## Very relevant to demonstrate that HDRF actually prefers to replicate elements with high degree

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# Todo: extend it to produce the plot for multiple graphs at once

hgraphs_folder = '../resources/'
hgraph_files = ['sparsine.mtx.hgr','2cubes_sphere.mtx.hgr']
experiment_prefix = '../test_default_'
partition = 'zoltanVertex'
num_processes = 12

storePlot = False
image_format = 'pdf'
image_names = ["a" + str(i) for i,v in enumerate(hgraph_files)]

for i,hgraph in enumerate(hgraph_files):
    hgraph_file = hgraphs_folder + hgraph

    #load partitioning scheme
    partitioning = np.genfromtxt(experiment_prefix + hgraph + '_' + partition + '_partitioning__' + str(num_processes),skip_header=0,delimiter=",")
    partitioning = [int(i)  for i in partitioning]

    # load pin degrees from hgraph file
    with open(hgraph_file) as f:
        f.readline() #skip header info
        elements = [x.rstrip().split(" ") for x in f]
        elements = [list( map(int,i) ) for i in elements]
        pin_degrees = [len(x) for x in elements]

    # create pandas frame with two columns:
    #   pin degree
    #   whether it is replicated
    replicated = []
    for element in elements:
        parts = set()
        for pins in element:
            parts.add(partitioning[pins-1])
        replicated.append(len(parts) > 1)


    df = pd.DataFrame({'Pin degrees' : pin_degrees, 'Replicated' : replicated})

    sns.catplot( x="Replicated", y="Pin degrees", data=df, legend=False, kind='violin')

    if storePlot:
        plt.savefig(image_names[i]+ image_format,format=image_format,dpi=1000)
    plt.show()