## Plot element degree (vertex degree for edge partitioning, hedge degree for vertex partitioning) vs whether or not the element has been replicated (i.e. partitioned)
## Very relevant to demonstrate that HDRF actually prefers to replicate elements with high degree

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# Todo: extend it to produce the plot for multiple graphs at once

hgraph_file = '../resources/sparsine.mtx.hgr'
partitioning_file = '../test_default_sparsine.mtx.hgr_zoltanVertex_partitioning__12'

storePlot = False
image_format = 'pdf'

#load partitioning scheme
partitioning = np.genfromtxt(partitioning_file,skip_header=0,delimiter=",")
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
    plt.savefig("a"+ image_format,format=image_format,dpi=1000)
plt.show()