## Plot to show hgraph metrics as the streaming partitioning progresses (from history files)
# These plots can help investigate how quality metrics progress as the partitioning goes
# At the moment it can track cummulative hyperedge cut and SOED

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


hgraph_folder = '../resources/synthetic_hgraphs/'
hgraph_file = "small_dense_powerlaw.hgr"
experiment_folder = '../results/ncom_bal/'
experiment_prefix = 'ncom_bal'
experiment_name = ['staggered_overlap_lambda10_parallelVertex','hyperPraw_bandwidth_lambda01','hyperPraw_bandwidth_lambda1']
partitioning = ['parallelVertex','hyperPrawVertex','hyperPrawVertex']
colours = ['blue', 'red','green']
legend_labels = ['baseline','lambda01','lambda1']

num_streams = 1
num_partitions = 96

useStepRelativeMeasurements = True
image_format = 'pdf'
image_names = ["a_" + str(i) for i,v in enumerate(partitioning)]

cuts = []
soeds = []
steps = []
for i,experiment in enumerate(experiment_name):
    # load stream of elements and partitioning
    print("Processing {}...".format(experiment_name[i]))
    hgraph = hgraph_folder + hgraph_file
    part_allocation = np.genfromtxt(experiment_folder + experiment_prefix + '_' + experiment + '_' + str(num_streams) + '_' + 
                                hgraph_file + '_' + partitioning[i] + '_partitioning__' + str(num_partitions),
                                skip_header=0,delimiter=",")

    # load pin degrees from hgraph file
    # each element is a hyperedge, and contains all vertices belonging to it
    with open(hgraph) as f:
        f.readline() #skip header info
        streamed_elements = [x.rstrip().split(" ") for x in f]
        streamed_elements = [list( map(int,i) ) for i in streamed_elements]

    # go through the stream of elements in the same order as when partitioning (based on the number of streams)
    num_elements = len(part_allocation)
    stream_steps = num_elements // num_streams
    if num_elements % num_streams > 0:
        stream_steps += 1
    steps.append(list(range(stream_steps)))
    cuts_per_step = [0 for _ in range(stream_steps)]
    soed_per_step = [0 for _ in range(stream_steps)]
    cummulative_elements = 0
    cummulative_max_soed = 0
    current_cut = 0
    current_soed = 0
    for step in range(stream_steps):
        for st in range(num_streams):
            cummulative_elements += 1
            element_id = step * num_streams + st
            if element_id >= num_elements:
                continue
            # check if the element (hyperedge) is cut
            # it is cut if at least two vertices in it are assigned to different partitions
            parts = set()
            for pin in streamed_elements[element_id]:
                parts.add(part_allocation[pin-1])
            if len(parts) > 1:
                current_cut += 1
            current_soed += len(parts) - 1
            cummulative_max_soed += min(num_partitions,len(streamed_elements[element_id])) - 1
        if useStepRelativeMeasurements:
            cuts_per_step[step] = current_cut/cummulative_elements
            soed_per_step[step] = current_soed/cummulative_max_soed
        else:
            cuts_per_step[step] = current_cut
            soed_per_step[step] = current_soed
                    
    # normalise cut values
    if not useStepRelativeMeasurements:
        cuts_per_step = [c / num_elements for c in cuts_per_step]
    soeds.append(soed_per_step)
    cuts.append(cuts_per_step)



# plot results of all cuts (one per partitioning)
for i in range(len(experiment_name)):
    df = pd.DataFrame({'steps' : steps[i], 'cut': cuts[i]})
    sns.lineplot(x='steps',y='cut',data=df,color=colours[i],label=legend_labels[i])

plt.show()

# plot soeds
for i in range(len(experiment_name)):
    df = pd.DataFrame({'steps' : steps[i], 'soed': soeds[i]})
    sns.lineplot(x='steps',y='soed',data=df,color=colours[i],label=legend_labels[i])

plt.show()       
    