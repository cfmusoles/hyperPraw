# Utility script to show the memory limits of hypergraph partitioning algorithms:
#   - Zoltan (global)
#   - HyperPRAW (streaming)

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Parameters
# CM: 80000 neurons, 300000000 synapses, MVC: 4130000 neurons, 24200000000 synapses
num_vertex = 4130000
num_hedges = 4130000
num_pins = 24200000000
num_partitions = [768, 1536, 3072, 6144, 12288, 12288*2]
# hyperPRAW parameters
num_streams = [768, 1536, 3072, 6144, 12288, 12288*2]
avg_hedge_replication_factor = [0.01, 0.05, 0.1, 0.5] # fraction of total partitions
max_hedge_replication_factor = num_pins / num_vertex
shared_mem_node_size = [768, 1536, 3072, 6144, 12288, 12288*2]
############

# graph details
as_bar_plot = False
geometric_scaling = True
show_title = True
zoltan_colour = "blue"
zoltan_linestyle = "-"
zoltan_legend = "zoltan"
colours = ["springgreen","mediumseagreen","green", "darkgreen"]
linestyles = ["--","--","-.","-"]
legend_labels = ['HyperPRAW (0.01)','HyperPRAW (0.05)','HyperPRAW (0.1)','HyperPRAW (0.5)']
plot_title = "Memory requirements to partition MVC"
plot_xlabel = "Number of partitions"
plot_ylabel = "Memory (GB)"
image_format = 'pdf'
plot_name = "partitioning_mem_requirements"

bar_plot_size = 0.8 / len(num_partitions)

# CONSTANTS
mem_scale_factor = 1e-9
double_size = 8 # size of doubles in memory
single_size = 4 # size of ints in memory
half_size = 2 # size of short in memory


def get_zoltan_limit(num_vertex, num_pins, num_partitions):
    '''
    Calculates the memory requirements to run a specific hypergraph through zoltan (Vertex compression)

    Zoltan implemented in ZoltanCompressedVertexPartitioning.h

    Assumes basic memory optimisations (seen_pins table is shared amongst processes belonging to the same node)
    This can be further improved by sharing once amongst all processes with RDMA.
    '''
    # global memory
    # double type structures
    global_limit = double_size * (num_vertex * 2 + num_pins)
    # single type structures
    global_limit += single_size * num_vertex

    # per process memory
    per_process_limit = single_size * 3 + double_size * num_vertex / num_partitions

    return global_limit + num_partitions * per_process_limit


def get_hyperPRAW_limit(num_vertex, num_hedges, num_streams, num_partitions, avg_hedge_replication_factor, node_size):
    '''
    Calculates the memory requirements to run a specific hypergraph through hyperPRAW

    HyperPRAW implemented in PRAW.h:parallelHDRF

    Assumes basic memory optimisations. The following data structures are shared amongst processes belonging to same node
    
    - seen_pins table (done)
    - comm cost table 
    - partitioning
    '''
    # cap max replication factor (to the max pins per vertex)
    replication = min(avg_hedge_replication_factor * num_partitions, max_hedge_replication_factor)
    replication_factor = replication / num_partitions
    # global memory
    # dominant at low scales (low number of streams / partitions)
    seen_table = num_hedges * (single_size + half_size * (num_partitions * replication_factor * 2))
    # apply memory optimisation to share central datastructure
    global_limit = seen_table * max(num_streams / node_size, 1)

    # per stream memory
    # this becomes dominant at large scales
    # the storage of bandwidth cost becomes dominant
    bandwidth_cost_matrix = (single_size * num_partitions * num_partitions) / 2
    stream_datastructures = single_size * (num_vertex + num_hedges) + double_size * num_partitions + double_size * num_vertex + bandwidth_cost_matrix + single_size * 3 * num_streams
    # apply moemory optimisations to share partitioning datastructures
    streams_limit = stream_datastructures * max(num_streams / node_size, 1)

    return global_limit + streams_limit

def plot(x,y,title,xlabel,ylabel,name,colour,legend,show,global_counter,linestyle='-'):
	if as_bar_plot:
		rx = x + global_counter * bar_plot_size*np.array(x)
		plt.bar(rx,y,width=bar_plot_size*np.array(x),color=colour,label=legend)
		
	else:
		plt.errorbar(x, y,linewidth=1,color=colour,label=legend,marker='s',markersize=5,ls=linestyle)
	if geometric_scaling:
		#plt.yscale("log",basey=10)
		plt.yscale("linear")
		plt.xscale("log",basex=10)
	else:
	#	plt.yscale("linear")
		plt.xscale("linear")
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	if show_title:
		plt.title(title)
	plt.tick_params(axis='x',which='minor',bottom=False,labelbottom=False)
	plt.xticks(num_partitions,num_partitions)
	#plt.tight_layout()
	plt.gcf().subplots_adjust(left=0.17)
	plt.legend(loc='best')
	if show:
		plt.savefig(name + "." + image_format,format=image_format,dpi=1000)
		plt.show()


# plot zoltan vs hyperPRAW memory requirements as problem scales:
#   - increased parallelism (greater number of processes and streams)
#   - increased hgraph size
# use best, average and worse case for HyperPRAW (controlling average hedge replication factor)


plt.figure()
#zoltan
zoltans = []
for i, partitions in enumerate(num_partitions):
    zoltan_mem_req = get_zoltan_limit(num_vertex, num_pins, partitions) * mem_scale_factor
    zoltans.append(zoltan_mem_req)
plot(num_partitions,zoltans,plot_title,plot_xlabel,plot_ylabel,plot_name,zoltan_colour,zoltan_legend, False,0,linestyle=zoltan_linestyle)

#hyperPRAWs
for j, replication_factor in enumerate(avg_hedge_replication_factor):
    hyperPRAWs = []
    for i, partitions in enumerate(num_partitions):
        hyperPRAW_mem_req = get_hyperPRAW_limit(num_vertex, num_hedges, num_streams[i], partitions, replication_factor, shared_mem_node_size[i]) * mem_scale_factor
        hyperPRAWs.append(hyperPRAW_mem_req)
    plot(num_partitions,hyperPRAWs,plot_title,plot_xlabel,plot_ylabel,plot_name,colours[j],legend_labels[j], j == (len(avg_hedge_replication_factor)-1),1+j,linestyle=linestyles[j])



#print("{0:.0f}: Zoltan memory requirements".format(get_zoltan_limit(num_vertex, num_pins, num_partitions)))

#print("{0:.0f}: HyperPRAW memory requirements".format(get_hyperPRAW_limit(num_vertex, num_hedges, num_streams, num_partitions, avg_hedge_replication_factor, shared_mem_node_size)))

