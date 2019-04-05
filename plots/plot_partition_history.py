## Python helper file to plot and store graphs to represent partition history resulting from praw partitioning
# Graphs produced: as many as columns that want to be compared

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

num_processes = 144

as_bar_plot = False

folder = "../results/stop_cnd/"
experiments = ["stop_cnd_refine_1700","stop_cnd_refine_1000","stop_cnd_refine_950","stop_cnd_hard"] 
graph_name = "sat14_itox_vc1130.cnf.dual.hgr"
# each element on the following arrays corresponds to an experiment run (collection of files)
colours = ["green","orange","blue","red"] # as many as the number of experiments included
legend_labels = ['Refinement 1.7','Refinement 1.0','Refinement 0.95','Imbalance reached']

# Each element on the following arrays corresponds to a column in columns_to_plot
columns_to_plot = [1,2,5]#,0]
plot_title = ["Hedge cut","Edge cut","Communication cost","Workload imbalance"]
plot_xlabel = ["Iteration","Iteration","Iteration","Iteration"]
plot_ylabel = ["Hyperedgecut ratio","Edgecut ratio","Edge comm cost","Imbalance"]
image_format = 'pdf'
plot_name = ["a_" + str(x) for x in range(len(columns_to_plot))] #["a1","a2","a3","a4","a5","a6","a7"]

bar_plot_size = 0.5 / len(experiments)

# general plot settings
plt.rcParams['figure.facecolor'] = 'white'
fig_settings = {  
        'lines.linewidth': 0.5,
        'axes.linewidth': 0.5,
        'axes.labelsize': 'small',
        'legend.fontsize': 'small',
        'font.size': 14,
        'savefig.dpi': 200,
}
plt.rcParams.update(fig_settings)

def get_data_from_csv(filename):
	data = np.genfromtxt(filename,skip_header=1,delimiter=",")
	return data

def plot(x,y, title,xlabel,ylabel,name,colour,legend,show,global_counter):
	if as_bar_plot:
		rx = x + global_counter * bar_plot_size*np.array(x)
		plt.bar(rx,y,width=bar_plot_size*np.array(x),color=colour,label=legend)
	else:	
		plt.errorbar(x, y,linestyle="-",linewidth=1,color=colour,label=legend,marker='s',markersize=2)

	#plt.yscale("linear")
	#plt.xscale("linear")
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)
	plt.tick_params(axis='x',which='minor',bottom=False,labelbottom=False)
	plt.xticks(x,x)
	#plt.tight_layout()
	plt.gcf().subplots_adjust(left=0.17)
	if len(experiments) > 1:
		plt.legend(loc='best')
	if show:
		plt.savefig(name + "." + image_format,format=image_format,dpi=1000)
		plt.show()


# one plot per column
for i in range(len(columns_to_plot)):
	# creating figure for the column plot
	plt.figure()
	for j in range(len(experiments)):
		file_name = experiments[j] + "_" + graph_name + "_partition_history___" + str(num_processes)
		data = get_data_from_csv(folder + file_name)
		y = data[:,columns_to_plot[i]]
		x = range(1,len(y)+1)
		
		#plot each column separately
		plot(x,y,plot_title[i],plot_xlabel[i],plot_ylabel[i],plot_name[i],colours[j],legend_labels[j], j == (len(experiments)-1),j)	


    
