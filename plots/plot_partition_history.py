## Python helper file to plot and store graphs to represent partition history resulting from praw partitioning
# Graphs produced: as many as columns that want to be compared

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

num_processes = 12

as_bar_plot = False

folder = "../"
experiments = ["test_default","test_bandwidth"] 
graph_name = "venkat01.mtx.hgr"
# each element on the following arrays corresponds to an experiment run (collection of files)
colours = ["red","green"] # as many as the number of experiments included
legend_labels = ['Default PRAW','PRAW-bandwidth']

# Each element on the following arrays corresponds to a column in columns_to_plot
columns_to_plot = [0,1,2]
plot_title = ["Workload imbalance","Hedge cut","Edge cut"]
plot_xlabel = ["Iteration","Iteration","Iteration"]
plot_ylabel = ["Imbalance","Hyperedgecut ratio","Edgecut ratio"]
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
		plt.errorbar(x, y,linewidth=1,color=colour,label=legend,marker='s',markersize=5)

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


    
