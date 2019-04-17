## Python helper file to plot and store graphs to represent partition history resulting from praw partitioning
# Graphs produced: as many as columns that want to be compared

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

num_processes = 144

as_bar_plot = False
show_title = False

# "sat14_itox_vc1130.cnf.dual.hgr" $SEED 2 #Y for esim
# "2cubes_sphere.mtx.hgr" $SEED 3 #Y for esim
# "ABACUS_shell_hd.mtx.hgr" $SEED 40 #Y
# "sparsine.mtx.hgr" $SEED 2 


folder = "../results/stop_cnd/"
experiments = ["stop_cnd_refine_1000","stop_cnd_refine_950","stop_cnd_hard"] 
graph_name = "sparsine.mtx.hgr"
# each element on the following arrays corresponds to an experiment run (collection of files)
colours = ["red","blue","black"] # as many as the number of experiments included
linestyles = ["--",":","-"]
legend_labels = ['Refinement 1.0','Refinement 0.95','No refinement']

# Each element on the following arrays corresponds to a column in columns_to_plot
columns_to_plot = [5]#,1,2,0]
scale = [1e6,1,1,1]
plot_title = ["Partitioning comm cost","Hedge cut","Edge cut","Workload imbalance"]
plot_xlabel = ["Iteration","Iteration","Iteration","Iteration"]
plot_ylabel = ["Partitioning comm cost (millions)","Hyperedgecut ratio","Edgecut ratio","Imbalance"]
image_format = 'pdf'
plot_name = ["a_" + str(x) for x in range(len(columns_to_plot))] #["a1","a2","a3","a4","a5","a6","a7"]

bar_plot_size = 0.5 / len(experiments)

max_x_value = 0

# general plot settings
plt.rcParams['figure.facecolor'] = 'white'
fig_settings = {  
        'lines.linewidth': 0.5,
        'axes.linewidth': 0.5,
        'axes.labelsize': 'medium',
        'legend.fontsize': 'medium',
        'font.size': 16,
        'savefig.dpi': 200,
}
plt.rcParams.update(fig_settings)

def get_data_from_csv(filename):
	data = np.genfromtxt(filename,skip_header=1,delimiter=",")
	return data

def plot(x,y, title,xlabel,ylabel,name,colour,linestyle,legend,show,global_counter):
	global max_x_value
	if as_bar_plot:
		rx = x + global_counter * bar_plot_size*np.array(x)
		plt.bar(rx,y,width=bar_plot_size*np.array(x),color=colour,label=legend)
	else:	
		plt.errorbar(x, y,linestyle=linestyle,linewidth=1,color=colour,label=legend,marker='s',markersize=2)

	#plt.yscale("linear")
	#plt.xscale("linear")
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	if show_title:
		plt.title(title)
	plt.tick_params(axis='x',which='minor',bottom=False,labelbottom=False)
	if len(x) > max_x_value:
		max_x_value = len(x)
		xtick = [n if n % 10 == 0 else '' for n in x]
		plt.xticks(x,xtick)
	#plt.tight_layout()
	plt.gcf().subplots_adjust(left=0.17,bottom=0.15)
	if len(experiments) > 1:
		plt.legend(loc='best')
	if show:
		plt.savefig(name + "." + image_format,format=image_format,dpi=200)
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
		y = y/scale[i]
		#plot each column separately
		plot(x,y,plot_title[i],plot_xlabel[i],plot_ylabel[i],plot_name[i],colours[j],linestyles[j],legend_labels[j], j == (len(experiments)-1),j)	


    
