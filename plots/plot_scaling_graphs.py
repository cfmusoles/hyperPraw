## Python helper file to plot and store graphs to represent scaling in distributed simulations
# Graphs produced: as many as columns that want to be compared

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

geometric_scaling = True
min_num_processes = 144
# for linear scaling of processors
max_num_processes = 288
process_step = 32
#for geometric scaling of processors
num_experiments = 1
geometric_step = 2

show_error = True
as_bar_plot = True

# "sat14_E02F20.cnf.hgr" $SEED 8 #Y
# "sat14_itox_vc1130.cnf.dual.hgr" $SEED 2 #Y for esim
# "2cubes_sphere.mtx.hgr" $SEED 3 #Y for esim
# "ABACUS_shell_hd.mtx.hgr" $SEED 40 #Y
# "sparsine.mtx.hgr" $SEED 2 
# "venkat01.mtx.hgr" $SEED 4 #Y

# "pdb1HYS.mtx.hgr" $SEED 1 20 #Y # hedge sim is too short
# "parabolic_fem.mtx.hgr" $SEED 4 1 #N 
# "sat14_10pipe_q0_k.cnf.primal.hgr" $SEED 1 1 #Y
# "sat14_E02F22.cnf.hgr" $SEED 3 1 #Y
# "sat14_openstacks-p30_3.085-SAT.cnf.dual.hgr" $SEED 1 1 #Y
# "webbase-1M.mtx.hgr" $SEED 1 1 #Y
# "sat14_dated-10-17-u.cnf.dual.hgr" $SEED 4 1 #~
# "ship_001.mtx.hgr" $SEED 1 30 #Y # hedge sim is too short

folder = "../results/runtime/"
experiment_name = "runtime"
graph_name = "sat14_E02F20.cnf.hgr"
# each element on the following arrays corresponds to an experiment run (collection of files)
experiments = [experiment_name +  "_zoltan_" + graph_name + "_zoltan",experiment_name + "_default_" + graph_name + "_prawS",experiment_name + "_bandwidth_" + graph_name + "_prawS",experiment_name + "_refinement_" + graph_name + "_prawSref"]
colours = ["red","green","blue","orange"] # as many as the number of experiments included
legend_labels = ['Zoltan','PRAW','PRAW-arc-aware','PRAW-refinement']

# Each element on the following arrays corresponds to a column in columns_to_plot
columns_to_plot = [1,2,4,3,5]#,8,9,10,11]
reference_values = [0,2,1,6,7,8,3,1,1] # used to take values on each column divided by these
use_ref_values = False
scale_plots = [1,1,1,1e-3,1,1,1,1]
plot_title = ["EdgeSim time","HedgeSim time","Edge cut","Hyperedge cut","SOED","Edge comm cost","Hedge comm cost","Messages sent (edge)","Messages sent (hedge)"]
plot_xlabel = ["Number of processes","Number of processes","Number of processes","Number of processes","Number of processes","Number of processes","Number of processes","Number of processes","Number of processes"]
plot_ylabel = ["Time(s)","Time(s)","Cut ratio","Cut ratio","SOED (thousands)","Cost","Cost","Messages sent","Messages sent"]
image_format = 'pdf'
plot_name = ["a_" + str(x) for x in range(len(columns_to_plot))] #["a1","a2","a3","a4","a5","a6","a7"]

bar_plot_size = 0.4 / len(experiments)

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

if geometric_scaling:
	w = num_experiments
	experiment_range = [min_num_processes * geometric_step ** (n-1) for n in range (1, num_experiments+1)]
else:
	w = (max_num_processes - min_num_processes+1) / process_step
	experiment_range = range(min_num_processes, max_num_processes+1,process_step)
	
def get_data_from_csv(filename):
	data = np.genfromtxt(filename,skip_header=1,delimiter=",")
	return data

def plot(x,y, error,title,xlabel,ylabel,name,colour,legend,show,global_counter):
	if as_bar_plot:
		rx = x + global_counter * bar_plot_size*np.array(x)
		plt.bar(rx,y,width=bar_plot_size*np.array(x),color=colour,label=legend)
		
	else:
		if show_error:
			plt.errorbar(x, y, error,linewidth=1,color=colour,label=legend,marker='s',markersize=5)
		else:
			plt.errorbar(x, y,linewidth=1,color=colour,label=legend,marker='s',markersize=5)
	if geometric_scaling:
	#	plt.yscale("log",basey=10)
		plt.xscale("log",basex=10)
	else:
	#	plt.yscale("linear")
		plt.xscale("linear")
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)
	plt.tick_params(axis='x',which='minor',bottom=False,labelbottom=False)
	plt.xticks(experiment_range,experiment_range)
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
		h = len(columns_to_plot)
		means = [[] for y in range(h)]
		stdevs = [[] for y in range(h)]
		
		for p in experiment_range:
			data = get_data_from_csv(folder + experiments[j] + "__" + str(p))
			if data.ndim <= 1:
				data = data.reshape(1,len(data))
			if use_ref_values:
				ref = [row[columns_to_plot[i]] / row[reference_values[i]] for row in data]
				means[i].append(np.mean(ref))
				stdevs[i].append(np.std(ref))
			else:
				means[i].append(np.mean(data[:,columns_to_plot[i]]))
				stdevs[i].append(np.std(data[:,columns_to_plot[i]]))
		# means contains a list of arrays (one per columns analysed) with the average values per column across files
		# stdevs the same but with the corresponding st deviation

		#plot each column separately
		y = [x * scale_plots[i] for x in means[i]]
		error = [x * scale_plots[i] for x in stdevs[i]]
		plot(experiment_range,y,error,plot_title[i],plot_xlabel[i],plot_ylabel[i],plot_name[i],colours[j],legend_labels[j], j == (len(experiments)-1),j)	


    
