## Python helper file to plot and store graphs to represent scaling in distributed simulations
# Graphs produced: as many as columns that want to be compared

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

num_processes = 576

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
# "webbase-1M.mtx.hgr" $SEED 1 1 #Y
# "sat14_dated-10-17-u.cnf.dual.hgr" $SEED 4 1 #~
# "ship_001.mtx.hgr" $SEED 1 30 #Y # hedge sim is too short

folder = "../results/final/"
experiment_name = "final"
graph_name = "ship_001.mtx.hgr"
columns = [1,3,4,5,6,7]
# each element on the following arrays corresponds to an experiment run (collection of files)
experiments = [experiment_name +  "_zoltan_" + graph_name + "_zoltan",experiment_name + "_bandwidth_" + graph_name + "_prawS"]
colours = ["red","blue","orange"] # as many as the number of experiments included
legend_labels = ['Zoltan','PRAW-arc-aware','PRAW-refinement']

# Each element on the following arrays corresponds to a column in columns_to_plot

	
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

for e in experiments:
	experiment = folder + e + "__" + str(num_processes)
	data = get_data_from_csv(experiment)
	print(e)
	for c in columns:
		print(np.mean(data[:,c]))
	print("\n")

