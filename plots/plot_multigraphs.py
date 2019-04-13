## Python helper file to plot and store graphs to represent scaling in distributed simulations
# Graphs produced: as many as columns that want to be compared

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

num_processes = 144

show_error = True
as_bar_plot = True
log_scale = False
show_title = False

# "sat14_itox_vc1130.cnf.dual.hgr" $SEED 2 #Y for esim
# "2cubes_sphere.mtx.hgr" $SEED 3 #Y for esim
# "ABACUS_shell_hd.mtx.hgr" $SEED 40 #Y
# "sparsine.mtx.hgr" $SEED 2 

# "pdb1HYS.mtx.hgr" $SEED 1 20 #Y # hedge sim is too short
# "sat14_atco_enc1_opt1_05_21.cnf.dual.hgr" $SEED 4 1 #N 
# "sat14_10pipe_q0_k.cnf.primal.hgr" $SEED 1 1 #Y
# "sat14_E02F22.cnf.hgr" $SEED 3 1 #Y
# "webbase-1M.mtx.hgr" $SEED 1 1 #Y
# "ship_001.mtx.hgr" $SEED 1 30 #Y # hedge sim is too short

folder = "../results/runtime/"
experiment_name = "runtime"
graphs = ["sat14_itox_vc1130.cnf.dual.hgr","2cubes_sphere.mtx.hgr","ABACUS_shell_hd.mtx.hgr","sparsine.mtx.hgr","pdb1HYS.mtx.hgr","sat14_atco_enc1_opt1_05_21.cnf.dual.hgr","sat14_10pipe_q0_k.cnf.primal.hgr","sat14_E02F22.cnf.hgr","webbase-1M.mtx.hgr","ship_001.mtx.hgr"]
graph_names = ["sat14 itox","2cubes","ABACUS","sparsine","pdb1HYS","sat14 atco dual","sat14 10pipe primal","sat14 E02F22","webbase-1M","ship 001"]
# each element on the following arrays corresponds to an experiment run (collection of files)
#experiments_name = [experiment_name +  "_zoltan_" + graph_name + "_zoltan",experiment_name + "_default_" + graph_name + "_prawS",experiment_name + "_bandwidth_" + graph_name + "_prawS",experiment_name + "_refinement_" + graph_name + "_prawSref"]
experiments_name = ["zoltan","default","bandwidth","refinement"]
experiments_partitioning = ["zoltan","prawS","prawS","prawSref"]
colours = ["black","tomato","seagreen","yellow"] # as many as the number of experiments included
patterns = ["///" , "|||" , "---" , "xxx" ]
legend_labels = ['Zoltan','PRAW','PRAW-arc-aware','PRAW-refinement']

# Each element on the following arrays corresponds to a column in columns_to_plot
columns_to_plot = [1,4,3,5,8]#,9,10,11]
reference_values = [0,2,1,6,7,8,3,1,1] # used to take values on each column divided by these
use_ref_values = False
scale_plots = [1,1,1,1e-3,1,1,1,1]
plot_title = ["EdgeSim time","Edge cut","Hyperedge cut","SOED","Edge comm cost","Hedge comm cost","Messages sent (edge)","Messages sent (hedge)"]
plot_ylabel = ["Time(s)","Cut ratio","Cut ratio","SOED (thousands)","Cost","Cost","Messages sent","Messages sent"]
image_format = 'pdf'
plot_name = ["a_" + str(x) for x in range(len(columns_to_plot))] #["a1","a2","a3","a4","a5","a6","a7"]

bar_plot_size = 0.7 / len(experiments_name)

# general plot settings
plt.rcParams['figure.facecolor'] = 'white'
fig_settings = {  
        'lines.linewidth': 0.5,
        'axes.linewidth': 0.5,
        'axes.labelsize': 'small',
        'legend.fontsize': 'small',
        'font.size': 14,
        'savefig.dpi': 200,
		'figure.figsize': {14.6,6}
}
plt.rcParams.update(fig_settings)

# required workaround to be able to export hatch patterns
def savepdfviasvg( name, **kwargs):
    import subprocess
    plt.savefig(name+".svg", format="svg", **kwargs)
    incmd = ["inkscape", name+".svg", "--export-pdf={}.pdf".format(name),
             "--export-pdf-version=1.5"] #"--export-ignore-filters",
    subprocess.check_output(incmd)


def get_data_from_csv(filename):
	data = np.genfromtxt(filename,skip_header=1,delimiter=",")
	return data

def plot(x,y, error,title,ylabel,name,colour,pattern,legend,show,global_counter):
	if as_bar_plot:
		rx = x + global_counter * bar_plot_size*np.ones(len(x))
		plt.bar(rx,y,width=bar_plot_size,color=colour,label=legend,hatch=pattern)
		
	else:
		if show_error:
			plt.errorbar(x, y, error,linewidth=1,color=colour,label=legend,marker='s',markersize=5)
		else:
			plt.errorbar(x, y,linewidth=1,color=colour,label=legend,marker='s',markersize=5)
	
	if log_scale:
		plt.yscale("log")
		#plt.ylim(1,3000)
	plt.ylabel(ylabel)
	if show_title:
		plt.title(title)
	plt.tick_params(axis='x',which='minor',bottom=False,labelbottom=False)
	#tick_names = [x.rsplit('.')[0] for x in graphs]
	plt.xticks(x,graph_names,rotation=25)
	#plt.tight_layout()
	plt.gcf().subplots_adjust(left=0.17,bottom=0.25)
	if len(graphs) > 1:
		plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
          ncol=len(experiments_name), fancybox=True)
	if show:
		savepdfviasvg(name)
		#plt.savefig(name + "." + image_format,format=image_format,dpi=1000)
		plt.show()


#####################
# one plot per column
#####################
for i in range(len(columns_to_plot)):
	# creating figure for the column plot
	plt.figure()
	for j in range(len(experiments_name)):
		h = len(columns_to_plot)
		means = [[] for y in range(h)]
		stdevs = [[] for y in range(h)]
		
		for exp in graphs:
			filename = experiment_name + "_" + experiments_name[j] + "_" + exp + "_" + experiments_partitioning[j] + "__" + str(num_processes)
			data = get_data_from_csv(folder + filename)
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
		plot(range(len(graphs)),y,error,plot_title[i],plot_ylabel[i],plot_name[i],colours[j],patterns[j],legend_labels[j], j == (len(experiments_name)-1),j)