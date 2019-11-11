## Python helper file to plot and store graphs to represent network bandwidth and MPI performance

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from pylab import *

num_processes = 96
plot_bandwidth = False			# plot network bandwidth data
plot_sent_data = True			# plot application sent data
plot_comm_cost = False			# plot combined comm cost
storeResults = False
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
#
# msdoor.mtx.hgr
# xenon2.mtx.hgr
# tmt_unsym.mtx.hgr
# BenElechi1.mtx.hgr
# xenon2.mtx.hgr
# IMDB.mtx.hgr

# atmosmodj.mtx.hgr
# kkt_power.mtx.hgr 

# small_dense_uniform.hgr
# small_dense_powerlaw.hgr
# large_sparse_uniform.hgr
# large_sparse_powerlaw.hgr

folder = "../results/norm_cb/"
bandwidth_send_experiment_name = 'results_mpi_send_bandwidth_1_' + str(num_processes)
graph_name = "small_dense_powerlaw.hgr"
partitioning = 'hyperPrawVertex'
test_name = 'norm_cb_hyperPraw_bandwidth_lambda1_16'

sim_sent_experiment = test_name + '_' + graph_name + '_' + partitioning + '_hedgeSim_comm_cost__' + str(num_processes)
#sim_sent_experiment = test_name + '_' + partitioning + '_comm_matrix_' + str(num_processes)

xlabel = "Process"
ylabel = "Process"
zlabels = ["MB/sec","Bytes"]
titles = ["P2P Bandwidth (praw-noB)","Sent data (praw-withB)","Cost of communication (praw-withB)"]
filenames = ["a1","a2","a3"]
image_format = 'pdf'

# general plot settings
plt.rcParams['figure.facecolor'] = 'white'
fig_settings = {  
        'lines.linewidth': 0.5,
        'axes.linewidth': 0.5,
        'axes.labelsize': 'medium',
        'legend.fontsize': 'medium',
        'font.size': 14,
        'savefig.dpi': 200,
}
plt.rcParams.update(fig_settings)

def get_data_from_csv(filename,delimiter=","):
	data = np.genfromtxt(filename,skip_header=0,delimiter=delimiter)
	return data

def plot_3dgraph(data,title,zlabel,filename):
	fig = plt.figure()
	ax = Axes3D(fig)

	lx= len(data[0])            # Work out matrix dimensions
	ly= len(data[:,0])
	xpos = np.arange(0,lx,1)    # Set up a mesh of positions
	ypos = np.arange(0,ly,1)
	xpos, ypos = np.meshgrid(xpos+0.5, ypos+0.5)

	xpos = xpos.flatten()   # Convert positions to 1D array
	ypos = ypos.flatten()
	zpos = np.ones(lx*ly)*1e-10

	dx = 1. * np.ones_like(zpos)
	dy = dx.copy()
	dz = data.flatten()

	colors = plt.cm.jet(data.flatten()/float(data.max()))
	ax.bar3d(xpos,ypos,zpos, dx, dy, dz,  color=colors, alpha=1., zsort='max',linewidth=0)
	
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	ax.set_zlabel(zlabel)
	
	if show_title:
		plt.title(title)
	if storeResults:
		plt.savefig(filename + "." + image_format,format=image_format)
	plt.show()

def plot_2dgraph(data,title,label,fname):
	#transform 0 values
	data[data == 0] = np.min(data[np.nonzero(data)])
	# transform to log scale
	if log_scale:
		data = np.log(data)
	# colour maps http://scipy-cookbook.readthedocs.io/items/Matplotlib_Show_colormaps.html
	fig, axes = plt.subplots(1,1)
	# plot bandwidth
	c = axes.pcolor(data,cmap=get_cmap("jet"),rasterized=True)
	axes.set_xlabel(xlabel)
	axes.set_ylabel(ylabel)
	if show_title:
		axes.set_title(title)
	cbar = fig.colorbar(c,ax=axes)
	cbar.ax.set_ylabel(label, rotation=90)
	
	fig.tight_layout()
	#fig.set_size_inches(11, 5)

	if storeResults:
		plt.savefig(fname + "." + image_format,format=image_format)
	plt.show()


## creating figure for bandwidth estimation (send performance)
if plot_bandwidth:

	#data = get_data_from_csv(folder + bandwidth_send_experiment_name,'\t')
	#plot_3dgraph(data,titles[0],zlabels[0],filenames[0])

	## Plot bandwidth and data sent separately in a double figure
	bandwidth_data = get_data_from_csv(folder + bandwidth_send_experiment_name,'\t')
	plot_2dgraph(bandwidth_data,"P2P Bandwidth",'log (MB/s)',filenames[1])

	

## creating combined figure for cost of communication
## each value corresponds to a process-process pair: total data sent / bandwidth 
if plot_sent_data:
	## IN 3D##
	#data = get_data_from_csv(folder + sim_sent_experiment,' ')
	#plot_3dgraph(data,titles[1],zlabels[1],filenames[1])
	####

	## Plot bandwidth and data sent separately in a double figure
	comm_data = get_data_from_csv(folder + sim_sent_experiment,' ')
	plot_2dgraph(comm_data,"Actual data sent",'log (Bytes sent)',filenames[2])
	

if plot_comm_cost:
	## Plot bandwidth and data sent separately in a double figure
	bandwidth_data = get_data_from_csv(folder + bandwidth_send_experiment_name,'\t')
	comm_data = get_data_from_csv(folder + sim_sent_experiment,' ')
	#transform 0 values
	bandwidth_data[bandwidth_data == 0] = np.max(bandwidth_data)
	comm_data[comm_data == 0] = 0.1
	# calculate cost of comm (theoretical)
	cost_data = np.array([c/b for b, c in zip(bandwidth_data, comm_data)])
	cost_data = np.nan_to_num(cost_data)
	print('Total cost: ' + str(cost_data.sum()))
	

	## Compound graph##
	#bandwidth_data = get_data_from_csv(folder + bandwidth_send_experiment_name,'\t')
	#comm_data = get_data_from_csv(folder + sim_sent_experiment,' ')
	#cost_data = np.array([c/b for b, c in zip(bandwidth_data, comm_data)])
	#cost_data = np.nan_to_num(cost_data)
	#print('Total cost: ' + str(cost_data.sum()))
	# colour maps http://scipy-cookbook.readthedocs.io/items/Matplotlib_Show_colormaps.html
	#pcolor(cost_data,cmap=get_cmap("binary"))
	#cbar = colorbar()
	#plt.title(titles[2])
	#plt.xlabel(xlabel)
	#plt.ylabel(ylabel)
	#cbar.ax.set_ylabel("Cost (time)",rotation=90)
	#if storeResults:
	#	plt.savefig(filenames[2] + "." + image_format,format=image_format,dpi=300)
	#plt.show()

	####

