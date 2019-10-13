## Plot to show the relative order of a hypergraph stream
# The graph is processed in stream order and two lines are created:
#   cummulative number of pins processed
#   cummulative number of unique pins processed
# The difference between the two and their growth gives an impression of how ordered the stream is

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


hgraphs_folder = '../resources/'
#hgraph_files = ["atmosmodj.mtx.hgr","kkt_power.mtx.hgr","sat14_velev-vliw-uns-2.0-uq5.cnf.dual.hgr"]
hgraph_files = ["shuffled_ABACUS_shell_hd.mtx.hgr"]#"sat14_itox_vc1130.cnf.dual.hgr","2cubes_sphere.mtx.hgr","ABACUS_shell_hd.mtx.hgr","sparsine.mtx.hgr","pdb1HYS.mtx.hgr","sat14_10pipe_q0_k.cnf.primal.hgr","sat14_E02F22.cnf.hgr","webbase-1M.mtx.hgr"]
total_colour = 'blue'
unique_colour = 'orange'
pins_scale = 10e-4

storePlot = False
image_format = 'pdf'
image_names = ["streamOrder_" + str(i) for i,v in enumerate(hgraph_files)]

for i,hgraph in enumerate(hgraph_files):
    hgraph_file = hgraphs_folder + hgraph

    # load pin degrees from hgraph file
    with open(hgraph_file) as f:
        f.readline() #skip header info
        elements = [x.rstrip().split(" ") for x in f]
        elements = [list( map(int,i) ) for i in elements]

    cummulative_pins_visited = []
    cummulative_unique_pins_visited = []

    unique_pins_visited = set()
    total_pins = 0

    for pins in elements:
        total_pins += len(pins)
        cummulative_pins_visited.append(total_pins)
        unique_pins_visited.update(pins)
        cummulative_unique_pins_visited.append(len(unique_pins_visited))
    
    print("{} total pins, {} unique pins".format(total_pins,len(cummulative_unique_pins_visited)))

    #scale
    cummulative_pins_visited = [pin * pins_scale for pin in cummulative_pins_visited]
    cummulative_unique_pins_visited = [pin * pins_scale for pin in cummulative_unique_pins_visited]

    # create pandas frame 
    df = pd.DataFrame({'Total pins' : cummulative_pins_visited, 'Unique pins' : cummulative_unique_pins_visited})

    fig, ax = plt.subplots() # plot in the same subfigure

    g = sns.lineplot(data=df['Total pins'],ax=ax,label='Total pins',color=total_colour) # plot total pins
    g = sns.lineplot(data=df['Unique pins'],ax=ax,label='Unique pins',color=unique_colour) # plot total pins

    ax.legend()
    ax.set(xlabel='Streamed elements', ylabel='Number of pins (x' + str(int(1/pins_scale)) + ')')

    # fill area under the curve
    l1 = ax.lines[0]
    l2 = ax.lines[1]
    x1 = l1.get_xydata()[:,0]
    y1 = l1.get_xydata()[:,1]
    x2 = l2.get_xydata()[:,0]
    y2 = l2.get_xydata()[:,1]
    ax.fill_between(x1,y1, color=total_colour, alpha=0.5)
    ax.fill_between(x2,y2, color=unique_colour, alpha=1.0)

    #ax = fig.get_axes()[0].set_yscale('log')

    if storePlot:
        plt.savefig(image_names[i]+ image_format,format=image_format,dpi=1000)
    plt.show()
        

    