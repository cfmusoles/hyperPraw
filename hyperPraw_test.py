# Script to test hyperPraw partitioning

from mpi4py import MPI
import sys
import numpy as np
import random
import hyperPraw

comm = MPI.COMM_WORLD
num_processes = comm.Get_size()
process_id = comm.Get_rank()

if len(sys.argv) < 4:
    print('Usage: python hyperPraw_test.py [hypergraph_filename] [max_iterations] [imbalance tolerance]')
    exit()

## INPUT: hypergraph_input_file + max_iterations + imbalance_tolerance
try:
    hgraph_file = sys.argv[1]
    max_iterations = int(sys.argv[2])
    imbalance_tolerance = float(sys.argv[3])
except:
    print('Error in parsing parameters')
    print('Usage: python hyperPraw_test.py [hypergraph_filename] [max_iterations] [imbalance tolerance]')
    exit()

# Read hypergraph and get hyperedges
hyperedges = []
hedge_ptr = []
num_vertices = 0
num_hyperedges = 0
with open(hgraph_file,'r') as f:
    #read header
    header = f.readline().split(' ')
    num_hyperedges = int(header[0])
    num_vertices = int(header[1])
    hedge_ptr = [[] for x in range(num_vertices)]
    #process each hypergraph
    counter = 0
    for line in f.readlines():
        ids = [int(x)-1 for x in line.split(' ')]
        hyperedges.append(ids)
        for id in ids:
            hedge_ptr[id].append(counter)
    counter += 1

# Broadcast random seed
rand_seed = comm.bcast(int(random.random() * 2**32), root=0)
np.random.seed(rand_seed)

# create (or load) communication cost matrix
intra_node_speed = 1
comm_cost_matrix = np.ones(shape=(num_processes,num_processes))
np.fill_diagonal(comm_cost_matrix,intra_node_speed)

# initialise vertex weights (can be 1)
vertex_weight = np.ones(num_vertices)

# Obtain partitioning calling hyperPraw.ParallelStreamingPartitioning()
partitioning = np.random.randint(0,num_processes,size=num_vertices)
partitioning = hyperPraw.ParallelStreamingPartitioning(partitioning,num_vertices,comm_cost_matrix,hyperedges,hedge_ptr,vertex_weight,max_iterations,imbalance_tolerance)

