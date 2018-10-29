# Script to test hyperPraw partitioning

from mpi4py import MPI

comm = MPI.COMM_WORLD
num_processes = comm.Get_size()
process_id = comm.Get_rank()

## INPUT: hypergraph_input_file + max_iterations + imbalance_tolerance

# Read hypergraph and get hyperedges

# Broadcast random seed

# create (or load) communication cost matrix

# initialise vertex weights (can be 1)

# Obtain partitioning calling hyperPraw.ParallelStreamingPartitioning()