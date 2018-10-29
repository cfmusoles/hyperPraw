# Python implementation of HyperPRAW, a hypergraph parallel restreaming partitioning algorithm

from mpi4py import MPI

comm = MPI.COMM_WORLD
num_processes = comm.Get_size()
process_id = comm.Get_rank()