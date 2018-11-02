# Python implementation of HyperPRAW, a hypergraph parallel restreaming partitioning algorithm

from mpi4py import MPI
import math
import numpy as np
import sys

# Calculate current imbalance in partitioning
def calculateCurrentImbalance(partitioning,num_processes):
    # initiate workload balance per partition and local stream
    part_load = np.zeros(num_processes)
    for p in range(num_processes):
        part_load[p] = (partitioning == p).sum()

    return part_load.max() / part_load.mean()

# Function to get the cost of a partition

# Function for partitioning stats

# SequentialStreamingPartitioning

# ParallelStreamingPartitioning
def ParallelStreamingPartitioning(partitioning,num_vertices,comm_cost_matrix,hyperedges,hedge_ptr,vertex_weight,max_iterations,imbalance_tolerance):
    # Parameters from Bataglino 2015
    g = 1.5
    a = math.sqrt(2) * len(hyperedges) / pow(num_vertices,g)
    ta = 1.7

    comm = MPI.COMM_WORLD
    num_processes = comm.Get_size()
    process_id = comm.Get_rank()

    for iter in range(max_iterations):
        # initiate workload balance per partition and local stream
        part_load = np.zeros(num_processes)
        for p in range(num_processes):
            part_load[p] = (partitioning == p).sum()
        local_stream_partitioning = np.zeros(num_vertices,dtype=np.int)
        # process each stream, reassigning local vertices where required
        for vid in range(num_vertices):
            if partitioning[vid] != process_id:
                continue
            # reevaluate objective function per partition
            # |P^t_i union N(v)| = number of vertices in partition i that are neighbours of vertex v 
            # where are neighbours located
            # new communication cost incurred
            current_neighbours_in_partition = np.zeros(num_processes)
            comm_cost_per_partition = np.zeros(num_processes)
            for he_id in hedge_ptr[vid]:
                for neighbour in hyperedges[he_id]:
                    if neighbour == vid:
                        continue
                    part_dest = partitioning[neighbour]
                    current_neighbours_in_partition[part_dest] += 1
                    # commCost(v,Pi) = forall edge in edges(Pi) cost += w(e) * c(Pi,Pj) where i != j
                    for fp in range(num_processes):
                        comm_cost_per_partition[fp] += 1 * comm_cost_matrix[fp][part_dest]
            # find partition with better value cost function
            current_values = [current_neighbours_in_partition[pp] - comm_cost_per_partition[pp]  - a * g/2 * math.pow(part_load[pp],g-1) \
                                for pp in range(num_processes)]
            best_partition = current_values.index(max(current_values))
            # assign new partition to local stream
            local_stream_partitioning[vid] = best_partition
            # update intermediate workload and assignment values
            part_load[best_partition] += vertex_weight[vid]
            part_load[partitioning[vid]] -= vertex_weight[vid]
            partitioning[vid] = best_partition

        # share new partitioning with other streams
        comm.Allreduce(local_stream_partitioning,partitioning,op=MPI.MAX)

        # check if desired imbalance has been reached
        imbalance = calculateCurrentImbalance(partitioning,num_processes)
        if imbalance <= imbalance_tolerance:
            break
        
        #print stats
        if process_id == 0:
            print('Iteration {}: {}'.format(iter,imbalance))

        # update parameters
        a *= ta

    return partitioning

