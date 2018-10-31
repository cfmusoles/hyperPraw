# Python implementation of HyperPRAW, a hypergraph parallel restreaming partitioning algorithm

from mpi4py import MPI
import math
import numpy as np


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
        for partition in partitioning:
            part_load[partition] += 1
        local_stream_partitioning = np.zeros(num_vertices)
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
            

    return []

'''
int ParallelStreamingPartitioning(idx_t* partitioning, float** comm_cost_matrix, int num_vertices, std::vector<std::vector<int> >* hyperedges, std::vector<std::vector<int> >* hedge_ptr, int* vtx_wgt, int max_iterations, float imbalance_tolerance) {
        //PARAMETERS: From Battaglino 2015 //
        float g = 1.5;
        float a = sqrt(2) * hyperedges->size() / pow(num_vertices,g);
        float ta = 1.7;
        //float tta = 0.98;
        ///////////////
        
        int process_id;
        MPI_Comm_rank(MPI_COMM_WORLD,&process_id);
        int num_processes;
        MPI_Comm_size(MPI_COMM_WORLD,&num_processes);

        //printPartitionStats(partitioning,num_processes,num_vertices,hyperedges,hedge_ptr,vtx_wgt,comm_cost_matrix);
        
        int* part_load = (int*)calloc(num_processes,sizeof(int));
        idx_t* local_stream_partitioning = (idx_t*)malloc(num_vertices*sizeof(idx_t));
        float* comm_cost_per_partition = (float*)malloc(num_processes*sizeof(float));
        int* current_neighbours_in_partition = (int*)malloc(num_processes*sizeof(int));
                
        for(int iter=0; iter < max_iterations; iter++) {
            memset(local_stream_partitioning,0,num_vertices * sizeof(idx_t));
            memset(part_load,0,num_processes * sizeof(int));
            for(int ii=0; ii < num_vertices; ii++) {
                part_load[partitioning[ii]] += vtx_wgt[ii]; // workload for vertex
            }
            // go through own vertex list and reassign
            for(int ii=0; ii < num_vertices; ii++) {
                if(partitioning[ii] != process_id) continue;
                // reevaluate objective function per partition
                // |P^t_i union N(v)| = number of vertices in partition i that are neighbours of vertex v 
                // where are neighbours located
                // new communication cost incurred
                memset(current_neighbours_in_partition,0,num_processes * sizeof(int));
                memset(comm_cost_per_partition,0,num_processes * sizeof(float));
                for(int he = 0; he < hedge_ptr->at(ii).size(); he++) {
                    int he_id = hedge_ptr->at(ii)[he];
                    for(int vt = 0; vt < hyperedges->at(he_id).size(); vt++) {
                        int dest_vertex = hyperedges->at(he_id)[vt];
                        if(dest_vertex == ii) continue;
                        int dest_part = partitioning[dest_vertex];
                        current_neighbours_in_partition[dest_part] += 1;
                        // recalculate comm cost for all possible partition assignments of ii
                        //  commCost(v,Pi) = forall edge in edges(Pi) cost += w(e) * c(Pi,Pj) where i != j
                        for(int fp=0; fp < num_processes; fp++) {
                            comm_cost_per_partition[fp] += 1 * comm_cost_matrix[fp][dest_part];
                        }
                    }
                }
                float max_value = std::numeric_limits<float>::lowest();
                int best_partition = partitioning[ii];
                for(int pp=0; pp < num_processes; pp++) {
                    // objective function is a mix of Battaglino 2015 (second part) and Zheng 2016 (communication cost part)
                    // (|P^t_i union N(v)| - commCost(v,Pi) - a * g/2 * |B|^(g-1))
                    float current_value = (float)current_neighbours_in_partition[pp] - comm_cost_per_partition[pp]  - a * g/2 * pow(part_load[pp],g-1);
                    if(current_value > max_value) {
                        max_value = current_value;
                        best_partition = pp;
                    }
                }
                    
                local_stream_partitioning[ii] = best_partition;
                // update intermediate workload and assignment values
                part_load[best_partition] += vtx_wgt[ii];
                part_load[partitioning[ii]] -= vtx_wgt[ii];
                partitioning[ii] = best_partition;
            }
            
            // share new partitioning with other streams
            MPI_Allreduce(local_stream_partitioning,partitioning,num_vertices,MPI_LONG,MPI_MAX,MPI_COMM_WORLD);
            
            // check if desired imbalance has been reached
            float imbalance = calculateImbalance(partitioning,num_processes,num_vertices,vtx_wgt);
            printf("%i: %f (%f | %f)\n",iter,imbalance,a,ta);
            if(imbalance < imbalance_tolerance) break;

            //printf("Iteration %i\n",iter);
            //printPartitionStats(partitioning,num_processes,num_vertices,hyperedges,hedge_ptr,vtx_wgt,comm_cost_matrix);

            // update parameters
            a *= ta;
            //if(ta > 1.05f) ta *= tta;
        }
        printPartitionStats(partitioning,num_processes,num_vertices,hyperedges,hedge_ptr,vtx_wgt,comm_cost_matrix);

        // clean up
        free(part_load);
        free(local_stream_partitioning);
        free(comm_cost_per_partition);
        free(current_neighbours_in_partition);

        // return successfully
        return 0;
    }
'''