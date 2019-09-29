#ifndef SIMPLE_PARALLEL_VERTEX_PARTITION_PARTITIONING
#define SIMPLE_PARALLEL_VERTEX_PARTITION_PARTITIONING

#include <vector>
#include <algorithm>
#include <set>
#include <cstring>
#include <cstdio>
#include <stdlib.h>
#include <string>
#include "Partitioning.h"
#include "PRAW.h"

/* DESCRIPTION OF THE STREAMING ALGORITHM
Parallel version of vertex partitioning. 
Based on Alistairh's streaming partitioning.
Stream is given to each process as a list of the hyperedge_id each vertex belongs to.
Requires processing whole file for each vertex, but it only stores the data relevant for local stream (scalable)
 */

class SimpleParallelVertexPartitioning : public Partitioning {
public:
	
	SimpleParallelVertexPartitioning(char* experimentName, char* graph_file, float imbalance_tolerance, int iterations, char* comm_bandwidth_file, bool useBandwidth, bool proportionalCommCost, bool saveHistory ) : Partitioning(graph_file,imbalance_tolerance) {
		experiment_name = experimentName;
        comm_bandwidth_filename = comm_bandwidth_file;
        use_bandwidth_file = useBandwidth;
        max_iterations = iterations;
        proportional_comm_cost = proportionalCommCost;
        save_partitioning_history = saveHistory;
	}
	virtual ~SimpleParallelVertexPartitioning() {}
	
	virtual void perform_partitioning(int num_processes,int process_id, int* iterations) {
		if(num_processes <= 1) {
			PRINTF("Partitioning not required\n");
			return;
		}

        // initialise comm cost matrix
        double** comm_cost_matrix = (double**)malloc(sizeof(double*) * num_processes);
        for(int ii=0; ii < num_processes; ii++) {
            comm_cost_matrix[ii] = (double*)calloc(num_processes,sizeof(double));
        }
        
        if(use_bandwidth_file)
            PRAW::get_comm_cost_matrix_from_bandwidth(comm_bandwidth_filename,comm_cost_matrix,num_processes,proportional_comm_cost);
        else 
            PRAW::get_comm_cost_matrix_from_bandwidth(NULL,comm_cost_matrix,num_processes,proportional_comm_cost);
        
        // initialise vertex weight values
        int* vtx_wgt = (int*)calloc(num_vertices,sizeof(int));
        for(int ii =0; ii < num_vertices; ii++) {
            vtx_wgt[ii] = 1;
        }


        // USED FOR ParallelHyperedgePartitioning_he_stream
        std::vector<std::vector<int> > hedge_ptr;
        PRAW::load_hedge_ptr_from_file_dist_CSR(hgraph_name, &hedge_ptr, process_id, num_processes, partitioning);

        *iterations = PRAW::ParallelHyperedgePartitioning_he_stream(experiment_name,partitioning,comm_cost_matrix, hgraph_name, num_vertices,&hedge_ptr,vtx_wgt,max_iterations, imbalance_tolerance,save_partitioning_history);

        // clean up operations
        for(int ii=0; ii < num_processes; ii++) {
            free(comm_cost_matrix[ii]);
        }
        free(comm_cost_matrix);

        // clean up operations
        free(vtx_wgt);
	}

private:
    char* experiment_name = NULL;
    char* comm_bandwidth_filename = NULL;
    bool use_bandwidth_file = false;
    int max_iterations;
    bool proportional_comm_cost = false;
    bool save_partitioning_history; 
};


#endif
