#ifndef PRAW_PARTITION_PARTITIONING
#define PRAW_PARTITION_PARTITIONING

#include <vector>
#include <algorithm>
#include <set>
#include <cstring>
#include <cstdio>
#include <stdlib.h>
#include <string>
#include "Partitioning.h"
#include "PRAW.h"

class HyperPRAWPartitioning : public Partitioning {
public:
	
	HyperPRAWPartitioning(char* graph_file, float imbalance_tolerance, int iterations, char* comm_bandwidth_file, bool parallel, bool useBandwidth, bool resetPartitioning) : Partitioning(graph_file,imbalance_tolerance) {
		comm_bandwidth_filename = comm_bandwidth_file;
        isParallel = parallel;
        use_bandwidth_file = useBandwidth;
        max_iterations = iterations;
        reset_partitioning = resetPartitioning;
	}
	virtual ~HyperPRAWPartitioning() {}
	
	virtual void perform_partitioning(int num_processes,int process_id) {
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
            PRAW::get_comm_cost_matrix_from_bandwidth(comm_bandwidth_filename,comm_cost_matrix,num_processes);
        else 
            PRAW::get_comm_cost_matrix_from_bandwidth(NULL,comm_cost_matrix,num_processes);
        
        // initialise vertex weight values
        int* vtx_wgt = (int*)calloc(num_vertices,sizeof(int));
        for(int ii =0; ii < num_vertices; ii++) {
            vtx_wgt[ii] = 1;
        }

        if(isParallel) {
            int max_outer_iters = 1;
            for(int outer_iter=0; outer_iter < max_outer_iters; outer_iter++) {
                //PRAW::ParallelIndependentRestreamingPartitioning(partitioning, comm_cost_matrix, hgraph_name, vtx_wgt, max_iterations, imbalance_tolerance, outer_iter == 0);
            }
            PRAW::ParallelIndependentRestreamingPartitioning(partitioning, comm_cost_matrix, hgraph_name, vtx_wgt, max_iterations, imbalance_tolerance, reset_partitioning);
        } else {
            if(process_id == 0) {
                PRAW::SequentialStreamingPartitioning(partitioning, num_processes, comm_cost_matrix, hgraph_name, vtx_wgt, max_iterations, imbalance_tolerance,reset_partitioning);
            } 
            MPI_Barrier(MPI_COMM_WORLD);
            // share new partitioning with other processes
            MPI_Bcast(partitioning,num_vertices,MPI_LONG,0,MPI_COMM_WORLD);
        }

        // clean up operations
        for(int ii=0; ii < num_processes; ii++) {
            free(comm_cost_matrix[ii]);
        }
        free(comm_cost_matrix);

        // clean up operations
        free(vtx_wgt);
	}

private:
    char* comm_bandwidth_filename = NULL;
    bool isParallel = false;
    bool use_bandwidth_file = false;
    int max_iterations;
    bool reset_partitioning = false;
};


#endif
