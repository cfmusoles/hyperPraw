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
	
	HyperPRAWPartitioning(char* graph_file, float imbalance_tolerance, int iterations, char* comm_bandwidth_file, bool parallel, bool useBandwidth) : Partitioning(graph_file,imbalance_tolerance) {
		comm_bandwidth_filename = comm_bandwidth_file;
        isParallel = parallel;
        use_bandwidth_file = useBandwidth;
        max_iterations = iterations;
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
            PRAW::ParallelIndependentRestreamingPartitioning(partitioning, comm_cost_matrix, hgraph_name, vtx_wgt, max_iterations, imbalance_tolerance);
        } else {
            if(process_id == 0) {
                // load complete model
                std::vector<std::vector<int> > hyperedges;
                std::vector<std::vector<int> > hedge_ptr;
                PRAW::load_hypergraph_from_file(hgraph_name, &hyperedges, &hedge_ptr);
                PRAW::SequentialStreamingPartitioning(partitioning,num_processes,comm_cost_matrix, num_vertices,&hyperedges,&hedge_ptr,vtx_wgt,max_iterations, imbalance_tolerance);
            }
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
};


#endif
