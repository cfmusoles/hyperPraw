#ifndef HYPEREDGE_PARTITION_PARTITIONING
#define HYPEREDGE_PARTITION_PARTITIONING

#include <vector>
#include <algorithm>
#include <set>
#include <cstring>
#include <cstdio>
#include <stdlib.h>
#include <string>
#include "Partitioning.h"
#include "PRAW.h"

class HyperedgePartitioning : public Partitioning {
public:
	
	HyperedgePartitioning(char* experimentName, char* graph_file, int iterations, float imbalance_tolerance, char* comm_bandwidth_file, bool parallel, bool useBandwidth, bool saveHistory) : Partitioning(graph_file,imbalance_tolerance,false) {
		experiment_name = experimentName;
        comm_bandwidth_filename = comm_bandwidth_file;
        isParallel = parallel;
        use_bandwidth_file = useBandwidth;
        save_partitioning_history = saveHistory;
        max_iterations = iterations;
	}
	virtual ~HyperedgePartitioning() {}
	
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
            PRAW::get_comm_cost_matrix_from_bandwidth(comm_bandwidth_filename,comm_cost_matrix,num_processes,false);
        else 
            PRAW::get_comm_cost_matrix_from_bandwidth(NULL,comm_cost_matrix,num_processes,false);
        
        // initialise vertex weight values
        int* he_wgt = (int*)calloc(num_hyperedges,sizeof(int));
        for(int ii =0; ii < num_hyperedges; ii++) {
            he_wgt[ii] = 1;
        }

        PRAW::ParallelHyperedgePartitioning(experiment_name,partitioning, comm_cost_matrix, hgraph_name, he_wgt, max_iterations, imbalance_tolerance, save_partitioning_history);
        
        // clean up operations
        for(int ii=0; ii < num_processes; ii++) {
            free(comm_cost_matrix[ii]);
        }
        free(comm_cost_matrix);

        // clean up operations
        free(he_wgt);
	}

private:
    char* experiment_name = NULL;
    char* comm_bandwidth_filename = NULL;
    bool isParallel = false;
    bool use_bandwidth_file = false;
    bool save_partitioning_history; 
    int max_iterations;
};


#endif
