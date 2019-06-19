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
	
	HyperedgePartitioning(char* experimentName, char* graph_file, float imbalance_tolerance, float ta_ref, int iterations, char* comm_bandwidth_file, bool parallel, bool useBandwidth, bool resetPartitioning, int stoppingCondition, bool proportionalCommCost, bool saveHistory) : Partitioning(graph_file,imbalance_tolerance) {
		experiment_name = experimentName;
        comm_bandwidth_filename = comm_bandwidth_file;
        isParallel = parallel;
        use_bandwidth_file = useBandwidth;
        max_iterations = iterations;
        reset_partitioning = resetPartitioning;
        stopping_condition = stoppingCondition;
        proportional_comm_cost = proportionalCommCost;
        ta_refinement = ta_ref;
        save_partitioning_history = saveHistory;
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
            PRAW::get_comm_cost_matrix_from_bandwidth(comm_bandwidth_filename,comm_cost_matrix,num_processes,proportional_comm_cost);
        else 
            PRAW::get_comm_cost_matrix_from_bandwidth(NULL,comm_cost_matrix,num_processes,proportional_comm_cost);
        
        // initialise vertex weight values
        int* vtx_wgt = (int*)calloc(num_vertices,sizeof(int));
        for(int ii =0; ii < num_vertices; ii++) {
            vtx_wgt[ii] = 1;
        }

        if(isParallel) {

            //PRAW::ParallelHyperedgePartitioning(experiment_name,partitioning, comm_cost_matrix, hgraph_name, vtx_wgt, max_iterations, imbalance_tolerance, ta_refinement, reset_partitioning,stopping_condition,save_partitioning_history);
            // alternative based on Alistairh minmax streaming
            std::vector<std::vector<int> > hedge_ptr;
            PRAW::load_hedge_ptr_from_file_dist_CSR(hgraph_name, &hedge_ptr, process_id, num_processes, partitioning);
            PRAW::ParallelHyperedgePartitioning_he_stream(experiment_name,partitioning,comm_cost_matrix, hgraph_name, num_vertices,&hedge_ptr,vtx_wgt,max_iterations, imbalance_tolerance,save_partitioning_history);
        } else {
            if(process_id == 0) {
                PRAW::SequentialHyperedgePartitioning(experiment_name,partitioning, num_processes, comm_cost_matrix, hgraph_name, vtx_wgt, max_iterations, imbalance_tolerance,ta_refinement,reset_partitioning,stopping_condition,save_partitioning_history);
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
    char* experiment_name = NULL;
    char* comm_bandwidth_filename = NULL;
    bool isParallel = false;
    bool use_bandwidth_file = false;
    int max_iterations;
    int stopping_condition;
    bool reset_partitioning = false;
    bool proportional_comm_cost = false;
    float ta_refinement;
    bool save_partitioning_history; 
};


#endif
