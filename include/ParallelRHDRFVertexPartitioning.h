#ifndef PARALLEL_R_HDFR_VERTEX_PARTITION_PARTITIONING
#define PARALLEL_R_HDFR_VERTEX_PARTITION_PARTITIONING

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
Parallel version of reverse HDRF partitioning 
HDRF is for hedge partitioning, this version transforms the stream to make it a vertex partitioning. 
Stream is given to each process as a separate file 
Each line represents a vertex and the hyperedge_id it  belongs to.
 */

class ParallelRHDRFVertexPartitioning : public Partitioning {
public:
	
	ParallelRHDRFVertexPartitioning(char* experimentName, char* graph_file, float imbalance_tolerance, int iterations, char* comm_bandwidth_file, bool useBandwidth, bool proportionalCommCost, bool saveHistory, int syncBatchSize) : Partitioning(graph_file,imbalance_tolerance) {
		experiment_name = experimentName;
        comm_bandwidth_filename = comm_bandwidth_file;
        use_bandwidth_file = useBandwidth;
        max_iterations = iterations;
        proportional_comm_cost = proportionalCommCost;
        save_partitioning_history = saveHistory;
        sync_batch_size = syncBatchSize;
	}
	virtual ~ParallelRHDRFVertexPartitioning() {}
	
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

        // Transform hgraph format
        // From: each line contains the vertices belonging to a hyperedge
        // To: each line contains the hyperedges a vertex belongs to
        std::vector<std::vector<int> > hedge_ptr;
        PRAW::load_hedge_ptr_from_file_dist_CSR(hgraph_name, &hedge_ptr, process_id, num_processes, partitioning);

        std::string hgraph_file = hgraph_name;
        hgraph_file += "_";
        char str_int[16];
        sprintf(str_int,"%i",num_processes);
        hgraph_file += str_int;
        hgraph_file += "_";
        sprintf(str_int,"%i",process_id);
        hgraph_file += str_int;
        hgraph_file += ".hgr";

        PRINTF("%i: Storing model in file %s\n",process_id,hgraph_file.c_str());
        FILE *fp = fopen(hgraph_file.c_str(), "w+");
        
        // write header: NUM_HYPEREDGES NUM_VERTICES
        fprintf(fp,"%i %i",num_vertices,num_hyperedges); // this needs reversing because each line represents a vertex, not a hyperedge
        fprintf(fp,"\n");

        // write reminder of hyperedges per vertex
        for(int ii=0; ii < hedge_ptr.size(); ii++) {
            for(int he=0; he < hedge_ptr[ii].size(); he++) {
                fprintf(fp,"%i ",hedge_ptr[ii][he]);
            }
            fprintf(fp,"\n");
            
        }
        fclose(fp);
        hedge_ptr.clear();
        hedge_ptr.swap(hedge_ptr);
        ///////////////////////////

        *iterations = PRAW::ParallelHDRF(experiment_name,partitioning, comm_cost_matrix, hgraph_file, vtx_wgt, max_iterations, imbalance_tolerance, save_partitioning_history,false,sync_batch_size);
        
        // clean up operations
        for(int ii=0; ii < num_processes; ii++) {
            free(comm_cost_matrix[ii]);
        }
        free(comm_cost_matrix);
        free(vtx_wgt);

        // remove graph file
        if( remove(hgraph_file.c_str()) != 0 )
            printf( "Error deleting temporary hgraph file %s\n",hgraph_file.c_str() );
	}

private:
    char* experiment_name = NULL;
    char* comm_bandwidth_filename = NULL;
    bool use_bandwidth_file = false;
    int max_iterations;
    bool proportional_comm_cost = false;
    bool save_partitioning_history; 
    int sync_batch_size = 1;
};


#endif
