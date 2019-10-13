#ifndef BASELINE_SEQUENTIAL_VERTEX_PARTITION_PARTITIONING
#define BASELINE_SEQUENTIAL_VERTEX_PARTITION_PARTITIONING

#include <vector>
#include <algorithm>
#include <set>
#include <cstring>
#include <cstdio>
#include <stdlib.h>
#include <string>
#include "Partitioning.h"
#include "PRAW.h"

class BaselineSequentialVertexPartitioning : public Partitioning {
public:
	
	BaselineSequentialVertexPartitioning(char* experimentName, char* graph_file, float imbalance_tolerance) : Partitioning(graph_file,imbalance_tolerance) {
		experiment_name = experimentName;
	}
	virtual ~BaselineSequentialVertexPartitioning() {}
	
	virtual void perform_partitioning(int num_processes,int process_id, int* iterations) {

        // Uses Alistairh's basic hypergraph partitioning to assign vertices to processes
		if(num_processes <= 1) {
			PRINTF("Partitioning not required\n");
			return;
		}

        // initialise vertex weight values
        int* vtx_wgt = (int*)calloc(num_vertices,sizeof(int));
        for(int ii =0; ii < num_vertices; ii++) {
            vtx_wgt[ii] = 1;
        }
        
        if(process_id == 0) {
            *iterations = PRAW::BaselineSequentialVertexPartitioning(experiment_name,partitioning, num_processes, hgraph_name, vtx_wgt, imbalance_tolerance);
        } 
        MPI_Barrier(MPI_COMM_WORLD);
        // share new partitioning with other processes
        MPI_Bcast(partitioning,num_vertices,MPI_LONG,0,MPI_COMM_WORLD);
        

        // clean up operations
        
        free(vtx_wgt);
	}

private:
    char* experiment_name = NULL;
};


#endif
