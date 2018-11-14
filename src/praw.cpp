// Test harness for SPAAW (Streaming parallel Partitioning Architecture AWare)
#define VERBOSE

#include <mpi.h>
#include <cstdio>
#include <stdlib.h>
#include <vector>
#include "PRAW.h"
#include <iterator>
#include <numeric>

int main(int argc, char** argv) {

    if(argc <= 3) {
        printf("Error: usage praw [graph filename] [iterations] [max_imbalance]\n");    
        return -1;
    }
    
    MPI_Init(&argc,&argv);
    int process_id;
    int num_processes;
    MPI_Comm_rank(MPI_COMM_WORLD,&process_id);
    MPI_Comm_size(MPI_COMM_WORLD,&num_processes);

    std::string filename = argv[1];
    int iterations = atoi(argv[2]);
    float imbalance_tolerance = atof(argv[3]);

    std::vector<std::vector<int> > hyperedges;
    std::vector<std::vector<int> > hedge_ptr;

    if(PRAW::load_hypergraph_from_file(filename,&hyperedges,&hedge_ptr) != 0) {
        printf("Error, could not find hypergraph file\n");
        MPI_Finalize();
        return 0;
    }
    int num_vertices = hedge_ptr.size();
    idx_t* partitioning = (idx_t*)calloc(num_vertices,sizeof(idx_t));
    // Assign unique vertices to partitions (common seed)
    int rand_seed = time(NULL);
    MPI_Bcast(&rand_seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    srand(rand_seed);
    
    for(int ii=0; ii < num_vertices; ii++) {
        partitioning[ii] = (double)rand() / (double)RAND_MAX * num_processes;
        if(partitioning[ii] == num_processes) partitioning[ii] -= 1;
    }
    // initialise comm cost matrix
    float** comm_cost_matrix = (float**)malloc(sizeof(float*) * num_processes);
    for(int ii=0; ii < num_processes; ii++) {
        comm_cost_matrix[ii] = (float*)calloc(num_processes,sizeof(float));
    }
    for(int ii = 0; ii < num_processes;ii++) {
        for(int jj=0; jj < num_processes; jj++) {
            comm_cost_matrix[ii][jj] = ii==jj ? 1 : 1;
        }
    }
    // initialise vertex weight values
    int* vtx_wgt = (int*)calloc(num_vertices,sizeof(int));
    for(int ii =0; ii < num_vertices; ii++) {
        vtx_wgt[ii] = 1;
    }
    PRAW::ParallelStreamingPartitioning(partitioning,comm_cost_matrix, num_vertices,&hyperedges,&hedge_ptr,vtx_wgt,iterations, imbalance_tolerance);

    // clear memory
    free(partitioning);
    free(vtx_wgt);
    for(int ii=0; ii < num_processes; ii++) {
        free(comm_cost_matrix[ii]);
    }
    free(comm_cost_matrix);

    MPI_Finalize();
    return 0;
}