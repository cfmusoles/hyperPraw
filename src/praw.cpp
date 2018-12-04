// Test harness for SPAAW (Streaming parallel Partitioning Architecture AWare)
#define VERBOSE                 // extra debug info printed out during runtime
#define EVALUATE_PARTITIONING   // evaluate resulting partitioning with predicted null compute traffic

#include <mpi.h>
#include <cstdio>
#include <stdlib.h>
#include <vector>
#include "PRAW.h"
#include <iterator>
#include <numeric>

int main(int argc, char** argv) {

    if(argc <= 3) {
        printf("Error: usage praw graph_filename iterations max_imbalance [bandwidth_file]\n");    
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
    char* bandwidth_file = NULL;
    if(argc >= 5) {
        bandwidth_file = argv[4];
    }

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
    
    PRAW::get_comm_cost_matrix_from_bandwidth(bandwidth_file,comm_cost_matrix,num_processes);

    // initialise vertex weight values
    int* vtx_wgt = (int*)calloc(num_vertices,sizeof(int));
    for(int ii =0; ii < num_vertices; ii++) {
        vtx_wgt[ii] = 1;
    }
    PRAW::SequentialStreamingPartitioning(partitioning,comm_cost_matrix, num_vertices,&hyperedges,&hedge_ptr,vtx_wgt,iterations, imbalance_tolerance);

    if(process_id == 0) 
        PRAW::storePartitionStats(filename,partitioning,num_processes,num_vertices,&hyperedges,&hedge_ptr,vtx_wgt,comm_cost_matrix);

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