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

    // DEFAULT PARAMETERS
    char* graph_file = NULL;
    int iterations = 1;
    float imbalance_tolerance = 1.05f;
    char* bandwidth_file = NULL;
    bool use_bandwidth_in_partitioning = false;
    int rand_seed = time(NULL);

	// getting command line parameters
    extern char *optarg;
	extern int optind, opterr, optopt;
	int c;
	while( (c = getopt(argc,argv,"h:i:m:b:Ws:")) != -1 ) {
		switch(c) {
			case 'h': // hypergraph filename
				graph_file = optarg;
				break;
			case 'i': // max iterations
				iterations = atoi(optarg);
				break;
            case 's': // random seed
				rand_seed = atoi(optarg);
				break;
            case 'm': // max imbalance (in thousands)
				imbalance_tolerance = atoi(optarg) * 0.001f;
				break;
			case 'b':
				bandwidth_file = optarg;
				break;
			case 'W': // processes per node
				use_bandwidth_in_partitioning = true;
				break;
		}
	}

    if(graph_file == NULL || (use_bandwidth_in_partitioning && bandwidth_file == NULL)) {
        printf("Error: usage praw -h graph_filename -i iterations -m max_imbalance -b bandwidth_file [-W]\n");    
        return -1;
    }

    MPI_Init(&argc,&argv);
    int process_id;
    int num_processes;
    MPI_Comm_rank(MPI_COMM_WORLD,&process_id);
    MPI_Comm_size(MPI_COMM_WORLD,&num_processes);

    if(process_id == 0) {
        PRINTF("\nPARAMETER LIST:\n%s%s\n%s%i\n%s%f\n%s%s\n%s%i\n%s%i\n\n",
            "graph file: ",graph_file,
            "iterations: ",iterations,
            "imbalance tolerance: ",imbalance_tolerance,
            "bandwidth file: ",bandwidth_file,
            "use bandwidth file in partitionin: ",use_bandwidth_in_partitioning,
            "random seed: ",rand_seed);
    }

    std::string filename = graph_file;

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
    MPI_Bcast(&rand_seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    srand(rand_seed);
    
    for(int ii=0; ii < num_vertices; ii++) {
        partitioning[ii] = (double)rand() / (double)RAND_MAX * num_processes;
        if(partitioning[ii] == num_processes) partitioning[ii] -= 1;
    }
    // initialise comm cost matrix
    double** comm_cost_matrix = (double**)malloc(sizeof(double*) * num_processes);
    for(int ii=0; ii < num_processes; ii++) {
        comm_cost_matrix[ii] = (double*)calloc(num_processes,sizeof(double));
    }
    
    if(use_bandwidth_in_partitioning)
        PRAW::get_comm_cost_matrix_from_bandwidth(bandwidth_file,comm_cost_matrix,num_processes);
    else 
        PRAW::get_comm_cost_matrix_from_bandwidth(NULL,comm_cost_matrix,num_processes);
    
    // initialise vertex weight values
    int* vtx_wgt = (int*)calloc(num_vertices,sizeof(int));
    for(int ii =0; ii < num_vertices; ii++) {
        vtx_wgt[ii] = 1;
    }
    PRAW::SequentialStreamingPartitioning(partitioning,comm_cost_matrix, num_vertices,&hyperedges,&hedge_ptr,vtx_wgt,iterations, imbalance_tolerance);

    if(process_id == 0) {
        // if bandwidth file was not used in partitioning but was provided, use it in evaluation
        if(!use_bandwidth_in_partitioning && bandwidth_file != NULL) {
            PRAW::get_comm_cost_matrix_from_bandwidth(bandwidth_file,comm_cost_matrix,num_processes);
        }
        PRAW::storePartitionStats(filename,partitioning,num_processes,num_vertices,&hyperedges,&hedge_ptr,vtx_wgt,comm_cost_matrix);
    }
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