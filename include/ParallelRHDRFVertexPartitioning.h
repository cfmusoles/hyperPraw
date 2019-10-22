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
	
	ParallelRHDRFVertexPartitioning(char* experimentName, char* graph_file, char* streamFile, int max_processes, float imbalance_tolerance, int iterations, char* comm_bandwidth_file, bool useBandwidth, bool proportionalCommCost, int syncBatchSize, bool use_expected_workload, bool input_order) : Partitioning(graph_file,imbalance_tolerance) {
		experiment_name = experimentName;
        comm_bandwidth_filename = comm_bandwidth_file;
        use_bandwidth_file = useBandwidth;
        max_iterations = iterations;
        proportional_comm_cost = proportionalCommCost;
        sync_batch_size = syncBatchSize;
        stream_file = streamFile;
        input_order_round_robin = input_order;
        max_num_processes = max_processes;        

        split_and_configure_stream();

	}
	virtual ~ParallelRHDRFVertexPartitioning() {
        MPI_Comm_free(&partitioning_comm);
    }

    void split_and_configure_stream() {

        int process_id;
        MPI_Comm_rank(MPI_COMM_WORLD,&process_id);
        
        // Create MPI communicator group (of only max_number_processes)
        // Only those processes will work towards partitioning
        MPI_Comm_split(MPI_COMM_WORLD,
                        process_id < max_num_processes ? 1 : 0, // this is the colour, is used to match the processes to the new group
                        process_id,                                         // this is the key, used to generate new rank
                        &partitioning_comm);                               // the subgroup generated
        
        int num_processes;
        MPI_Comm_size(partitioning_comm,&num_processes); 

        is_active = process_id < max_num_processes;

        if(!is_active) return;

        // divide original hMetis file into as many streams as processes
        std::string hgraph_file = stream_file;
        hgraph_file += "_";
        char str_int[16];
        sprintf(str_int,"%i",num_processes);
        hgraph_file += str_int;
        hgraph_file += "_";
        sprintf(str_int,"%i",process_id);
        hgraph_file += str_int;
        hgraph_file += ".hgr";
        hgraph_part_file = hgraph_file;

        PRINTF("%i: Storing model in file %s\n",process_id,hgraph_file.c_str());
        FILE *fp = fopen(hgraph_file.c_str(), "w+");
        
        // load and parse full graph
        std::ifstream istream(stream_file);
        
        if(!istream) {
            printf("Error while opening hMETIS file %s\n",stream_file);
            return;
        }

        std::string line;
        // ignore header
        std::getline(istream,line);
        // write header: NUM_HYPEREDGES NUM_VERTICES (inverted because this is rHDRF)
        fprintf(fp,"%i %i",num_vertices,num_hyperedges);
        fprintf(fp,"\n");

        // read reminder of file (one line per vertex)
        int local_stream_size = num_vertices / num_processes + ((num_vertices % num_processes > process_id) ? 1 : 0);
        int first_element_in_stream = (num_vertices / num_processes * process_id + std::min(num_vertices % num_processes,process_id));
        int counter = 0;
        while(std::getline(istream,line)) {
            // Choice between round robin or bulk first //
            bool own_element;
            if(input_order_round_robin) {
                own_element = counter % num_processes == process_id;
            } else {
                own_element = counter < (first_element_in_stream + local_stream_size) && counter >= first_element_in_stream;
            }
            if(own_element) {
                char str[line.length() + 1]; 
                strcpy(str, line.c_str()); 
                char* token = strtok(str, " "); 
                while (token != NULL) { 
                    fprintf(fp,"%i ",atoi(token));
                    token = strtok(NULL, " "); 
                }
                fprintf(fp,"\n");
            }
            counter++;
            
        }
        istream.close();
        fclose(fp);
    }
	
	virtual void perform_partitioning(int num_partitions,int process_id, int* iterations) {
		if(num_partitions <= 1) {
			PRINTF("Partitioning not required\n");
			return;
		}

        if(is_active) {
            // initialise comm cost matrix
            double** comm_cost_matrix = (double**)malloc(sizeof(double*) * num_partitions);
            for(int ii=0; ii < num_partitions; ii++) {
                comm_cost_matrix[ii] = (double*)calloc(num_partitions,sizeof(double));
            }
            
            if(use_bandwidth_file)
                PRAW::get_comm_cost_matrix_from_bandwidth(comm_bandwidth_filename,comm_cost_matrix,num_partitions,proportional_comm_cost);
            else 
                PRAW::get_comm_cost_matrix_from_bandwidth(NULL,comm_cost_matrix,num_partitions,proportional_comm_cost);
            
            // initialise vertex weight values
            int* vtx_wgt = (int*)calloc(num_vertices,sizeof(int));
            for(int ii =0; ii < num_vertices; ii++) {
                vtx_wgt[ii] = 1;
            }

            *iterations = PRAW::ParallelHDRF(experiment_name,partitioning, num_partitions, partitioning_comm, comm_cost_matrix, hgraph_part_file.c_str(), vtx_wgt, max_iterations, imbalance_tolerance,false,sync_batch_size,input_order_round_robin);

            // clean up operations
            for(int ii=0; ii < num_partitions; ii++) {
                free(comm_cost_matrix[ii]);
            }
            free(comm_cost_matrix);
            free(vtx_wgt);

            // remove graph file
            if(remove(hgraph_part_file.c_str()) != 0 )
                printf( "Error deleting temporary hgraph file %s\n",hgraph_part_file.c_str() );
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // broadcast results to not participating processes
        MPI_Bcast(partitioning, num_vertices, MPI_LONG, 0,MPI_COMM_WORLD);

        
	}

private:
    char* experiment_name = NULL;
    char* comm_bandwidth_filename = NULL;
    bool use_bandwidth_file = false;
    int max_iterations;
    bool proportional_comm_cost = false;
    int sync_batch_size = 1;
    char* stream_file = NULL;
    std::string hgraph_part_file;
    bool input_order_round_robin = true;
    MPI_Comm partitioning_comm;
    int max_num_processes;
    bool is_active = true;
};


#endif
