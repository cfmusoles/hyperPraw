#ifndef PARALLEL_VERTEX_PARTITION_PARTITIONING
#define PARALLEL_VERTEX_PARTITION_PARTITIONING

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
Parallel version of Alistairh streaming partitioning 
Stream is given to each process as a separate file 
Two versions where the elements can represent a hyperedge or a vertex
Determined by parameter element_is_vertex
 */

class ParallelStreamingPartitioning : public Partitioning {
public:
	
	ParallelStreamingPartitioning(char* experimentName, char* graph_file, char* streamFile, int max_processes, float imbalance_tolerance, int windowSize, bool input_order,bool elementIsVertex, bool useHDRF, bool staggeredStreams, bool useBalanceCost, float l) : Partitioning(graph_file,imbalance_tolerance,elementIsVertex) {
		experiment_name = experimentName;
        stream_window_size = windowSize;
        stream_file = streamFile;
        input_order_round_robin = input_order;
        max_num_processes = max_processes;  
        element_is_vertex = elementIsVertex; 
        use_hdrf = useHDRF;   
        staggered_streams = staggeredStreams; 
        use_balance_cost = useBalanceCost;
        lambda = l;

        split_and_configure_stream();

	}
	virtual ~ParallelStreamingPartitioning() {
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

        // TODO: should allow for elements to be vertices or hyperedges, based on parameter element_is_vertex
        int num_elements = element_is_vertex ? num_vertices : num_hyperedges;
        int num_pins = element_is_vertex ? num_hyperedges : num_vertices;

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
        // write header: NUM_HYPEREDGES NUM_VERTICES
        fprintf(fp,"%i %i",num_elements,num_pins);
        fprintf(fp,"\n");

        // read reminder of file (one line per vertex)
        int local_stream_size = num_elements / num_processes + ((num_elements % num_processes > process_id) ? 1 : 0);
        int first_element_in_stream = (num_elements / num_processes * process_id + std::min(num_elements % num_processes,process_id));
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
        int num_elements = element_is_vertex ? num_vertices : num_hyperedges;

        if(is_active) {
            
            // initialise element weight values
            int* element_wgt = (int*)calloc(num_elements,sizeof(int));
            for(int ii =0; ii < num_elements; ii++) {
                element_wgt[ii] = 1;
            }
            bool local_replica_degree_updates_only = false;
            *iterations = PRAW::ParallelStreaming(experiment_name,partitioning, num_partitions, partitioning_comm, hgraph_part_file.c_str(), element_wgt, imbalance_tolerance,local_replica_degree_updates_only,stream_window_size,input_order_round_robin, use_hdrf,staggered_streams,use_balance_cost,lambda);

            // clean up operations
            free(element_wgt);

            // remove graph file
            if(remove(hgraph_part_file.c_str()) != 0 )
                printf( "Error deleting temporary hgraph file %s\n",hgraph_part_file.c_str() );
        }
        MPI_Barrier(MPI_COMM_WORLD);


        // broadcast results to not participating processes
        MPI_Bcast(partitioning, num_elements, MPI_LONG, 0,MPI_COMM_WORLD);

        
	}

private:
    char* experiment_name = NULL;
    int stream_window_size = 1;
    char* stream_file = NULL;
    std::string hgraph_part_file;
    bool input_order_round_robin = true;
    MPI_Comm partitioning_comm;
    int max_num_processes;
    bool is_active = true;
    bool element_is_vertex;
    bool use_hdrf = false;
    bool staggered_streams = true;
    bool use_balance_cost = false;
    float lambda = 1.0f;
};


#endif
