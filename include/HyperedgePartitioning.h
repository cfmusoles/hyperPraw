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
	
	HyperedgePartitioning(char* experimentName, char* graph_file, int iterations, float imbalance_tolerance, char* comm_bandwidth_file, bool useBandwidth, bool saveHistory) : Partitioning(graph_file,imbalance_tolerance,false) {
		experiment_name = experimentName;
        comm_bandwidth_filename = comm_bandwidth_file;
        use_bandwidth_file = useBandwidth;
        save_partitioning_history = saveHistory;
        max_iterations = iterations;
	}
	virtual ~HyperedgePartitioning() {}
	
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
            PRAW::get_comm_cost_matrix_from_bandwidth(comm_bandwidth_filename,comm_cost_matrix,num_processes,false);
        else 
            PRAW::get_comm_cost_matrix_from_bandwidth(NULL,comm_cost_matrix,num_processes,false);
        
        // initialise vertex weight values
        int* he_wgt = (int*)calloc(num_hyperedges,sizeof(int));
        for(int ii =0; ii < num_hyperedges; ii++) {
            he_wgt[ii] = 1;
        }

        // divide original hMetis file into as many streams as processes
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
        
        // load and parse full graph
        std::ifstream istream(hgraph_name);
        
        if(!istream) {
            printf("Error while opening hMETIS file %s\n",hgraph_name);
            return;
        }

        std::string line;
        // ignore header
        std::getline(istream,line);
        // write header: NUM_HYPEREDGES NUM_VERTICES
        //  num hyperedges == number of vertices, since each hyperedge represents a presynaptic neuron and all its connecting post synaptic neighbours
        fprintf(fp,"%i %i",num_hyperedges,num_vertices);
        fprintf(fp,"\n");

        // read reminder of file (one line per hyperedge)
        int counter = 0;
        while(std::getline(istream,line)) {
            if(counter % num_processes == process_id) {
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


        *iterations = PRAW::ParallelHDRF(experiment_name,partitioning, comm_cost_matrix, hgraph_file, he_wgt, max_iterations, imbalance_tolerance, save_partitioning_history);

        /*std::vector<std::vector<int> > hyperedges;
        std::vector<std::vector<int> > hedge_ptr;
        std::string filename = hgraph_name;
        PRAW::load_hypergraph_from_file(filename, &hyperedges, &hedge_ptr);

        *iterations = PRAW::ParallelVertexPartitioning(experiment_name,partitioning, comm_cost_matrix, hedge_ptr.size(), hyperedges.size(), &hyperedges, he_wgt, max_iterations, imbalance_tolerance, save_partitioning_history);
        */

        // remove graph file
        if( remove(hgraph_file.c_str()) != 0 )
            printf( "Error deleting temporary hgraph file %s\n",hgraph_file.c_str() );
        
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
    bool use_bandwidth_file = false;
    bool save_partitioning_history; 
    int max_iterations;
};


#endif



