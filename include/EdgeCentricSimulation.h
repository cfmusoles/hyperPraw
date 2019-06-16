// Hyperedge centric simulations to evaluate hypergraph partitioning
// Hyperedges are assigned to partitions
// Communication is proportional to vertex replica factor

#ifndef EDGECENTRICSIMULATION__H
#define EDGECENTRICSIMULATION__H

#include <mpi.h>
#include <cstdio>
#include <stdlib.h>
#include <vector>
#include <set>
#include "PRAW.h"
#include <iterator>
#include <numeric>
#include "SimulationUtils.h"

namespace EdgeCentricSimulation {

    void runSimulation(char* experiment_name, char* graph_file, char* part_method, char* bandwidth_file,  idx_t* partitioning, float partition_timer, int num_vertices, int simulation_iterations, int sim_steps, int hedge_sim_steps_multiplier, int message_size, bool proportional_comm_cost) {

        int process_id;
        int num_processes;
        MPI_Comm_rank(MPI_COMM_WORLD,&process_id);
        MPI_Comm_size(MPI_COMM_WORLD,&num_processes);

        // load graph and get vertex replicas
        std::ifstream istream(graph_file);
            
        if(!istream) {
            printf("Error while opening hMETIS file %s\n",graph_file);
        }
        std::string line;
        // process header
        std::getline(istream,line);
        char str[line.length() + 1]; 
        strcpy(str, line.c_str()); 
        char* token = strtok(str, " "); 
        std::vector<int> tokens;
        while (token != NULL) { 
            tokens.push_back(atoi(token)); 
            token = strtok(NULL, " "); 
        } 
        num_vertices = tokens[1];
        int num_hyperedges = tokens[0];

        std::vector<std::set<int> > vertex_replicas(num_vertices);
        int current_hyperedge = 0;

        while(std::getline(istream,line)) {
            char str[line.length() + 1]; 
            strcpy(str, line.c_str()); 
            char* token = strtok(str, " "); 
            std::vector<int> vertices;
            while (token != NULL) { 
                int vid = atoi(token) - 1;
                vertex_replicas[vid].insert(partitioning[current_hyperedge]);
                token = strtok(NULL, " "); 
            } 
            current_hyperedge++;
        }


        ////////////////////////////////////////////////////////
        // repeat simulation with same partitioning distribution
        ////////////////////////////////////////////////////////
        for(int iteration=0; iteration < simulation_iterations; iteration++) {

            MPI_Barrier(MPI_COMM_WORLD);

            //////////////////////////////////////////////////////////////////////////////////
            // Simulation
            //////////////////////////////////////////////////////////////////////////////////

            // Parallel communication should increase with
            //        Vertex replica
            int* buffer = (int*)malloc(sizeof(int)*message_size);
            double timer = MPI_Wtime();
            long int messages_sent = 0;
        #ifdef SAVE_COMM_COST
            int* sent_communication = (int*)calloc(num_processes,sizeof(int));
        #endif   
            
            ////////////////////////
            // SIMULATION CODE
            ////////////////////////
            for(int tt = 0; tt < sim_steps; tt++) {
                for(int vid=0; vid < num_vertices; vid++) {
                    std::set<int>::iterator sit;
                    for(sit =vertex_replicas[vid].begin(); sit != vertex_replicas[vid].end(); sit++) {
                        int sender_id = *sit;
                        if(sender_id == process_id) {
                            // send data
                            std::set<int>::iterator dit;
                            for(dit =vertex_replicas[vid].begin(); dit != vertex_replicas[vid].end(); dit++) {
                                int dest_part = *dit;
                                if(dest_part == sender_id) continue;
                                messages_sent++;
        #ifdef SAVE_COMM_COST
                                sent_communication[dest_part] += 1;
        #endif
                                MPI_Send(buffer,message_size,MPI_INT,dest_part,vid,MPI_COMM_WORLD);
                            }
                        } else {
                            // receive data if on the list
                            std::set<int>::iterator rit;
                            for(rit =vertex_replicas[vid].begin(); rit != vertex_replicas[vid].end(); rit++) {
                                int receiver_id = *rit;
                                if(receiver_id == process_id && receiver_id != sender_id) {
                                    MPI_Recv(buffer,message_size,MPI_INT,sender_id,vid,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                                }
                            }
                        }
                    }
                }
            }

            // wait for all processes to finish
            MPI_Barrier(MPI_COMM_WORLD);

            double total_edge_sim_time = MPI_Wtime() - timer;
            //total number of messages exchanged
            long int total_edge_messages_sent;
            MPI_Allreduce(&messages_sent, &total_edge_messages_sent, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

        #ifdef SAVE_COMM_COST
            SimulationUtils::storeSimCommunication(sent_communication,process_id,num_processes,0,experiment_name,graph_file,part_method);
        #endif

            //////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////

            //Store metrics
            //    simulation time
            //    communication time
            //    partitioning stats (hedge cut, SOED, absorption)
            
            if(process_id == 0) {
                printf("%i: simulation time (%i steps): %f secs\n",process_id,sim_steps,total_edge_sim_time);
                // used to calculate the theoretical cost of communication
                // if bandwidth file is not provided, then assumes all costs are equal
                // initialise comm cost matrix (for theoretical cost analysis)
                double** comm_cost_matrix = (double**)malloc(sizeof(double*) * num_processes);
                for(int ii=0; ii < num_processes; ii++) {
                    comm_cost_matrix[ii] = (double*)calloc(num_processes,sizeof(double));
                }
                PRAW::get_comm_cost_matrix_from_bandwidth(bandwidth_file,comm_cost_matrix,num_processes,proportional_comm_cost);
                    
                // calculate partitioning stats
                float vertex_replication_factor;
                float max_hedge_imbalance;
                float hedgecut;
                
                std::string filename = graph_file;

                int* he_wgt = (int*)calloc(num_hyperedges,sizeof(int));
                for(int ii =0; ii < num_hyperedges; ii++) {
                    he_wgt[ii] = 1;
                }
                
                PRAW::getEdgeCentricPartitionStatsFromFile(partitioning, num_processes, filename, he_wgt,comm_cost_matrix,
                                        &vertex_replication_factor, &max_hedge_imbalance, &hedgecut);
                
                printf("Partition time %.2fs\nVertex replication factor %.3f, %.3f (hedge imbalance), %.3f (hyperedge cut)\nMessages sent %li\n",partition_timer,vertex_replication_factor,max_hedge_imbalance,hedgecut,total_edge_messages_sent);
                
                // store stats in file
                filename = experiment_name;
                filename += "_";
                std::string graph_string = graph_file;
                filename += PRAW::getFileName(graph_string);
                filename += "_";
                filename += part_method;
                char str_int[16];
                sprintf(str_int,"%i",num_processes);
                filename += "__";
                filename +=  str_int;
                bool fileexists = access(filename.c_str(), F_OK) != -1;
                FILE *fp = fopen(filename.c_str(), "ab+");
                if(fp == NULL) {
                    printf("Error when storing results into file\n");
                } else {
                    if(!fileexists) // file does not exist, add header
                        fprintf(fp,"%s,%s,%s,%s,%s,%s\n","Partition time","Sim time","Vertex replication factor","Hedge imbalance","Hedgecut","Sim Messages sent");
                    fprintf(fp,"%.3f,%.3f,%.3f,%.3f,%.0f,%li\n",partition_timer,total_edge_sim_time,vertex_replication_factor,max_hedge_imbalance,hedgecut,total_edge_messages_sent);
                }
                fclose(fp);

                // clean up operations
                for(int ii=0; ii < num_processes; ii++) {
                    free(comm_cost_matrix[ii]);
                }
                free(comm_cost_matrix);
                free(he_wgt);
            }

            // clean up
            free(buffer);
        #ifdef SAVE_COMM_COST
            free(sent_communication);
        #endif


        }
    }

}

#endif