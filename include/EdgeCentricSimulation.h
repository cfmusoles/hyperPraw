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
            // TODO SIMULATION CODE
            ////////////////////////

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
                float max_vertex_imbalance;
                float max_hedge_imbalance;
                double total_sim_comm_cost;
                
                std::string filename = graph_file;
                
                PRAW::getEdgeCentricPartitionStatsFromFile(partitioning, num_processes, filename, NULL,comm_cost_matrix,
                                        &vertex_replication_factor, &max_vertex_imbalance, &max_hedge_imbalance, &total_sim_comm_cost);
                
                printf("Partition time %.2fs\nVertex replication factor %.3f, %.3f (vertex imbalance), %.3f (hedge imbalance), %.0f (sim comm cost), \nMessages sent %li\n",partition_timer,vertex_replication_factor,max_vertex_imbalance,max_hedge_imbalance,total_sim_comm_cost,total_edge_messages_sent);
                
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
                        fprintf(fp,"%s,%s,%s,%s,%s,%s,%s\n","Partition time","Sim time","Vertex replication factor","Vertex imbalance","Hedge imbalance","Sim Comm cost","Sim Messages sent");
                    fprintf(fp,"%.3f,%.3f,%.3f,%.3f,%.3f,%.0f,%li\n",partition_timer,total_edge_sim_time,vertex_replication_factor,max_vertex_imbalance,max_hedge_imbalance,total_sim_comm_cost,total_edge_messages_sent);
                }
                fclose(fp);

                // clean up operations
                for(int ii=0; ii < num_processes; ii++) {
                    free(comm_cost_matrix[ii]);
                }
                free(comm_cost_matrix);
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