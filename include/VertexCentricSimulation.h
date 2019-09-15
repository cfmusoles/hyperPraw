// Vertex centric simulations to evaluate hypergraph partitioning
// Vertices are assigned to partitions
// Communication is proportional to edge cut and hyperedge cut

#ifndef VERTEXCENTRICSIMULATION__H
#define VERTEXCENTRICSIMULATION__H

#include <mpi.h>
#include <cstdio>
#include <stdlib.h>
#include <vector>
#include <chrono>
#include <thread>
#include <set>
#include <random>
#include "PRAW.h"
#include <iterator>
#include <numeric>
#include "SimulationUtils.h"

namespace VertexCentricSimulation {

    
    void runSimulation(char* experiment_name, char* graph_file, char* part_method, char* bandwidth_file,  idx_t* partitioning, float partition_timer, int num_vertices, int simulation_iterations, int edge_sim_steps, int hedge_sim_steps, int fake_compute_time, float fake_compute_std, int message_size, bool proportional_comm_cost) {

        int process_id;
        int num_processes;
        MPI_Comm_rank(MPI_COMM_WORLD,&process_id);
        MPI_Comm_size(MPI_COMM_WORLD,&num_processes);

        // set simulation to test hypergraph partitioning
        // load model (only local hyperedges loaded)
        std::vector<std::vector<int> > hyperedges;
        std::vector<std::vector<int> > hedge_ptr;
        PRAW::load_hypergraph_from_file_dist_CSR(graph_file, &hyperedges, &hedge_ptr, process_id, partitioning);


        ////////////////////////////////////////////////////////
        // repeat simulation with same partitioning distribution
        ////////////////////////////////////////////////////////
        for(int iteration=0; iteration < simulation_iterations; iteration++) {

            MPI_Barrier(MPI_COMM_WORLD);

            //////////////////////////////////////////////////////////////////////////////////
            // Edge sim
            //////////////////////////////////////////////////////////////////////////////////

            // Parallel communication should increase with
            //        Hedge cut
            //        SOED
            //    Parallel communication should decrease with
            //        Absorption
            int* buffer = (int*)malloc(sizeof(int)*message_size);
            double timer = MPI_Wtime();
            long int messages_sent = 0;
        #ifdef SAVE_COMM_COST
            int* sent_communication = (int*)calloc(num_processes,sizeof(int));
        #endif   

            // random generator used for fake computation
            std::default_random_engine generator;
            std::normal_distribution<double> distribution(fake_compute_time,fake_compute_std);
            
            for(int tt = 0; tt < edge_sim_steps; tt++) {
                
                // fake compute based on stochastic time sleep
                if(fake_compute_time > 0) {
                    // do compute for every local vertex
                    for(int vid=0; vid < num_vertices; vid++) {
                        if(partitioning[vid] == process_id) {
                            int t = distribution(generator);
                            t = std::max(0,t);
                            std::this_thread::sleep_for(std::chrono::microseconds(t));
                        }
                    }
                }

                // for each local hyperedge
                //      if any vertex is not local, add destination to target list
                //      send messages all to all for processes in target list plus local
                //          use hedge id as flag for the messages
                //          send messages in a ring order
                for(int he_id = 0; he_id < hyperedges.size(); he_id++) {
                    // communication is proportional to edgecut
                    for(int vid = 0; vid < hyperedges[he_id].size(); vid++) {
                        int origin_vertex = hyperedges[he_id][vid];
                        int origin_part = partitioning[origin_vertex];
                        for(int did = 0; did < hyperedges[he_id].size(); did++) {
                            int dest_vertex = hyperedges[he_id][did];
                            int dest_part = partitioning[dest_vertex];
                            if (origin_vertex == dest_vertex || origin_part == dest_part) continue;
                            if(origin_part == process_id ) {
                                // send
                                messages_sent++;
        #ifdef SAVE_COMM_COST
                                sent_communication[dest_part] += 1;
        #endif
                                MPI_Send(buffer,message_size,MPI_INT,dest_part,he_id,MPI_COMM_WORLD);
                            } else {
                                if(dest_part == process_id) {
                                    // receive
                                    MPI_Recv(buffer,message_size,MPI_INT,origin_part,he_id,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                                }
                            }
                        }
                    }

                    // why is this needed?
                    //MPI_Barrier(MPI_COMM_WORLD);
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
            // Hedge sim
            //////////////////////////////////////////////////////////////////////////////////

            MPI_Barrier(MPI_COMM_WORLD);
            messages_sent = 0;
        #ifdef SAVE_COMM_COST
            memset(sent_communication,0,num_processes * sizeof(int));
        #endif
            int* neighbouring_partitions = (int*)calloc(num_processes,sizeof(int));
            timer = MPI_Wtime();
            
            for(int tt = 0; tt < hedge_sim_steps; tt++) {
                // for each local hyperedge
                //      if any vertex is not local, add destination to target list
                //      send messages all to all for processes in target list plus local
                //          use hedge id as flag for the messages
                //          send messages in a ring order
                for(int he_id = 0; he_id < hyperedges.size(); he_id++) {
                    // communication is proportional to hedge cut alone
                    memset(neighbouring_partitions,0,num_processes * sizeof(int));
                    std::set<int> partitions;
                    for(int vid = 0; vid < hyperedges[he_id].size(); vid++) {
                        int dest_vertex = hyperedges[he_id][vid];
                        partitions.insert(partitioning[dest_vertex]);
                        neighbouring_partitions[partitioning[dest_vertex]] += 1;
                    }
                    if(partitions.size() > 1) {
                        for (std::set<int>::iterator it=partitions.begin(); it!=partitions.end(); ++it) {
                            int sender_id = *it;
                            if(sender_id == process_id) { // turn to send messages
                                // send as many messages as targets (not counting local)
                                for (std::set<int>::iterator receiver=partitions.begin(); receiver!=partitions.end(); ++receiver) {
                                    int receiver_id = *receiver;
                                    if(sender_id == receiver_id) continue;
                                    messages_sent++;
                                    int m_size = message_size * neighbouring_partitions[receiver_id];
                                    int* bf = (int*)malloc(m_size * sizeof(int));
        #ifdef SAVE_COMM_COST
                                    sent_communication[receiver_id] += 1;
        #endif
                                    //MPI_Send(buffer,message_size,MPI_INT,receiver_id,he_id,MPI_COMM_WORLD);
                                    MPI_Send(bf,m_size,MPI_INT,receiver_id,he_id,MPI_COMM_WORLD);
                                    free(bf);
                                }
                            } else { // turn to receive messages
                                // receive one message from sender id
                                //MPI_Recv(buffer,message_size,MPI_INT,sender_id,he_id,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                                int m_size = message_size * neighbouring_partitions[process_id];
                                int* bf = (int*)malloc(m_size * sizeof(int));
                                MPI_Recv(bf,m_size,MPI_INT,sender_id,he_id,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                                free(bf);
                            }
                        }
                    }

                    // why is this needed?
                    //MPI_Barrier(MPI_COMM_WORLD);
                }
            }
            // wait for all processes to finish
            MPI_Barrier(MPI_COMM_WORLD);

            double total_hedge_sim_time = MPI_Wtime() - timer;
            //total number of messages exchanged
            long int total_hedge_messages_sent;
            MPI_Allreduce(&messages_sent, &total_hedge_messages_sent, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

#ifdef SAVE_COMM_COST
            SimulationUtils::storeSimCommunication(sent_communication,process_id,num_processes,1,experiment_name,graph_file,part_method);
#endif
            //////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////

            //Store metrics
            //    simulation time
            //    communication time
            //    partitioning stats (hedge cut, SOED, absorption)
            
            if(process_id == 0) {
                printf("%i: Edge simulation time (%i steps): %f secs, Hedge simulation time: %f secs\n",process_id,edge_sim_steps,total_edge_sim_time,total_hedge_sim_time);
                // used to calculate the theoretical cost of communication
                // if bandwidth file is not provided, then assumes all costs are equal
                // initialise comm cost matrix (for theoretical cost analysis)
                double** comm_cost_matrix = (double**)malloc(sizeof(double*) * num_processes);
                for(int ii=0; ii < num_processes; ii++) {
                    comm_cost_matrix[ii] = (double*)calloc(num_processes,sizeof(double));
                }
                PRAW::get_comm_cost_matrix_from_bandwidth(bandwidth_file,comm_cost_matrix,num_processes,proportional_comm_cost);
                    
                // calculate partitioning stats
                float hyperedges_cut_ratio;
                float edges_cut_ratio;
                int soed;
                float absorption;
                float max_imbalance;
                double total_edge_comm_cost;
                double total_hedge_comm_cost;
                float vertex_replication_factor;
                
                std::string filename = graph_file;
                
                PRAW::getVertexCentricPartitionStatsFromFile(partitioning, num_processes, filename, NULL,comm_cost_matrix,
                                        &hyperedges_cut_ratio, &edges_cut_ratio, &soed, &absorption, &max_imbalance, &total_edge_comm_cost,&total_hedge_comm_cost,&vertex_replication_factor);
                
                printf("Partition time %.2fs\nHedgecut, %.3f, %.3f (cut net), %i (SOED), %.1f (absorption) %.3f (max imbalance), %.0f (edge comm cost), %.0f (hedge comm cost)\nEdgesim messages sent %li,Hedgesim messages sent %li, Vertex replication factor %.3f\n",partition_timer,hyperedges_cut_ratio,edges_cut_ratio,soed,absorption,max_imbalance,total_edge_comm_cost,total_hedge_comm_cost,total_edge_messages_sent,total_hedge_messages_sent,vertex_replication_factor);
                
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
                        fprintf(fp,"%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n","Partition time","Edge Sim time","Hedge sim time","Hedge cut ratio","Cut net","SOED","Absorption","Max imbalance","Edge Comm cost","Hedge comm cost","Edgesim Messages sent","Hedgesim Messages sent","Vertex replication factor");
                    fprintf(fp,"%.3f,%.3f,%.3f,%.3f,%.3f,%i,%.1f,%.3f,%.0f,%.0f,%li,%li,%.3f\n",partition_timer,total_edge_sim_time,total_hedge_sim_time,hyperedges_cut_ratio,edges_cut_ratio,soed,absorption,max_imbalance,total_edge_comm_cost,total_hedge_comm_cost,total_edge_messages_sent,total_hedge_messages_sent,vertex_replication_factor);
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