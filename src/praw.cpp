// Test harness for SPAAW (Streaming parallel Partitioning Architecture AWare)
//#define VERBOSE                 // extra debug info printed out during runtime
#define SAVE_COMM_COST      // store actual p2p communication based on partitioning

#include <mpi.h>
#include <cstdio>
#include <stdlib.h>
#include <vector>
#include <set>
#include "PRAW.h"
#include "Partitioning.h"
#include "RandomPartitioning.h"
#include "ZoltanPartitioning.h"
#include "HyperPRAWPartitioning.h"
#include <iterator>
#include <numeric>

void storeSimCommunication(int* sent_communication,int process_id, int num_processes, int mode, char* experiment_name, char* graph_file, char* part_method) {
    // gather all results in node 0
    if(process_id == 0) {
        int** comm = (int**)malloc(num_processes * sizeof(int*));
        for(int ii=0; ii < num_processes; ii++) {
            comm[ii] = (int*)calloc(num_processes,sizeof(int));
        }
        memcpy(comm[0],sent_communication,sizeof(int) * num_processes);
        for(int ii=1; ii < num_processes; ii++) {
            MPI_Recv(comm[ii],num_processes,MPI_INT,ii,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        // store matrix in a file
        std::string filename = experiment_name;
        filename += "_";
        std::string graph_string = graph_file;
        filename += PRAW::getFileName(graph_string);
        filename += "_";
        filename += part_method;
        if(mode == 0) {
            filename += "_edgeSim";
        } else if (mode == 1) {
            filename += "_hedgeSim";
        }
        filename += "_comm_cost";
        char str_int[16];
        sprintf(str_int,"%i",num_processes);
        filename += "__";
        filename +=  str_int;
        bool fileexists = access(filename.c_str(), F_OK) != -1;
        FILE *fp = fopen(filename.c_str(), "w");
        if(fp == NULL) {
            printf("Error when storing comm cost into file\n");
        } else {
            for(int jj=0; jj < num_processes; jj++) {
                for(int ii=0; ii < num_processes; ii++) {
                    fprintf(fp,"%i ",comm[jj][ii]);
                }
                fprintf(fp,"\n");
            }
        }
        fclose(fp);
    } else {
        MPI_Send(sent_communication,num_processes,MPI_INT,0,0,MPI_COMM_WORLD);
    }
    
    
}


int main(int argc, char** argv) {

    // initialise MPI
    MPI_Init(&argc,&argv);
    int process_id;
    int num_processes;
    MPI_Comm_rank(MPI_COMM_WORLD,&process_id);
    MPI_Comm_size(MPI_COMM_WORLD,&num_processes);

    // DEFAULT PARAMETERS
    char* experiment_name = NULL;
    char* graph_file = NULL;
    int iterations = 1;
    float imbalance_tolerance = 1.05f;
    char* bandwidth_file = NULL;
    bool use_bandwidth_in_partitioning = false;
    int rand_seed = time(NULL);
    char* part_method = NULL;
    int sim_steps  = 100;
    int message_size = 1;
    int stopping_condition = 0;
    bool proportional_comm_cost = false;
    float ta_refinement = 1.0f;
    bool save_partitioning_history = false;
    int simulation_iterations = 1;

    // getting command line parameters
    extern char *optarg;
	extern int optind, opterr, optopt;
	int c;
	while( (c = getopt(argc,argv,"n:h:i:m:b:Ws:p:t:k:o:c:r:Hq:")) != -1 ) {
		switch(c) {
			case 'n': // test name
				experiment_name = optarg;
				break;
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
			case 'b': // network bandwidth file
				bandwidth_file = optarg;
				break;
			case 'W': // use bandwith file in partitioning
				use_bandwidth_in_partitioning = true;
				break;
            case 'p': // partitioning method
				part_method = optarg;
				break;
            case 't': // simulated steps
				sim_steps = atoi(optarg);
				break;
            case 'k': // message size during simulation
				message_size = atoi(optarg);
				break;
            case 'o': // type of stopping condition
				stopping_condition = atoi(optarg);
				break;
            case 'c': // type of comm cost mapping
				proportional_comm_cost = atoi(optarg) == 1;
				break;
            case 'r': // tempering alpha when within imbalance tolerance
				ta_refinement = atoi(optarg) * 0.001f;
				break;
			case 'H': // save streaming partitioning history
				save_partitioning_history = true;
				break;
            case 'q': // simulation iterations
				simulation_iterations = atoi(optarg);
				break;
		}
	}

    // set and propagate random seed
    MPI_Bcast(&rand_seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    srand(rand_seed);

    double partition_timer = MPI_Wtime();

    // Create partition object to hold connections and partitioning information across processes
	Partitioning* partition;
	if(strcmp(part_method,"zoltan") == 0) {
		PRINTF("%i: Partitioning: zoltan\n",process_id);
		partition = new ZoltanPartitioning(graph_file,imbalance_tolerance);
	} else if(strcmp(part_method,"prawP") == 0) {  
		PRINTF("%i: Partitioning: parallel hyperPRAW\n",process_id);
        //Partitioning* p1 = new ZoltanPartitioning(graph_file,imbalance_tolerance);
        //Partitioning* p1 = new HyperPRAWPartitioning(graph_file,imbalance_tolerance,iterations,NULL,true,false,false);
        //p1->perform_partitioning(num_processes,process_id);
		partition = new HyperPRAWPartitioning(experiment_name,graph_file,imbalance_tolerance,ta_refinement,iterations,bandwidth_file,true,use_bandwidth_in_partitioning,false,stopping_condition,proportional_comm_cost,save_partitioning_history);
        //memcpy(partition->partitioning,p1->partitioning,partition->num_vertices * sizeof(idx_t));
        //free(p1);
    } else if(strcmp(part_method,"prawSref") == 0) {  
		PRINTF("%i: Partitioning: sequential refinement hyperPRAW\n",process_id);
        Partitioning* p1 = new ZoltanPartitioning(graph_file,imbalance_tolerance);
        p1->perform_partitioning(num_processes,process_id);
		partition = new HyperPRAWPartitioning(experiment_name,graph_file,imbalance_tolerance,ta_refinement,iterations,bandwidth_file,false,use_bandwidth_in_partitioning,false,stopping_condition,proportional_comm_cost,save_partitioning_history);
        memcpy(partition->partitioning,p1->partitioning,partition->num_vertices * sizeof(idx_t));
        free(p1);
    } else if(strcmp(part_method,"prawS") == 0) {  
		PRINTF("%i: Partitioning: sequential hyperPRAW\n",process_id);
		partition = new HyperPRAWPartitioning(experiment_name,graph_file,imbalance_tolerance,ta_refinement,iterations,bandwidth_file,false,use_bandwidth_in_partitioning,true,stopping_condition,proportional_comm_cost,save_partitioning_history);
	} else { // default is random
		PRINTF("%i: Partitioning: random\n",process_id);
		partition = new RandomPartitioning(graph_file,imbalance_tolerance);
	}
    srand(rand_seed);
	partition->perform_partitioning(num_processes,process_id);

    partition_timer = MPI_Wtime() - partition_timer;

    // set simulation to test hypergraph partitioning
    // load model (only local hyperedges loaded)
    std::vector<std::vector<int> > hyperedges;
    std::vector<std::vector<int> > hedge_ptr;
    PRAW::load_hypergraph_from_file_dist_CSR(graph_file, &hyperedges, &hedge_ptr, process_id, partition->partitioning);


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
        
        for(int tt = 0; tt < sim_steps; tt++) {
            // for each local hyperedge
            //      if any vertex is not local, add destination to target list
            //      send messages all to all for processes in target list plus local
            //          use hedge id as flag for the messages
            //          send messages in a ring order
            for(int he_id = 0; he_id < hyperedges.size(); he_id++) {
                // communication is proportional to edgecut
                for(int vid = 0; vid < hyperedges[he_id].size(); vid++) {
                    int origin_vertex = hyperedges[he_id][vid];
                    int origin_part = partition->partitioning[origin_vertex];
                    for(int did = 0; did < hyperedges[he_id].size(); did++) {
                        int dest_vertex = hyperedges[he_id][did];
                        int dest_part = partition->partitioning[dest_vertex];
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
                MPI_Barrier(MPI_COMM_WORLD);
            }
        }
        // wait for all processes to finish
        MPI_Barrier(MPI_COMM_WORLD);

        double total_edge_sim_time = MPI_Wtime() - timer;
        //total number of messages exchanged
        long int total_edge_messages_sent;
        MPI_Allreduce(&messages_sent, &total_edge_messages_sent, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    #ifdef SAVE_COMM_COST
        storeSimCommunication(sent_communication,process_id,num_processes,0,experiment_name,graph_file,part_method);
    #endif

        //////////////////////////////////////////////////////////////////////////////////
        // Hedge sim
        //////////////////////////////////////////////////////////////////////////////////

        MPI_Barrier(MPI_COMM_WORLD);
        timer = MPI_Wtime();
        messages_sent = 0;
    #ifdef SAVE_COMM_COST
        memset(sent_communication,0,num_processes * sizeof(int));
    #endif
        
        for(int tt = 0; tt < sim_steps; tt++) {
            // for each local hyperedge
            //      if any vertex is not local, add destination to target list
            //      send messages all to all for processes in target list plus local
            //          use hedge id as flag for the messages
            //          send messages in a ring order
            for(int he_id = 0; he_id < hyperedges.size(); he_id++) {
                // communication is proportional to hedge cut alone
                std::set<int> partitions;
                for(int vid = 0; vid < hyperedges[he_id].size(); vid++) {
                    int dest_vertex = hyperedges[he_id][vid];
                    partitions.insert(partition->partitioning[dest_vertex]);
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
    #ifdef SAVE_COMM_COST
                                sent_communication[receiver_id] += 1;
    #endif
                                MPI_Send(buffer,message_size,MPI_INT,receiver_id,he_id,MPI_COMM_WORLD);
                            }
                        } else { // turn to receive messages
                            // receive one message from sender id
                            MPI_Recv(buffer,message_size,MPI_INT,sender_id,he_id,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                        }
                    }
                }

                // why is this needed?
                MPI_Barrier(MPI_COMM_WORLD);
            }
        }
        // wait for all processes to finish
        MPI_Barrier(MPI_COMM_WORLD);

        double total_hedge_sim_time = MPI_Wtime() - timer;
        //total number of messages exchanged
        long int total_hedge_messages_sent;
        MPI_Allreduce(&messages_sent, &total_hedge_messages_sent, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

        storeSimCommunication(sent_communication,process_id,num_processes,1,experiment_name,graph_file,part_method);

        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////

        //Store metrics
        //    simulation time
        //    communication time
        //    partitioning stats (hedge cut, SOED, absorption)
        
        if(process_id == 0) {
            printf("%i: Edge simulation time (%i steps): %f secs, Hedge simulation time: %f secs\n",process_id,sim_steps,total_edge_sim_time,total_hedge_sim_time);
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
            
            /*hyperedges.clear();
            hyperedges.swap(hyperedges);
            hedge_ptr.clear();
            hedge_ptr.swap(hedge_ptr);*/
            std::string filename = graph_file;
            
            PRAW::getPartitionStatsFromFile(partition->partitioning, num_processes, partition->num_vertices, filename, NULL,comm_cost_matrix,
                                    &hyperedges_cut_ratio, &edges_cut_ratio, &soed, &absorption, &max_imbalance, &total_edge_comm_cost,&total_hedge_comm_cost);
            
            printf("Partition time %.2fs\nHedgecut, %.3f, %.3f (cut net), %i (SOED), %.1f (absorption) %.3f (max imbalance), %.0f (edge comm cost), %.0f (hedge comm cost)\nEdgesim messages sent %li,Hedgesim messages sent %li\n",partition_timer,hyperedges_cut_ratio,edges_cut_ratio,soed,absorption,max_imbalance,total_edge_comm_cost,total_hedge_comm_cost,total_edge_messages_sent,total_hedge_messages_sent);
            
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
                    fprintf(fp,"%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n","Partition time","Edge Sim time","Hedge sim time","Hedge cut ratio","Cut net","SOED","Absorption","Max imbalance","Edge Comm cost","Hedge comm cost","Edgesim Messages sent","Hedgesim Messages sent");
                fprintf(fp,"%.3f,%.3f,%.3f,%.3f,%.3f,%i,%.1f,%.3f,%.0f,%.0f,%li,%li\n",partition_timer,total_edge_sim_time,total_hedge_sim_time,hyperedges_cut_ratio,edges_cut_ratio,soed,absorption,max_imbalance,total_edge_comm_cost,total_hedge_comm_cost,total_edge_messages_sent,total_hedge_messages_sent);
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
    

    free(partition);

    // finalise MPI and application
    MPI_Finalize();
    return 0;
}
