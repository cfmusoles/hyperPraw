// Parallel Restreaming Architecture aWare partitioning algorithm

#ifndef PRAW__H
#define PRAW__H

#define SAVE_HISTORY        // stores partitioning iteration history in file
#define SAVE_COMM_COST      // store theoretical p2p communication based on partitioning

#include <vector>
#include <time.h>
#include <mpi.h>
#include <cmath>
#include <cstring>
#include <limits>
#include <set>
#include <cstdio>
#include <string>
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <unistd.h>
#include <libgen.h>

#ifdef VERBOSE
#define PRINTF(...) printf(__VA_ARGS__)
#else
#define PRINTF(...) 
#endif

const int MASTER_NODE = 0;

typedef int64_t idx_t; // needs to probably match the type in METIS

namespace PRAW {
    std::string getFileName(std::string filePath)
    {
        int n = filePath.length();  
      
        // declaring character array 
        char char_array[n+1];  
        
        // copying the contents of the  
        // string to char array 
        strcpy(char_array, filePath.c_str());  
        return std::string(basename(char_array));
    }

    float calculateImbalance(idx_t* partitioning, int num_processes, int num_vertices, int* vtx_wgt) {
        long int* workload = (long int*)calloc(num_processes,sizeof(long int));
        long int total_workload = 0;
        for(int vid=0; vid < num_vertices; vid++) {
            workload[partitioning[vid]] += vtx_wgt[vid];
            total_workload += vtx_wgt[vid];
        }
        float max_imbalance = 0;
        float expected_work = total_workload / num_processes;
        for(int ii=0; ii < num_processes; ii++) {
            PRINTF("%i: %li\n",ii,workload[ii]);
            float imbalance = (workload[ii] / expected_work);
            if(imbalance > max_imbalance)
                max_imbalance = imbalance;
        }

        free(workload);

        return  max_imbalance;
    }

    void getPartitionStats(idx_t* partitioning, int num_processes, int num_vertices, std::vector<std::vector<int> >* hyperedges, std::vector<std::vector<int> >* hedge_ptr, int* vtx_wgt,double** comm_cost_matrix, // input
                            float* hyperedges_cut_ratio, float* edges_cut_ratio, int* soed, float* absorption, float* max_imbalance, double* total_comm_cost) { // output
        
        // For hypergraph partitioning, two common metrics for objective functions Heuer 2017 (connectivity vs hyperedge cut):
        // communication cost estimate (given network bandwidth estimates)
        // quality measurements from hMETIS: Sum of External Degrees, Absorption
        *hyperedges_cut_ratio=0;
        *edges_cut_ratio=0;
        *soed=0;
        *absorption=0;
        *max_imbalance=0;
        *total_comm_cost=0;

        int* workload = (int*)calloc(num_processes,sizeof(int));
        int total_workload = 0;
        long int edgecut = 0;
        int hyperedges_cut = 0;
        long int total_edges = 0;
        int* vertices_in_partition = (int*)calloc(num_processes,sizeof(int)); // used for absorption
        
        for(int vid=0; vid < num_vertices; vid++) {
            if(vtx_wgt != NULL) {
                workload[partitioning[vid]] += vtx_wgt[vid];
                total_workload += vtx_wgt[vid];
            } else {
                workload[partitioning[vid]] += 1;
                total_workload += 1;
            }
        }
        for(int ii=0; ii < hyperedges->size(); ii++){
            std::set<int> connectivity;
            memset(vertices_in_partition,0,sizeof(int) * num_processes);
            for(int ff = 0; ff < hyperedges->at(ii).size(); ff++) {
                int from = hyperedges->at(ii)[ff];
                connectivity.insert(partitioning[from]);
                vertices_in_partition[partitioning[from]]++; 

                for(int tt=0; tt < hyperedges->at(ii).size(); tt++) {
                    if(tt==ff) continue;
                    int to = hyperedges->at(ii)[tt];
                    total_edges++;
                    int to_part = partitioning[to];
                    //vertices_in_partition[to_part]++;
                    if(to_part != partitioning[from]) {
                        edgecut++;
                        connectivity.insert(to_part);
                    }
                    // this cost measures mainly edge cut
                    if(comm_cost_matrix != NULL) {
                        *total_comm_cost += comm_cost_matrix[partitioning[from]][to_part];
                    }
                }
            }
            // metrics per hyperedge
            if(connectivity.size() > 1) {
                *soed += connectivity.size(); // counts as 1 external degree per partition participating
                hyperedges_cut++;
                // communication cost based on hyperedge cut
                /*for (std::set<int>::iterator sender=connectivity.begin(); sender!=connectivity.end(); ++sender) {
                    int sender_id = *sender;
                    for (std::set<int>::iterator receiver=connectivity.begin(); receiver!=connectivity.end(); ++receiver) {
                        int receiver_id = *receiver;    
                        if (sender_id != receiver_id)
                            *total_comm_cost += comm_cost_matrix[sender_id][receiver_id];
                    }
                }*/
            }
            if(hyperedges->at(ii).size() > 1) {
                for(int pp = 0; pp < num_processes; pp++) {
                    if(vertices_in_partition[pp] == 0) continue;
                    *absorption += (float)(vertices_in_partition[pp]-1) / (float)(hyperedges->at(ii).size()-1);
                } 
            }
        }

        float expected_work = total_workload / num_processes;
        for(int ii=0; ii < num_processes; ii++) {
            float imbalance = (workload[ii] / expected_work);
            if(imbalance  > *max_imbalance) {
                *max_imbalance = imbalance; 
            }
        }

        *hyperedges_cut_ratio=(float)hyperedges_cut/hyperedges->size();
        *edges_cut_ratio=(float)edgecut/total_edges;
        

        PRINTF("Quality: %i (hedgecut, %.3f total) %.3f (cut net), %i (SOED), %.1f (absorption) %.3f (max imbalance), %f (comm cost)\n",hyperedges_cut,*hyperedges_cut_ratio,*edges_cut_ratio,*soed,*absorption,*max_imbalance,*total_comm_cost);
        
        // clean up
        free(workload);
        free(vertices_in_partition);
    }

    
    void storePartitionStats(std::string experiment_name, idx_t* partitioning, int num_processes, int num_vertices, std::vector<std::vector<int> >* hyperedges, std::vector<std::vector<int> >* hedge_ptr, int* vtx_wgt,double** comm_cost_matrix) {
        float hyperedges_cut_ratio;
        float edges_cut_ratio;
        int soed;
        float absorption;
        float max_imbalance;
        double total_comm_cost;
        getPartitionStats(partitioning, num_processes, num_vertices, hyperedges, hedge_ptr, vtx_wgt,comm_cost_matrix,
                            &hyperedges_cut_ratio, &edges_cut_ratio, &soed, &absorption, &max_imbalance, &total_comm_cost);
        
        printf("Hedgecut, %.3f, %.3f (cut net), %i (SOED), %.1f (absorption) %.3f (max imbalance), %f (comm cost)\n",hyperedges_cut_ratio,edges_cut_ratio,soed,absorption,max_imbalance,total_comm_cost);
        // store stats in a file
        experiment_name = getFileName(experiment_name);
        experiment_name += "_part_stats_";
		char str_int[16];
        sprintf(str_int,"%i",num_processes);
		experiment_name += "__";
		experiment_name +=  str_int;
		bool fileexists = access(experiment_name.c_str(), F_OK) != -1;
		FILE *fp = fopen(experiment_name.c_str(), "ab+");
		if(fp == NULL) {
		    printf("Error when storing results into file\n");
		} else {
			if(!fileexists) // file does not exist, add header
			    fprintf(fp,"%s,%s,%s,%s,%s,%s\n","Hedge cut ratio","Cut net","SOED","Absorption","Max imbalance","Comm cost");
			fprintf(fp,"%f,%f,%i,%f,%f,%f\n",hyperedges_cut_ratio,edges_cut_ratio,soed,absorption,max_imbalance,total_comm_cost);
        }
        fclose(fp);
    }

    // get number of vertices in hMETIS file
    void get_hypergraph_file_header(std::string hypergraph_file, int* num_vertices, int* num_hyperedges) {
        std::ifstream istream(hypergraph_file.c_str());
        
        if(!istream) {
            printf("Error while opening hMETIS file %s\n",hypergraph_file.c_str());
        }
        std::string line;
        // process header
        std::getline(istream,line);
        std::istringstream buf(line);
        std::istream_iterator<int> beg(buf), end;
        std::vector<int> tokens(beg, end);
        *num_vertices = tokens[1];
        *num_hyperedges = tokens[0];

        istream.close();

    }
    // Load hMETIS file (hypergraph)	
    int load_hypergraph_from_file(std::string filename, std::vector<std::vector<int> >* hyperedges, std::vector<std::vector<int> >* hedge_ptr) {
        
        // get header info
        int total_vertices;;
        int total_hyperedges;
        get_hypergraph_file_header(filename,&total_vertices,&total_hyperedges);

        hyperedges->resize(total_hyperedges);
        hedge_ptr->resize(total_vertices);
        
        PRINTF("Loaded from file: Vertices: %i; hyperedges %i:\n",total_vertices,total_hyperedges);
        
        std::ifstream istream(filename.c_str());
        
        if(!istream) {
            printf("Error while opening hMETIS file %s\n",filename.c_str());
            return -1;
        }
        std::string line;
        // skip header
        std::getline(istream,line);
        // read reminder of file (one line per hyperedge)
        int counter = 0;
        while(std::getline(istream,line)) {
            std::istringstream buf(line);
            std::istream_iterator<int> beg(buf), end;
            std::vector<int> tokens(beg, end);
            for(int ii=0; ii < tokens.size(); ii++) {
                int vertex_id = tokens[ii]-1;
                hedge_ptr->at(vertex_id).push_back(counter);
                hyperedges->at(counter).push_back(vertex_id);
            }
            counter++;
            
        }
        istream.close();
        return 0;
        
    }

    // load to distributed CSR format
    int load_hypergraph_from_file_dist_CSR(std::string filename, std::vector<std::vector<int> >* hyperedges, std::vector<std::vector<int> >* hedge_ptr, int process_id, idx_t* partitioning) {
        
        // get header info
        int total_vertices;;
        int total_hyperedges;
        get_hypergraph_file_header(filename,&total_vertices,&total_hyperedges);

        hyperedges->resize(total_hyperedges);
        hedge_ptr->resize(total_vertices);
        
        PRINTF("Loaded from file: Vertices: %i; hyperedges %i:\n",total_vertices,total_hyperedges);
        
        std::ifstream istream(filename.c_str());
        
        if(!istream) {
            printf("Error while opening hMETIS file %s\n",filename.c_str());
            return -1;
        }
        
        std::string line;
        // skip header
        std::getline(istream,line);
        // read reminder of file (one line per hyperedge)
        int counter = 0;
        int nonlocal = 0;
        while(std::getline(istream,line)) {
            std::istringstream buf(line);
            std::istream_iterator<int> beg(buf), end;
            std::vector<int> tokens(beg, end);
            // check if any vertex is local
            bool local = false;
            for(int ii=0; ii < tokens.size(); ii++) {
                int vertex_id = tokens[ii]-1;
                if(partitioning[vertex_id] == process_id) {
                    local = true;
                    break;
                }
            }
            // add to memory only if local
            if(local) {
                for(int ii=0; ii < tokens.size(); ii++) {
                    int vertex_id = tokens[ii]-1;
                    if(partitioning[vertex_id] == process_id) {
                        hedge_ptr->at(vertex_id).push_back(counter);
                    }
                    hyperedges->at(counter).push_back(vertex_id);
                }
            } else nonlocal++;
            counter++;            
        }
        istream.close();
        PRINTF("%i: non local hyperedges: %i\n",process_id,nonlocal);
        return 0;
        
    }

    void getPartitionStatsFromFile(idx_t* partitioning, int num_processes, int num_vertices, std::string hgraph_filename, int* vtx_wgt,double** comm_cost_matrix, // input
                            float* hyperedges_cut_ratio, float* edges_cut_ratio, int* soed, float* absorption, float* max_imbalance, double* total_comm_cost, bool save_theoretical_comm) { // output

        // take partitioning statistics directly from reading a graph file
        // loading the entire graph on a process should be avoided for scalability
        *hyperedges_cut_ratio=0;
        *edges_cut_ratio=0;
        *soed=0;
        *absorption=0;
        *max_imbalance=0;
        *total_comm_cost=0;

        int* workload = (int*)calloc(num_processes,sizeof(int));
        int total_workload = 0;
        long int edgecut = 0;
        int hyperedges_cut = 0;
        long int total_edges = 0;
        int* vertices_in_partition = (int*)calloc(num_processes,sizeof(int)); // used for absorption
        
        // get header info
        int total_vertices;;
        int total_hyperedges;
        get_hypergraph_file_header(hgraph_filename,&total_vertices,&total_hyperedges);

        std::vector<std::vector<int> > hyperedges(total_hyperedges);
        std::vector<std::vector<int> > hedge_ptr(total_vertices);
        
        std::ifstream istream(hgraph_filename.c_str());

#ifdef SAVE_COMM_COST
        int** theoretical_comm = (int**)malloc(num_processes * sizeof(int*));
        for(int ii=0; ii < num_processes; ii++) {
            theoretical_comm[ii] = (int*)calloc(num_processes,sizeof(int));
        }
#endif
        
        if(!istream) {
            printf("Error while opening hMETIS file %s\n",hgraph_filename.c_str());
            return;
        }
        std::string line;

        // compute workload for all vertices and partitions
        for(int vid=0; vid < total_vertices; vid++) {
            if(vtx_wgt != NULL) {
                workload[partitioning[vid]] += vtx_wgt[vid];
                total_workload += vtx_wgt[vid];
            } else {
                workload[partitioning[vid]] += 1;
                total_workload += 1;
            }
        }

        // skip header
        std::getline(istream,line);
        // read reminder of file (one line per hyperedge)
        while(std::getline(istream,line)) {
            // each line corresponds to a hyperedge
            std::istringstream buf(line);
            std::istream_iterator<int> beg(buf), end;
            std::vector<int> tokens(beg, end);

            std::set<int> connectivity;
            memset(vertices_in_partition,0,sizeof(int) * num_processes);
            for(int ff=0; ff < tokens.size(); ff++) {
                // each vertex in the hyperedge
                int from = tokens[ff]-1;
                connectivity.insert(partitioning[from]);
                vertices_in_partition[partitioning[from]]++;
                for(int tt=0; tt < tokens.size(); tt++) {
                    if(tt==ff) continue;
                    int to = tokens[tt]-1;
                    total_edges++;
                    int to_part = partitioning[to];
                    //vertices_in_partition[to_part]++;
                    if(to_part != partitioning[from]) {
                        edgecut++;
                        connectivity.insert(to_part);
                    }
                    // this cost measures mainly edge cut
                    if(comm_cost_matrix != NULL) {
                        *total_comm_cost += comm_cost_matrix[partitioning[from]][to_part];
#ifdef SAVE_COMM_COST
                        if(save_theoretical_comm) theoretical_comm[partitioning[from]][to_part] += 1;
#endif
                    }
                }
            }
            // metrics per hyperedge
            if(connectivity.size() > 1) {
                *soed += connectivity.size(); // counts as 1 external degree per partition participating
                hyperedges_cut++;
                // communication cost based on hyperedge cut
                /*for (std::set<int>::iterator sender=connectivity.begin(); sender!=connectivity.end(); ++sender) {
                    int sender_id = *sender;
                    for (std::set<int>::iterator receiver=connectivity.begin(); receiver!=connectivity.end(); ++receiver) {
                        int receiver_id = *receiver;    
                        if (sender_id != receiver_id)
                            *total_comm_cost += comm_cost_matrix[sender_id][receiver_id];
                    }
                }*/
            }
            if(tokens.size() > 1) {
                for(int pp = 0; pp < num_processes; pp++) {
                    if(vertices_in_partition[pp] == 0) continue;
                    *absorption += (float)(vertices_in_partition[pp]-1) / (float)(tokens.size()-1);
                } 
            }
            
        }
        istream.close();

        float expected_work = total_workload / num_processes;
        for(int ii=0; ii < num_processes; ii++) {
            float imbalance = (workload[ii] / expected_work);
            if(imbalance  > *max_imbalance) {
                *max_imbalance = imbalance; 
            }
        }

        *hyperedges_cut_ratio=(float)hyperedges_cut/total_hyperedges;
        *edges_cut_ratio=(float)edgecut/total_edges;

#ifdef SAVE_COMM_COST
        if(save_theoretical_comm) {
            std::string comm_cost_file = getFileName(hgraph_filename);
            comm_cost_file += "_theoretical_comm";
            char str_int[16];
            sprintf(str_int,"%i",num_processes);
            comm_cost_file += "__";
            comm_cost_file +=  str_int;
            // remove  file if exists
            FILE *fp = fopen(comm_cost_file.c_str(), "w");
            if(fp == NULL) {
                printf("Error when storing partitioning history into file\n");
            } else {
                for(int ii=0; ii < num_processes; ii++) {
                    for(int jj=0; jj < num_processes; jj++) {
                        fprintf(fp,"%i ",theoretical_comm[ii][jj]);
                    }
                    fprintf(fp,"\n");
                }
            }
            fclose(fp);
        }
        
            
        for(int ii=0; ii < num_processes; ii++) {
            free(theoretical_comm[ii]);
        }
        free(theoretical_comm);
#endif
        

        PRINTF("Quality: %i (hedgecut, %.3f total) %.3f (cut net), %i (SOED), %.1f (absorption) %.3f (max imbalance), %f (comm cost)\n",hyperedges_cut,*hyperedges_cut_ratio,*edges_cut_ratio,*soed,*absorption,*max_imbalance,*total_comm_cost);
        
        // clean up
        free(workload);
        free(vertices_in_partition);
    }

    void get_comm_cost_matrix_from_bandwidth(char* comm_bandwidth_filename, double** comm_cost_matrix, int partitions) {
        std::ifstream input_stream(comm_bandwidth_filename);
        if(input_stream) {
            std::string line;
            int current_process = 0;
            while(std::getline(input_stream,line)) {
                std::istringstream buf(line);
                std::istream_iterator<double> beg(buf), end;
                std::vector<double> tokens(beg, end);
                for(int ii=0; ii < tokens.size(); ii++) {
                    comm_cost_matrix[current_process][ii] = tokens[ii];
                }
                current_process++;
            }
            // standardise bandwidth matrix and transform to cost matrix
            
            // transform bandwidth to cost (relative to min bandwidth)
            float min_bandwidth = std::numeric_limits<float>::max();
            float max_bandwidth = 0;
            for(int ii = 0; ii < partitions; ii++) {
                for(int jj = 0; jj < partitions; jj++) {
                    if(ii == jj) continue;
                    if(comm_cost_matrix[ii][jj] < min_bandwidth)
                        min_bandwidth = comm_cost_matrix[ii][jj];
                    if(comm_cost_matrix[ii][jj] > max_bandwidth)
                        max_bandwidth = comm_cost_matrix[ii][jj];
                }
            }
            for(int ii = 0; ii < partitions;ii++) {
                std::transform(comm_cost_matrix[ii],comm_cost_matrix[ii]+partitions,comm_cost_matrix[ii],
                            [min_bandwidth,max_bandwidth] (double value) {  
                                float ratio = max_bandwidth / min_bandwidth;
                                return value <= std::numeric_limits<float>::epsilon() ? 0 : (1-(value-min_bandwidth)/(max_bandwidth-min_bandwidth)) * ratio + 1;
                                //return value <= std::numeric_limits<float>::epsilon() ? 0 : 2 - ( (value-min_bandwidth)/(max_bandwidth-min_bandwidth) );
                            }   
                );
            }
        } else {
            // file not found, default values
            for(int ii = 0; ii < partitions;ii++) {
                for(int jj=0; jj < partitions; jj++) {
                    comm_cost_matrix[ii][jj] = ii == jj ? 0 : 1;
                }
            }
        }
    }

    int SequentialStreamingPartitioning(idx_t* partitioning, int num_processes, double** comm_cost_matrix, std::string hypergraph_filename, int* vtx_wgt, int iterations, float imbalance_tolerance, bool reset_partitioning) {
        
        // get meta info (num vertices and hyperedges)
        int num_vertices, num_hyperedges;
        get_hypergraph_file_header(hypergraph_filename, &num_vertices, &num_hyperedges);
        
        //PARAMETERS: From Battaglino 2015 //
        // https://github.com/cjbattagl/GraSP
        // g and a determine load balance importance in cost function; was 1.5
        double g = 1.5; // same as FENNEL Tsourakakis 2012
        // battaglino's initial alpha, was sqrt(2) * num_hyperedges / pow(num_vertices,g);
        double a = sqrt(num_processes) * num_hyperedges / pow(num_vertices,g); // same as FENNEL Tsourakakis 2012
        // ta is the update rate of parameter a; was 1.7
        double ta_start = 1.7; // used when imbalance is far from imbalance_tolerance
        double ta_refine = 1.3; // used when imbalance is close to imbalance_tolerance
        // minimum number of iterations run (not checking imbalance threshold)
        // removed whilst we are using hyperPraw as refinement algorithm
        //      hence, if balanced is kept after first iteration, that's good enough
        int frozen_iters = 0;
        ///////////////
        
        // algorithm from GraSP (Battaglino 2016)
        // 1 - Distributed vertices over partitions (partition = vertex_id % num_partitions)
        // needs to load num_vertices from file
        if(reset_partitioning) {
            for (int vid=0; vid < num_vertices; vid++) {
                partitioning[vid] = MASTER_NODE;
            }
        }
        
        // 2 - Divide the graph in a distributed CSR format (like ParMETIS)
        //  compressed vertex or compressed hedge format? --> see zoltan
        //  for each local vertex, store the list of vertices adjacent to it (belonging to same hedges)
        std::vector<std::vector<int> > hyperedges;
        std::vector<std::vector<int> > hedge_ptr;
        load_hypergraph_from_file_dist_CSR(hypergraph_filename, &hyperedges, &hedge_ptr, MASTER_NODE, partitioning);
        

        // each process must read from file only the info relevant to its data
        // 3 - Initiate N number of iterations on each process:
        //      a - one vertex at a time, assign to best partition (based on eval function)
        //      b - update tempering parameters
        //      c - share with all new partition assignments
        
#ifdef SAVE_HISTORY
        std::string history_file = getFileName(hypergraph_filename);
        history_file += "_partition_history_parallel_";
        char str_int[16];
        sprintf(str_int,"%i",num_processes);
        history_file += "__";
        history_file +=  str_int;
        // remove history file if exists
        FILE *fp = fopen(history_file.c_str(), "w");
        if(fp == NULL) {
            printf("Error when storing partitioning history into file\n");
        } else {
            fprintf(fp,"%s\n","Imbalance");
        }
        fclose(fp);
        
            
#endif
        long int* part_load = (long int*)calloc(num_processes,sizeof(long int));
        double* comm_cost_per_partition = (double*)malloc(num_processes*sizeof(double));
        int* current_neighbours_in_partition = (int*)malloc(num_processes*sizeof(int));
        // overfit variables
        bool check_overfit = false;
        idx_t* last_partitioning = NULL;
        float last_cut_metric;
        bool rollback = false;
        float last_imbalance = num_processes;
        //double timing = 0;
        //double ttt;
        for(int iter=0; iter < iterations; iter++) {
            //timing = 0;
            memset(part_load,0,num_processes * sizeof(long int));
            double total_workload = 0;
            for(int ii=0; ii < num_vertices; ii++) {
                part_load[partitioning[ii]] += vtx_wgt[ii]; // workload for vertex
                total_workload += vtx_wgt[ii];
            }
            double expected_workload = total_workload / num_processes;
            
            // go through own vertex list and reassign
            for(int vid=0; vid < num_vertices; vid++) {
                memset(current_neighbours_in_partition,0,num_processes * sizeof(int));
                memset(comm_cost_per_partition,0,num_processes * sizeof(double));

                int total_neighbours = 1;
                double max_comm_cost = 0;
                // reevaluate objective function per partition
                // |P^t_i union N(v)| = number of vertices in partition i that are neighbours of vertex v 
                // where are neighbours located
                // new communication cost incurred
                
                // does not double count vertices that are present in multiple hyperedges
                // communication cost should be based on hedge cut?
                for(int he = 0; he < hedge_ptr[vid].size(); he++) {
                    //bool* visited = (bool*)calloc(num_vertices,sizeof(bool));
                    int he_id = hedge_ptr[vid][he];
                    for(int vt = 0; vt < hyperedges[he_id].size(); vt++) {
                        int dest_vertex = hyperedges[he_id][vt];
                        if(dest_vertex == vid) continue;
                        total_neighbours++;
                        int dest_part = partitioning[dest_vertex];
                        //if(!visited[dest_vertex]) 
                            current_neighbours_in_partition[dest_part] += 1;
                        // recalculate comm cost for all possible partition assignments of vid
                        //  commCost(v,Pi) = forall edge in edges(Pi) cost += w(e) * c(Pi,Pj) where i != j
                        for(int fp=0; fp < num_processes; fp++) {
                            comm_cost_per_partition[fp] += 1 * comm_cost_matrix[fp][dest_part];
                            if(comm_cost_per_partition[fp] > max_comm_cost)
                                max_comm_cost = comm_cost_per_partition[fp];
                        }
                        //visited[dest_vertex] = true;
                    }
                    //free(visited);
                }
                

                if(max_comm_cost < std::numeric_limits<double>::epsilon()) max_comm_cost = 1;
                
                // allocate vertex (for local heuristically, for non local speculatively)
                double max_value = std::numeric_limits<double>::lowest();
                int best_partition = partitioning[vid];
                //std::vector<int> best_parts;
                for(int pp=0; pp < num_processes; pp++) {
                    // total cost of communication (edgecuts * number of participating partitions)
                    long int total_comm_cost = 0;
                    for(int jj=0; jj < num_processes; jj++) {
                        if(pp != jj)
                            total_comm_cost += current_neighbours_in_partition[jj] > 0 ? 1 : 0;
                    }
                    
                    // TODO: Is total_comm_cost making things worse when using bandwidth info??
                    double current_value = current_neighbours_in_partition[pp]/(double)total_neighbours -comm_cost_per_partition[pp] - a * (part_load[pp]/expected_workload);
                    //double current_value = current_neighbours_in_partition[pp]/(double)total_neighbours -(double)total_comm_cost / (double)num_processes * comm_cost_per_partition[pp] / max_comm_cost - a * (part_load[pp]/expected_workload);
                    // double current_value  = current_neighbours_in_partition[pp] -(double)total_comm_cost * comm_cost_per_partition[pp] - a * g/2 * pow(part_load[pp],g-1);
                    
                    // lesson learned, global hygergraph partitioners use connectivity metric as cost function
                    // try lotfifar 2015
                    //  cost function is connectivity degree * weight for each he
                    //  balance is bounded on both sides, +and- imbalance tolerance
                    // can we try coarsening the hypergraph to speed up partitioning?
                    //  can use zoltan's / patoh inner product matching / heavy connectivity matching
                    //  other similarity metrics such as Jaccard Index or Cosine measure (lotfifar 2015)
                    // we are not measuring migration costs (cataluyrek 2007 models it well)
                    if(current_value > max_value) {
                        max_value = current_value;
                        best_partition = pp;
                        //best_parts.clear();
                        //best_parts.push_back(pp);
                    } /*else if(fabs(current_value-max_value) <= std::numeric_limits<double>::epsilon()*2) {
                        best_parts.push_back(pp);
                    }*/
                }
                
                //best_partition = best_parts[(int)(best_parts.size() * (double)rand() / (double)RAND_MAX)];
                
                
                // update intermediate workload and assignment values
                part_load[best_partition] += vtx_wgt[vid];
                part_load[partitioning[vid]] -= vtx_wgt[vid];
                 
                // update partitioning assignment
                partitioning[vid] = best_partition;
                
                
            }
            
            // check if desired imbalance has been reached
            float imbalance = calculateImbalance(partitioning,num_processes,num_vertices,vtx_wgt);
            PRINTF("%i: %f (%f | %f)\n",iter,imbalance,a,ta_start);

#ifdef SAVE_HISTORY
            FILE *fp = fopen(history_file.c_str(), "ab+");
            if(fp == NULL) {
                printf("Error when storing partitioning history into file\n");
            } else {
                fprintf(fp,"%.3f\n",imbalance);
            }
            fclose(fp);
            
            
#endif
            // stop the process in the following conditions
            //  1. imbalance tolerance has been reached
            //      record current cut metric and partitioning and do one more iteration
            //      if imbalance is still ok 
            //          metric has not been improved, take recorded partitioning and stop
            //          metric has been improved, store partitioning and do one more iteration
            // ALL PROCESS MUST STOP to check if 0 has broken out of the loop
            if(frozen_iters <= iter) {
                // problem! once check_overfit is set, must check if partitioning result was better before
                // reproduce issue with venkat01 at 24 processes with update part load at 2500
                if (imbalance < imbalance_tolerance) {
                    // get cut metric
                    float hyperedges_cut_ratio;
                    float edges_cut_ratio;
                    int soed;
                    float absorption;
                    float max_imbalance;
                    double total_comm_cost;
                    PRAW::getPartitionStatsFromFile(partitioning, num_processes, num_vertices, hypergraph_filename, NULL,comm_cost_matrix,
                                &hyperedges_cut_ratio, &edges_cut_ratio, &soed, &absorption, &max_imbalance, &total_comm_cost,false);
                    double cut_metric = hyperedges_cut_ratio + edges_cut_ratio;//hyperedges_cut_ratio;

                    if(!check_overfit) {
                        // record partitioning and cut metric
                        last_cut_metric = cut_metric;
                        if(last_partitioning == NULL) {
                            last_partitioning = (idx_t*)malloc(num_vertices*sizeof(idx_t));
                        }
                        memcpy(last_partitioning,partitioning,num_vertices * sizeof(idx_t));
                        check_overfit = true;
                    } else {
                        // check if cut metric has improved
                        if(cut_metric >= last_cut_metric) {
                            // send signal to stop
                            rollback = true;
                            break;
                        } else {
                            last_cut_metric = cut_metric;
                            memcpy(last_partitioning,partitioning,num_vertices * sizeof(idx_t));
                        }
                    }
                    
                } else {
                    /*if(check_overfit) {
                        rollback = true;
                        break;
                    }*/
                    check_overfit = false;
                }  
            }
            //if(frozen_iters <= iter && imbalance < imbalance_tolerance) break;

            // update parameters
            if(imbalance > imbalance_tolerance) {
                if(imbalance > 1.2f * imbalance_tolerance) {
                    a *= ta_start;
                } else {
                    a *= ta_refine;
                }
                
            }
            last_imbalance = imbalance;
        }

        if(rollback) {
            // share last partitioning with all
            memcpy(partitioning,last_partitioning,num_vertices * sizeof(idx_t));
            free(last_partitioning);
            
        }

        // clean up
        free(part_load);
        free(comm_cost_per_partition);
        free(current_neighbours_in_partition);

        return 0;
    }

    int ParallelIndependentRestreamingPartitioning(idx_t* partitioning, double** comm_cost_matrix, std::string hypergraph_filename, int* vtx_wgt, int iterations, float imbalance_tolerance, bool reset_partitioning) {
        
        int process_id;
        MPI_Comm_rank(MPI_COMM_WORLD,&process_id);
        int num_processes;
        MPI_Comm_size(MPI_COMM_WORLD,&num_processes);

        // get meta info (num vertices and hyperedges)
        int num_vertices, num_hyperedges;
        get_hypergraph_file_header(hypergraph_filename, &num_vertices, &num_hyperedges);
        
        //PARAMETERS: From Battaglino 2015 //
        // https://github.com/cjbattagl/GraSP
        // g and a determine load balance importance in cost function; was 1.5
        double g = 1.5; // same as FENNEL Tsourakakis 2012
        // battaglino's initial alpha, was sqrt(2) * num_hyperedges / pow(num_vertices,g);
        double a = sqrt(num_processes) * num_hyperedges / pow(num_vertices,g); // same as FENNEL Tsourakakis 2012
        // ta is the update rate of parameter a; was 1.7
        double ta_start = 1.7; // used when imbalance is far from imbalance_tolerance
        double ta_refine = 1.3; // used when imbalance is close to imbalance_tolerance
        // after how many vertices checked in the stream the partitio load is sync across processes
        int part_load_update_after_vertices = 4000;//sqrt(num_processes) * 300; // in the paper it is 4096
        // minimum number of iterations run (not checking imbalance threshold)
        // removed whilst we are using hyperPraw as refinement algorithm
        //      hence, if balanced is kept after first iteration, that's good enough
        int frozen_iters = 0;
        ///////////////
        
        // algorithm from GraSP (Battaglino 2016)
        // 1 - Distributed vertices over partitions (partition = vertex_id % num_partitions)
        // needs to load num_vertices from file
        if (reset_partitioning) {
            for (int vid=0; vid < num_vertices; vid++) {
                partitioning[vid] = vid % num_processes;
            }
            //frozen_iters = ceil(0.1f * iterations);
        }
        
        // 2 - Divide the graph in a distributed CSR format (like ParMETIS)
        //  compressed vertex or compressed hedge format? --> see zoltan
        //  for each local vertex, store the list of vertices adjacent to it (belonging to same hedges)
        std::vector<std::vector<int> > hyperedges;
        std::vector<std::vector<int> > hedge_ptr;
        load_hypergraph_from_file_dist_CSR(hypergraph_filename, &hyperedges, &hedge_ptr, process_id, partitioning);
        

        // each process must read from file only the info relevant to its data
        // 3 - Initiate N number of iterations on each process:
        //      a - one vertex at a time, assign to best partition (based on eval function)
        //      b - update tempering parameters
        //      c - share with all new partition assignments
        
#ifdef SAVE_HISTORY
        std::string history_file = getFileName(hypergraph_filename);
        history_file += "_partition_history_parallel_";
        char str_int[16];
        sprintf(str_int,"%i",num_processes);
        history_file += "__";
        history_file +=  str_int;
        // remove history file if exists
        if(process_id == 0) {
            FILE *fp = fopen(history_file.c_str(), "w");
            if(fp == NULL) {
                printf("Error when storing partitioning history into file\n");
            } else {
                fprintf(fp,"%s\n","Imbalance");
            }
            fclose(fp);
        }
            
#endif
        long int* part_load = (long int*)calloc(num_processes,sizeof(long int));
        idx_t* local_stream_partitioning = (idx_t*)malloc(num_vertices*sizeof(idx_t));
        double* comm_cost_per_partition = (double*)malloc(num_processes*sizeof(double));
        int* current_neighbours_in_partition = (int*)malloc(num_processes*sizeof(int));
        long int* part_load_update = (long int*)calloc(num_processes,sizeof(long int));
        long int* part_load_speculative_update = (long int*)calloc(num_processes,sizeof(long int));
        // overfit variables
        bool check_overfit = false;
        idx_t* last_partitioning = NULL;
        float last_cut_metric;
        bool rollback = false;
        float last_imbalance = num_processes;
        //double timing = 0;
        //double ttt;
        for(int iter=0; iter < iterations; iter++) {
            //timing = 0;
            memset(local_stream_partitioning,0,num_vertices * sizeof(idx_t));
            memset(part_load,0,num_processes * sizeof(long int));
            memset(part_load_update,0,num_processes * sizeof(long int));
            memset(part_load_speculative_update,0,num_processes * sizeof(long int));
            double total_workload = 0;
            for(int ii=0; ii < num_vertices; ii++) {
                part_load[partitioning[ii]] += vtx_wgt[ii]; // workload for vertex
                total_workload += vtx_wgt[ii];
            }
            double expected_workload = total_workload / num_processes;
            
            // go through own vertex list and reassign
            for(int vid=0; vid < num_vertices; vid++) {
                // share updated partition loads after constant number of iterations
                if(vid % part_load_update_after_vertices == 0) {
                    for(int pl=0; pl < num_processes; pl++) {
                        // we need to discount the local part update from the total since the total includes it
                        // heuristic local update
                        part_load[pl] -= part_load_update[pl];
                        // speculative local update
                        part_load[pl] -= part_load_speculative_update[pl];
                    }
                    MPI_Allreduce(MPI_IN_PLACE,part_load_update,num_processes,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
                    for(int pl=0; pl < num_processes; pl++) {
                        part_load[pl] += part_load_update[pl];
                    }
                    memset(part_load_update,0,num_processes * sizeof(long int));
                    memset(part_load_speculative_update,0,num_processes * sizeof(long int));
                }
                memset(current_neighbours_in_partition,0,num_processes * sizeof(int));
                memset(comm_cost_per_partition,0,num_processes * sizeof(double));

                bool isLocal = hedge_ptr[vid].size() > 0;

                int total_neighbours = 1;
                double max_comm_cost = 0;
                // if local vertex, calculate full heuristic (cost of communication...)
                // if non local vertex, speculatively place it based on current partitioning load balance
                // this alleviates the problems of parallel streams maintaining workload balance when 
                // alpha parameter is high and partition load update is sparse
                if(isLocal) {
                    // reevaluate objective function per partition
                    // |P^t_i union N(v)| = number of vertices in partition i that are neighbours of vertex v 
                    // where are neighbours located
                    // new communication cost incurred
                    //ttt = MPI_Wtime();
                    
                    // does not double count vertices that are present in multiple hyperedges
                    // communication cost should be based on hedge cut?
                    //bool* visited = (bool*)calloc(num_vertices,sizeof(bool));
                    for(int he = 0; he < hedge_ptr[vid].size(); he++) {
                        int he_id = hedge_ptr[vid][he];
                        for(int vt = 0; vt < hyperedges[he_id].size(); vt++) {
                            int dest_vertex = hyperedges[he_id][vt];
                            if(dest_vertex == vid) continue;
                            total_neighbours++;
                            int dest_part = partitioning[dest_vertex];
                            //if(!visited[dest_vertex]) 
                                current_neighbours_in_partition[dest_part] += 1;
                            // recalculate comm cost for all possible partition assignments of vid
                            //  commCost(v,Pi) = forall edge in edges(Pi) cost += w(e) * c(Pi,Pj) where i != j
                            for(int fp=0; fp < num_processes; fp++) {
                                comm_cost_per_partition[fp] += 1 * comm_cost_matrix[fp][dest_part];
                                if(comm_cost_per_partition[fp] > max_comm_cost)
                                    max_comm_cost = comm_cost_per_partition[fp];
                            }
                            //visited[dest_vertex] = true;
                        }
                    }
                    //free(visited);
                    //timing += MPI_Wtime() - ttt;
                } 

                if(max_comm_cost < std::numeric_limits<double>::epsilon()) max_comm_cost = 1;
                
                // allocate vertex (for local heuristically, for non local speculatively)
                double max_value = std::numeric_limits<double>::lowest();
                int best_partition = partitioning[vid];
                //std::vector<int> best_parts;
                for(int pp=0; pp < num_processes; pp++) {
                    // total cost of communication (edgecuts * number of participating partitions)
                    long int total_comm_cost = 0;
                    if(isLocal) {
                        for(int jj=0; jj < num_processes; jj++) {
                            if(pp != jj)
                                total_comm_cost += current_neighbours_in_partition[jj] > 0 ? 1 : 0;
                        }
                    } 

                    float random_factor = 0;
                    if (!isLocal) {
                        random_factor =  0.2f * ((float)rand() / (float)RAND_MAX - 0.5f);
                    }

                    double current_value = random_factor + current_neighbours_in_partition[pp]/(double)total_neighbours -(double)total_comm_cost / (double)num_processes * comm_cost_per_partition[pp] / max_comm_cost - a * (part_load[pp]/expected_workload);
                    //double current_value =  (float)current_neighbours_in_partition[pp]/(float)total_neighbours - (double)total_comm_cost/(double)num_processes * comm_cost_per_partition[pp] - a * (part_load[pp]/expected_workload);
                    // double current_value  = current_neighbours_in_partition[pp] -(double)total_comm_cost * comm_cost_per_partition[pp] - a * g/2 * pow(part_load[pp],g-1);
                    
                    // lesson learned, global hygergraph partitioners use connectivity metric as cost function
                    // try lotfifar 2015
                    //  cost function is connectivity degree * weight for each he
                    //  balance is bounded on both sides, +and- imbalance tolerance
                    // can we try coarsening the hypergraph to speed up partitioning?
                    //  can use zoltan's / patoh inner product matching / heavy connectivity matching
                    //  other similarity metrics such as Jaccard Index or Cosine measure (lotfifar 2015)
                    // we are not measuring migration costs (cataluyrek 2007 models it well)
                    if(current_value > max_value) {
                        max_value = current_value;
                        best_partition = pp;
                        //best_parts.clear();
                        //best_parts.push_back(pp);
                    } /*else if(fabs(current_value-max_value) <= std::numeric_limits<double>::epsilon()*2) {
                        best_parts.push_back(pp);
                    }*/
                }
                
                //best_partition = best_parts[(int)(best_parts.size() * (double)rand() / (double)RAND_MAX)];
                
                
                // update intermediate workload and assignment values
                part_load[best_partition] += vtx_wgt[vid];
                part_load[partitioning[vid]] -= vtx_wgt[vid];
                 
                if(isLocal) {
                    // update local changes counter
                    part_load_update[partitioning[vid]] -= vtx_wgt[vid];
                    part_load_update[best_partition] += vtx_wgt[vid];
                    // update partitioning assignment
                    partitioning[vid] = best_partition;
                    local_stream_partitioning[vid] = best_partition;
                } else {
                    // keep a record of speculative load update (does not need to be propagated later)
                    part_load_speculative_update[partitioning[vid]] -= vtx_wgt[vid];
                    part_load_speculative_update[best_partition] += vtx_wgt[vid];                
                }
                
            }
            //printf("%f\n",timing);
            
            // share new partitioning with other streams
            MPI_Allreduce(local_stream_partitioning,partitioning,num_vertices,MPI_LONG,MPI_MAX,MPI_COMM_WORLD);

            // check if desired imbalance has been reached
            float imbalance = calculateImbalance(partitioning,num_processes,num_vertices,vtx_wgt);
            PRINTF("%i: %f (%f | %f)\n",iter,imbalance,a,ta_start);

#ifdef SAVE_HISTORY
            if(process_id == MASTER_NODE) {
                FILE *fp = fopen(history_file.c_str(), "ab+");
                if(fp == NULL) {
                    printf("Error when storing partitioning history into file\n");
                } else {
                    fprintf(fp,"%.3f\n",imbalance);
                }
                fclose(fp);
            }
            
#endif
            // stop the process in the following conditions
            //  1. imbalance tolerance has been reached
            //      record current cut metric and partitioning and do one more iteration
            //      if imbalance is still ok 
            //          metric has not been improved, take recorded partitioning and stop
            //          metric has been improved, store partitioning and do one more iteration
            // ALL PROCESS MUST STOP to check if 0 has broken out of the loop
            if(frozen_iters <= iter) {
                // problem! once check_overfit is set, must check if partitioning result was better before
                // reproduce issue with venkat01 at 24 processes with update part load at 2500
                if (imbalance < imbalance_tolerance) {
                    if(process_id == MASTER_NODE) {
                        // get cut metric
                        float hyperedges_cut_ratio;
                        float edges_cut_ratio;
                        int soed;
                        float absorption;
                        float max_imbalance;
                        double total_comm_cost;
                        PRAW::getPartitionStatsFromFile(partitioning, num_processes, num_vertices, hypergraph_filename, NULL,comm_cost_matrix,
                                    &hyperedges_cut_ratio, &edges_cut_ratio, &soed, &absorption, &max_imbalance, &total_comm_cost,false);
                        double cut_metric = hyperedges_cut_ratio + edges_cut_ratio;//hyperedges_cut_ratio;

                        if(!check_overfit) {
                            // record partitioning and cut metric
                            last_cut_metric = cut_metric;
                            if(last_partitioning == NULL) {
                                last_partitioning = (idx_t*)malloc(num_vertices*sizeof(idx_t));
                            }
                            memcpy(last_partitioning,partitioning,num_vertices * sizeof(idx_t));
                            check_overfit = true;
                        } else {
                            // check if cut metric has improved
                            if(cut_metric >= last_cut_metric) {
                                // send signal to stop
                                int message = 0;
                                for(int dest=0; dest < num_processes; dest++) {
                                    if(dest == MASTER_NODE) continue;
                                    MPI_Send(&message,1,MPI_INT,dest,0,MPI_COMM_WORLD);
                                }
                                rollback = true;
                                break;
                            } else {
                                last_cut_metric = cut_metric;
                                memcpy(last_partitioning,partitioning,num_vertices * sizeof(idx_t));
                            }
                        }
                        // send signal to continue
                        int message = 1;
                        for(int dest=0; dest < num_processes; dest++) {
                            if(dest == MASTER_NODE) continue;
                            MPI_Send(&message,1,MPI_INT,dest,0,MPI_COMM_WORLD);
                        }
                    } else {
                        // all other processes must wait for the signal
                        int message;
                        MPI_Recv(&message,1,MPI_INT,MASTER_NODE,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                        if(message == 0) {
                            rollback = true;
                            break;
                        }
                        check_overfit = true;
                    }
                } else {
                    /*if(check_overfit) {
                        rollback = true;
                        break;
                    }*/
                    check_overfit = false;
                }  
            }
            //if(frozen_iters <= iter && imbalance < imbalance_tolerance) break;

            // update parameters
            if(imbalance > imbalance_tolerance) {
                if(imbalance > 1.2f * imbalance_tolerance) {
                    a *= ta_start;
                } else {
                    a *= ta_refine;
                }
                
            }
            last_imbalance = imbalance;
        }

        if(rollback) {
            if(process_id == MASTER_NODE) {
                // share last partitioning with all
                memcpy(partitioning,last_partitioning,num_vertices * sizeof(idx_t));
                free(last_partitioning);
                for(int dest=0; dest < num_processes; dest++) {
                    if(dest == MASTER_NODE) continue;
                    MPI_Send(partitioning,num_vertices,MPI_LONG,dest,0,MPI_COMM_WORLD);
                }
            } else {
                // update partitioning from 0
                MPI_Recv(partitioning,num_vertices,MPI_LONG,MASTER_NODE,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }

        // clean up
        free(local_stream_partitioning);
        free(part_load);
        free(comm_cost_per_partition);
        free(current_neighbours_in_partition);
        free(part_load_update);
        free(part_load_speculative_update);

        return 0;
    }
    

    
}

#endif