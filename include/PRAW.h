// Parallel Restreaming Architecture aWare partitioning algorithm

#ifndef PRAW__H
#define PRAW__H

#define SAVE_HISTORY        // stores partitioning iteration history in file

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
            if(imbalance > *max_imbalance)
                *max_imbalance = imbalance;
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
            }
            counter++;            
        }
        istream.close();
        return 0;
        
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
            /*double* totals = (double*)calloc(partitions,sizeof(double));
            for(int ii = 0; ii < partitions;ii++) {
                totals[ii] = std::accumulate(comm_cost_matrix[ii],comm_cost_matrix[ii]+partitions,0.0);
            }
            double total_bandwidth = std::accumulate(totals,totals+partitions,0.0);
            for(int ii = 0; ii < partitions;ii++) {
                std::transform(comm_cost_matrix[ii],comm_cost_matrix[ii]+partitions,comm_cost_matrix[ii],[total_bandwidth] (double value) { return value == 0 ? 0 : 1-(value / total_bandwidth); });
            }
            free(totals);*/
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
                std::transform(comm_cost_matrix[ii],comm_cost_matrix[ii]+partitions,comm_cost_matrix[ii],[min_bandwidth,max_bandwidth] (double value) { return value == 0 ? 0 : 2 - ( (value-min_bandwidth)/(max_bandwidth-min_bandwidth) ); });
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

    int SequentialStreamingPartitioning(idx_t* partitioning, int num_processes, double** comm_cost_matrix, int num_vertices, std::vector<std::vector<int> >* hyperedges, std::vector<std::vector<int> >* hedge_ptr, int* vtx_wgt, int max_iterations, float imbalance_tolerance) {
        
        //PARAMETERS: From Battaglino 2015 //
        float g = 1.5;
        float a = sqrt(2) * hyperedges->size() / pow(num_vertices,g);
        float ta = 1.7;
        ///////////////
        
        // Do partitioning only on master node
        double expected_workload = 0;
        int* part_load = (int*)calloc(num_processes,sizeof(int));
        for(int ii=0; ii < num_vertices; ii++) {
            expected_workload += vtx_wgt[ii];
            part_load[partitioning[ii]] += vtx_wgt[ii]; // workload for vertex
        }
        expected_workload /= num_processes;
        double* comm_cost_per_partition = (double*)malloc(num_processes*sizeof(double));
        int* current_neighbours_in_partition = (int*)malloc(num_processes*sizeof(int));
        int* current_neighbours_elsewhere = (int*)malloc(num_processes*sizeof(int));
        for(int iter=0; iter < max_iterations; iter++) {
            // go through own vertex list and reassign
            for(int vid=0; vid < num_vertices; vid++) {
                // reevaluate objective function per partition
                // |P^t_i union N(v)| = number of vertices in partition i that are neighbours of vertex v 
                // where are neighbours located
                // new communication cost incurred
                memset(current_neighbours_in_partition,0,num_processes * sizeof(int));
                memset(current_neighbours_elsewhere,0,num_processes * sizeof(int));
                memset(comm_cost_per_partition,0,num_processes * sizeof(double));
                for(int he = 0; he < hedge_ptr->at(vid).size(); he++) {
                    int he_id = hedge_ptr->at(vid)[he];
                    for(int vt = 0; vt < hyperedges->at(he_id).size(); vt++) {
                        int dest_vertex = hyperedges->at(he_id)[vt];
                        if(dest_vertex == vid) continue;
                        int dest_part = partitioning[dest_vertex];
                        current_neighbours_in_partition[dest_part] += 1; //  we may be counting twice a dest_node that appears in more than one hedge
                        // recalculate comm cost for all possible partition assignments of ii
                        //  commCost(v,Pi) = forall edge in edges(Pi) cost += w(e) * c(Pi,Pj) where i != j
                        for(int fp=0; fp < num_processes; fp++) {
                            comm_cost_per_partition[fp] += 1 * comm_cost_matrix[fp][dest_part];
                            if(fp != dest_part) {
                                current_neighbours_elsewhere[fp] += 1;
                            }
                        }
                    }
                }
                float max_value = std::numeric_limits<float>::lowest();
                int best_partition = partitioning[vid];
                for(int pp=0; pp < num_processes; pp++) {

                    // total cost of communication (edgecuts * number of participating partitions)
                    int total_comm_cost = 0;
                    for(int jj=0; jj < num_processes; jj++) {
                        if(pp != jj)
                            total_comm_cost += current_neighbours_in_partition[jj] > 0 ? 1 : 0;
                    }
                        
                    // testing objective function
                    //float current_value = current_neighbours_in_partition[pp] - current_neighbours_elsewhere[pp] - comm_cost_per_partition[pp]  - a * g/2 * pow(part_load[pp],g-1);
                    float current_value = -total_comm_cost * comm_cost_per_partition[pp]  - a * g/2 * pow(part_load[pp],g-1);
                        
                    // objective function is a mix of Battaglino 2015 (second part) and Zheng 2016 (communication cost part)
                    // (|P^t_i union N(v)| - commCost(v,Pi) - a * g/2 * |B|^(g-1))
                    //float current_value = current_neighbours_in_partition[pp] - comm_cost_per_partition[pp]  - a * g/2 * pow(part_load[pp],g-1);
                                        
                    // alternative from Battaglino 2015
                    //float current_value = current_neighbours_in_partition[pp] - a * g/2 * pow(part_load[pp],g-1);
                    // alternative from ARGO
                    //float current_value = (1.0f/(comm_cost_per_partition[pp]+1)) * (1-part_load[pp]/expected_workload);
                    // alternative from Alistarh 2015
                    //if(part_load[pp] >= expected_workload) continue;
                    //float current_value = current_neighbours_in_partition[pp];
                        
                    if(current_value > max_value) {
                        max_value = current_value;
                        best_partition = pp;
                    }
                }
                    
                // update intermediate workload and assignment values
                part_load[partitioning[vid]] -= vtx_wgt[vid];
                partitioning[vid] = best_partition;
                part_load[best_partition] += vtx_wgt[vid];
            }
                
            // check if desired imbalance has been reached
            float imbalance = calculateImbalance(partitioning,num_processes,num_vertices,vtx_wgt);
            PRINTF("%i: %f (%f | %f)\n",iter,imbalance,a,ta);
            if(imbalance < imbalance_tolerance) break;

            // update parameters
            a *= ta;
        }
        // clean up
        free(part_load);
        free(current_neighbours_in_partition);
    
        // return successfully
        return 0;
    }

    int ParallelIndependentRestreamingPartitioning(idx_t* partitioning, double** comm_cost_matrix, std::string hypergraph_filename, int* vtx_wgt, int iterations, float imbalance_tolerance) {
        
        int process_id;
        MPI_Comm_rank(MPI_COMM_WORLD,&process_id);
        int num_processes;
        MPI_Comm_size(MPI_COMM_WORLD,&num_processes);

        // get meta info (num vertices and hyperedges)
        int num_vertices, num_hyperedges;
        get_hypergraph_file_header(hypergraph_filename, &num_vertices, &num_hyperedges);
        
        //PARAMETERS: From Battaglino 2015 //
        // g and a determine load balance importance in cost function
        double g = 1.5;
        // battaglino's a
        //double a = sqrt(2) * num_hyperedges / pow(num_vertices,g);
        // our proposed a
        double a = sqrt(num_processes) * num_hyperedges / pow(num_vertices,g);
        // ta is the update rate of parameter a
        double ta = 1.7;
        // after how many vertices checked in the stream the partitio load is sync across processes
        int part_load_update_after_vertices = 100; // in the paper it is 4096
        // minimum number of iterations run (not checking imbalance threshold)
        // removed whilst we are using hyperPraw as refinement algorithm
        //      hence, if balanced is kept after first iteration, that's good enough
        int frozen_iters = 0;//ceil(0.1f * iterations);
        ///////////////
        
        // algorithm from GraSP (Battaglino 2016)
        // 1 - Distributed vertices over partitions (partition = vertex_id % num_partitions)
        // needs to load num_vertices from file
        for (int vid=0; vid < num_vertices; vid++) {
            //partitioning[vid] = vid % num_processes;
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
        int* part_load_update = (int*)calloc(num_processes,sizeof(int));
        int* part_load_update_recv = (int*)malloc(num_processes*sizeof(int));
        bool* communicating_with = (bool*)malloc(num_processes*sizeof(bool));
        //double timing = 0;
        //double ttt;
        for(int iter=0; iter < iterations; iter++) {
            //timing = 0;
            memset(local_stream_partitioning,0,num_vertices * sizeof(idx_t));
            memset(part_load,0,num_processes * sizeof(long int));
            memset(part_load_update,0,num_processes * sizeof(int));
            for(int ii=0; ii < num_vertices; ii++) {
                part_load[partitioning[ii]] += vtx_wgt[ii]; // workload for vertex
            }
            // go through own vertex list and reassign
            for(int vid=0; vid < num_vertices; vid++) {
                // share updated partition loads after constant number of iterations
                if(vid % part_load_update_after_vertices == 0) {
                    MPI_Allreduce(part_load_update,part_load_update_recv,num_processes,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
                    for(int pl=0; pl < num_processes; pl++) {
                        if(pl == process_id) continue;
                        part_load[pl] += part_load_update_recv[pl];
                    }
                    memset(part_load_update,0,num_processes * sizeof(int));
                }
                // always iterate through the same local list of vertices
                //if(vid % num_processes != process_id) continue;
                if(hedge_ptr[vid].size() == 0) continue;
                // reevaluate objective function per partition
                // |P^t_i union N(v)| = number of vertices in partition i that are neighbours of vertex v 
                // where are neighbours located
                // new communication cost incurred
                //ttt = MPI_Wtime();
                memset(current_neighbours_in_partition,0,num_processes * sizeof(int));
                memset(comm_cost_per_partition,0,num_processes * sizeof(double));
                memset(communicating_with,0,num_processes * sizeof(bool));
                for(int he = 0; he < hedge_ptr[vid].size(); he++) {
                    int he_id = hedge_ptr[vid][he];
                    for(int vt = 0; vt < hyperedges[he_id].size(); vt++) {
                        int dest_vertex = hyperedges[he_id][vt];
                        if(dest_vertex == vid) continue;
                        int dest_part = partitioning[dest_vertex];
                        current_neighbours_in_partition[dest_part] += 1;
                        // recalculate comm cost for all possible partition assignments of vid
                        //  commCost(v,Pi) = forall edge in edges(Pi) cost += w(e) * c(Pi,Pj) where i != j
                        for(int fp=0; fp < num_processes; fp++) {
                            comm_cost_per_partition[fp] += 1 * comm_cost_matrix[fp][dest_part];
                        }
                    }
                }
                //timing += MPI_Wtime() - ttt;
                double max_value = std::numeric_limits<double>::lowest();
                int best_partition = partitioning[vid];
                std::vector<int> best_parts;
                best_parts.push_back(partitioning[vid]);
                for(int pp=0; pp < num_processes; pp++) {
                    // total cost of communication (edgecuts * number of participating partitions)
                    long int total_comm_cost = 0;
                    for(int jj=0; jj < num_processes; jj++) {
                        if(pp != jj)
                            total_comm_cost += current_neighbours_in_partition[jj] > 0 ? 1 : 0;
                    }
                    
                    double current_value = -(double)total_comm_cost/(double)num_processes * comm_cost_per_partition[pp]  - a * g/2 * pow(part_load[pp],g-1);
                    
                    
                    //float current_value = current_neighbours_in_partition[pp] - comm_cost_per_partition[pp]  - a * g/2 * pow(part_load[pp],g-1);
                    
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
                        //best_partition = pp;
                        best_parts.clear();
                        best_parts.push_back(pp);
                    } else if(fabs(current_value-max_value) <= std::numeric_limits<double>::epsilon()*2) {
                        best_parts.push_back(pp);
                    }
                }
                
                best_partition = best_parts[(int)(best_parts.size() * (double)rand() / (double)RAND_MAX)];
                
                local_stream_partitioning[vid] = best_partition;
                // update intermediate workload and assignment values
                // Battaglino does not update this array mid stream
                part_load[best_partition] += vtx_wgt[vid];
                part_load[partitioning[vid]] -= vtx_wgt[vid];
                // update local changes counter
                part_load_update[partitioning[vid]] -= vtx_wgt[vid];
                part_load_update[best_partition] += vtx_wgt[vid];
                // update partitioning assignment
                // Battaglino does not update this array mid stream
                partitioning[vid] = best_partition;
                
            }
            //printf("%f\n",timing);
            
            // share new partitioning with other streams
            MPI_Allreduce(local_stream_partitioning,partitioning,num_vertices,MPI_LONG,MPI_MAX,MPI_COMM_WORLD);

            // check if desired imbalance has been reached
            float imbalance = calculateImbalance(partitioning,num_processes,num_vertices,vtx_wgt);
            PRINTF("%i: %f (%f | %f)\n",iter,imbalance,a,ta);

#ifdef SAVE_HISTORY
            if(process_id == 0) {
                FILE *fp = fopen(history_file.c_str(), "ab+");
                if(fp == NULL) {
                    printf("Error when storing partitioning history into file\n");
                } else {
                    fprintf(fp,"%.3f\n",imbalance);
                }
                fclose(fp);
            }
            
#endif
            if(frozen_iters <= iter && imbalance < imbalance_tolerance) break;

            // update parameters
            a *= ta;
        }

        // clean up
        free(local_stream_partitioning);
        free(part_load);
        free(comm_cost_per_partition);
        free(current_neighbours_in_partition);
        free(part_load_update_recv);
        
        return 0;
    }
    

    
}

#endif