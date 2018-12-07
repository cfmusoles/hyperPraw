// Parallel Restreaming Architecture aWare partitioning algorithm

#ifndef PRAW__H
#define PRAW__H

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
        int* workload = (int*)calloc(num_processes,sizeof(int));
        int total_workload = 0;
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

    void storePartitionStats(std::string experiment_name, idx_t* partitioning, int num_processes, int num_vertices, std::vector<std::vector<int> >* hyperedges, std::vector<std::vector<int> >* hedge_ptr, int* vtx_wgt,double** comm_cost_matrix) {
        // For hypergraph partitioning, two common metrics for objective functions Heuer 2017 (connectivity vs hyperedge cut):
        // communication cost estimate (given network bandwidth estimates)
        // quality measurements from hMETIS: Sum of External Degrees, Absorption
        int* workload = (int*)calloc(num_processes,sizeof(int));
        int total_workload = 0;
        int edgecut = 0;
        int hyperedges_cut = 0;
        int total_edges = 0;
        int soed = 0; // Sum of External degrees
        int connectivity_metric = 0; // hyperedges cut weighted by the number of participating partitions - 1
        int* vertices_in_partition = (int*)calloc(num_processes,sizeof(int)); // used for absorption
        float absorption = 0;
        std::set<int> hyperedges_visited;
        double total_comm_cost = 0;

#ifdef EVALUATE_PARTITIONING
        // Evaluate partitioning by predicting amount of data transferred in null compute sim
        int** part_connectivity = (int**)malloc(sizeof(int*) * num_processes);
        for(int ii=0; ii < num_processes; ii++) {
            part_connectivity[ii] = (int*)calloc(num_processes,sizeof(int));
        }
#endif

        for(int vid=0; vid < num_vertices; vid++) {
            if(vtx_wgt != NULL) {
                workload[partitioning[vid]] += vtx_wgt[vid];
                total_workload += vtx_wgt[vid];
            } else {
                workload[partitioning[vid]] += 1;
                total_workload += 1;
            }
            for(int he = 0; he < hedge_ptr->at(vid).size(); he++) {
                int he_id = hedge_ptr->at(vid)[he];
                bool visited = hyperedges_visited.count(he_id) > 0;
                total_edges++;
                std::set<int> connectivity;
                connectivity.insert(partitioning[vid]);
                memset(vertices_in_partition,0,sizeof(int) * num_processes);
                vertices_in_partition[partitioning[vid]]++;
                for(int vt = 0; vt < hyperedges->at(he_id).size(); vt++) {
                    int dest_vertex = hyperedges->at(he_id)[vt];
                    if(dest_vertex == vid) continue;
                    int dest_part = partitioning[dest_vertex];
                    vertices_in_partition[dest_part]++;
                    if(dest_part != partitioning[vid]) {
                        edgecut++;
                        connectivity.insert(dest_part);
#ifdef EVALUATE_PARTITIONING
                        part_connectivity[partitioning[vid]][dest_part] += 1;
#endif
                    }
                    if(comm_cost_matrix != NULL) {
                        total_comm_cost += comm_cost_matrix[partitioning[vid]][dest_part];
                    }
                    
                    total_edges++;
                }
                if(!visited) {
                    connectivity_metric += connectivity.size()-1;
                    if(connectivity.size() > 1) {
                        soed += connectivity.size(); // counts as 1 external degree per partition participating
                        hyperedges_cut++;
                    }
                    hyperedges_visited.insert(he_id);
                    for(int ii = 0; ii < num_processes; ii++) {
                        if(vertices_in_partition[ii] == 0) continue;
                        absorption += (float)(vertices_in_partition[ii]-1) / (float)(hyperedges->at(he_id).size()-1);
                    } 
                }
            }
        }
        float max_imbalance = 0;
        float expected_work = total_workload / num_processes;
        for(int ii=0; ii < num_processes; ii++) {
            PRINTF("%i workload: %i\n",ii,workload[ii]);
            float imbalance = (workload[ii] / expected_work);
            if(imbalance > max_imbalance)
                max_imbalance = imbalance;
        }

#ifdef EVALUATE_PARTITIONING
        // store partitioning connectivity in file
        std::string fname = "prediction_comm_matrix_";
        char str_to_int[16];
        sprintf(str_to_int,"%i",num_processes);
		fname +=  str_to_int;
        FILE *comm_matrix_file = fopen(fname.c_str(), "w");
        if(comm_matrix_file == NULL) {
		    printf("Error when storing results into file\n");
		} else {
            for(int from =0; from < num_processes; from++) {
                for(int to =0; to < num_processes; to++) {
                    fprintf(comm_matrix_file,"%i ",part_connectivity[from][to]);
                }
                fprintf(comm_matrix_file,"\n");
            }
        }
        fclose(comm_matrix_file);
        // clear memory
        for(int ii=0; ii < num_processes; ii++) {
            free(part_connectivity[ii]);
        }
        free(part_connectivity);
#endif

        printf("Quality: %i (hedgecut, %.3f total) %.3f (cut net), %i (connectivity metric), %i (SOED), %.1f (absorption) %.3f (max imbalance), %f (comm cost)\n",hyperedges_cut,(float)hyperedges_cut/hyperedges->size(),(float)edgecut/total_edges,connectivity_metric,soed,absorption,max_imbalance,total_comm_cost);
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
			    fprintf(fp,"%s,%s,%s,%s,%s,%s,%s,%s\n","Quality (hedges cut)","Hedge cut ratio","Cut net","Connectivity metric","SOED","Absorption","Max imbalance","Comm cost");
			fprintf(fp,"%i,%f,%f,%i,%i,%f,%f,%f\n",hyperedges_cut,(float)hyperedges_cut/hyperedges->size(),(float)edgecut/total_edges,connectivity_metric,soed,absorption,max_imbalance,total_comm_cost);
        }
        fclose(fp);

        // clean up
        free(workload);
        free(vertices_in_partition);
    }

    // Load hMETIS file (hypergraph)	
    int load_hypergraph_from_file(std::string filename, std::vector<std::vector<int> >* hyperedges, std::vector<std::vector<int> >* hedge_ptr) {
        
        std::ifstream istream(filename.c_str());
        
        if(!istream) {
            printf("Error while opening hMETIS file %s\n",filename.c_str());
            return -1;
        }

        std::string line;
        // process header
        std::getline(istream,line);
        std::istringstream buf(line);
        std::istream_iterator<int> beg(buf), end;
        std::vector<int> tokens(beg, end);
        int total_vertices = tokens[1];
        int total_hyperedges = tokens[0];
        hyperedges->resize(total_hyperedges);
        hedge_ptr->resize(total_vertices);
        
        PRINTF("Loaded from file: Vertices: %i; hyperedges %i:\n",total_vertices,total_hyperedges);
        
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

    int ParallelStreamingPartitioning(idx_t* partitioning, double** comm_cost_matrix, int num_vertices, std::vector<std::vector<int> >* hyperedges, std::vector<std::vector<int> >* hedge_ptr, int* vtx_wgt, int max_iterations, float imbalance_tolerance) {
        //PARAMETERS: From Battaglino 2015 //
        float g = 1.5;
        float a = sqrt(2) * hyperedges->size() / pow(num_vertices,g);
        float ta = 1.7;
        //float tta = 0.98;
        ///////////////
        
        int process_id;
        MPI_Comm_rank(MPI_COMM_WORLD,&process_id);
        int num_processes;
        MPI_Comm_size(MPI_COMM_WORLD,&num_processes);

        //printPartitionStats(partitioning,num_processes,num_vertices,hyperedges,hedge_ptr,vtx_wgt,comm_cost_matrix);
        
        int* part_load = (int*)calloc(num_processes,sizeof(int));
        idx_t* local_stream_partitioning = (idx_t*)malloc(num_vertices*sizeof(idx_t));
        double* comm_cost_per_partition = (double*)malloc(num_processes*sizeof(double));
        int* current_neighbours_in_partition = (int*)malloc(num_processes*sizeof(int));
                
        for(int iter=0; iter < max_iterations; iter++) {
            memset(local_stream_partitioning,0,num_vertices * sizeof(idx_t));
            memset(part_load,0,num_processes * sizeof(int));
            for(int ii=0; ii < num_vertices; ii++) {
                part_load[partitioning[ii]] += vtx_wgt[ii]; // workload for vertex
            }
            // go through own vertex list and reassign
            for(int ii=0; ii < num_vertices; ii++) {
                if(partitioning[ii] != process_id) continue;
                // reevaluate objective function per partition
                // |P^t_i union N(v)| = number of vertices in partition i that are neighbours of vertex v 
                // where are neighbours located
                // new communication cost incurred
                memset(current_neighbours_in_partition,0,num_processes * sizeof(int));
                memset(comm_cost_per_partition,0,num_processes * sizeof(double));
                for(int he = 0; he < hedge_ptr->at(ii).size(); he++) {
                    int he_id = hedge_ptr->at(ii)[he];
                    for(int vt = 0; vt < hyperedges->at(he_id).size(); vt++) {
                        int dest_vertex = hyperedges->at(he_id)[vt];
                        if(dest_vertex == ii) continue;
                        int dest_part = partitioning[dest_vertex];
                        current_neighbours_in_partition[dest_part] += 1;
                        // recalculate comm cost for all possible partition assignments of ii
                        //  commCost(v,Pi) = forall edge in edges(Pi) cost += w(e) * c(Pi,Pj) where i != j
                        for(int fp=0; fp < num_processes; fp++) {
                            comm_cost_per_partition[fp] += 1 * comm_cost_matrix[fp][dest_part];
                        }
                    }
                }
                float max_value = std::numeric_limits<float>::lowest();
                int best_partition = partitioning[ii];
                for(int pp=0; pp < num_processes; pp++) {
                    // objective function is a mix of Battaglino 2015 (second part) and Zheng 2016 (communication cost part)
                    // (|P^t_i union N(v)| - commCost(v,Pi) - a * g/2 * |B|^(g-1))
                    float current_value = (float)current_neighbours_in_partition[pp] - comm_cost_per_partition[pp]  - a * g/2 * pow(part_load[pp],g-1);
                    if(current_value > max_value) {
                        max_value = current_value;
                        best_partition = pp;
                    }
                }
                    
                local_stream_partitioning[ii] = best_partition;
                // update intermediate workload and assignment values
                part_load[best_partition] += vtx_wgt[ii];
                part_load[partitioning[ii]] -= vtx_wgt[ii];
                partitioning[ii] = best_partition;
            }
            
            // share new partitioning with other streams
            MPI_Allreduce(local_stream_partitioning,partitioning,num_vertices,MPI_LONG,MPI_MAX,MPI_COMM_WORLD);
            
            // check if desired imbalance has been reached
            float imbalance = calculateImbalance(partitioning,num_processes,num_vertices,vtx_wgt);
            PRINTF("%i: %f (%f | %f)\n",iter,imbalance,a,ta);
            if(imbalance < imbalance_tolerance) break;

            //printf("Iteration %i\n",iter);
            //printPartitionStats(partitioning,num_processes,num_vertices,hyperedges,hedge_ptr,vtx_wgt,comm_cost_matrix);

            // update parameters
            a *= ta;
            //if(ta > 1.05f) ta *= tta;
        }
        
        // clean up
        free(part_load);
        free(local_stream_partitioning);
        free(comm_cost_per_partition);
        free(current_neighbours_in_partition);

        // return successfully
        return 0;
    }

    int SequentialStreamingPartitioning(idx_t* partitioning, double** comm_cost_matrix, int num_vertices, std::vector<std::vector<int> >* hyperedges, std::vector<std::vector<int> >* hedge_ptr, int* vtx_wgt, int max_iterations, float imbalance_tolerance) {
        
        //PARAMETERS: From Battaglino 2015 //
        float g = 1.5;
        float a = sqrt(2) * hyperedges->size() / pow(num_vertices,g);
        float ta = 1.7;
        ///////////////
        
        int process_id;
        MPI_Comm_rank(MPI_COMM_WORLD,&process_id);
        int num_processes;
        MPI_Comm_size(MPI_COMM_WORLD,&num_processes);

        if(process_id == 0) {
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
                        // testing objective function
                        float current_value = current_neighbours_in_partition[pp] - current_neighbours_elsewhere[pp] - comm_cost_per_partition[pp]  - a * g/2 * pow(part_load[pp],g-1);
                        
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
        }

        // share new partitioning with other processes
        MPI_Bcast(partitioning,num_vertices,MPI_LONG,0,MPI_COMM_WORLD);

        //printPartitionStats(partitioning,num_processes,num_vertices,hyperedges,hedge_ptr,vtx_wgt,comm_cost_matrix);
        
        // return successfully
        return 0;
    }

    
}

#endif