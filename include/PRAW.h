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
#include <unordered_map>

#ifdef VERBOSE
#define PRINTF(...) printf(__VA_ARGS__)
#else
#define PRINTF(...) 
#endif


typedef int64_t idx_t; // needs to probably match the type in METIS

namespace PRAW {

    const int MASTER_NODE = 0;

    // data structure used by the HDFR algorithm
    struct pin_data {
        int partial_degree = 0; // how many times has the pin  appeared in the stream so far
        std::set<int> A;  // in what partitions it has been replicated so far
    };

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
            //workload[partitioning[vid]] += vtx_wgt[vid];
            //total_workload += vtx_wgt[vid];
            if(vtx_wgt != NULL) {
                workload[partitioning[vid]] += vtx_wgt[vid];
                total_workload += vtx_wgt[vid];
            } else {
                workload[partitioning[vid]] += 1;
                total_workload += 1;
            }
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
                            float* hyperedges_cut_ratio, float* edges_cut_ratio, int* soed, float* absorption, float* max_imbalance, double* total_edge_comm_cost, double* total_hedge_comm_cost) { // output
        
        // For hypergraph partitioning, two common metrics for objective functions Heuer 2017 (connectivity vs hyperedge cut):
        // communication cost estimate (given network bandwidth estimates)
        // quality measurements from hMETIS: Sum of External Degrees, Absorption
        *hyperedges_cut_ratio=0;
        *edges_cut_ratio=0;
        *soed=0;
        *absorption=0;
        *max_imbalance=0;
        *total_edge_comm_cost=0;
        *total_hedge_comm_cost=0;

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
                        *total_edge_comm_cost += comm_cost_matrix[partitioning[from]][to_part];
                    }
                }
            }
            // metrics per hyperedge
            if(connectivity.size() > 1) {
                *soed += connectivity.size(); // counts as 1 external degree per partition participating
                hyperedges_cut++;
                // communication cost based on hyperedge cut
                for (std::set<int>::iterator sender=connectivity.begin(); sender!=connectivity.end(); ++sender) {
                    int sender_id = *sender;
                    for (std::set<int>::iterator receiver=connectivity.begin(); receiver!=connectivity.end(); ++receiver) {
                        int receiver_id = *receiver;    
                        if (sender_id != receiver_id)
                            *total_hedge_comm_cost += comm_cost_matrix[sender_id][receiver_id];
                    }
                }
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
        

        PRINTF("Quality: %i (hedgecut, %.3f total) %.3f (cut net), %i (SOED), %.1f (absorption) %.3f (max imbalance), %f (edge cost), %f (hedge cost)\n",hyperedges_cut,*hyperedges_cut_ratio,*edges_cut_ratio,*soed,*absorption,*max_imbalance,*total_edge_comm_cost,*total_hedge_comm_cost);
        
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
        double total_edge_comm_cost;
        double total_hedge_comm_cost;
        getPartitionStats(partitioning, num_processes, num_vertices, hyperedges, hedge_ptr, vtx_wgt,comm_cost_matrix,
                            &hyperedges_cut_ratio, &edges_cut_ratio, &soed, &absorption, &max_imbalance, &total_edge_comm_cost,&total_hedge_comm_cost);
        
        printf("Hedgecut, %.3f, %.3f (cut net), %i (SOED), %.1f (absorption) %.3f (max imbalance), %f (edge comm cost), %f (hedge cost)\n",hyperedges_cut_ratio,edges_cut_ratio,soed,absorption,max_imbalance,total_edge_comm_cost,total_hedge_comm_cost);
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
			    fprintf(fp,"%s,%s,%s,%s,%s,%s,%s\n","Hedge cut ratio","Cut net","SOED","Absorption","Max imbalance","Edge comm cost","Hedge comm cost");
			fprintf(fp,"%f,%f,%i,%f,%f,%f,%f\n",hyperedges_cut_ratio,edges_cut_ratio,soed,absorption,max_imbalance,total_edge_comm_cost,total_hedge_comm_cost);
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

        if(hyperedges != NULL) hyperedges->resize(total_hyperedges);
        if(hedge_ptr != NULL) hedge_ptr->resize(total_vertices);
        
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
            /*std::istringstream buf(line);
            std::istream_iterator<int> beg(buf), end;
            std::vector<int> tokens(beg, end);*/

            char str[line.length() + 1]; 
            strcpy(str, line.c_str()); 
            char* token = strtok(str, " "); 
            std::vector<int> tokens;
            while (token != NULL) { 
                tokens.push_back(atoi(token)); 
                token = strtok(NULL, " "); 
            } 

            for(int ii=0; ii < tokens.size(); ii++) {
                int vertex_id = tokens[ii]-1;
                if(hedge_ptr != NULL) hedge_ptr->at(vertex_id).push_back(counter);
                if(hyperedges != NULL) hyperedges->at(counter).push_back(vertex_id);
            }
            counter++;
            
        }
        istream.close();
        return 0;
        
    }

    // load to distributed CSR format
    int load_hypergraph_from_file_dist_CSR(std::string filename, std::vector<std::vector<int> >* hyperedges, std::vector<std::vector<int> >* hedge_ptr, int process_id, idx_t* partitioning, bool load_all = false) {
        
        // get header info
        int total_vertices;;
        int total_hyperedges;
        get_hypergraph_file_header(filename,&total_vertices,&total_hyperedges);

        if(hyperedges != NULL) hyperedges->resize(total_hyperedges);
        if(hedge_ptr != NULL) hedge_ptr->resize(total_vertices);
        
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
            /*std::istringstream buf(line);
            std::istream_iterator<int> beg(buf), end;
            std::vector<int> tokens(beg, end);*/

            char str[line.length() + 1]; 
            strcpy(str, line.c_str()); 
            char* token = strtok(str, " "); 
            std::vector<int> tokens;
            while (token != NULL) { 
                tokens.push_back(atoi(token)); 
                token = strtok(NULL, " "); 
            } 

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
            if(local || load_all) {
                for(int ii=0; ii < tokens.size(); ii++) {
                    int vertex_id = tokens[ii]-1;
                    if(partitioning[vertex_id] == process_id || load_all) {
                        if(hedge_ptr != NULL) hedge_ptr->at(vertex_id).push_back(counter);
                    }
                    if(hyperedges != NULL) hyperedges->at(counter).push_back(vertex_id);
                }
            } else nonlocal++;
            counter++;            
        }
        istream.close();
        PRINTF("%i: non local hyperedges: %i\n",process_id,nonlocal);
        return 0;
        
    }

    // load hyperedge ids per vertex (only local vertex stream)
    int load_hedge_ptr_from_file_dist_CSR(std::string filename, std::vector<std::vector<int> >* hedge_ptr, int process_id, int num_processes, idx_t* partitioning) {
        
        // get header info
        int total_vertices;
        int total_hyperedges;
        get_hypergraph_file_header(filename,&total_vertices,&total_hyperedges);

        hedge_ptr->resize(ceil((float)total_vertices / num_processes));
        
        std::ifstream istream(filename.c_str());
        
        if(!istream) {
            printf("Error while opening hMETIS file %s\n",filename.c_str());
            return -1;
        }
        
        std::string line;
        // skip header
        std::getline(istream,line);
        // read reminder of file (one line per hyperedge)
        int counter = 1;
        while(std::getline(istream,line)) {
            char str[line.length() + 1]; 
            strcpy(str, line.c_str()); 
            char* token = strtok(str, " "); 
            std::vector<int> tokens;
            while (token != NULL) { 
                tokens.push_back(atoi(token)); 
                token = strtok(NULL, " "); 
            } 

            for(int ii=0; ii < tokens.size(); ii++) {
                int vertex_id = tokens[ii]-1;
                if(vertex_id % num_processes == process_id) {
                    hedge_ptr->at(vertex_id / num_processes).push_back(counter);
                }
            }
            counter++;            
        }
        istream.close();

        return 0;
        
    }

    // vertex centric partitioning stats for hypergraphs
    // SOED
    // Hyperedge cut
    // Absorption
    // imbalance
    // sim comm cost
    void getVertexCentricPartitionStatsFromFile(idx_t* partitioning, int num_processes, std::string hgraph_filename, int* vtx_wgt,double** comm_cost_matrix, // input
                            float* hyperedges_cut_ratio, float* edges_cut_ratio, int* soed, float* absorption, float* max_imbalance, double* total_edge_comm_cost, double* total_hedge_comm_cost,float* vertex_replication_factor) { // output

        // take partitioning statistics directly from reading a graph file
        // loading the entire graph on a process should be avoided for scalability
        *hyperedges_cut_ratio=0;
        *edges_cut_ratio=0;
        *soed=0;
        *absorption=0;
        *max_imbalance=0;
        *total_edge_comm_cost=0;
        *vertex_replication_factor=0;
        if(total_hedge_comm_cost != NULL) *total_hedge_comm_cost=0;
        
        //int* workload = (int*)calloc(num_processes,sizeof(int));
        int total_workload = 0;
        long int edgecut = 0;
        int hyperedges_cut = 0;
        long int total_edges = 0;
        int* vertices_in_partition = (int*)calloc(num_processes,sizeof(int)); // used for absorption
        
        // get header info
        int total_vertices;
        int total_hyperedges;
        get_hypergraph_file_header(hgraph_filename,&total_vertices,&total_hyperedges);
        
        std::ifstream istream(hgraph_filename.c_str());

        
        if(!istream) {
            printf("Error while opening hMETIS file %s\n",hgraph_filename.c_str());
            return;
        }
        std::string line;

        // compute workload for all vertices and partitions
        *max_imbalance = calculateImbalance(partitioning, num_processes, total_vertices, vtx_wgt);

        // keep track of vertex replicas on each partition with a set per partition
        std::vector<std::set<int> > vertex_replicas(num_processes);

        // skip header
        std::getline(istream,line);
        // read reminder of file (one line per hyperedge)
        while(std::getline(istream,line)) {
            // each line corresponds to a hyperedge

            char str[line.length() + 1]; 
            strcpy(str, line.c_str()); 
            char* token = strtok(str, " "); 
            std::vector<int> tokens;
            while (token != NULL) { 
                tokens.push_back(atoi(token)); 
                token = strtok(NULL, " "); 
            } 

            std::set<int> connectivity;
            memset(vertices_in_partition,0,sizeof(int) * num_processes);
            for(int ff=0; ff < tokens.size(); ff++) {
                // each vertex in the hyperedge
                int from = tokens[ff]-1;
                connectivity.insert(partitioning[from]);
                vertices_in_partition[partitioning[from]]++;
                vertex_replicas[partitioning[from]].insert(from);
                for(int tt=0; tt < tokens.size(); tt++) {
                    if(tt==ff) continue;
                    int to = tokens[tt]-1;
                    total_edges++;
                    int to_part = partitioning[to];
                    //vertices_in_partition[to_part]++;
                    vertex_replicas[to_part].insert(to);
                    if(to_part != partitioning[from]) {
                        edgecut++;
                        connectivity.insert(to_part);
                        vertex_replicas[to_part].insert(from);
                        vertex_replicas[partitioning[from]].insert(to);
                    }
                    // this cost measures mainly edge cut
                    if(comm_cost_matrix != NULL) {
                        *total_edge_comm_cost += comm_cost_matrix[partitioning[from]][to_part];
                    }
                }
            }
            // metrics per hyperedge
            if(connectivity.size() > 1) {
                *soed += connectivity.size(); // counts as 1 external degree per partition participating
                hyperedges_cut++;
                if(total_hedge_comm_cost != NULL) {
                // communication cost based on hyperedge cut
                    for (std::set<int>::iterator sender=connectivity.begin(); sender!=connectivity.end(); ++sender) {
                        int sender_id = *sender;
                        for (std::set<int>::iterator receiver=connectivity.begin(); receiver!=connectivity.end(); ++receiver) {
                            int receiver_id = *receiver;    
                            if (sender_id != receiver_id)
                                *total_hedge_comm_cost += comm_cost_matrix[sender_id][receiver_id];
                        }
                    }
                }
            }
            if(tokens.size() > 1) {
                for(int pp = 0; pp < num_processes; pp++) {
                    if(vertices_in_partition[pp] == 0) continue;
                    *absorption += (float)(vertices_in_partition[pp]-1) / (float)(tokens.size()-1);
                } 
            }
            
        }
        istream.close();

        *hyperedges_cut_ratio=(float)hyperedges_cut/total_hyperedges;
        *edges_cut_ratio=(float)edgecut/total_edges;

        // total vertex replica factor is the sum of all vertex lists (of each partition) over the number vertices
        for(int ii=0; ii < num_processes; ii++) {
            *vertex_replication_factor += vertex_replicas[ii].size();
        }
        *vertex_replication_factor /= total_vertices;
        

        PRINTF("Quality: %i (hedgecut, %.3f total) %.3f (cut net), %i (SOED), %.1f (absorption) %.3f (max imbalance), %f (edge comm cost),, %f (hedge comm cost)\n",hyperedges_cut,*hyperedges_cut_ratio,*edges_cut_ratio,*soed,*absorption,*max_imbalance,*total_edge_comm_cost,total_hedge_comm_cost == NULL ? 0 : *total_hedge_comm_cost);
        
        // clean up
        free(vertices_in_partition);
    }

    // parallel version of the stats getter
    void getVertexCentricPartitionStatsFromFile_parallel(idx_t* partitioning, int num_processes, int process_id, std::string hgraph_filename, int* vtx_wgt,double** comm_cost_matrix, // input
                            float* hyperedges_cut_ratio, float* edges_cut_ratio, int* soed, float* absorption, float* max_imbalance, double* total_edge_comm_cost, double* total_hedge_comm_cost) { // output

        // take partitioning statistics directly from reading a graph file
        // loading the entire graph on a process should be avoided for scalability
        *hyperedges_cut_ratio=0;
        *edges_cut_ratio=0;
        *soed=0;
        *absorption=0;
        *max_imbalance=0;
        *total_edge_comm_cost=0;
        if(total_hedge_comm_cost != NULL) *total_hedge_comm_cost=0;
        
        //int* workload = (int*)calloc(num_processes,sizeof(int));
        int total_workload = 0;
        long int edgecut = 0;
        int hyperedges_cut = 0;
        long int total_edges = 0;
        int* vertices_in_partition = (int*)calloc(num_processes,sizeof(int)); // used for absorption
        
        // get header info
        int total_vertices;
        int total_hyperedges;
        get_hypergraph_file_header(hgraph_filename,&total_vertices,&total_hyperedges);
        
        std::ifstream istream(hgraph_filename.c_str());

        
        if(!istream) {
            printf("Error while opening hMETIS file %s\n",hgraph_filename.c_str());
            return;
        }
        std::string line;

        // compute workload for all vertices and partitions
        *max_imbalance = calculateImbalance(partitioning, num_processes, total_vertices, vtx_wgt);

        // skip header
        std::getline(istream,line);
        // read reminder of file (one line per hyperedge)
        int current_line = 1;
        while(std::getline(istream,line)) {
            // each line corresponds to a hyperedge
            current_line++;
            if(current_line % num_processes != process_id) continue;

            char str[line.length() + 1]; 
            strcpy(str, line.c_str()); 
            char* token = strtok(str, " "); 
            std::vector<int> tokens;
            while (token != NULL) { 
                tokens.push_back(atoi(token)); 
                token = strtok(NULL, " "); 
            } 

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
                        *total_edge_comm_cost += comm_cost_matrix[partitioning[from]][to_part];
                    }
                }
            }
            // metrics per hyperedge
            if(connectivity.size() > 1) {
                *soed += connectivity.size(); // counts as 1 external degree per partition participating
                hyperedges_cut++;
                if(total_hedge_comm_cost != NULL) {
                // communication cost based on hyperedge cut
                    for (std::set<int>::iterator sender=connectivity.begin(); sender!=connectivity.end(); ++sender) {
                        int sender_id = *sender;
                        for (std::set<int>::iterator receiver=connectivity.begin(); receiver!=connectivity.end(); ++receiver) {
                            int receiver_id = *receiver;    
                            if (sender_id != receiver_id)
                                *total_hedge_comm_cost += comm_cost_matrix[sender_id][receiver_id];
                        }
                    }
                }
            }
            if(tokens.size() > 1) {
                for(int pp = 0; pp < num_processes; pp++) {
                    if(vertices_in_partition[pp] == 0) continue;
                    *absorption += (float)(vertices_in_partition[pp]-1) / (float)(tokens.size()-1);
                } 
            }
            
        }
        istream.close();

        // gather parallel results
        int32_t total_hedgecut = 0;
        MPI_Allreduce(&hyperedges_cut,&total_hedgecut,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        hyperedges_cut = total_hedgecut;
        long int total_edgecut = 0;
        MPI_Allreduce(&edgecut,&total_edgecut,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
        edgecut = total_edgecut;
        long int sum_edges = 0;
        MPI_Allreduce(&total_edges,&sum_edges,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
        total_edges = sum_edges;
        int total_soed = 0;
        MPI_Allreduce(soed,&total_soed,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        *soed = total_soed;
        float total_absorption = 0;
        MPI_Allreduce(absorption,&total_absorption,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
        *absorption = total_absorption;
        double total_edgecost = 0;
        MPI_Allreduce(total_edge_comm_cost,&total_edgecost,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        *total_edge_comm_cost = total_edgecost;
        if(total_hedge_comm_cost != NULL) {
            double total_hedgecost = 0;
            MPI_Allreduce(total_hedge_comm_cost,&total_hedgecost,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            *total_hedge_comm_cost = total_hedgecost;
        }


        *hyperedges_cut_ratio=(float)hyperedges_cut/total_hyperedges;
        *edges_cut_ratio=(float)edgecut/total_edges;
        

        if(process_id == MASTER_NODE)
            PRINTF("Quality: %i (hedgecut, %.3f total) %.3f (cut net), %i (SOED), %.1f (absorption) %.3f (max imbalance), %f (edge comm cost),, %f (hedge comm cost)\n",hyperedges_cut,*hyperedges_cut_ratio,*edges_cut_ratio,*soed,*absorption,*max_imbalance,*total_edge_comm_cost,total_hedge_comm_cost == NULL ? 0 : *total_hedge_comm_cost);
        
        // clean up
        free(vertices_in_partition);
    }

    // Get total edge comm cost for a vertex centric partitioning
    void getVertexCentricEdgeCommCostFromFile(idx_t* partitioning, int num_processes, int process_id, std::string hgraph_filename,double** comm_cost_matrix, // input
                                                double* total_edge_comm_cost) { // output

        // take partitioning statistics directly from reading a graph file
        // loading the entire graph on a process should be avoided for scalability
        *total_edge_comm_cost=0;
        
        std::ifstream istream(hgraph_filename.c_str());

        if(!istream) {
            printf("Error while opening hMETIS file %s\n",hgraph_filename.c_str());
            return;
        }
        std::string line;

        // skip header
        std::getline(istream,line);
        // read reminder of file (one line per hyperedge)
        int current_line = 0;
        while(std::getline(istream,line)) {
            current_line++;
            if(current_line % num_processes != process_id) continue;
            // each line corresponds to a hyperedge

            char str[line.length() + 1]; 
            strcpy(str, line.c_str()); 
            char* token = strtok(str, " "); 
            std::vector<int> tokens;
            while (token != NULL) { 
                tokens.push_back(atoi(token)); 
                token = strtok(NULL, " "); 
            } 

            for(int ff=0; ff < tokens.size(); ff++) {
                // each vertex in the hyperedge
                int from = tokens[ff]-1;
                for(int tt=0; tt < tokens.size(); tt++) {
                    if(tt==ff) continue;
                    int to = tokens[tt]-1;
                    int to_part = partitioning[to];
                    // this cost measures mainly edge cut
                    if(comm_cost_matrix != NULL) {
                        *total_edge_comm_cost += comm_cost_matrix[partitioning[from]][to_part];
                    }
                }
            }            
        }
        istream.close();
        
        // gather all parallel counts
        double total_cost = 0;
        MPI_Allreduce(total_edge_comm_cost,&total_cost,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        *total_edge_comm_cost = total_cost;
        
    }

    // hypergraph partition stats for edge partitioning
    // vertex replica factor
    // sim comm cost
    // imbalance
    void getEdgeCentricPartitionStatsFromFile(idx_t* partitioning, int num_processes, std::string hgraph_filename, int* he_wgt,double** comm_cost_matrix, // input
                            float* vertex_replication_factor, float* max_hedge_imbalance, float* hedgecut) { // output

        // take partitioning statistics directly from reading a graph file
        // loading the entire graph on a process should be avoided for scalability
        *vertex_replication_factor=0;
        *max_hedge_imbalance=0;
        *hedgecut=0;
        
        long int total_workload = 0;
        
        // get header info
        int total_vertices;
        int total_hyperedges;
        get_hypergraph_file_header(hgraph_filename,&total_vertices,&total_hyperedges);
        
        std::ifstream istream(hgraph_filename.c_str());

        
        if(!istream) {
            printf("Error while opening hMETIS file %s\n",hgraph_filename.c_str());
            return;
        }
        std::string line;

        // keep track of vertex replicas on each partition with a set per partition
        std::vector<std::set<int> > vertex_list(num_processes);
        // keep track of number of hyperedges per partition
        int* hedges_size = (int*)calloc(num_processes,sizeof(int));
        // skip header
        std::getline(istream,line);
        // read reminder of file (one line per hyperedge)
        int he_id = 0;
        while(std::getline(istream,line)) {
            // each line corresponds to a hyperedge
            char str[line.length() + 1]; 
            strcpy(str, line.c_str()); 
            char* token = strtok(str, " "); 
            std::vector<int> tokens;
            while (token != NULL) { 
                tokens.push_back(atoi(token)); 
                token = strtok(NULL, " "); 
            } 

            int dest_partition = partitioning[he_id];
            // Petroni uses number of edges on a partition to measure load balance
            // should we use number of hyperedges, or the number of vertices on a partition?
            // It depends on the modelled application
            // If the computation is proportional to number of hedges, then count number of hedges
            // Counting number of hedges:
            hedges_size[dest_partition] += he_wgt[he_id];
            total_workload += he_wgt[he_id];

            // add all vertices (replicas) to partition
            for(int ii=0; ii < tokens.size(); ii++) {
                int vertex_id = tokens[ii] - 1;
                // Counting the number of vertices per partition:
                vertex_list[dest_partition].insert(vertex_id);
            }
            he_id++;            
        }
        istream.close();

        // total vertex replica factor is the sum of all vertex lists (of each partition) over the number vertices
        for(int ii=0; ii < num_processes; ii++) {
            *vertex_replication_factor += vertex_list[ii].size();
            if(hedges_size[ii] > *max_hedge_imbalance) {
                *max_hedge_imbalance = hedges_size[ii];
            }
            
        }        
        *max_hedge_imbalance = *max_hedge_imbalance / ((float)total_workload / num_processes);
        *vertex_replication_factor = *vertex_replication_factor / total_vertices;

        // calculate hedgecut
        // only applicable if num vertices == num hyperedges
        if(total_vertices == total_hyperedges) {
            // a hyperedge is cut if at least one of its vertices reside in a different partition than the host (where the hyperedge id has been assigned)
            std::ifstream istream(hgraph_filename.c_str());
            // skip header
            std::getline(istream,line);
            // read reminder of file (one line per hyperedge)
            int hedges_cut = 0;
            int he_id = 0;
            while(std::getline(istream,line)) {
                // each line corresponds to a hyperedge
                char str[line.length() + 1]; 
                strcpy(str, line.c_str()); 
                char* token = strtok(str, " "); 
                while (token != NULL) { 
                    int vertex_id = atoi(token) - 1; 
                    if(partitioning[he_id] != partitioning[vertex_id]) {
                        hedges_cut++;
                        break;
                    }
                    token = strtok(NULL, " "); 
                } 
                he_id++;
            }
            istream.close();  

            *hedgecut = (float)hedges_cut / total_hyperedges;      
        } else {
            *hedgecut = -1;
        }


        PRINTF("Vertex replication factor: %.3f, %.3f (hedge imbalance) %.3f (hyperedge cut)\n",*vertex_replication_factor,*max_hedge_imbalance,*hedgecut);
        
        // clean up
        free(hedges_size);
        
    }

    // parallel version of edge centric stats
    void getEdgeCentricReplicationFactor(std::unordered_map<int,pin_data>* seen_vertices, int num_vertices, // input
                            float* vertex_replication_factor) { // output

        *vertex_replication_factor = 0;
        for(int vid = 0; vid < num_vertices; vid++) {
            std::unordered_map<int,pin_data>::const_iterator it = seen_vertices->find (vid);
            if(it == seen_vertices->end()) continue;
            *vertex_replication_factor += seen_vertices->at(vid).A.size();
        }
        *vertex_replication_factor /= num_vertices;
        
    }

    void get_comm_cost_matrix_from_bandwidth(char* comm_bandwidth_filename, double** comm_cost_matrix, int partitions, bool proportional_comm_cost) {
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
                            [proportional_comm_cost,min_bandwidth,max_bandwidth] (double value) {  
                                if(proportional_comm_cost) {
                                    float ratio = max_bandwidth / min_bandwidth;
                                    return value <= std::numeric_limits<float>::epsilon() ? 0 : (1-(value-min_bandwidth)/(max_bandwidth-min_bandwidth)) * ratio + 1;
                                } else {
                                    return value <= std::numeric_limits<float>::epsilon() ? 0 : 2 - ( (value-min_bandwidth)/(max_bandwidth-min_bandwidth) );
                                }
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

    int SequentialVertexPartitioning(char* experiment_name, idx_t* partitioning, int num_processes, double** comm_cost_matrix, std::string hypergraph_filename, int* vtx_wgt, int iterations, float imbalance_tolerance, float ta_refine, bool reset_partitioning, int stopping_condition, bool save_partitioning_history) {
        
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
        double ta_start = 1.7; // used when imbalance is above imbalance_tolerance
        //double ta_refine = 0.95; // used when imbalance is below imbalance_tolerance
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
                partitioning[vid] = vid % num_processes;//MASTER_NODE;
            }
        }
        
        // 2 - Divide the graph in a distributed CSR format (like ParMETIS)
        //  compressed vertex or compressed hedge format? --> see zoltan
        //  for each local vertex, store the list of vertices adjacent to it (belonging to same hedges)
        std::vector<std::vector<int> > hyperedges;
        std::vector<std::vector<int> > hedge_ptr;
        load_hypergraph_from_file_dist_CSR(hypergraph_filename, &hyperedges, &hedge_ptr, MASTER_NODE, partitioning, true);
        

        // each process must read from file only the info relevant to its data
        // 3 - Initiate N number of iterations on each process:
        //      a - one vertex at a time, assign to best partition (based on eval function)
        //      b - update tempering parameters
        //      c - share with all new partition assignments
        
        std::string history_file = experiment_name;
        
        if(save_partitioning_history) {
            history_file += "_";
            history_file += getFileName(hypergraph_filename);
            history_file += "_partition_history_";
            char str_int[16];
            sprintf(str_int,"%i",num_processes);
            history_file += "__";
            history_file +=  str_int;
            // remove history file if exists
            FILE *fp = fopen(history_file.c_str(), "w");
            if(fp == NULL) {
                printf("Error when storing partitioning history into file\n");
            } else {
                fprintf(fp,"%s\n","Imbalance, Hedges cut, Edges cut, SOED, Absorption, Edge sim comm cost, Hedge sim comm cost");
            }
            fclose(fp);
        }
        

        long int* part_load = (long int*)calloc(num_processes,sizeof(long int));
        int* current_neighbours_in_partition = (int*)malloc(num_processes*sizeof(int));

        // communication balance variables
        /*double* sum_all_bandwidth = (double*)calloc(num_processes,sizeof(double));
        for (int ii=0; ii < num_processes ; ii++) {
            for (int jj=0; jj < num_processes ; jj++) {
                sum_all_bandwidth[ii] += 1.0/comm_cost_matrix[ii][jj];
            }
        }*/

        // overfit variables
        bool check_overfit = false;
        idx_t* last_partitioning = NULL;
        float last_cut_metric;
        bool rollback = false;
        float last_imbalance = num_processes;
        //double timing = 0;
        //double ttt;
        memset(part_load,0,num_processes * sizeof(long int));
        double total_workload = 0;
        for(int ii=0; ii < num_vertices; ii++) {
            part_load[partitioning[ii]] += vtx_wgt[ii]; // workload for vertex
            total_workload += vtx_wgt[ii];
        }
        double expected_workload = total_workload / num_processes;
        for(int iter=0; iter < iterations; iter++) {
            //timing = 0;
            // go through own vertex list and reassign
            for(int vid=0; vid < num_vertices; vid++) {
                memset(current_neighbours_in_partition,0,num_processes * sizeof(int));
                //long int total_edges = 0;

                //int total_neighbours = 1;
                // where are neighbours located
                // new communication cost incurred
                
                // does not double count vertices that are present in multiple hyperedges
                // communication cost should be based on hedge cut?
                for(int he = 0; he < hedge_ptr[vid].size(); he++) {
                    int he_id = hedge_ptr[vid][he];
                    for(int vt = 0; vt < hyperedges[he_id].size(); vt++) {
                        int dest_vertex = hyperedges[he_id][vt];
                        if(dest_vertex == vid) continue;
                        int dest_part = partitioning[dest_vertex];
                        current_neighbours_in_partition[dest_part] += 1;
                        //total_neighbours++;
                        //total_edges++;
                    }
                }
                
                // allocate vertex
                double max_value = std::numeric_limits<double>::lowest();
                int best_partition = partitioning[vid];
                for(int pp=0; pp < num_processes; pp++) {
                    double connectivity_imbalance = 0;
                    double total_comm_cost = 0;
                    int neighbouring_partitions = 0;
                    //long int outer_edges = total_edges - current_neighbours_in_partition[pp];
                    for(int jj=0; jj < num_processes; jj++) {
                        if(pp != jj) {
                            // total cost of communication (edgecuts * comm cost * number of participating partitions)
                            total_comm_cost += current_neighbours_in_partition[jj] * comm_cost_matrix[pp][jj];
                            neighbouring_partitions += current_neighbours_in_partition[jj] > 0 ? 1 : 0;
                            // load starting connectivity
                            // inspired in Xue, et al 2015 (architecture aware streaming graph partitioning)
                            // not only keep workload balance, but p2p comms balanced (according to bandwidth)
                            // see page 6 of their paper for more info
                            // expected num edge between each i j partition = total edges * bandwidth(i,j) / sum ( all i,j bandwidths)
                            // edge imbalance of partition i = sum for each j ( expected edges between i,j - actual edges between i,j - i neighbours in j )
                            // the highest value means vid needs to be placed in i
                            // do not count edges to partition that will be local
                            //connectivity_imbalance += (outer_edges * (1.0/comm_cost_matrix[pp][jj]) / sum_all_bandwidth[pp] - current_neighbours_in_partition[jj]);
                        }
                        

                    }
                    
                    //double current_value =  -(double)neighbouring_partitions/(double)num_processes * total_comm_cost - a * (part_load[pp]/expected_workload) - 1.5*connectivity_imbalance;
                    double current_value =  -(double)neighbouring_partitions/(double)num_processes * total_comm_cost - a * (part_load[pp]/expected_workload);
                    //double current_value = current_neighbours_in_partition[pp]/(double)total_neighbours -(double)total_comm_cost / (double)num_processes * comm_cost_per_partition[pp] - a * (part_load[pp]/expected_workload);
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
                    } 
                }
                
                // update intermediate workload and assignment values
                part_load[best_partition] += vtx_wgt[vid];
                part_load[partitioning[vid]] -= vtx_wgt[vid];
                 
                // update partitioning assignment
                partitioning[vid] = best_partition;
                
                
            }
            
            // check if desired imbalance has been reached
            float imbalance = calculateImbalance(partitioning,num_processes,num_vertices,vtx_wgt);
            PRINTF("%i: %f (%f | %f)\n",iter,imbalance,a,ta_start);

            // get cut metric
            float hyperedges_cut_ratio;
            float edges_cut_ratio;
            int soed;
            float absorption;
            float max_imbalance;
            double total_edge_comm_cost;
            double total_hedge_comm_cost;
            float vertex_replicastion_factor;
            if(save_partitioning_history) {
                PRAW::getVertexCentricPartitionStatsFromFile(partitioning, num_processes, hypergraph_filename, NULL,comm_cost_matrix,
                            &hyperedges_cut_ratio, &edges_cut_ratio, &soed, &absorption, &max_imbalance, 
                            &total_edge_comm_cost,
                            stopping_condition == 3 ? &total_hedge_comm_cost : NULL,// only check hedge cost if it's going to be used
                            &vertex_replicastion_factor); 
                FILE *fp = fopen(history_file.c_str(), "ab+");
                if(fp == NULL) {
                    printf("Error when storing partitioning history into file\n");
                } else {
                    fprintf(fp,"%.3f,%.3f,%.3f,%i,%.3f,%.3f,%.3f\n",imbalance,hyperedges_cut_ratio,edges_cut_ratio,soed,absorption,total_edge_comm_cost,total_hedge_comm_cost);
                }
                fclose(fp);
            }
            
            
            // stop the process in the following conditions based on stopping_condition parameter
            // 0 = stop as soon as imbalance tolerance has been reached
            // > 0 -> 1 = hedge + edge cut metric, 2 = edge comm cost metric, 3 = hedge comm cost metric
                //  1. imbalance tolerance has been reached
                //      record current cut metric and partitioning and do one more iteration
                //      if imbalance is still ok 
                //          metric has not been improved, take recorded partitioning and stop
                //          metric has been improved, store partitioning and do one more iteration
            if(frozen_iters <= iter) {
                if(stopping_condition == 0) {
                    if (imbalance < imbalance_tolerance) {
                        break;
                    }
                } else {
                    if (imbalance < imbalance_tolerance) {
                        // get cut metric
                        if(!save_partitioning_history) {
                            // getting all metrics (way too expensive)
                            PRAW::getVertexCentricPartitionStatsFromFile(partitioning, num_processes,hypergraph_filename, NULL,comm_cost_matrix,
                                        &hyperedges_cut_ratio, &edges_cut_ratio, &soed, &absorption, &max_imbalance, 
                                        &total_edge_comm_cost,
                                        stopping_condition == 3 ? &total_hedge_comm_cost : NULL,// only check hedge cost if it's going to be used
                                        &vertex_replicastion_factor); 
                        }
                        double cut_metric;
                        if(stopping_condition == 1) cut_metric = hyperedges_cut_ratio + edges_cut_ratio;//hyperedges_cut_ratio + edges_cut_ratio;
                        if(stopping_condition == 2) cut_metric = ceil(total_edge_comm_cost);//hyperedges_cut_ratio;
                        if(stopping_condition == 3) cut_metric = ceil(total_hedge_comm_cost);//hyperedges_cut_ratio;

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
                
            }
            //if(frozen_iters <= iter && imbalance < imbalance_tolerance) break;

            if(false && imbalance > 10) {
                // if the assignment is too imbalanced, discard it and reset
                // avoids the issue of overloading most of the communication to and from one single process
                for (int vid=0; vid < num_vertices; vid++) {
                    partitioning[vid] = vid % num_processes;
                }
                memset(part_load,0,num_processes * sizeof(long int));
                for(int ii=0; ii < num_vertices; ii++) {
                    part_load[partitioning[ii]] += vtx_wgt[ii]; // workload for vertex
                }
            }
            // update parameters
            if(imbalance > imbalance_tolerance) {
                if(imbalance > ta_start) {
                    a *= imbalance;
                } else {
                    a *= ta_start;
                }
                
            } else {
                a *= ta_refine;
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
        free(current_neighbours_in_partition);
        //free(sum_all_bandwidth);

        return 1;
    }

    int ParallelHyperedgePartitioning(char* experiment_name, idx_t* partitioning, double** comm_cost_matrix, std::string hypergraph_filename, int* vtx_wgt, int iterations, float imbalance_tolerance, float ta_refine, bool reset_partitioning, int stopping_condition, bool save_partitioning_history) {
        
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
        // minimum number of iterations run (not checking imbalance threshold)
        // removed whilst we are using hyperPraw as refinement algorithm
        //      hence, if balanced is kept after first iteration, that's good enough
        int frozen_iters = 0;
        ///////////////
        
        // algorithm from GraSP (Battaglino 2016)
        // 1 - Distributed vertices over partitions (partition = vertex_id % num_partitions)
        // needs to load num_vertices from file
        for (int vid=0; vid < num_vertices; vid++) {
            partitioning[vid] = vid % num_processes;
        }
        //frozen_iters = ceil(0.1f * iterations);
        
        
        // 2 - Divide the graph in a distributed CSR format (like ParMETIS)
        //  compressed vertex or compressed hedge format? --> see zoltan
        //  for each local vertex, store the list of vertices adjacent to it (belonging to same hedges)
        std::vector<std::vector<int> > hyperedges;
        std::vector<std::vector<int> > hedge_ptr;
        load_hypergraph_from_file_dist_CSR(hypergraph_filename, &hyperedges, &hedge_ptr, process_id, partitioning);

        // Check balance guarantee 
        // The parallel algorithm is guaranteed to reach load imbalance tolerance if the hypergraphs safisfies:
        //      num_processes <= floor(sqrt(imbalance_tolerance * num_vertices - num_vertices))
        if(num_processes > floor(sqrt(imbalance_tolerance * num_vertices - num_vertices))) {
            int p = num_processes;
            int v = num_vertices;
            float i = imbalance_tolerance;
            int max_processes_for_guarantee = floor(sqrt(v * i - v));
            int min_hgraph_size = pow(p,2) / (i - 1);
            printf("WARNING: Current run is not guaranteed to reach load imbalance tolerance. Decrease the number of processes to %i.\nWith %i processes, %i vertices are required for guarantee\n",
                            max_processes_for_guarantee,p,min_hgraph_size);
        }        

        // each process must read from file only the info relevant to its data
        // 3 - Initiate N number of iterations on each process:
        //      a - one vertex at a time, assign to best partition (based on eval function)
        //      b - update tempering parameters
        //      c - share with all new partition assignments
        
        std::string history_file = experiment_name;
        
        if(save_partitioning_history) {
            history_file += "_";
            history_file += getFileName(hypergraph_filename);
            history_file += "_partition_history_";
            char str_int[16];
            sprintf(str_int,"%i",num_processes);
            history_file += "__";
            history_file +=  str_int;
            // remove history file if exists
            FILE *fp = fopen(history_file.c_str(), "w");
            if(fp == NULL) {
                printf("Error when storing partitioning history into file\n");
            } else {
                fprintf(fp,"%s\n","Imbalance, Hedges cut, Edges cut, SOED, Absorption, Edge sim comm cost, Hedge sim comm cost");
            }
            fclose(fp);
        }

        long int* part_load = (long int*)calloc(num_processes,sizeof(long int));
        int* current_neighbours_in_partition = (int*)malloc(num_processes*sizeof(int));
        // overfit variables
        bool check_overfit = false;
        idx_t* last_partitioning = NULL;
        float last_cut_metric;
        bool rollback = false;
        float last_imbalance = num_processes;
        
        memset(part_load,0,num_processes * sizeof(long int));
        double total_workload = 0;
        for(int ii=0; ii < num_vertices; ii++) {
            part_load[partitioning[ii]] += vtx_wgt[ii]; // workload for vertex
            total_workload += vtx_wgt[ii];
            
        }
        double expected_workload = total_workload / num_processes;
        /*double maxload = part_load[0];
        double minload = part_load[0];
        for(int ll=0; ll < num_processes; ll++) {
            if(part_load[ll] > maxload) maxload = part_load[ll];
            if(part_load[ll] < minload) minload = part_load[ll];
        }*/
        
        int iter = 0;
        for(iter=0; iter < iterations; iter++) {
            
            int best_partition = 0;
            int last_partition_update = 0;

            // go through own vertex list and reassign
            for(int vid=0; vid < num_vertices; vid++) {
                
                bool isLocal = vid % num_processes == process_id; //hedge_ptr[vid].size() > 0; // always go through the same list of vertices per process
                // if local vertex, calculate full heuristic (cost of communication...)
                // if non local vertex, speculatively place it based on current partitioning load balance
                // this alleviates the problems of parallel streams maintaining workload balance when 
                // alpha parameter is high and partition load update is sparse
                if(isLocal) {
                    memset(current_neighbours_in_partition,0,num_processes * sizeof(int));

                    // reevaluate objective function per partition
                    // |P^t_i union N(v)| = number of vertices in partition i that are neighbours of vertex v 
                    // where are neighbours located
                    // new communication cost incurred
                    
                    // does double count vertices that are present in multiple hyperedges
                    for(int he = 0; he < hedge_ptr[vid].size(); he++) {
                        int he_id = hedge_ptr[vid][he];
                        int cardinality = hyperedges[he_id].size();
                        for(int vt = 0; vt < cardinality; vt++) {
                            /*int dest_vertex = hyperedges[he_id][vt];
                            if(dest_vertex == vid) continue;
                            int dest_part = partitioning[dest_vertex];
                            current_neighbours_in_partition[dest_part] += 1;*/
                            if(hyperedges[he_id][vt] == vid) continue;
                            current_neighbours_in_partition[partitioning[hyperedges[he_id][vt]]] += 1;
                        }
                    }
                
                    // allocate vertex (for local heuristically, for non local speculatively)
                    double max_value = std::numeric_limits<double>::lowest();
                    best_partition = partitioning[vid];
                    for(int pp=0; pp < num_processes; pp++) {
                        // total cost of communication (edgecuts * number of participating partitions)
                        long int total_comm_cost = 0;
                        int neighbouring_partitions = 0;
                        
                        for(int jj=0; jj < num_processes; jj++) {
                            if(pp != jj) {
                                total_comm_cost += current_neighbours_in_partition[jj] * comm_cost_matrix[pp][jj];
                                neighbouring_partitions += current_neighbours_in_partition[jj] > 0 ? 1 : 0;
                            }
                        } 

                        double current_value =  -(double)neighbouring_partitions/(double)num_processes * total_comm_cost - a * (part_load[pp]/expected_workload);
                        //double current_value =  -(double)neighbouring_partitions/(double)num_processes * total_comm_cost + a * (maxload - part_load[pp]) / (maxload - minload);
                        
                        if(current_value > max_value) {
                            max_value = current_value;
                            best_partition = pp;
                        }
                    }
                }

                last_partition_update++;
                // share allocation decisions amongst processes
                // Do it once each processor has had a chance to place one vertex (every num_processes)
                if((vid+1) % num_processes == 0 || vid == num_vertices-1) {
                    // each process sends the partition allocation for its local vertex
                    // each process receives the partition allocation for all other vertices
                    int* v_mappings = (int*)malloc(num_processes * sizeof(int));
                    MPI_Allgather(&best_partition,1,MPI_INT,v_mappings,1,MPI_INT,MPI_COMM_WORLD);
                    // update local part_load
                    for(int ii=0; ii < last_partition_update; ii++) {
                        int id = vid + 1 - last_partition_update + ii;
                        int dest_partition = v_mappings[ii];
                        if(partitioning[id] != dest_partition) {
                            part_load[partitioning[id]] -= vtx_wgt[id];
                            part_load[dest_partition] += vtx_wgt[id];
                            //if(part_load[dest_partition] > maxload) maxload = part_load[dest_partition];
                            //if(part_load[partitioning[id]] < minload) minload = part_load[partitioning[id]];
                            partitioning[id] = dest_partition;
                        }
                    }
                    free(v_mappings);
                    last_partition_update = 0;
                }
                
                
            }

            // check if desired imbalance has been reached
            float imbalance = calculateImbalance(partitioning,num_processes,num_vertices,vtx_wgt);
            PRINTF("%i: %f (%f | %f)\n",iter,imbalance,a,ta_start);

        
            // get cut metric
            float hyperedges_cut_ratio;
            float edges_cut_ratio;
            int soed;
            float absorption;
            float max_imbalance;
            double total_edge_comm_cost;
            double total_hedge_comm_cost;
            if(save_partitioning_history) {
                PRAW::getVertexCentricPartitionStatsFromFile_parallel(partitioning, num_processes,process_id,hypergraph_filename, NULL,comm_cost_matrix,
                            &hyperedges_cut_ratio, &edges_cut_ratio, &soed, &absorption, &max_imbalance, 
                            &total_edge_comm_cost,
                            stopping_condition == 3 ? &total_hedge_comm_cost : NULL); // only check hedge cost if it's going to be used
                if(process_id == MASTER_NODE) {
                    FILE *fp = fopen(history_file.c_str(), "ab+");
                    if(fp == NULL) {
                        printf("Error when storing partitioning history into file\n");
                    } else {
                        fprintf(fp,"%.3f,%.3f,%.3f,%i,%.3f,%.3f,%.3f\n",imbalance,hyperedges_cut_ratio,edges_cut_ratio,soed,absorption,total_edge_comm_cost,total_hedge_comm_cost);
                    }
                    fclose(fp);
                }
            }

            // stop the process in the following conditions
            //  1. imbalance tolerance has been reached
            //      record current cut metric and partitioning and do one more iteration
            //      if imbalance is still ok 
            //          metric has not been improved, take recorded partitioning and stop
            //          metric has been improved, store partitioning and do one more iteration
            // ALL PROCESS MUST STOP to check if 0 has broken out of the loop
            if(frozen_iters <= iter) {
                if(stopping_condition == 0) {
                    if (imbalance < imbalance_tolerance) {
                        break;
                    }
                } else {
                    if (imbalance < imbalance_tolerance) {
                        // checking the quality of partitioning in parallel is faster than leaving it to the master node alone
                        if(!save_partitioning_history) {
                            PRAW::getVertexCentricEdgeCommCostFromFile(partitioning, num_processes, process_id, hypergraph_filename,comm_cost_matrix,
                                                                   &total_edge_comm_cost);
                            /*PRAW::getVertexCentricPartitionStatsFromFile_parallel(partitioning, num_processes,process_id,hypergraph_filename, NULL,comm_cost_matrix,
                                                    &hyperedges_cut_ratio, &edges_cut_ratio, &soed, &absorption, &max_imbalance, 
                                                    &total_edge_comm_cost,
                                                    stopping_condition == 3 ? &total_hedge_comm_cost : NULL); // only check hedge cost if it's going to be used
                            */
                        }
                        if(process_id == MASTER_NODE) {
                            double cut_metric = ceil(total_edge_comm_cost);//hyperedges_cut_ratio;

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
                        check_overfit = false;
                    }  
                }
                
            }

            // update parameters
            if(imbalance > imbalance_tolerance) {
                if(imbalance > ta_start) {
                    a *= imbalance;
                } else {
                    a *= ta_start;
                }
                
            } else {
                a *= ta_refine;
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
        free(part_load);
        free(current_neighbours_in_partition);

        return iter+1;
    }

    // Stream from multiple files / streams
    int ParallelHyperedgePartitioning_he_stream(char* experiment_name, idx_t* partitioning, double** comm_cost_matrix, char* hypergraph_filename, int num_vertices, std::vector<std::vector<int> >* he_stream, int* vtx_wgt, int max_iterations, float imbalance_tolerance, bool save_partitioning_history) {
        
        ////////////////////////
        // Parallel streaming algorithm based on streaming proposed by Alistairh (minmax hypergraph streaming)
        // Single pass, a decision is made looking only at the info for one vertex (vertex id plus all hyperedges it belongs to)
        // Necessary to store, for each partition, the list of hyperedges currently present in them
        // Novelty: 
        //      Parallel streams
        //      Inclusion of communication costs
        //      Inclusion of load balance parameter
        ////////////////////////


        // Parameters (from HDRF, Petroni 2015)
        float lambda = 0.0f;
        // own parameters
        float lambda_update = 0.1f;

        int process_id;
        MPI_Comm_rank(MPI_COMM_WORLD,&process_id);
        int num_processes;
        MPI_Comm_size(MPI_COMM_WORLD,&num_processes); 

        // create shared data structures (partitions workload, list of replica destinations for each vertex, partial degree for each vertex)
        long int* part_load = (long int*)calloc(num_processes, sizeof(long int));
        float max_comm_cost = 0;
        for(int ff=0; ff < num_processes; ff++) {
            float current_max_cost = 0;
            for(int tt=0; tt < num_processes; tt++) {
                current_max_cost += comm_cost_matrix[ff][tt];
            }
            if(current_max_cost > max_comm_cost) {
                max_comm_cost = current_max_cost;
            }
        }

        std::string history_file = experiment_name;
        
        if(save_partitioning_history && process_id == 0 && hypergraph_filename != NULL) {
            history_file += "_partition_history__";
            char str_int[16];
            sprintf(str_int,"%i",num_processes);
            history_file +=  str_int;
            // remove history file if exists
            FILE *fp = fopen(history_file.c_str(), "w");
            if(fp == NULL) {
                printf("Error when storing partitioning history into file\n");
            } else {
                fprintf(fp,"%s\n","Imbalance, Hedges cut, Edges cut, SOED, Absorption, Edge sim comm cost, Hedge sim comm cost");
            }
            fclose(fp);
        }

        int iter =0;
        for(iter=0; iter < max_iterations; iter++) {

            //int num_vertices = he_stream->size();
            //int num_hyperedges = tokens[0];
            
            PRINTF("Found in stream: Vertices: %i\n",num_vertices);

            // initialise workload
            // workload starts empty across all partitions
            long int total_workload = 0;
            memset(part_load,0,num_processes * sizeof(long int));
            
            long int maxsize = part_load[0]; // used for c_bal
            long int minsize = part_load[0]; // used for c_bal

            // Check balance guarantee 
            // The parallel algorithm is guaranteed to reach load imbalance tolerance if the hypergraphs safisfies:
             //      num_processes < floor(sqrt(num_vertices * imbalance_tolerance - num_vertices))
            if(num_processes >= floor(sqrt(num_vertices * imbalance_tolerance - num_vertices))) {
                int p = num_processes;
                float v = num_vertices;
                float i = imbalance_tolerance;
                int max_processes_for_guarantee = floor(sqrt(v * i - v));
                int min_graph_size = pow(p,2) / (i - 1);
                printf("WARNING: Current run is not guaranteed to reach load imbalance tolerance. Decrease the number of processes to %i.\nWith %i processes, %i vertices are required for guarantee\n",
                                max_processes_for_guarantee,p,min_graph_size);
            }      

            // HOW DO WE STORE AND COORDINATE THE DATASTRUCTURES? Balance between memory and communication
            // store for each hyperedge, store list of partitions
            std::unordered_map<int,std::set<int> > seen_hyperedges;

            // read reminder of file (one line per hyperedge)
            int vid = 0;
            int current_stream = 0; // current local position on local stream
            int v_mapping = -1; // mapping of current he (to partition)
            while(vid < num_vertices) {
                if(current_stream < he_stream->size()) {
                    // prevents out of bounds in last round of streaming for last processes
                    // process hyperedges to which vertex_id belongs and assign to partition
                    float max_value = 0;
                    int best_partition = 0;
                    int num_hedges = he_stream->at(current_stream).size();
                    for(int pp=0; pp < num_processes; pp++) {
                        float c_overlap = 0;
                        float c_comm = 0;
                        int top_connections;
                        for(int hid=0; hid < num_hedges; hid++) {
                            int he_id = he_stream->at(current_stream)[hid];
                            top_connections += seen_hyperedges[he_id].size();
                            std::set<int>::iterator it;
                            for(it = seen_hyperedges[he_id].begin(); it != seen_hyperedges[he_id].end(); it++) {
                                int dest = *it;
                                // calculate communication cost
                                c_comm += comm_cost_matrix[pp][dest];
                                // Calculate overlap of heges on current considered partition
                                if(dest == pp) c_overlap += 1;
                            }
                        }
                        c_overlap /= num_hedges;
                        c_comm = 1 - c_comm / (max_comm_cost * num_hedges);
                        
                        // calculate load balance component
                        float c_bal = lambda * (maxsize - part_load[pp]) / (0.1f + maxsize - minsize);
                        

                        // assign to partition that maximises C_rep + C_bal + C_comm
                        float current_value = c_overlap + c_comm + c_bal;
                        if(current_value > max_value) {
                            max_value = current_value;
                            best_partition = pp;
                            //printf("%f, %f, %f\n",c_overlap,c_bal,c_comm);
                        }

                    }
                    v_mapping = best_partition;
                }
                
                
                // synchronise data amongst parallel streams
                std::vector<int> new_he_presence;
                new_he_presence.push_back(v_mapping); // add in front the partition selected for vertex_id
                for(int ii=0; ii < he_stream->at(current_stream).size(); ii++) {
                    int he_id = he_stream->at(current_stream)[ii];
                    // only send it over to other processes if it has not been recorded before
                    if(seen_hyperedges[he_id].find(v_mapping) == seen_hyperedges[he_id].end()) {
                        // if the vertex has not been seen before
                        new_he_presence.push_back(he_id);
                    }
                }

                current_stream++;
                vid += num_processes; // since we are doing parallel streams, one per process

                // share send buffer size with other processes
                // size = number of new replicas + 1 (the partition selected)
                int* recvcounts = (int*)malloc(num_processes * sizeof(int));
                int send_size = new_he_presence.size();
                MPI_Allgather(&send_size,1,MPI_INT,recvcounts,1,MPI_INT,MPI_COMM_WORLD);

                // share partition selected and new replicas list to all
                int counter = 0;
                int* displs = (int*)malloc(sizeof(int) * num_processes);
                for(int ii=0; ii < num_processes; ii++) {
                    displs[ii] = counter;
                    counter += recvcounts[ii]; 
                }
                int* recvbuffer = (int*)malloc(sizeof(int) * counter);
                MPI_Allgatherv(&new_he_presence[0],send_size,MPI_INT,recvbuffer,recvcounts,displs,MPI_INT,MPI_COMM_WORLD);

                // process new he presence and add them to local datastructures
                for(int ii=0; ii < num_processes; ii++)  {
                    int dest_partition = recvbuffer[displs[ii]];
                    int vertex_id = vid - num_processes + ii;
                    if(vertex_id >= num_vertices) break;
                    part_load[dest_partition] += vtx_wgt[vertex_id];
                    total_workload += vtx_wgt[vertex_id];
                    partitioning[vertex_id] = dest_partition;
                    for(int jj=1; jj < recvcounts[ii]; jj++) {
                        int he_id = recvbuffer[displs[ii] + jj];
                        // update hedge info             
                        seen_hyperedges[he_id].insert(dest_partition);                
                    }
                    
                }

                free(recvcounts);
                free(displs);
                free(recvbuffer);

                // update workload limits
                long int max = part_load[0];
                long int min = part_load[0];
                for(int pp=0; pp < num_processes; pp++) {
                    if(part_load[pp] > max) max = part_load[pp];
                    if(part_load[pp] < min) min = part_load[pp];
                }
                minsize = min;
                maxsize = max;
                
            }

            // check for termination condition (tolerance imbalance reached)   
            float max_imbalance = ((float)maxsize) / ((float)total_workload/num_processes);
            PRINTF("***Hedge imbalance: %.3f\n",max_imbalance);

            if(save_partitioning_history && hypergraph_filename != NULL) {
                float hyperedges_cut_ratio;
                float edges_cut_ratio;
                int soed;
                float absorption;
                float max_imbalance;
                double total_edge_comm_cost;
                double total_hedge_comm_cost;
                // store partition history
                PRAW::getVertexCentricPartitionStatsFromFile_parallel(partitioning, num_processes,process_id,hypergraph_filename, vtx_wgt,comm_cost_matrix,
                            &hyperedges_cut_ratio, &edges_cut_ratio, &soed, &absorption, &max_imbalance, 
                            &total_edge_comm_cost,NULL);
                
                if(process_id == MASTER_NODE) {
                    FILE *fp = fopen(history_file.c_str(), "ab+");
                    if(fp == NULL) {
                        printf("Error when storing partitioning history into file\n");
                    } else {
                        fprintf(fp,"%.3f,%.3f,%.3f,%i,%.3f,%.3f,%.3f\n",max_imbalance,hyperedges_cut_ratio,edges_cut_ratio,soed,absorption,total_edge_comm_cost,total_hedge_comm_cost);
                    }
                    fclose(fp);
                }
            }

            if(max_imbalance <= imbalance_tolerance) break;

            // update lambda (importance of load balancing)
            lambda += lambda_update;
        }

        // clean up
        free(part_load);

        return iter+1;

    }


    // Stream from multiple files / streams
    int ParallelHDRF(char* experiment_name, idx_t* partitioning, double** comm_cost_matrix, std::string hypergraph_filename, int* element_wgt, int max_iterations, float imbalance_tolerance, bool save_partitioning_history, bool local_replica_degree_updates_only = false, int sync_batch_size = 1, bool use_max_expected_workload = true) {
        // Parallel Hyperedge Partitioning based algorithm
        // Because it can be applied to both vertex and hyperedge partitionings, we adopt the following nomenclature:
        //      element: what each line in the stream represent
        //          for vertex partitioning, each line is a hyperedge
        //          for hedge partitioning, each line is a vertex
        //      pin: each of the elements presented at a time by the stream
        //          for vertex partitioning, they are vertices that the hyperedge contain
        //          for hedge partitioning, they are hyperedges the vertex belongs to

        // The goal is to assign hyperedges to partitions   
        // Minimisation goal: 
        //      replication factor (sum of vertex replicas)
        //      reducing communication cost (taking comm cost matrix into account)
        // Constraint: keep load balance amongst partitions (weighed sum of hyperedges per partition = sum of vertices per partition?)
        // It can be a restreaming approach, in which in the first stream we use partial info (partial degree) and full info afterwards
        //      particularly when incorporating comm cost into the objective function

        // Input: 
        //      hypergraph, in a stream (one hyperedge at a time)
        //      number of partitions
        //      communication cost matrix
        // Output: hyperedge assignment to partitions

        // Algorithm employed: a variety of HDRF (High Degree Replicated First, Petroni). 
        // standalone implementation of HDRF https://github.com/fabiopetroni/VGP
        // It is justified if the hgraphs are highly skewed (few vertices with high degree) power law graphs
        // we can use lambda to modify the importance of workload (adaptive?) or use another parameter for the communication load part
        // Required information per process
        //      Partial vertex degree (number of times a vertex has appeared in the stream)
        //      List of neighbouring partitions per local vertex --> list of partitions that contain a replica of the vertex
        //      Partitions workload assignment (max and min)

        // Implementation questios
        //      How does the algorithm keep A(v)? It builds it up as it goes along, then shares it amongs parallel processes
        //      How often are A(v) and workload updated across partitions?
        //      Is the partial vertex degree shared? Does not seem so, but would it improve the results?

        // Pseudocode
        // 1 Decide which hyperedges will be local to each process (he_id % num_processes == process_id)
        // 2 Open stream and read one by one only local hyperedges
        // 3 For each he and each vertex in it, fetch and update
        //      A(v)
        //      Partial degree of v
        // 4 Calculate cost function for each partition
        //      Fetch workload for partition
        // 5 Assign he to best partition
        //      Update workload for partition
        // 6 Synchronise data amongst processes (how often? variable to optimise and experiment with)
        //      A(v)
        //      workload (experiment with and without)
        //      Partial degree? (experiment with and without)

        // Parameters (from HDRF, Petroni 2015)
        float lambda = use_max_expected_workload ? 0.1f : 1.0f;
        // own parameters
        float lambda_update = 1.1f;
        float lambda_refinement = 0.95f;
            
        
        int process_id;
        MPI_Comm_rank(MPI_COMM_WORLD,&process_id);
        int num_processes;
        MPI_Comm_size(MPI_COMM_WORLD,&num_processes); 

        // create shared data structures (partitions workload, list of replica destinations for each vertex, partial degree for each vertex)
        long int* part_load = (long int*)calloc(num_processes, sizeof(long int));
        /*double max_comm_cost = 0;
        for(int ff=0; ff < num_processes; ff++) {
            double current_max_cost = 0;
            for(int tt=0; tt < num_processes; tt++) {
                current_max_cost += comm_cost_matrix[ff][tt];
            }
            if(current_max_cost > max_comm_cost) {
                max_comm_cost = current_max_cost;
            }
        }*/

        std::string history_file = experiment_name;
        
        if(save_partitioning_history && process_id == 0) {
            history_file += "_";
            history_file += getFileName(hypergraph_filename);
            history_file += "_partition_history__";
            char str_int[16];
            sprintf(str_int,"%i",num_processes);
            history_file +=  str_int;
            // remove history file if exists
            FILE *fp = fopen(history_file.c_str(), "w");
            if(fp == NULL) {
                printf("Error when storing partitioning history into file\n");
            } else {
                fprintf(fp,"%s\n","Iteration, Lambda, Vertex replication factor, Hedge imbalance");
            }
            fclose(fp);
        }

        int num_pins;
        int num_elements;

        // avoid overfitting variables
        bool check_overfit = false;
        bool rollback = false;
        idx_t* last_partitioning = NULL;

        int iter = 0;
        for(iter=0; iter < max_iterations; iter++) {

            // Open stream
            std::ifstream istream(hypergraph_filename.c_str());
            
            if(!istream) {
                printf("Error while opening hMETIS file %s\n",hypergraph_filename.c_str());
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
            num_pins = tokens[1];
            num_elements = tokens[0];

            // used to force all solutions to be within imbalance tolerance
            // currently uses full graph knowledge (just number of elements)
            // assumes all elements have same workload (1)
            long int max_expected_workload = num_elements / num_processes * imbalance_tolerance - num_processes;
            if(use_max_expected_workload && max_expected_workload < 1) {
                PRINTF("Graph is too small! Too many processes (max expected workload limit is %li)\n",max_expected_workload);
                return 0;
            }
       
            
            PRINTF("Found in file: Pins: %i; elements %i:\n",num_pins,num_elements);

            // initialise workload
            // workload starts empty across all partitions
            long int total_workload = 0;
            memset(part_load,0,num_processes * sizeof(long int));
            
            long int maxsize = part_load[0]; // used for c_bal
            long int minsize = part_load[0]; // used for c_bal

            // Check balance guarantee 
            // The parallel algorithm is guaranteed to reach load imbalance tolerance if the hypergraphs safisfies:
            //      total_workload * imbalance_tolerance > total_workload +  average_cardinality * p
            //      that translates to:
            //      num_hyperedges * (average_cardinality / num_processes) * imbalance_tolerance > num_hyperedges * average_cardinality / num_processes + average_cardinality * num_processes
            //      which simplifies to, when solved for num_processes to: (average_cardinality is cancelled out)
            //      num_processes < floor(sqrt(num_hyperedges * imbalance_tolerance - num_hyperedges))
            if(num_processes >= floor(sqrt(num_elements * imbalance_tolerance - num_elements))) {
                int p = num_processes;
                float h = num_elements;
                float i = imbalance_tolerance;
                int max_processes_for_guarantee = floor(sqrt(h * i - h));
                int min_graph_size = pow(p,2) / (i - 1);
                printf("WARNING: Current run is not guaranteed to reach load imbalance tolerance. Decrease the number of processes to %i.\nWith %i processes, %i elements are required for guarantee\n",
                                max_processes_for_guarantee,p,min_graph_size);
            }      

            // HOW DO WE STORE AND COORDINATE THE DATASTRUCTURES? Balance between memory and communication
            // create an object Vertex that contains two variables
            //      int partial_degree --> store the current partial degree
            //      std::vector<int> A --> store the list of partitions that have a replica of the vertex
            // store the Vertex in a std::unordered_map (hashmap) of Vertex* of length num_vertices (vertex id is given by the index)
            std::unordered_map<int,pin_data> seen_pins;

            // read reminder of file (one line per elemennt)
            int element_id = 0; // current elemennt
            int element_mapping = -1; // mapping of current he (to partition)
            std::vector<int> local_pins; // list of hyperedges and vertices seen in the current batch that belong to process

            // data structures for batch synchronisation
            std::vector<int> local_pins_size;
            std::vector<int> new_replicas;

            while(element_id < num_elements) {
                if(std::getline(istream,line)) {
                    char str[line.length() + 1]; 
                    strcpy(str, line.c_str()); 
                    char* token = strtok(str, " "); 
                    local_pins.clear();
                    while (token != NULL) { 
                        local_pins.push_back(atoi(token)-1); 
                        token = strtok(NULL, " "); 
                    } 
                    
                    // local hyperedge, process and assign it
                    // calculate norm_part_degree for each vertex
                    //std::vector<double> normalised_part_degrees(vertices.size());
                    double normalised_part_degrees[local_pins.size()];
                    long int total_degrees = 0;
                    for(int ii=0; ii < local_pins.size(); ii++) {
                        int pin_id = local_pins[ii];
                        normalised_part_degrees[ii] = std::max(seen_pins[pin_id].partial_degree,1); // if vertex is newly seen, it will be counted in the next sync. But count it here too
                        total_degrees += normalised_part_degrees[ii];
                    }
                    std::transform(normalised_part_degrees,normalised_part_degrees+local_pins.size(),normalised_part_degrees,
                        [total_degrees] (double value) {  
                            return value / total_degrees;
                        }  
                    );

                    // calculate C_rep(he) per partition per vertex
                    //      sum 1 + (1-norm_part_degree(v)) if p exists in A(v)
                    // calculate C_bal(he) per partition
                    //      lambda * (maxsize - |p|) / (e + maxsize - minsize)
                    //      lambda --> > 1
                    //      e --> required to avoid dividing by 0 if maxsize == minsize
                    // TODO: can we avoid having to do two passes across all processes?
                    // can we avoid having these two datastructures?
                    float* c_total = (float*)calloc(num_processes,sizeof(float));
                    double* c_comms = (double*)calloc(num_processes,sizeof(double));   
                    float comm_min = std::numeric_limits<float>::max();
                    float comm_max = 0;
                    for(int pp=0; pp < num_processes; pp++) {
                        if(use_max_expected_workload && part_load[pp] >= max_expected_workload) {
                            continue;
                        }
                        double c_rep = 0;
                        double c_comm = 0;
                        for(int vv=0; vv < local_pins.size(); vv++) {
                            int pin_id = local_pins[vv];
                            bool present_in_partition = false;
                            std::set<int>::iterator it;
                            for (it = seen_pins[pin_id].A.begin(); it != seen_pins[pin_id].A.end(); ++it)
                            {
                                int part = *it;
                                present_in_partition |= part == pp;
                                // communication should be proportional to the duplication of pins
                                // if a pin is duplicated in two partitions, then communication will happen across those partitions
                                c_comm += comm_cost_matrix[pp][part];
                            }
                            c_rep += present_in_partition ? 1 + (1 - normalised_part_degrees[vv]) : 0;
                            
                        }
                        // c_rep and c_comm must be normalised to avoid them dominating the final equation
                        c_comms[pp] = c_comm;
                        if(c_comm > comm_max) comm_max = c_comm;
                        if(c_comm < comm_min) comm_min = c_comm;

                        c_rep = c_rep/(local_pins.size()+1);

                        float c_bal = lambda * (maxsize - part_load[pp]) / (0.1 + maxsize - minsize);

                        c_total[pp] = c_bal+ c_rep;

                    }

                    float max_value = 0;
                    int best_partition = 0;
                    for(int pp=0; pp < num_processes; pp++) {   
                        if(use_max_expected_workload && part_load[pp] >= max_expected_workload) {
                            continue;
                        } 
                        // normalise c_comms
                        c_comms[pp] = (comm_max - c_comms[pp]) / (0.1 + comm_max - comm_min);
                        
                        // assign to partition that maximises C_rep + C_bal + C_comm
                        float current_value = c_total[pp] + c_comms[pp];
                        if(current_value > max_value) {
                            max_value = current_value;
                            best_partition = pp;
                            //printf("%f, %f, %f\n",c_rep,c_bal,c_comm);
                        }
                    }
                    element_mapping = best_partition;

                    free(c_comms);
                    free(c_total);
                }

                element_id += num_processes;

                // synchronise data
                // ***** 
                //  TODO: test just synchronising partition sizes, not vertex degree
                //  assumption: if pins are sufficiently shuffled, partial local degree may be enough and saves data shared
                // *****
                // ***** 
                //  TODO: Batch synchronisation requres testing. Are we saving time? Are results correct?
                //  issues: tradeoff between less comm overhead and quality of partition (potentially higher graph size requirements as partition load is not updated often)
                // *****
                new_replicas.push_back(element_mapping); // add in front the partition selected
                
                int new_pins = 0;
                for(int ii=0; ii < local_pins.size(); ii++) {
                    int pin_id = local_pins[ii];
                    if(!local_replica_degree_updates_only) {
                        new_replicas.push_back(pin_id);
                        new_pins++;
                    } else {
                        // TODO: test performance degradation when only updating partial degree with local info
                        if(seen_pins[pin_id].partial_degree == 0) {
                            // if the pin has not been seen before
                            new_replicas.push_back(pin_id);
                            new_pins++;
                        } else if(seen_pins[pin_id].A.find(element_mapping) == seen_pins[pin_id].A.end()) {
                            // if it has been seen but it's the first replica on new partition
                            new_replicas.push_back(pin_id);
                            new_pins++;
                        } else {
                            seen_pins[pin_id].partial_degree += 1;
                        }
                    }
                }
                local_pins_size.push_back(new_pins);

                // batch synchronisation
                if((element_id / num_processes) % sync_batch_size == 0) {
                    // synchronise length of pins list to be sent
                    int total_pins_size = sync_batch_size * num_processes;
                    int* remote_pins_size = (int*)malloc(sizeof(int) * total_pins_size);

                    MPI_Allgather(&local_pins_size[0],sync_batch_size,MPI_INT,remote_pins_size,sync_batch_size,MPI_INT,MPI_COMM_WORLD);
                    
                    // synchronise list of pins
                    // share send buffer size with other processes
                    // size = number of new replicas + 1 (the partition selected)
                    int* recvcounts = (int*)malloc(num_processes * sizeof(int));
                    int send_size = new_replicas.size();
                    
                    MPI_Allgather(&send_size,1,MPI_INT,recvcounts,1,MPI_INT,MPI_COMM_WORLD);
                    
                    // share partition selected and new replicas list to all
                    int counter = 0;
                    int* displs = (int*)malloc(sizeof(int) * num_processes);
                    for(int ii=0; ii < num_processes; ii++) {
                        displs[ii] = counter;
                        counter += recvcounts[ii]; 
                    }
                    int* recvbuffer = (int*)malloc(sizeof(int) * counter);
                    MPI_Allgatherv(&new_replicas[0],send_size,MPI_INT,recvbuffer,recvcounts,displs,MPI_INT,MPI_COMM_WORLD);
                    
                    // process new replicas and add them to local datastructures
                    for(int ii=0; ii < num_processes; ii++)  {
                        int current_element_index = 0;
                        for(int el_order=0; el_order < sync_batch_size; el_order++) {
                            int current_pin_length = remote_pins_size[ii*sync_batch_size+el_order];
                            int dest_partition = recvbuffer[displs[ii] + current_element_index];
                            int current_element = element_id - (num_processes * (sync_batch_size - el_order)) + ii;
                            if(current_element >= num_elements) break;
                            part_load[dest_partition] += element_wgt[current_element];
                            total_workload += element_wgt[current_element];
                            partitioning[current_element] = dest_partition;
                            for(int jj=0; jj < current_pin_length; jj++) {
                                int pin_id = recvbuffer[displs[ii] + current_element_index + 1 + jj];
                                // update vertex info
                                seen_pins[pin_id].partial_degree += 1;
                                seen_pins[pin_id].A.insert(dest_partition);                  
                            }
                            current_element_index += current_pin_length+1;
                        }

                    }
                    free(recvcounts);
                    free(displs);
                    free(recvbuffer);
                    free(remote_pins_size);

                    // clear intermediate data structures
                    new_replicas.clear();
                    local_pins_size.clear();  

                    // update workload limits
                    /*long int max = part_load[0];
                    long int min = part_load[0];
                    for(int pp=0; pp < num_processes; pp++) {
                        if(part_load[pp] > max) max = part_load[pp];
                        if(part_load[pp] < min) min = part_load[pp];
                    }
                    minsize = min;
                    maxsize = max;*/
                    minsize = *std::min_element(part_load, part_load + num_processes);
                    maxsize = *std::max_element(part_load, part_load + num_processes);
                }

            }
            istream.close();

            // check for termination condition (tolerance imbalance reached)   
            float max_imbalance = ((float)maxsize) / ((float)total_workload/num_processes);
            PRINTF("***Imbalance: %.3f, current lambda: %f\n",max_imbalance,lambda);                            

            if(save_partitioning_history && process_id == MASTER_NODE) {
                // store partition history
                float pin_replication_factor;
                PRAW::getEdgeCentricReplicationFactor(&seen_pins,num_pins,
                                        &pin_replication_factor);

                FILE *fp = fopen(history_file.c_str(), "ab+");
                if(fp == NULL) {
                    printf("Error when storing partitioning history into file\n");
                } else {
                    fprintf(fp,"%i,%.2f, %.3f,%.3f\n",iter,lambda,pin_replication_factor,max_imbalance);
                }
                fclose(fp);
            }

            //if(max_imbalance <= imbalance_tolerance) break;

            // update lambda (importance of load balancing)
            // keep searching until the algorithm has dipped into tolerance imbalance and gone out once
            // then return last partitioning that was within tolerance imbalance
            if(max_imbalance > imbalance_tolerance) {
                if(check_overfit) {
                    // exit and return previous partitioning
                    rollback = true;
                    break;
                }
                lambda *= lambda_update;
            } else {
                if(use_max_expected_workload) break;
                lambda *= lambda_refinement;
                check_overfit = true;
                if(process_id == MASTER_NODE) {
                    if(last_partitioning == NULL) {
                        last_partitioning = (idx_t*)malloc(num_elements*sizeof(idx_t));
                    }
                    memcpy(last_partitioning,partitioning,num_elements * sizeof(idx_t));
                }
            }
        }

        // roll back to last partitioning that was inside imbalance tolerance
        if(rollback) {
            if(process_id == MASTER_NODE) {
                // share last partitioning with all
                memcpy(partitioning,last_partitioning,num_elements * sizeof(idx_t));
                free(last_partitioning);
                for(int dest=0; dest < num_processes; dest++) {
                    if(dest == MASTER_NODE) continue;
                    MPI_Send(partitioning,num_elements,MPI_LONG,dest,0,MPI_COMM_WORLD);
                }
            } else {
                // update partitioning from 0
                MPI_Recv(partitioning,num_elements,MPI_LONG,MASTER_NODE,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }

        // clean up
        free(part_load);

        return iter+1;

    }
    

    
}

#endif