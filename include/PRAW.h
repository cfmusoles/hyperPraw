// Parallel Restreaming Architecture aWare partitioning algorithm

#ifndef PRAW__H
#define PRAW__H

#include <vector>
#include <deque>
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
        std::unordered_map<unsigned short,unsigned short> P;
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
            //PRINTF("%i: %li\n",ii,workload[ii]);
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
                            float* hyperedges_cut_ratio, float* edges_cut_ratio, int* soed, float* absorption, float* max_imbalance, double* total_edge_comm_cost, double* total_hedge_comm_cost,double* vertex_replication_factor) { // output

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
                            float* hyperedges_cut_ratio, float* edges_cut_ratio, int* soed, float* absorption, float* max_imbalance, double* total_edge_comm_cost, double* total_hedge_comm_cost,double* vertex_replication_factor) { // output
        // take partitioning statistics directly from reading a graph file
        // loading the entire graph on a process should be avoided for scalability
        *hyperedges_cut_ratio=0;
        *edges_cut_ratio=0;
        *soed=0;
        *absorption=0;
        *max_imbalance=0;
        *total_edge_comm_cost=0;
        if(vertex_replication_factor != NULL) *vertex_replication_factor=0;
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

        // keep track of vertex replicas on each partition with a set per partition
        std::vector<std::set<int> > vertex_replicas(num_processes);
        
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
                if(vertex_replication_factor != NULL) vertex_replicas[partitioning[from]].insert(from);
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
                        if(vertex_replication_factor != NULL) vertex_replicas[to_part].insert(from);
                        if(vertex_replication_factor != NULL) vertex_replicas[partitioning[from]].insert(to);
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
        
        if(vertex_replication_factor != NULL) {
            for(int pin=0; pin < total_vertices; pin++) {
                for(int pp=0; pp < num_processes; pp++) {
                    int present = vertex_replicas[pp].find(pin) != vertex_replicas[pp].end() ? 1 : 0;
                    int result;
                    MPI_Allreduce(&present,&result,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
                    if(result > 0)*vertex_replication_factor += 1;
                }
            }
            *vertex_replication_factor /= total_vertices;
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
            double vertex_replicastion_factor;
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

    int BaselineSequentialVertexPartitioning(char* experiment_name, idx_t* partitioning, int num_processes, std::string hypergraph_filename, int* vtx_wgt, float imbalance_tolerance) {
        // This algorithm corresponds to a basic implementation of Alistairh's hypergraph partitioning
        // From partitions that are not overloaded, select the one for which there is higher topic intersection (hyperedge intersection)

        // get meta info (num vertices and hyperedges)
        int num_vertices, num_hyperedges;
        get_hypergraph_file_header(hypergraph_filename, &num_vertices, &num_hyperedges);
        
        // 2 - Divide the graph in a distributed CSR format (like ParMETIS)
        //  compressed vertex or compressed hedge format? --> see zoltan
        //  for each local vertex, store the list of vertices adjacent to it (belonging to same hedges)
        std::vector<std::vector<int> > hyperedges;
        std::vector<std::vector<int> > hedge_ptr;
        load_hypergraph_from_file_dist_CSR(hypergraph_filename, &hyperedges, &hedge_ptr, MASTER_NODE, partitioning, true);        

        long int* part_load = (long int*)calloc(num_processes,sizeof(long int));

        memset(part_load,0,num_processes * sizeof(long int));
        /*double total_workload = 0;
        for(int ii=0; ii < num_vertices; ii++) {
            total_workload += vtx_wgt[ii];
        }*/
        double max_workload = (double)num_vertices / num_processes * imbalance_tolerance;
        
        std::vector<std::set<int> > hyperedges_in_partition(num_processes);

#ifdef DEBUG
        long int first_time_pins = 0;
        long int total_pins = 0;
        long int partition_filled = 0;
        bool filled_parts[num_processes];
        for(int ii=0; ii<num_processes; ii++) {
            filled_parts[ii] = false;
        }  
#endif

        // go through own vertex list and reassign
        for(int vid=0; vid < num_vertices; vid++) {
            // calculate hyperedge overlap for current vertex in all partitions
            std::vector<int> overlap(num_processes,0);
            for(int he = 0; he < hedge_ptr[vid].size(); he++) {
                int he_id = hedge_ptr[vid][he];
                bool seen = false;
                for(int pp=0; pp < num_processes; pp++) {
                    if(hyperedges_in_partition[pp].find(he_id) != hyperedges_in_partition[pp].end()) {
                        overlap[pp]++;
                        seen = true;
                    }
                }
#ifdef DEBUG
                if(!seen) first_time_pins++;
                total_pins++;
#endif
            }

            // allocate vertex to underloaded partitions
            int max_value = 0;
            int best_partition = 0;
            for(int pp=0; pp < num_processes; pp++) {
                if(part_load[pp] >= max_workload) {
#ifdef DEBUG
                    partition_filled++;
                    filled_parts[pp] = true;
#endif
                    continue;
                }

                int current_value = overlap[pp];
                
                if(current_value > max_value ||                                                     // best overlap value
                        (current_value == max_value && part_load[best_partition] > part_load[pp])) {  // equal overlap but less loaded
                    max_value = current_value;
                    best_partition = pp;
                } 
            }
            
            // update intermediate workload and assignment values
            part_load[best_partition] += vtx_wgt[vid];
            // update partitioning assignment
            partitioning[vid] = best_partition;
            for(int he = 0; he < hedge_ptr[vid].size(); he++) {
                int he_id = hedge_ptr[vid][he];
                hyperedges_in_partition[best_partition].insert(he_id);
            }
        }
        
        // clean up
        free(part_load);

#ifdef DEBUG

        PRINTF("Pins seen for first time: %li (%.2f %%)\n",first_time_pins,first_time_pins*1.0f/total_pins*100);

        int parts_full = 0;
        for(int ii=0; ii<num_processes; ii++) {
            if(filled_parts[ii]) parts_full++;
        }  

        PRINTF("Partitions filled [%i]: %li (%.2f)\n",parts_full,partition_filled,partition_filled*1.0f/num_processes);
#endif

        return 1;
    }


    // Parallel version of Alistairh hypergraph streaming partitioning
    int ParallelStreaming(char* experiment_name, idx_t* partitioning, int num_partitions, MPI_Comm partitioning_comm, std::string hypergraph_filename, int* element_wgt, float imbalance_tolerance, bool local_replica_degree_updates_only = false, int window_size = 1, bool input_order_round_robin = true, bool use_hdrf = false, bool staggered_streams=true, bool use_balance_cost = false, float lambda = 1.0f) {
        // Parallel Hyperedge Partitioning based algorithm
        // Because it can be applied to both vertex and hyperedge partitionings, we adopt the following nomenclature:
        //      element: what each line in the stream represent
        //          for vertex partitioning, each line is a hyperedge
        //          for hedge partitioning, each line is a vertex
        //      pin: each of the elements presented at a time by the stream
        //          for vertex partitioning, they are vertices that the hyperedge contain
        //          for hedge partitioning, they are hyperedges the vertex belongs to

        // The goal is to assign elements to partitions based on the pin overlap 
        //      element is assigned to partition that maximises the pin overlap
        // Constraint: keep max partition load based on imbalance factor
        
        // Input: 
        //      hypergraph, in a stream (one element at a time)
        //      number of partitions
        // Output: element assignment to partitions
        
        int process_id;
        MPI_Comm_rank(partitioning_comm,&process_id);
        int num_processes;
        MPI_Comm_size(partitioning_comm,&num_processes);
        
        // create shared data structures (partitions workload, list of replica destinations for each vertex, partial degree for each vertex)
        long int* part_load = (long int*)calloc(num_partitions, sizeof(long int));
        
        int num_pins;
        int num_elements;
    
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
        // TODO@ needs to account for batch sync update
        double max_expected_workload = (double)num_elements / num_partitions * imbalance_tolerance - (num_processes-1) * (double)window_size;
        if(max_expected_workload < 1) {
            PRINTF("Graph is too small! Too many processes (max expected workload limit is %f)\n",max_expected_workload);
            return 0;
        }    
        // each process should start partition assignment eval from its process_id (to avoid initial cramming of elements on initial partitions)
        int start_process = (float)process_id / num_processes * num_partitions;                
        
        PRINTF("Found in file: Pins: %i; elements %i:\n",num_pins,num_elements);

        // initialise workload
        // workload starts empty across all partitions
        long int total_workload = 0;

#ifdef DEBUG
        // number of pins seen for the first time
        long int first_time_pins = 0;
        long int total_pins = 0;
        long int partition_filled = 0;
        bool filled_parts[num_partitions];
        for(int ii=0; ii < num_partitions; ii++) filled_parts[ii] = false;
#endif

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

        
        int max_stream_size = num_elements / num_processes + ((num_elements % num_processes > 0) ? 1 : 0);
        while(element_id < max_stream_size) {

            int actual_window_size = 0;
            
            // data structures for batch synchronisation
            std::vector<int> local_pins_size(window_size,0);
            std::vector<std::vector<int> > new_replicas(window_size);
            
            // Load a batch of elements and their pins
            std::vector<std::vector<int> > batch_elements;
            for(int ww=0; ww < window_size; ww++) {
                if(std::getline(istream,line)) {
                    char str[line.length() + 1]; 
                    strcpy(str, line.c_str()); 
                    char* token = strtok(str, " "); 
                    std::vector<int> current_line;
                    while (token != NULL) { 
                        int pin_id = atoi(token)-1;
                        current_line.push_back(pin_id);
#ifdef DEBUG
                        if(seen_pins[pin_id].partial_degree == 0) first_time_pins++;
                        total_pins++;
#endif
                        // add it to local knowledge so it can be used in windowed partitioning
                        seen_pins[pin_id].partial_degree += 1;
                        token = strtok(NULL, " "); 
                    }
                    actual_window_size++;
                    batch_elements.push_back(current_line);
                }
            }
            // prioritise elements based on count of already seen pins
            // similar to ADWISE (delay uninformed decisions)

            std::vector<int> element_priority(actual_window_size);
            for(int el=0; el < actual_window_size; el++) {
                int score = 0;
                for(int pin=1; pin < batch_elements[el].size(); pin++) {
                    int pin_id = batch_elements[el][pin];
                    score += seen_pins[pin_id].partial_degree;
                }
                element_priority[el] = score;
            }
            // sort batch_elements according to their score
            std::vector<int> index(actual_window_size, 0);
            for (int i = 0 ; i < actual_window_size ; i++) {
                index[i] = i;
            }
            // TODO: test both sortings (highest to lowest and viceversa)
            sort(index.begin(), index.end(),
                [&](const int& a, const int& b) {
                    return (element_priority[a] > element_priority[b]); // from highest to lowest
                    //return (element_priority[a] < element_priority[b]); // from lowest to highest
                }
            );

            // process each batched element in order
            for(int el=0; el < actual_window_size; el++) {
                int idx = index[el];
                int n_pins = batch_elements[idx].size();
                //printf("%i: (%i) %i\n",process_id,el,element_priority[idx]);
                ////// local_pins is batch_elements[idx]
                ////// element_id is element_id + idx

                // calculate norm_part_degree for each vertex ONLY FOR HDRF
                // similar to Petroni HDRF
                double normalised_part_degrees[n_pins];
                if(use_hdrf) {
                    long int total_degrees = 0;
                    for(int ii=0; ii < n_pins; ii++) {
                        int pin_id = batch_elements[idx][ii];
                        normalised_part_degrees[ii] = seen_pins[pin_id].partial_degree; // if vertex is newly seen, it will be counted in the next sync
                        total_degrees += normalised_part_degrees[ii];
                    }
                    if(total_degrees > 0) {
                        std::transform(normalised_part_degrees,normalised_part_degrees+n_pins,normalised_part_degrees,
                            [total_degrees] (double value) {  
                                return value / total_degrees;
                            }  
                        );
                    }


                }
            
                // calculate C_rep(he) per partition per vertex
                //      sum 1 + (1-norm_part_degree(v)) if p exists in A(v)
                double max_value = std::numeric_limits<double>::lowest();
                int best_partition = 0;
                for(int pp=0; pp < num_partitions; pp++) {
                    // choose which partition to evaluate (per stream)
                    int current_part;
                    if(staggered_streams) {
                        // each process should start from its process_id (to avoid initial cramming of elements on initial partitions)
                        current_part = start_process + pp;
                        if(current_part >= num_partitions) current_part -= num_partitions;

                    } else {
                        current_part = pp;
                    }
                    
                    if(part_load[current_part] >= max_expected_workload) {
#ifdef DEBUG
                        partition_filled++;
                        filled_parts[current_part] = true;
#endif
                        continue;
                    }
                    double c_rep = 0;
                    double c_comm = 0;
                    for(int vv=0; vv < n_pins; vv++) {
                        int pin_id = batch_elements[idx][vv];
                        bool present_in_partition = false;
                        std::set<int>::iterator it;
                        if(use_hdrf) {
                            for (it = seen_pins[pin_id].A.begin(); it != seen_pins[pin_id].A.end(); ++it)
                            {
                                int part = *it;
                                present_in_partition |= part == current_part;
                            }

                            c_rep += present_in_partition ? 1 + (1 - normalised_part_degrees[vv]) : 0;
                        } else {
                            for (it = seen_pins[pin_id].A.begin(); it != seen_pins[pin_id].A.end(); ++it)
                            {
                                int part = *it;
                                if(part == current_part) {
                                    c_rep += 1;
                                    break;
                                }
                            }
                        }
                    }

                    double current_value = c_rep;
                    if(use_balance_cost) {
                        float c_bal = 1 * pow(part_load[current_part],lambda);//lambda * pow(part_load[current_part],0.5f);
                        current_value -= c_bal;
                    }
                    
                    if(current_value > max_value ||                                                 
                                (current_value == max_value && part_load[best_partition] > part_load[current_part])) {
                        max_value = current_value;
                        best_partition = current_part;
                    }
                }
                element_mapping = best_partition;

                // synchronise data
                // ***** 
                //  TODO: Batch synchronisation requres testing. Are we saving time? Are results correct?
                //  issues: tradeoff between less comm overhead and quality of partition (potentially higher graph size requirements as partition load is not updated often)
                // *****
                
                new_replicas[idx].push_back(element_mapping); // add in front the partition selected
                part_load[element_mapping] += 1; //  TODO: not really the element id  !! should be using element_wgt[current_element_id]
                int new_pins = 0;
                for(int ii=0; ii < batch_elements[idx].size(); ii++) {
                    int pin_id = batch_elements[idx][ii];
                    if(!local_replica_degree_updates_only && use_hdrf) {
                        new_replicas[idx].push_back(pin_id);
                        new_pins++;
                    } else {
                        // TODO: test performance degradation when only updating partial degree with local info
                        if(seen_pins[pin_id].A.find(element_mapping) == seen_pins[pin_id].A.end()) {
                            // if it has been seen but it's the first replica on new partition
                            new_replicas[idx].push_back(pin_id);
                            new_pins++;
                        } else {
                            // we already account for this when reading the stream
                            //seen_pins[pin_id].partial_degree += 1;
                        }
                    }
                    // must update seen_pins data structure as it goes along
                    // then during remote sync avoid double counting local pins
                    // we already updated seen_pins[].part_degree when we read the stream
                    seen_pins[pin_id].A.insert(element_mapping);
                }
                local_pins_size[idx] = new_pins;                   
            }
            element_id += window_size;

            
            
            // batch synchronisation
            // synchronise length of pins list to be sent
            int* remote_pins_size = (int*)malloc(sizeof(int) * window_size * num_processes);

            MPI_Allgather(&local_pins_size[0],window_size,MPI_INT,remote_pins_size,window_size,MPI_INT,partitioning_comm);
            
            // synchronise list of pins
            // share send buffer size with other processes
            // size = number of new replicas + 1 (the partition selected)
            // flatten new_replicas first
            std::vector<int> new_replicas_sync;
            
            for(int idx = 0; idx < window_size; idx++) {
                if(idx >= actual_window_size) continue;
                new_replicas_sync.insert(new_replicas_sync.end(),new_replicas[idx].begin(),new_replicas[idx].end());
            }
            int* recvcounts = (int*)malloc(num_processes * sizeof(int));
            int send_size = new_replicas_sync.size();                
            
            MPI_Allgather(&send_size,1,MPI_INT,recvcounts,1,MPI_INT,partitioning_comm);
            
            // share partition selected and new replicas list to all
            int counter = 0;
            int* displs = (int*)malloc(sizeof(int) * num_processes);
            for(int ii=0; ii < num_processes; ii++) {
                displs[ii] = counter;
                counter += recvcounts[ii]; 
            }

            int* recvbuffer = (int*)malloc(sizeof(int) * counter);
            MPI_Allgatherv(&new_replicas_sync[0],send_size,MPI_INT,recvbuffer,recvcounts,displs,MPI_INT,partitioning_comm);
            
            // process new replicas and add them to local datastructures
            for(int ii=0; ii < num_processes; ii++)  {
                int current_element_index = 0;
                for(int el_order=0; el_order < window_size; el_order++) {
                    // Choice between round robin or bulk first //
                    int current_element;
                    if(input_order_round_robin) {
                        current_element = element_id * num_processes - (num_processes * (window_size - el_order)) + ii;
                        // check expected limits for global element id
                        if(current_element >= num_elements) break;
                    } else {
                        int first_element_in_stream = num_elements / num_processes * ii + std::min(num_elements % num_processes,ii);
                        current_element = first_element_in_stream + element_id - (window_size - el_order);
                        // check expected limits for current process id
                        int current_stream_size = num_elements / num_processes + ((num_elements % num_processes > ii) ? 1 : 0);
                        if(current_element >= first_element_in_stream + current_stream_size) break;
                    }
                    int current_pin_length = remote_pins_size[ii*window_size+el_order];
                    int dest_partition = recvbuffer[displs[ii] + current_element_index];
                    
                    // move current_element_index to start reading pins
                    current_element_index++;
                    
                    total_workload += element_wgt[current_element];
                    partitioning[current_element] = dest_partition;

                    if(ii != process_id) {
                        // only update pins from remote processes (they were accounted for during local streaming)
                        part_load[dest_partition] += element_wgt[current_element];
                        for(int jj=0; jj < current_pin_length; jj++) {
                            int pin_id = recvbuffer[displs[ii] + current_element_index + jj];
                            
                            // update remote vertex info
                            seen_pins[pin_id].partial_degree += 1;
                            seen_pins[pin_id].A.insert(dest_partition);  
                                    
                        }
                    }
                    // shift current_element_index to read the next element
                    current_element_index += current_pin_length;
                }

            }
            free(recvcounts);
            free(displs);
            free(recvbuffer);
            free(remote_pins_size);
        }
        istream.close();
        
        

#ifdef DEBUG 
        // check for termination condition (tolerance imbalance reached)   
        long int max_sz = *std::max_element(part_load, part_load + num_partitions);
        double max_imbalance = ((double)max_sz) / ((double)total_workload/num_partitions);
        
        PRINTF("***Imbalance: %.3f\n",max_imbalance); 

        PRINTF("%i: Pins seen for the first time: %li (%.2f %%)\n",process_id,first_time_pins,first_time_pins*1.0f/total_pins*100);

        int parts_full = 0;
        for(int ii=0; ii<num_partitions; ii++) {
            if(filled_parts[ii]) parts_full++;
        }  

        PRINTF("%i: partitions filled [%i]: %li (%li)\n",process_id,parts_full,partition_filled,partition_filled*num_partitions);

#endif

        // clean up
        free(part_load);

        return 1;

    }

    
    // Stream from multiple files / streams
    int ParallelHDRF(char* experiment_name, idx_t* partitioning, int num_partitions, MPI_Comm partitioning_comm, double** comm_cost_matrix, std::string hypergraph_filename, int* element_wgt, int max_iterations, float imbalance_tolerance, bool local_replica_degree_updates_only = false, int window_size = 1, bool input_order_round_robin = true, float lambda = 1.0f,bool save_partitioning_history = false) {
        // Parallel Hyperedge Partitioning based algorithm
        // Because it can be applied to both vertex and hyperedge partitionings, we adopt the following nomenclature:
        //      element: what each line in the stream represent
        //          for vertex partitioning, each line is a hyperedge
        //          for hedge partitioning, each line is a vertex
        //      pin: each of the elements presented at a time by the stream
        //          for vertex partitioning, they are vertices that the hyperedge contain
        //          for hedge partitioning, they are hyperedges the vertex belongs to

        // The goal is to assign elements to partitions   
        // Minimisation goal: 
        //      replication factor (sum of pin replicas)
        //      reducing communication cost (taking comm cost matrix into account)
        // Constraint: keep load balance amongst partitions (weighed sum of hyperedges per partition = sum of vertices per partition?)
        // It can be a restreaming approach, in which in the first stream we use partial info (partial degree) and full info afterwards
        //      particularly when incorporating comm cost into the objective function

        // Input: 
        //      hypergraph, in a stream (one element at a time)
        //      number of partitions
        //      communication cost matrix
        // Output: element assignment to partitions

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

        // own parameters
        float lambda_update = 0.05f; // must be greater than 1. Used when partitions are too imbalanced at the end of the pass
        float lambda_refine = -0.02f; // must be lower than 1. Used when partitions are within imbalance limits at the end of the pass

        int process_id;
        MPI_Comm_rank(partitioning_comm,&process_id);
        int num_processes;
        MPI_Comm_size(partitioning_comm,&num_processes);
        
        // create shared data structures (partitions workload, list of replica destinations for each vertex, partial degree for each vertex)
        long int* part_load = (long int*)calloc(num_partitions, sizeof(long int));

        int num_pins;
        int num_elements;

        // avoid overfitting variables
        bool check_overfit = false;
        bool rollback = false;
        idx_t* last_partitioning = NULL;

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
            // TODO@ needs to account for batch sync update
            double max_expected_workload = (double)num_elements / num_partitions * imbalance_tolerance - (num_processes-1) * (double)window_size;
            double average_expected_workload = (double)num_elements / num_partitions;
            if(max_expected_workload < 1) {
                PRINTF("Graph is too small! Too many processes (max expected workload limit is %f)\n",max_expected_workload);
                return 0;
            }    
            // each process should start partition assignment eval from its process_id (to avoid initial cramming of elements on initial partitions)
            int start_process = (float)process_id / num_processes * num_partitions;  
            
            PRINTF("Found in file: Pins: %i; elements %i:\n",num_pins,num_elements);

            // initialise workload
            memset(part_load,0,num_partitions * sizeof(long int));
            // workload starts empty across all partitions
            long int total_workload = 0;


    #ifdef DEBUG
            // number of pins seen for the first time
            long int first_time_pins = 0;
            long int total_pins = 0;
            long int partition_filled = 0;
            bool filled_parts[num_partitions];
            for(int ii=0; ii < num_partitions; ii++) filled_parts[ii] = false;
    #endif

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
            //      std::unordered_map<short,short> --> store the partitions that have a replica of the pin plus the pins and the number of replications on that partition
            // store the Vertex in a std::unordered_map (hashmap) of Vertex* of length num_vertices (vertex id is given by the index)
            std::unordered_map<int,pin_data> seen_pins;

            // read reminder of file (one line per elemennt)
            int element_id = 0; // current elemennt
            int element_mapping = -1; // mapping of current he (to partition)

            
            int max_stream_size = num_elements / num_processes + ((num_elements % num_processes > 0) ? 1 : 0);
            while(element_id < max_stream_size) {
                int actual_window_size = 0;
                
                // data structures for batch synchronisation
                std::vector<int> local_pins_size(window_size,0);
                std::vector<std::deque<int> > new_replicas(window_size);
                
                // Load a batch of elements and their pins
                std::vector<std::vector<int> > batch_elements;
                for(int ww=0; ww < window_size; ww++) {
                    if(std::getline(istream,line)) {
                        char str[line.length() + 1]; 
                        strcpy(str, line.c_str()); 
                        char* token = strtok(str, " "); 
                        std::vector<int> current_line;
                        while (token != NULL) { 
                            int pin_id = atoi(token)-1;
                            current_line.push_back(pin_id);
    #ifdef DEBUG
                            if(seen_pins[pin_id].partial_degree == 0) first_time_pins++;
                            total_pins++;
    #endif
                            token = strtok(NULL, " "); 
                        }
                        actual_window_size++;
                        batch_elements.push_back(current_line);
                    }
                }
                // prioritise elements based on count of already seen pins
                // similar to ADWISE (delay uninformed decisions)
                std::vector<int> element_priority(actual_window_size);
                for(int el=0; el < actual_window_size; el++) {
                    int score = 0;
                    for(int pin=1; pin < batch_elements[el].size(); pin++) {
                        int pin_id = batch_elements[el][pin];
                        score += seen_pins[pin_id].partial_degree;
                    }
                    element_priority[el] = score;
                }
                // sort batch_elements according to their score
                std::vector<int> index(actual_window_size, 0);
                for (int i = 0 ; i < actual_window_size ; i++) {
                    index[i] = i;
                }
                // TODO: test both sortings (highest to lowest and viceversa)
                sort(index.begin(), index.end(),
                    [&](const int& a, const int& b) {
                        return (element_priority[a] > element_priority[b]); // from highest to lowest
                    }
                );
                // process each batched element in order
                for(int el=0; el < actual_window_size; el++) {
                    int idx = index[el];
                    int num_pins = batch_elements[idx].size();
                    //printf("%i: (%i) %i\n",process_id,el,element_priority[idx]);
                    ////// local_pins is batch_elements[idx]
                    ////// element_id is element_id + idx
                
                
                    // local hyperedge, process and assign it
                    // calculate norm_part_degree for each vertex
                    // similar to Petroni HDRF
                    /*double normalised_part_degrees[num_pins];
                    long int total_degrees = 0;
                    for(int ii=0; ii < num_pins; ii++) {
                        int pin_id = batch_elements[idx][ii];
                        normalised_part_degrees[ii] = seen_pins[pin_id].partial_degree; // if vertex is newly seen, it will be counted in the next sync
                        total_degrees += normalised_part_degrees[ii];
                    }
                    if(total_degrees > 0) {
                        std::transform(normalised_part_degrees,normalised_part_degrees+num_pins,normalised_part_degrees,
                            [total_degrees] (double value) {  
                                return value / total_degrees;
                            }  
                        );
                    }*/

                    // instead of looping through the pins on every partition, loop once through pins
                    // create cost buckets for each partition and fill them once.
                    // Keep track of lowest cost partition
                    // still, staggered start when deciding what partition has highest score
                    // optimisation trick to avoid having to calculate communication cost to partitions that do not have any replicas in it
                    std::unordered_map<unsigned int,unsigned long> neighbours_in_partition;
                    double max_value = std::numeric_limits<double>::lowest();
                    int best_partition = 0;
                    int new_pins = 0;
                    for(int vv=0; vv < num_pins; vv++) {
                        int pin_id = batch_elements[idx][vv];

                        // update new replicas (for parallel synchronisation)
                        new_replicas[idx].push_back(pin_id);
                        new_pins++;

                        std::unordered_map<unsigned short,unsigned short>::iterator it;
                        for (it = seen_pins[pin_id].P.begin(); it != seen_pins[pin_id].P.end(); ++it)
                        {
                            unsigned short part = it->first;
                            unsigned short replicas = it->second;
                            //present_in_partition |= part == current_part;
                            // communication should be proportional to the duplication of pins
                            // if a pin is duplicated in two partitions, then communication will happen across those partitions
                            neighbours_in_partition[part] += replicas;
                            // if we were to use c_rep (overlap measure) we would need to do extra computation here
                        }
                    }
                    local_pins_size[idx] = new_pins; 

                    // calculate cost of communication for all partition candidates
                    // loop through the map only once
                    double* c_comms = (double*)calloc(num_partitions,sizeof(double));
                    std::unordered_map<unsigned int, unsigned long>::iterator it = neighbours_in_partition.begin();
                    long max_replicaset_load = 0;
                    long min_replicaset_load = std::numeric_limits<long>::max();
                    while(it != neighbours_in_partition.end()) {
                        unsigned int part = it->first;
                        unsigned long replicas = it->second;
                        for(int origin=0; origin < num_partitions; origin++) {
                            c_comms[origin] += comm_cost_matrix[origin][part] * replicas;
                        }
                        // we can calculate a upper load limit for acceptable partitions
                        // based on the minimum cost of communication (max replicas) amongst the partitions with current replicas
                        // i.e. the limit of c_bal a partition without replicas would have to be below to even compete with those with replicas
                        // this removes the need to check many partitions later
                        if(part_load[part] > max_replicaset_load) {
                            max_replicaset_load = part_load[part];
                        } else if(part_load[part] < min_replicaset_load) {
                            min_replicaset_load = part_load[part];
                        }
                        it++;
                    }
                    for(int pp=0; pp < num_partitions; pp++) {
                        // each process should start from its process_id (to avoid initial cramming of elements on initial partitions)
                        int current_part = start_process + pp;
                        if(current_part >= num_partitions) current_part -= num_partitions;
                        long current_load = part_load[current_part];
                        if(current_load >= max_expected_workload) { // partition full
    #ifdef DEBUG
                            partition_filled++;
                            filled_parts[current_part] = true;
    #endif
                            continue;
                        }
                        // check for part load upper and lower limits
                        // this avoids us having to calculate any further (avoids pow call)
                        if((max_replicaset_load > 0 && current_load > max_replicaset_load)
                                        || (neighbours_in_partition[current_part] == 0 && current_load >= min_replicaset_load)) {
                            continue;
                        }
                        
                        // calculate total cost
                        float c_bal = pow(current_load,lambda);
                        double current_value = - c_comms[current_part] - c_bal;
                        
                        if(current_value > max_value ||                                                 
                                    current_value == max_value && part_load[best_partition] > current_load) {
                            max_value = current_value;
                            best_partition = current_part;
                        }

                    }
                    free(c_comms);

                    element_mapping = best_partition;
                    new_replicas[idx].push_front(element_mapping); // add in front the partition selected
                    part_load[element_mapping] += element_wgt[element_id + idx];

                    // synchronise data
                    // ***** 
                    //  TODO: test just synchronising partition sizes, not vertex degree
                    //  assumption: if pins are sufficiently shuffled, partial local degree may be enough and saves data shared
                    // *****
                    // ***** 
                    //  TODO: Batch synchronisation requres testing. Are we saving time? Are results correct?
                    //  issues: tradeoff between less comm overhead and quality of partition (potentially higher graph size requirements as partition load is not updated often)
                    // *****
                    /*new_replicas[idx].push_back(element_mapping); // add in front the partition selected
                    part_load[element_mapping] += element_wgt[element_id + idx];
                    int new_pins = 0;
                    for(int ii=0; ii < batch_elements[idx].size(); ii++) {
                        int pin_id = batch_elements[idx][ii];
                        
                        if(!local_replica_degree_updates_only) {
                            new_replicas[idx].push_back(pin_id);
                            new_pins++;
                        } else {
                            // only update remote streams when new in - partitioning mappings are found
                            // TODO: test performance degradation when only updating partial degree with local info
                            if(seen_pins[pin_id].A.find(element_mapping) == seen_pins[pin_id].A.end()) {
                                // if it has been seen but it's the first replica on new partition
                                new_replicas[idx].push_back(pin_id);
                                new_pins++;
                            } else {
                                seen_pins[pin_id].partial_degree += 1;
                            }
                        }
                        // must update seen_pins data structure as it goes along
                        // then during remote sync avoid double counting local pins
                        // we already updated seen_pins[].part_degree when we read the stream
                        // TODO: we don't need both seen_pins.A and .P, whichever we don't use can be inhibited
                        // seen_pins.A is still used by getEdgeCentricReplicationFactor to calculate vertex replication
                        seen_pins[pin_id].A.insert(element_mapping);
                        if(seen_pins[pin_id].P.find(element_mapping) == seen_pins[pin_id].P.end()) {
                            seen_pins[pin_id].P[element_mapping] = 1;
                        } else {
                            seen_pins[pin_id].P[element_mapping] += 1;
                        }
                    }
                    local_pins_size[idx] = new_pins; 
                    */
                                    
                }

                element_id += window_size;


                // batch synchronisation
                // synchronise list of pins
                // share send buffer size with other processes
                // size = number of new replicas + 1 (the partition selected)
                // flatten new_replicas first
                // send partial pin counts (one per element) and total send size (all pins for all elements) simultaneously
                std::vector<int> new_replicas_sync;
                
                for(int idx = 0; idx < window_size; idx++) {
                    if(idx >= actual_window_size) continue;
                    new_replicas_sync.insert(new_replicas_sync.end(),new_replicas[idx].begin(),new_replicas[idx].end());
                }
                int send_size = new_replicas_sync.size();

                local_pins_size.push_back(send_size);
                int* remote_pins_size = (int*)malloc(sizeof(int) * (window_size+1) * num_processes);

                MPI_Allgather(&local_pins_size[0],window_size+1,MPI_INT,remote_pins_size,window_size+1,MPI_INT,partitioning_comm);
                
                // share partition selected and new replicas list to all
                int counter = 0;
                int* displs = (int*)malloc(sizeof(int) * num_processes);
                int* recvcounts = (int*)malloc(num_processes * sizeof(int));
                for(int ii=0; ii < num_processes; ii++) {
                    displs[ii] = counter;
                    recvcounts[ii] = remote_pins_size[ii*(window_size+1) + window_size];
                    counter += recvcounts[ii]; 
                }
                int* recvbuffer = (int*)malloc(sizeof(int) * counter);
                MPI_Allgatherv(&new_replicas_sync[0],send_size,MPI_INT,recvbuffer,recvcounts,displs,MPI_INT,partitioning_comm);
                
                // process new replicas and add them to local datastructures
                for(int ii=0; ii < num_processes; ii++)  {
                    int current_element_index = 0;
                    for(int el_order=0; el_order < window_size; el_order++) {
                        // Choice between round robin or bulk first //
                        int current_element;
                        if(input_order_round_robin) {
                            current_element = element_id * num_processes - (num_processes * (window_size - el_order)) + ii;
                            // check expected limits for global element id
                            if(current_element >= num_elements) break;
                        } else {
                            int first_element_in_stream = num_elements / num_processes * ii + std::min(num_elements % num_processes,ii);
                            current_element = first_element_in_stream + element_id - (window_size - el_order);
                            // check expected limits for current process id
                            int current_stream_size = num_elements / num_processes + ((num_elements % num_processes > ii) ? 1 : 0);
                            if(current_element >= first_element_in_stream + current_stream_size) break;
                        }
                        int current_pin_length = remote_pins_size[ii*(window_size+1)+el_order];
                        int dest_partition = recvbuffer[displs[ii] + current_element_index];
                        // move current_element_index to start reading pins
                        current_element_index++;
                        
                        total_workload += element_wgt[current_element];
                        partitioning[current_element] = dest_partition;
                        // update pins from remote processes (including local ones too)
                        if(ii != process_id) {
                            // part load for local elements has already been updated
                            part_load[dest_partition] += element_wgt[current_element];
                        }
                        for(int jj=0; jj < current_pin_length; jj++) {
                            int pin_id = recvbuffer[displs[ii] + current_element_index + jj];
                            // update remote vertex info
                            seen_pins[pin_id].partial_degree += 1;
                            seen_pins[pin_id].A.insert(dest_partition); 
                            if(seen_pins[pin_id].P.find(dest_partition) == seen_pins[pin_id].P.end()) {
                                seen_pins[pin_id].P[dest_partition] = 1;
                            } else {
                                seen_pins[pin_id].P[dest_partition] += 1;
                            } 
                                    
                        }
                        // shift current_element_index to read the next element
                        current_element_index += current_pin_length;
                    }

                }
                free(recvcounts);
                free(displs);
                free(recvbuffer);
                free(remote_pins_size);  

            }
            istream.close();

            // check for termination condition
            // update lambda (importance of load balancing)
            // keep searching until the last cut metric is not improved
            // this algorithm is always guaranteed to find balanced partitions
            
            // test load balance (give more importance to load balance if the ratio max/min is too high)
            float maxsize = *std::max_element(part_load, part_load + num_partitions);
            float minsize = *std::min_element(part_load, part_load + num_partitions);
            float max_imbalance = minsize <= 0 ? num_partitions : maxsize / minsize;
            PRINTF("***Max-min ratio: %.3f, current lambda: %f.\n",max_imbalance,lambda); 

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

            // check if within imbalance allowance (max over min)
            if(max_imbalance < (imbalance_tolerance / (1.0f/imbalance_tolerance) * imbalance_tolerance)) {
                check_overfit = true;
                if(process_id == MASTER_NODE) {
                    if(last_partitioning == NULL) {
                        last_partitioning = (idx_t*)malloc(num_elements*sizeof(idx_t));
                    }
                    memcpy(last_partitioning,partitioning,num_elements * sizeof(idx_t));
                }              
                lambda += lambda_refine;   
            } else {
                // still too imbalanced
                if(check_overfit) {
                    rollback = true;
                    break;
                }
                lambda += lambda_update;
                if (lambda <= 0) break;
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
                    MPI_Send(partitioning,num_elements,MPI_LONG,dest,0,partitioning_comm);
                }
            } else {
                // update partitioning from 0
                MPI_Recv(partitioning,num_elements,MPI_LONG,MASTER_NODE,0,partitioning_comm,MPI_STATUS_IGNORE);
            }
        }

        
        

#ifdef DEBUG 
        // check for termination condition (tolerance imbalance reached)   
        long int max_sz = *std::max_element(part_load, part_load + num_partitions);
        double max_imbalance = ((double)max_sz) / ((double)total_workload/num_partitions);
        
        PRINTF("***Imbalance: %.3f\n",max_imbalance); 

        PRINTF("%i: Pins seen for the first time: %li (%.2f %%)\n",process_id,first_time_pins,first_time_pins*1.0f/total_pins*100);

        int parts_full = 0;
        for(int ii=0; ii<num_partitions; ii++) {
            if(filled_parts[ii]) parts_full++;
        }  

        PRINTF("%i: partitions filled [%i]: %li (%li)\n",process_id,parts_full,partition_filled,partition_filled*num_partitions);

#endif

        // clean up
        free(part_load);

        return iter+1;

    }
    

    /*
    // Stream from multiple files / streams
    int ParallelHDRF(char* experiment_name, idx_t* partitioning, double** comm_cost_matrix, std::string hypergraph_filename, int* element_wgt, int max_iterations, float imbalance_tolerance, bool save_partitioning_history, bool local_replica_degree_updates_only = false, int sync_batch_size = 1, bool use_max_expected_workload = true, bool input_order_round_robin = true) {
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
        //double max_comm_cost = 0;
        //for(int ff=0; ff < num_processes; ff++) {
        //    double current_max_cost = 0;
        //    for(int tt=0; tt < num_processes; tt++) {
        //        current_max_cost += comm_cost_matrix[ff][tt];
        //    }
        //    if(current_max_cost > max_comm_cost) {
        //        max_comm_cost = current_max_cost;
        //    }
        //}

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
            // TODO@ needs to account for batch sync update
            long int max_expected_workload = num_elements / num_processes * imbalance_tolerance - num_processes * sync_batch_size;
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

            int local_stream_size = num_elements / num_processes + ((num_elements % num_processes > 0) ? 0 : 1);
            while(element_id < local_stream_size) {
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
                        }
                    }
                    element_mapping = best_partition;

                    free(c_comms);
                    free(c_total);
                }

                element_id += 1;//num_processes;

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
                //if((element_id / num_processes) % sync_batch_size == 0) {
                if((element_id ) % sync_batch_size == 0) {
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
                            // Choice between round robin or bulk first //
                            int current_element;
                            if(input_order_round_robin) {
                                current_element = element_id * num_processes - (num_processes * (sync_batch_size - el_order)) + ii;
                            } else {
                                int first_element_in_stream = (num_elements / num_processes + std::min(num_elements % num_processes,ii)) * ii;
                                current_element = first_element_in_stream + element_id - (sync_batch_size - el_order);
                            }
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
                    //long int max = part_load[0];
                    //long int min = part_load[0];
                    //for(int pp=0; pp < num_processes; pp++) {
                    //    if(part_load[pp] > max) max = part_load[pp];
                    //    if(part_load[pp] < min) min = part_load[pp];
                    //}
                    //minsize = min;
                    //maxsize = max;
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
    */
     
    

    
}

#endif