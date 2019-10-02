#ifndef ZOLTAN_COMPRESSED_VERTEX_PARTITIONING
#define ZOLTAN_COMPRESSED_VERTEX_PARTITIONING

#include "zoltan.h"
#include <vector>
#include "Partitioning.h"
#include "PRAW.h"
#include <cstring>
#include <cstdio>
#include <stdlib.h>
#include <string>

/* VERTEX PARTITIONING VERSION
TODO: use distributed CSR to avoid having to load the entire graph on all processes */

//////////////////////////////////////////////////
// ZOLTAN SPECIFIC DATA STRUCTURES AND FUNCTIONS
//////////////////////////////////////////////////

typedef struct{
	// Zoltan will partition vertices, while minimizing edge cuts 
	int numMyVertices;  // number of vertices that I own initially 
	ZOLTAN_ID_TYPE *vtxGID;        // global ID of these vertices 
	float* vtx_wgts;	// weights of each of the vertices 
	//int numMyHEdges;    // number of my hyperedges 
	//float* edge_wgts;	// weight of hyperedges 
	int npins; // number of vertices in my hyperedges 
	//ZOLTAN_ID_TYPE *edgeGID;       // global ID of each of my hyperedges 
	int *vtxedge_ptr;     // index into pin_GID array of edge's vertices  // vtxedge_ptr
	ZOLTAN_ID_TYPE *pin_GID;  // Edges of vertex vtxGID[i] begin at pin_GID[vtxedge_ptr[i]] // pin_GID
	int process_id;	// what process this graph resides in
	idx_t* partitioning;	// pointer to the partitioning array
} CV_HGRAPH_DATA;

// Application defined query functions. 
#pragma region Zoltan_app_functions
int cv_get_number_of_vertices(void *data, int *ierr) {
	CV_HGRAPH_DATA *hg = (CV_HGRAPH_DATA *)data;
	PRINTF("%i: get_number_of_vertices\n\n",hg->process_id);
	*ierr = ZOLTAN_OK;
	return hg->numMyVertices;
}
void cv_get_vertex_list(void *data, int num_gid_entries, int num_lid_entries,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr) {
	CV_HGRAPH_DATA *hg= (CV_HGRAPH_DATA *)data;
	PRINTF("%i: get_vertex_list\n\n",hg->process_id);

	*ierr = ZOLTAN_OK;
	// For weights: 
	// 	wgt_dim is the number of weights (set by parameter OBJ_WEIGHT_DIM)
	// 	obj_wgts: associated weights. Weights for object i are stored in obj_wgts[(i-1)*wgt_dim:i*wgt_dim-1]
	//
	for (int i=0; i<hg->numMyVertices; i++){
		globalID[i] = hg->vtxGID[i];
		localID[i] = i;
		if(wgt_dim > 0) obj_wgts[i] = hg->vtx_wgts[i];
	}

}
void cv_get_hypergraph_size(void *data, int *num_lists, int *num_pins,
                                int *format, int *ierr){
	CV_HGRAPH_DATA *hg = (CV_HGRAPH_DATA *)data;
	PRINTF("%i: get_hypergraph_size (%i, %i)\n\n",hg->process_id,hg->numMyVertices,hg->npins);
	*ierr = ZOLTAN_OK;
	*num_lists = hg->numMyVertices; // number of vertices
	*num_pins = hg->npins; // number of elements in pin_GID
	*format = ZOLTAN_COMPRESSED_VERTEX;
	
}
void cv_get_hypergraph(void *data, int num_gid_entries, int num_vtx_edge, int num_pins,
                           int format, ZOLTAN_ID_PTR vtxedge_GID, int *vtxedge_ptr,
                           ZOLTAN_ID_PTR pin_GIDs, int *ierr) {
	CV_HGRAPH_DATA *hg = (CV_HGRAPH_DATA *)data;
	PRINTF("%i: get_hypergraph\n\n",hg->process_id);
	*ierr = ZOLTAN_OK;

	if ( (num_vtx_edge != hg->numMyVertices) || (num_pins != hg->npins) ||
	   (format != ZOLTAN_COMPRESSED_VERTEX)) {
		*ierr = ZOLTAN_FATAL;
		return;
	}
	int i;
	for (i=0; i < num_vtx_edge; i++){
		vtxedge_GID[i] = hg->vtxGID[i];
		vtxedge_ptr[i] = hg->vtxedge_ptr[i];
	}
	if(num_vtx_edge > 0)
		vtxedge_ptr[num_vtx_edge] = hg->vtxedge_ptr[num_vtx_edge];
	for (i=0; i < num_pins; i++){
		pin_GIDs[i] = hg->pin_GID[i];
	}
}
void cv_user_migration_multi_object_size(void *data, int num_gid_entries, int num_lid_entries, int num_ids, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,int* sizes, int *ierr) {
	CV_HGRAPH_DATA *hg = (CV_HGRAPH_DATA *)data;
	PRINTF("%i: user_migration_multi_object_size\n\n",hg->process_id);
	// set sizes. We will transfer the pins (connections) for each vertex and the vertex weight (represented as int here)
	for(int ii=0; ii < num_ids; ii++) {
		ZOLTAN_ID_TYPE local_id = local_ids[ii];
		int pins_for_vertex = hg->vtxedge_ptr[local_id+1] - hg->vtxedge_ptr[local_id];
		sizes[ii] = pins_for_vertex * sizeof(ZOLTAN_ID_TYPE) + sizeof(int);
	}
	*ierr = ZOLTAN_OK;
	
}
void cv_user_pack_multi_obj(void *data, int num_gid_entries, int num_lid_entries, int num_ids, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int* dest, int* sizes, int* idx, char *buf, int *ierr) {
	
	CV_HGRAPH_DATA *hg = (CV_HGRAPH_DATA *)data;
	PRINTF("%i: user_pack_multi_obj\n\n",hg->process_id);
	bool pins_migrated = false;
	// transfer ownership
	std::vector<ZOLTAN_ID_TYPE> local_indices(num_ids);
	for(int ii=0; ii < num_ids; ii++) {
		local_indices[ii] = local_ids[ii*num_lid_entries];
	} 

	// DON'T DO MIGRATION (POSSIBLE OUT OF MEMORY IN ARCHER)
	for(int ii=0; ii < num_ids; ii++) {
		// TODO: needs to put data into buf, and index the data in idx array
		//idx_t* partitioning;	// pointer to the partitioning array
		hg->partitioning[global_ids[ii*num_gid_entries]] = dest[ii];
	}
}
void cv_user_unpack_multi_obj(void *data, int num_gid_entries, int num_ids, ZOLTAN_ID_PTR global_ids, int* sizes, int* idx, char *buf, int *ierr) {
	CV_HGRAPH_DATA *hg = (CV_HGRAPH_DATA *)data;
	PRINTF("%i: user_unpack_multi_obj\n\n",hg->process_id);
	if(num_ids <= 0 ) {
		*ierr = ZOLTAN_OK;
		return;
	}

	// DON'T DO MIGRATION (POSSIBLE OUT OF MEMORY IN ARCHER)
	for(int ii=0; ii < num_ids; ii++) {
		hg->partitioning[global_ids[ii*num_gid_entries]] = hg->process_id;
	}	
}
#pragma endregion


class ZoltanCompressedVertexPartitioning : public Partitioning {
public:
	
	ZoltanCompressedVertexPartitioning(char* graph_file, float imbalance_tolerance) : Partitioning(graph_file, imbalance_tolerance) {
		float ver;
		rc = Zoltan_Initialize(0, NULL, &ver);
	
		if (rc != ZOLTAN_OK){
			printf("Error initialising Zoltan...\n");
			MPI_Finalize();
			exit(0);
		}

		zz = Zoltan_Create(MPI_COMM_WORLD);

		setZoltanParams();
		// initialise hg struct
		hg.vtxGID = NULL;
		hg.vtx_wgts = NULL;
		hg.vtxedge_ptr = NULL;
		hg.pin_GID = NULL;
	}
	virtual ~ZoltanCompressedVertexPartitioning() {
		Zoltan_Destroy(&zz);
	}
	
	virtual void perform_partitioning(int num_processes,int process_id, int* iterations) {
		
		*iterations = 1;

		
		hg.process_id = process_id; 

		// Assign unique vertices to partitions
		for(int ii=0; ii < num_vertices; ii++) {
            partitioning[ii] = ii % num_processes;
        }

		hg.partitioning = partitioning;

		// 1 Process the hgraph file and build local datastructures
        std::ifstream istream(hgraph_name);
        
        if(!istream) {
            printf("Error while opening hMETIS file %s\n",hgraph_name);
            return;
        }

		// Assign local hyperedges to each partition
		// round robin allocation based on hedge id
		hg.numMyVertices = num_vertices / num_processes + (num_vertices % num_processes > process_id ? 1 : 0);
		
		std::vector<ZOLTAN_ID_TYPE> pin_GID;

        std::string line;
        // skip header
        std::getline(istream,line);

		std::vector<std::vector<int> > hyperedges(num_vertices);
			
        // read reminder of file (one line per hyperedge)
		int current_hedge_id = 0;
		while(std::getline(istream,line)) {
			// add it to local datastructure if hedge id is assigned to local partition
			char str[line.length() + 1]; 
			strcpy(str, line.c_str()); 
			char* token = strtok(str, " "); 
			while (token != NULL) { 
				int vtx = atoi(token) - 1;
				if(vtx % num_processes == process_id) {
					hyperedges[vtx].push_back(current_hedge_id);
				}
				token = strtok(NULL, " ");
			}
			current_hedge_id++;
		}
		
		// initialise zoltan datastructures
		hg.vtxGID = (ZOLTAN_ID_TYPE*)malloc(sizeof(ZOLTAN_ID_TYPE) * hg.numMyVertices);
		hg.vtxedge_ptr = (int*)malloc(sizeof(int) * (hg.numMyVertices+1));
		hg.vtx_wgts = (float*)calloc(hg.numMyVertices, sizeof(float));
		for(int ii=0; ii < hg.numMyVertices; ii++) {
			hg.vtx_wgts[ii] = 1;
		}
		
		// build structures
		int local_vtx_index = 0;
		hg.npins = 0;
		for(int vtx=0; vtx < num_vertices; vtx++) {
			if(vtx % num_processes == process_id) {
				hg.vtxGID[local_vtx_index] = vtx;
				hg.vtxedge_ptr[local_vtx_index] = hg.npins;
				for(int he=0; he < hyperedges[vtx].size(); he++) {
					int he_id = hyperedges[vtx][he];
					pin_GID.push_back(he_id);
					hg.npins++;
				}
				local_vtx_index++;
			}
		}
		hg.vtxedge_ptr[hg.numMyVertices] = hg.npins;	
		
		hg.pin_GID = &(pin_GID[0]);
		
		





		
		/* OLD METHOD, TOO HEAVY IN MEMORY */
		/*
		// Assign unique vertices to partitions
		hg.numMyVertices = 0;
		for(int ii=0; ii < num_vertices; ii++) {
            partitioning[ii] = ii % num_processes;
			if(partitioning[ii] == process_id) hg.numMyVertices++;
        }
        
        // 1 Load model from hmetis file hgraph_name
		std::vector<std::vector<int> > hyperedges;
        std::vector<std::vector<int> > hedge_ptr;
        PRAW::load_hypergraph_from_file(hgraph_name,&hyperedges,&hedge_ptr);

		// 2 Convert hmetis format to zoltan
		hg.process_id = process_id; 
		hg.partitioning = partitioning;

		//int* total_vtx_wgts_v = (int*)calloc(num_vertices,sizeof(int));
		hg.vtxGID = (ZOLTAN_ID_TYPE*)malloc(sizeof(ZOLTAN_ID_TYPE) * hg.numMyVertices);
		hg.vtxedge_ptr = (int*)malloc(sizeof(int) * (hg.numMyVertices+1));
		int current_local_vertex = 0;
		hg.npins = 0;
		for(int from=0; from < num_vertices; from++) {
			//total_vtx_wgts_v[from] += 1; // accounts for neuron update computation 
			if(partitioning[from] == process_id) {
				hg.vtxGID[current_local_vertex] = from;
				hg.vtxedge_ptr[current_local_vertex] = hg.npins;
				current_local_vertex++;
				for(int jj=0; jj < hedge_ptr[from].size(); jj++) {
					int he_id = hedge_ptr[from][jj];
					for(int pp=0; pp < hyperedges[he_id].size(); pp++) {
						int to = hyperedges[he_id][pp];
                        if(to == from) continue; //
						hg.npins++;
						// Take into consideration synaptic computational load to weight vertices load 
						// issue with hmetis hypergraph --> a connection between two neurons may exist more than once in the model->intereconnections array
						// incoming synaptic computational balance must only count each once
						//total_vtx_wgts_v[to] += hyperedges[he_id].size()-1;
					}
					
				}
				
			}
		}
		hg.vtxedge_ptr[hg.numMyVertices] = hg.npins;
		hg.pin_GID = (ZOLTAN_ID_TYPE*)malloc(sizeof(ZOLTAN_ID_TYPE) * hg.npins);
		int counter = 0;
		for(int from=0; from < num_vertices; from++) {
			if(partitioning[from] == process_id) {
				for(int jj=0; jj < hedge_ptr[from].size(); jj++) {
					int he_id = hedge_ptr[from][jj];
					for(int pp=0; pp < hyperedges[he_id].size(); pp++) {
						if(from == hyperedges[he_id][pp]) continue; //
						hg.pin_GID[counter] = hyperedges[he_id][pp];
						counter++;
					}
				}
				
			}
		}
		
		hg.vtx_wgts = (float*)malloc(sizeof(float) * hg.numMyVertices);
		current_local_vertex = 0;
		for(int ii=0; ii < num_vertices; ii++) {
			if(partitioning[ii] == process_id) {
				//hg.vtx_wgts[current_local_vertex] = 1;
				hg.vtx_wgts[current_local_vertex] = 1;//model->null_compute ? 1 : total_vtx_wgts_v[ii];
				current_local_vertex++;
			}
		}
		//free(total_vtx_wgts_v);
		*/
		
		printf("%i: %i\n",process_id,hg.npins);
		

		// 3 Do partitioning
		int changes, numGidEntries, numLidEntries, numImport, numExport;
		ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
		int *importProcs, *importToPart, *exportProcs, *exportToPart;
		
		PRINTF("%i: Calling Zoltan\n",process_id);
		// Zoltan_LB_Set_Part_Sizes to specify the size of the partitions http://www.cs.sandia.gov/zoltan/ug_html/ug_interface_lb.html#Zoltan_LB_Set_Part_Sizes
		// Zoltan_LB_Eval to evaluate quality of decomposition http://www.cs.sandia.gov/zoltan/ug_html/ug_interface_lb.html#Zoltan_LB_Eval
		//ZOLTAN_BALANCE_EVAL balance_eval;
		//ZOLTAN_GRAPH_EVAL graph_eval;
		//ZOLTAN_HG_EVAL hg_eval;
		//Zoltan_LB_Eval(zz,1,&balance_eval,&graph_eval,&hg_eval);
			
		rc = Zoltan_LB_Partition(zz, // input (all remaining fields are output) 
			&changes,        // 1 if partitioning was changed, 0 otherwise  
			&numGidEntries,  // Number of integers used for a global ID 
			&numLidEntries,  // Number of integers used for a local ID 
			&numImport,      // Number of vertices to be sent to me 
			&importGlobalGids,  // Global IDs of vertices to be sent to me 
			&importLocalGids,   // Local IDs of vertices to be sent to me 
			&importProcs,    // Process rank for source of each incoming vertex 
			&importToPart,   // New partition for each incoming vertex 
			&numExport,      // Number of vertices I must send to other processes
			&exportGlobalGids,  // Global IDs of the vertices I must send 
			&exportLocalGids,   // Local IDs of the vertices I must send 
			&exportProcs,    // Process to which I send each of the vertices 
			&exportToPart);  // Partition to which each vertex will belong 
		
		if (rc != ZOLTAN_OK){
			printf("Zoltan error after partitioning...\n");
			MPI_Finalize();
			Zoltan_Destroy(&zz);
			exit(0);
		}
		
		// We are doing AUTO_MIGRATE
		// http://www.cs.sandia.gov/zoltan/ug_html/ug_interface_mig.html
		// after migration, partitioning within each part needs to be purged for non-owned vertices 
		// as the information may no longer be accurate
		
		idx_t* partAssign = (idx_t*)calloc(num_vertices,sizeof(idx_t)); // calloc will init to 0 all
		for (int i=0; i < num_vertices; i++) {
			if(partitioning[i] == process_id)
				partAssign[i] = process_id;		
		}
		MPI_Allreduce(partAssign,partitioning,num_vertices,MPI_LONG,MPI_MAX,MPI_COMM_WORLD);
		free(partAssign);

		Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, 
						&importProcs, &importToPart);
		Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, 
						&exportProcs, &exportToPart);	
		if(hg.vtxGID != NULL) free(hg.vtxGID);
		//if(hg.pin_GID != NULL) free(hg.pin_GID);
		if(hg.vtxedge_ptr != NULL) free(hg.vtxedge_ptr);
		if(hg.vtx_wgts != NULL) free(hg.vtx_wgts);
	}


private:

	struct Zoltan_Struct *zz;	
	CV_HGRAPH_DATA hg;
	int rc;
	
	void setZoltanParams() {
		/* General parameters */
		std::string imbalance = std::to_string(imbalance_tolerance);

		Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
		Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH");   /* partitioning method */
		Zoltan_Set_Param(zz, "HYPERGRAPH_PACKAGE", "PHG"); /* version of method */
		Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");/* global IDs are integers */
		Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");/* local IDs are integers */
		Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL"); /* export AND import lists */
		Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1"); /* use Zoltan default vertex weights */
		Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "0");/* use Zoltan default hyperedge weights */
		Zoltan_Set_Param(zz, "IMBALANCE_TOL", imbalance.c_str());/* imbalance tolerance */

		// parameters for PHG: http://www.cs.sandia.gov/zoltan/ug_html/ug_alg_phg.html
		
		Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
		Zoltan_Set_Param(zz, "CHECK_HYPERGRAPH", "0");
		Zoltan_Set_Param(zz, "PHG_EDGE_WEIGHT_OPERATION", "MAX");
		Zoltan_Set_Param(zz, "PHG_EDGE_SIZE_THRESHOLD", "1.0");
		Zoltan_Set_Param(zz, "PHG_CUT_OBJECTIVE", "CONNECTIVITY"); // CONNECTIVITY, HYPEREDGES
		Zoltan_Set_Param(zz, "PHG_COARSENING_METHOD", "AGG");
		Zoltan_Set_Param(zz, "PHG_COARSEPARTITION_METHOD", "AUTO");
		Zoltan_Set_Param(zz, "PHG_REFINEMENT_METHOD", "FM");
		Zoltan_Set_Param(zz, "PHG_REFINEMENT_QUALITY", "1");
		Zoltan_Set_Param(zz, "AUTO_MIGRATE", "TRUE");  
		
		/* Application defined query functions */
		/* To set partitioning callbacks, see: http://www.cs.sandia.gov/zoltan/ug_html/ug_query_lb.html */
		Zoltan_Set_Num_Obj_Fn(zz, cv_get_number_of_vertices, &hg);	// number of objects owned by processor
		Zoltan_Set_Obj_List_Fn(zz, cv_get_vertex_list, &hg);			// list of weights for owned objects
		Zoltan_Set_HG_Size_CS_Fn(zz, cv_get_hypergraph_size, &hg);	// dimensionality of the problem
		Zoltan_Set_HG_CS_Fn(zz, cv_get_hypergraph, &hg);				//  coordinate system (or connectivity)
		//Zoltan_Set_HG_Edge_Wts_Fn(zz, get_hyperedges_weights, &hg);  // get hyperedges weights
		//Zoltan_Set_HG_Size_Edge_Wts_Fn(zz, get_hyperedges_size_weights, &hg); // get length of hyperedges weights
		/* for migration */
		Zoltan_Set_Obj_Size_Multi_Fn(zz, cv_user_migration_multi_object_size, &hg);
		Zoltan_Set_Pack_Obj_Multi_Fn(zz, cv_user_pack_multi_obj, &hg);
		Zoltan_Set_Unpack_Obj_Multi_Fn(zz, cv_user_unpack_multi_obj, &hg);
	}
};


#endif
