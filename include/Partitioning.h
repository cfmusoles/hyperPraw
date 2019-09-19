#ifndef PARTITIONING_H
#define PARTITIONING_H

#include <vector>
#include <cstdlib>	
#include "PRAW.h"

typedef int64_t idx_t; // needs to probably match the type in METIS

class Partitioning {
public:
	idx_t* partitioning;
	char* hgraph_name = NULL;
	int num_vertices;
	int num_hyperedges;

	Partitioning(char* graph_name, float imbalance, bool isVertexCentric = true) {
		hgraph_name = graph_name;
		imbalance_tolerance = imbalance;

		PRAW::get_hypergraph_file_header(hgraph_name,&num_vertices,&num_hyperedges);

		if(!isVertexCentric) {
			// hyperedge partitioning
			partitioning = (idx_t*)calloc(num_hyperedges,sizeof(idx_t));
		} else {
			// vertex partitioning
			partitioning = (idx_t*)calloc(num_vertices,sizeof(idx_t));
		}
	}
	virtual ~Partitioning() {
		free(partitioning);
	}
	
	virtual void perform_partitioning(int num_processes,int process_id, int* iterations) {}
	
protected:
	float imbalance_tolerance;
};

#endif
