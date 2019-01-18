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

	Partitioning(char* graph_name, float imbalance) {
		hgraph_name = graph_name;
		imbalance_tolerance = imbalance;

		int vertices;
		int hedges;
		PRAW::get_hypergraph_file_header(hgraph_name,&vertices,&hedges);

		partitioning = (idx_t*)calloc(vertices,sizeof(idx_t));
	}
	virtual ~Partitioning() {
		free(partitioning);
	}
	
	virtual void perform_partitioning(int num_processes,int process_id) {}
	
protected:
	float imbalance_tolerance;
};

#endif
