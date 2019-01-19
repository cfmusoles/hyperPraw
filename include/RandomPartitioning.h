#ifndef RANDOM_PARTITIONING
#define RANDOM_PARTITIONING

#include <vector>
#include <algorithm>
#include "Partitioning.h"

class RandomPartitioning :  public Partitioning {
public:
	
	RandomPartitioning(char* graph_name, float imbalance) : Partitioning(graph_name,imbalance) {}
	
	virtual ~RandomPartitioning() {}
	
	virtual void perform_partitioning(int num_processes,int process_id) {
		
		if(num_processes == 1) {
			PRINTF("Partitioning not required\n");
		} else {
			PRINTF("Random-balanced assignment of elements...\n");
			for(int ii=0; ii < num_vertices; ii++) {
				partitioning[ii] = (double)rand() / (double)RAND_MAX * num_processes;
				if(partitioning[ii] == num_processes) partitioning[ii] -= 1;
			}
		}
		
	}
};


#endif
