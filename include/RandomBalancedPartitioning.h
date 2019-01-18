#ifndef RANDOM_BALANCED_PARTITIONING
#define RANDOM_BALANCED_PARTITIONING

#include <vector>
#include <algorithm>
#include "Partitioning.h"

class RandomBalancedPartitioning :  public Partitioning {
public:
	
	RandomBalancedPartitioning(char* graph_name, float imbalance) : Partitioning(graph_name,imbalance) {}
	
	virtual ~RandomBalancedPartitioning() {}
	
	virtual void perform_partitioning(int num_processes,int process_id) {
		/*
		if(partitions == 1) {
			PRINTF("Partitioning not required\n");
			for(int ii=0; ii < model->population_size; ii++) 
				partitioning[ii] = 0;
		} else {
			PRINTF("Random-balanced assignment of elements...\n");
			// generate datastructure holding info of #incoming synapses for each neuron
			int* incoming_connections = (int*)calloc(model->population_size,sizeof(int));
			for(int ii=0; ii < model->population_size; ii++) {
				for(int jj=0; jj < model->interconnections_size[ii]; jj++) {
					incoming_connections[abs(model->interconnections[ii][jj])] += 1;
				}
			}
			// calculate expected weight per process (number of syns + neurons)
			float expected_weight = model->total_connections / partitions;
			// for each partition
			//	pick a neuron at random and add it to the partition
			//	substract its weight (number of synapses) from the expected weight 
			//	do until actual weight <= expected weight
			std::vector<int> actual_partition_weight(partitions,0);
			std::vector<int> unassigned_elements(model->population_size);
			for(int ii=0; ii < model->population_size; ii++) {
				unassigned_elements[ii] = ii;
			}
			std::random_shuffle(unassigned_elements.begin(),unassigned_elements.end());
			for(int ii=0; ii < partitions; ii++) {
				int next = 0;
				while(true) {
					if(next >= unassigned_elements.size()) break;
					int element = unassigned_elements[next];
					if(actual_partition_weight[ii] + incoming_connections[element] * 0.5f < expected_weight) {
						actual_partition_weight[ii] += incoming_connections[element];
						partitioning[element] = ii;
						unassigned_elements.erase(std::remove(unassigned_elements.begin(), unassigned_elements.end(), element), unassigned_elements.end());
					} else {
						if(actual_partition_weight[ii] > expected_weight) break;
						next++;
					}
				}
			}
			
			// if any leftover, assign randomly
			for(int ii=0; ii < unassigned_elements.size(); ii++) {
				partitioning[unassigned_elements[ii]] = floor(rand01() * partitions);
			}
			free(incoming_connections);
		}
		*/
	}
};


#endif
