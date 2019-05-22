// Common utility functions for all types of simulations
// Helps to evaluate hypergraph partitionings

#ifndef SIMULATIONUTILS__H
#define SIMULATIONUTILS__H

namespace SimulationUtils {

    void storeSimCommunication(int* sent_communication,int process_id, int num_processes, int mode, char* experiment_name, char* graph_file, char* part_method) {
        // gather all results in node 0
        if(process_id == 0) {
            int** comm = (int**)malloc(num_processes * sizeof(int*));
            for(int ii=0; ii < num_processes; ii++) {
                comm[ii] = (int*)calloc(num_processes,sizeof(int));
            }
            memcpy(comm[0],sent_communication,sizeof(int) * num_processes);
            for(int ii=1; ii < num_processes; ii++) {
                MPI_Recv(comm[ii],num_processes,MPI_INT,ii,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            // store matrix in a file
            std::string filename = experiment_name;
            filename += "_";
            std::string graph_string = graph_file;
            filename += PRAW::getFileName(graph_string);
            filename += "_";
            filename += part_method;
            if(mode == 0) {
                filename += "_edgeSim";
            } else if (mode == 1) {
                filename += "_hedgeSim";
            }
            filename += "_comm_cost";
            char str_int[16];
            sprintf(str_int,"%i",num_processes);
            filename += "__";
            filename +=  str_int;
            bool fileexists = access(filename.c_str(), F_OK) != -1;
            FILE *fp = fopen(filename.c_str(), "w");
            if(fp == NULL) {
                printf("Error when storing comm cost into file\n");
            } else {
                for(int jj=0; jj < num_processes; jj++) {
                    for(int ii=0; ii < num_processes; ii++) {
                        fprintf(fp,"%i ",comm[jj][ii]);
                    }
                    fprintf(fp,"\n");
                }
            }
            fclose(fp);
        } else {
            MPI_Send(sent_communication,num_processes,MPI_INT,0,0,MPI_COMM_WORLD);
        }
        
        
    }
}

#endif