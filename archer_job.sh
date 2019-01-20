#!/bin/bash --login

# name of the job
#PBS -N simulation_job  
# how many nodes
#PBS -l select=8
# walltime
#PBS -l walltime=0:20:0
# budget code
#PBS -A e582  

PROCESSES=192
TEST_REPETITIONS=2

# bandwidth probing parameters
SIZE=512
ITERATIONS=20
WINDOW=10
# comm benchmark parameters
SIM_TIME=10000
BYTES_PER_PROCESS=100
# simulation parameters
HYPERGRAPH_FILE="2D_54019_highK.mtx.hgr" #"sparsine.mtx.hgr"

# This shifts to the directory that you submitted the job from
cd $PBS_O_WORKDIR


# bandwidth matrix creation
BM_FILE="results_mpi_send_bandwidth_"$PROCESSES
aprun -n $PROCESSES mpi_perf $SIZE $ITERATIONS $WINDOW


# test bandwidth matrix
for i in $(seq 1 $TEST_REPETITIONS)
do
	SEED=$RANDOM
        #aprun -n $PROCESSES comm_benchmark $SIM_TIME $BYTES_PER_PROCESS benchmark_default
        #sleep 1
        #aprun -n $PROCESSES comm_benchmark $SIM_TIME $BYTES_PER_PROCESS benchmark_bm $BM_FILE
        #sleep 1	
	aprun -n $PROCESSES hyperPraw -n zoltan -h $HYPERGRAPH_FILE -i 100 -m 1100 -p zoltan -t 50 -s 111 -b $BM_FILE -W
	sleep 1
	aprun -n $PROCESSES hyperPraw -n praw_default -h $HYPERGRAPH_FILE -i 100 -m 1100 -p praw -t 50 -s 111 -b $BM_FILE
	sleep 1
	aprun -n $PROCESSES hyperPraw -n praw_bandwidth -h $HYPERGRAPH_FILE -i 100 -m 1100 -p praw -t 50 -s 111 -b $BM_FILE -W	
	sleep 1
	aprun -n $PROCESSES hyperPraw -n random -h $HYPERGRAPH_FILE -i 100 -m 1100 -p random -t 50 -s 111 -b $BM_FILE -W
	sleep 1
done




