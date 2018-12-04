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
TEST_REPETITIONS=3

# bandwidth probing parameters
SIZE=4096
ITERATIONS=20
WINDOW=10
# comm benchmark parameters
SIM_TIME=15000
BYTES_PER_PROCESS=100
# simulation parameters
HYPERGRAPH_FILE="2D_54019_highK.mtx.hgr"
COMM_PATTERN="nbx"

# This shifts to the directory that you submitted the job from
cd $PBS_O_WORKDIR


# bandwidth matrix creation
BM_FILE="results_mpi_send_bandwidth_"$PROCESSES
aprun -n $PROCESSES mpi_perf $SIZE $ITERATIONS $WINDOW


# test bandwidth matrix
for i in $(seq 1 $TEST_REPETITIONS)
do
	SEED=$RANDOM
        aprun -n $PROCESSES comm_benchmark $SIM_TIME $BYTES_PER_PROCESS benchmark_default
        sleep 1
        aprun -n $PROCESSES comm_benchmark $SIM_TIME $BYTES_PER_PROCESS benchmark_bm $BM_FILE
        sleep 1
	aprun -n $PROCESSES distSim -n nullcompute_zoltanFile -c $COMM_PATTERN -p "zoltanFile" -s $SEED -k 1000 -f 1000 -t 500 -h $HYPERGRAPH_FILE -i 24 -N
	aprun -n $PROCESSES distSim -n nullcompute_praw_default -c $COMM_PATTERN -p "prawS" -s $SEED -k 1000 -f 1000 -t 500 -h $HYPERGRAPH_FILE -i 24 -N
	aprun -n $PROCESSES distSim -n nullcompute_praw_bandwidth -c $COMM_PATTERN -p "prawS" -s $SEED -k 1000 -f 1000 -t 500 -h $HYPERGRAPH_FILE -i 24 -N -b $BM_FILE
	sleep 1
done




