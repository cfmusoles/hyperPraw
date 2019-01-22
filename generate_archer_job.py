# Create ARCHER job files based on parameters passed

import sys
import math

#job templates
template_1 = '''#!/bin/bash --login

# name of the job
#PBS -N '''
template_2 = '''
# how many nodes
#PBS -l select='''
template_3=''':bigmem='''
template_4='''
# walltime
#PBS -l walltime=1:00:0
# budget code
#PBS -A e582

TEST_REPETITIONS=2
PROCESSES='''
template_5='''
# bandwidth probing parameters
SIZE=512
ITERATIONS=20
WINDOW=10
# comm benchmark parameters
SIM_TIME=10000
BYTES_PER_PROCESS=100
# simulation parameters
HYPERGRAPH_FILE='''
template_6='''
SIM_STEPS='''
template_7='''
EXPERIMENT_NAME='''
template_8='''
# This shifts to the directory that you submitted the job from
cd $PBS_O_WORKDIR


# bandwidth matrix creation
BM_FILE="results_mpi_send_bandwidth_"$PROCESSES
aprun -n $PROCESSES mpi_perf $SIZE $ITERATIONS $WINDOW


# test bandwidth matrix
for i in $(seq 1 $TEST_REPETITIONS)
do
	SEED=$RANDOM
    aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME -h $HYPERGRAPH_FILE -i 100 -m 1100 -p zoltan -t $SIM_STEPS -s 111 -b $BM_FILE -W
	sleep 1
	aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"_default" -h $HYPERGRAPH_FILE -i 100 -m 1100 -p prawS -t $SIM_STEPS -s 111 -b $BM_FILE
	sleep 1
	aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"_bandwidth" -h $HYPERGRAPH_FILE -i 100 -m 1100 -p prawS -t $SIM_STEPS -s 111 -b $BM_FILE -W
done

'''


if len(sys.argv) < 8:
	print("Input error: usage -> python generate_archer_job.py jobName min_processes num_experiments geometric_step big_mem[true|false] hypergraph_file simulation_steps")
	exit()

test_name = sys.argv[1]
min_processes = int(sys.argv[2])
num_experiments = int(sys.argv[3])
geometric_step = int(sys.argv[4])
big_mem = (sys.argv[5] == "true" or sys.argv[5] == "True")
hgraph_file = sys.argv[5]
sim_steps = int(sys.argv[6])

process_counts = [min_processes * geometric_step ** (n-1) for n in range (1, num_experiments+1)]
print("Generating experiments")
print(process_counts)

for p in process_counts:
	nodes = max(int(math.ceil(p / 24)),1)
	writebuffer = open("archer_job_" + test_name + "_" + str(p) + ".sh",'w')
	writebuffer.write(template_1 + test_name)
	writebuffer.write(template_2 + str(nodes))
	writebuffer.write(template_3 + str(big_mem).lower())
	writebuffer.write(template_4 + str(p))
	writebuffer.write(template_5 + hgraph_file)
	writebuffer.write(template_6 + str(sim_steps))
	writebuffer.write(template_7 + test_name)
	writebuffer.write(template_8)




