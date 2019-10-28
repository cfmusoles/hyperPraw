# Create ARCHER job files based on parameters passed

# parallel: This experiment compares baseline sequential hypergraph partitioning (Alistairh) to a parallelised version with multiple substreams
	

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
#PBS -l walltime=6:00:0
# budget code
#PBS -A e582

PARTITIONS='''
template_5='''
EXPERIMENT_NAME='''
template_6='''
MESSAGE_SIZE='''
template_7='''
# This shifts to the directory that you submitted the job from
cd $PBS_O_WORKDIR

run_experiment() {
	HYPERGRAPH_FILE="$1"
	SEED="$2"
	MAX_PROCESSES="$3"
	PART="$4"
	WINDOW_SIZE="$5"
	E_SIM_STEPS=0
	H_SIM_STEPS=0
	GRAPH_STREAM="inverted_"$HYPERGRAPH_FILE

	aprun -n $PARTITIONS hyperPraw -n $EXPERIMENT_NAME"_"$PART"_"$MAX_PROCESSES"_"$WINDOW_SIZE -h $HYPERGRAPH_FILE -i 100 -m 1200 -p $PART -t $E_SIM_STEPS -x $H_SIM_STEPS -s $SEED -k $MESSAGE_SIZE -e $GRAPH_STREAM -P -K $MAX_PROCESSES -g $WINDOW_SIZE
	sleep 1
}

# baseline strategy only run once
SEED=$RANDOM
run_experiment "small.hgr" $SEED 1 "baselineSequential" 1
run_experiment "shuffled_2cubes_sphere.mtx.hgr" $SEED 1 "baselineSequential" 1
#run_experiment "small_dense_uniform.hgr" $SEED 1 "baselineSequential" 1
#run_experiment "small_dense_powerlaw.hgr" $SEED 1 "baselineSequential" 1
#run_experiment "large_sparse_uniform.hgr" $SEED 1 "baselineSequential" 1
#run_experiment "large_sparse_powerlaw.hgr" $SEED 1 "baselineSequential" 1
#run_experiment "2cubes_sphere.mtx.hgr" $SEED 1 "baselineSequential" 1
#run_experiment "ABACUS_shell_hd.mtx.hgr" $SEED 1 "baselineSequential" 1
#run_experiment "sparsine.mtx.hgr" $SEED 1 "baselineSequential" 1

# run parallel versions
NUM_PARALLEL_EXPERIMENTS=5
PROCESSES="3"
FACTOR="2"
for p in $(seq 1 $NUM_PARALLEL_EXPERIMENTS)
do
	SEED=$RANDOM
	#synthetic graphs
	#run_experiment "small_dense_uniform.hgr" $SEED $PROCESSES "parallelVertex" 1
	run_experiment "small.hgr" $SEED $PROCESSES "parallelVertex" 1
	run_experiment "shuffled_2cubes_sphere.mtx.hgr" $SEED $PROCESSES "parallelVertex" 1
	#run_experiment "small_dense_powerlaw.hgr" $SEED $PROCESSES "parallelVertex" 1
	#run_experiment "large_sparse_uniform.hgr" $SEED $PROCESSES "parallelVertex" 1
	#run_experiment "large_sparse_powerlaw.hgr" $SEED $PROCESSES "parallelVertex" 1
	#run_experiment "2cubes_sphere.mtx.hgr" $SEED $PROCESSES "parallelVertex" 1
	#run_experiment "ABACUS_shell_hd.mtx.hgr" $SEED $PROCESSES "parallelVertex" 1
	#run_experiment "sparsine.mtx.hgr" $SEED $PROCESSES "parallelVertex" 1

	PROCESSES=$(($PROCESSES * $FACTOR))
	
done



'''


if len(sys.argv) < 7:
	print("Input error: usage -> python generate_archer_job.py jobName min_processes num_experiments geometric_step big_mem[true|false] simulation_steps")
	exit()

test_name = sys.argv[1]
min_processes = int(sys.argv[2])
num_experiments = int(sys.argv[3])
geometric_step = int(sys.argv[4])
big_mem = (sys.argv[5] == "true" or sys.argv[5] == "True")
message_size = int(sys.argv[6])

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
	writebuffer.write(template_5 + test_name)
	writebuffer.write(template_6 + str(message_size))
	writebuffer.write(template_7)




