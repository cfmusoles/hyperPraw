# Create ARCHER job files based on parameters passed

# COMM_COST_EXPERIMENT: This experiment demonstrates the effectiveness of architecture aware partitioning, evaluating different bandwidth to comm cost mappings
# Strategies compared:
	# prawS default: not using bandwidth info (baseline)
	# prawS 0_1: uses bandwidth info as a 0 - 1 mapping to comm cost
	# prawS proportional: uses bandwidth info with a proportional mapping to comm cost (from 0 to ratio max / min bandwidth)
# stable parameters
	# hedgeEdge as stopping condition
	# imbalance tolerance 1.1
	# 100 max iterations
	
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
#PBS -l walltime=0:20:0
# budget code
#PBS -A e582
# bandwidth probing parameters
SIZE=512
ITERATIONS=20
WINDOW=10

TEST_REPETITIONS=2
PROCESSES='''
template_5='''
# simulation parameters
SIM_STEPS='''
template_6='''
EXPERIMENT_NAME='''
template_7='''
MESSAGE_SIZE='''
template_8='''
# This shifts to the directory that you submitted the job from
cd $PBS_O_WORKDIR

# bandwidth matrix creation
BM_FILE="results_mpi_send_bandwidth_"$PROCESSES
aprun -n $PROCESSES mpi_perf $SIZE $ITERATIONS $WINDOW

run_experiment() {
	HYPERGRAPH_FILE="$1"
	SEED="$2"
	# best default strategy from experiment 1
	aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"_default_hedgeEdge" -h $HYPERGRAPH_FILE -i 100 -m 1100 -p prawS -t $SIM_STEPS -s $SEED -k $MESSAGE_SIZE -o 1 -b $BM_FILE
	sleep 1
	# bandwidth mapped to 0 - 1 default stopping condition
	aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"_bandwidth_0_1_hedgeEdge" -h $HYPERGRAPH_FILE -i 100 -m 1100 -p prawS -t $SIM_STEPS -s $SEED -k $MESSAGE_SIZE -o 1 -b $BM_FILE -W -c 0
	sleep 1
	# bandwidth mapped to proportional default stopping condition
	aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"_bandwidth_proportional_hedgeEdge" -h $HYPERGRAPH_FILE -i 100 -m 1100 -p prawS -t $SIM_STEPS -s $SEED -k $MESSAGE_SIZE -o 1 -b $BM_FILE -W -c 1
	sleep 1
	#aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"_zoltan" -h $HYPERGRAPH_FILE -i 100 -m 1050 -p zoltan -t $SIM_STEPS -s $SEED -k $MESSAGE_SIZE -b $BM_FILE -c 1
	#sleep 1
}

for i in $(seq 1 $TEST_REPETITIONS)
do
	SEED=$RANDOM
	#run_experiment "sat14_E02F20.cnf.hgr" $SEED
	run_experiment "crashbasis.mtx.hgr" $SEED
	run_experiment "sat14_aaai10-planning-ipc5-pathways-17-step21.cnf.dual.hgr" $SEED
	run_experiment "sparsine.mtx.hgr" $SEED
	run_experiment "venkat01.mtx.hgr" $SEED
done

'''


if len(sys.argv) < 8:
	print("Input error: usage -> python generate_archer_job.py jobName min_processes num_experiments geometric_step big_mem[true|false] simulation_steps")
	exit()

test_name = sys.argv[1]
min_processes = int(sys.argv[2])
num_experiments = int(sys.argv[3])
geometric_step = int(sys.argv[4])
big_mem = (sys.argv[5] == "true" or sys.argv[5] == "True")
sim_steps = int(sys.argv[6])
message_size = int(sys.argv[7])

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
	writebuffer.write(template_5 + str(sim_steps))
	writebuffer.write(template_6 + test_name)
	writebuffer.write(template_7 + str(message_size))
	writebuffer.write(template_8)



