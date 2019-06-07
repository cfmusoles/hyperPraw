# Create ARCHER job files based on parameters passed

# imbalance: This experiment demonstrates the impact of imbalance tolerance factors over runtime to show the tradeoff between imbalance tolerance and runtime improvement
# Strategies compared:
	# prawP with different levels of imbalance
	# prawP with different amounts of fake compute time and std
# stable parameters
	# imbalance tolerance 1.1 (zoltan has 1.070 since streaming tends to reduce the imbalance significantly under the tolerance)
	# 100 max iterations
	# total edge cost communication as stopping condition
	# 0.95 tempering  refinement
	

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
#PBS -l walltime=4:00:0
# budget code
#PBS -A e582
# bandwidth probing parameters
SIZE=512
ITERATIONS=20
WINDOW=10

TEST_REPETITIONS=1
PROCESSES='''
template_5='''
EXPERIMENT_NAME='''
template_6='''
MESSAGE_SIZE='''
template_7='''
# This shifts to the directory that you submitted the job from
cd $PBS_O_WORKDIR

# bandwidth matrix creation
#renaming is necessary to avoid clashes between simultaneous jobs
ORIGINAL_BM_FILE="results_mpi_send_bandwidth_"$PROCESSES
aprun -n $PROCESSES mpi_perf $SIZE $ITERATIONS $WINDOW
for p in $(seq 1 10)
do
	FILENAME="results_mpi_send_bandwidth_"$p"_"$PROCESSES
	if [ ! -f $FILENAME ]; then
	    BM_FILE="results_mpi_send_bandwidth_"$p"_"$PROCESSES
	    break
	fi
done

mv $ORIGINAL_BM_FILE $BM_FILE

run_experiment() {
	HYPERGRAPH_FILE="$1"
	SEED="$2"
	E_SIM_STEPS="$3"
	IMBALANCE_TOLERANCE="$4"
	MEAN_WORKLOAD="$5"
	STD_WORKLOAD="$6"
	NAME=$EXPERIMENT_NAME"_imbalance_"$IMBALANCE_TOLERANCE"_std_"$STD
	aprun -n $PROCESSES hyperPraw -n $NAME -h $HYPERGRAPH_FILE -i 100 -m $IMBALANCE_TOLERANCE -p prawP -t $E_SIM_STEPS -s $SEED -k $MESSAGE_SIZE -o 2 -b $BM_FILE -W -r 950 -f $MEAN_WORKLOAD -u $STD_WORKLOAD
	sleep 1
}

for i in $(seq 1 $TEST_REPETITIONS)
do
	# average computation
	MEAN=1000
	for IMBALANCE in $(seq 1100 150 1550)
	do
		SEED=$RANDOM
		for STD in $(seq 0 200 1000)
		do
			#run_experiment "sat14_itox_vc1130.cnf.dual.hgr" $SEED 2 $IMBALANCE $MEAN $STD 
			run_experiment "2cubes_sphere.mtx.hgr" $SEED 3 $IMBALANCE $MEAN $STD 
			#run_experiment "ABACUS_shell_hd.mtx.hgr" $SEED 40 $IMBALANCE $MEAN $STD 
			#run_experiment "sparsine.mtx.hgr" $SEED 2 $IMBALANCE $MEAN $STD

			#large graphs
			#run_experiment "pdb1HYS.mtx.hgr" $SEED 1 $IMBALANCE $MEAN $STD
			#run_experiment "sat14_10pipe_q0_k.cnf.primal.hgr" $SEED 1 $IMBALANCE $MEAN $STD 
			#run_experiment "sat14_E02F22.cnf.hgr" $SEED 2 $IMBALANCE $MEAN $STD 
			#run_experiment "webbase-1M.mtx.hgr" $SEED 1 $IMBALANCE $MEAN $STD
			#run_experiment "ship_001.mtx.hgr" $SEED 1 $IMBALANCE $MEAN $STD 
			#run_experiment "sat14_atco_enc1_opt1_05_21.cnf.dual.hgr" $SEED 1 $IMBALANCE $MEAN $STD
		done

	done
	

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




