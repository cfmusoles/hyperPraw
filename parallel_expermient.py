# Create ARCHER job files based on parameters passed

# parallel: This experiment compares praw sequential and parallel implementations. Should be able to compare quality of partitions and runtime impact
# Strategies compared:
	# prawS: serial architecture aware
	# prawP: parallel architecture aware
	# zoltan: multilevel partitioning benchmark
# stable parameters
	# imbalance tolerance 1.2 (zoltan has 1.070 since streaming tends to reduce the imbalance significantly under the tolerance)
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
#PBS -l walltime=6:00:0
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
	H_SIM_STEPS="$4"
	GRAPH_STREAM="inverted_"$HYPERGRAPH_FILE

	#aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"_sequential" -h $HYPERGRAPH_FILE -i 100 -m 1200 -p sequential -t $E_SIM_STEPS -x $H_SIM_STEPS -s $SEED -k $MESSAGE_SIZE -o 2 -b $BM_FILE -W -c 0 -r 950
	#sleep 1
	#aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"__parallelHyperedge" -h $HYPERGRAPH_FILE -i 100 -m 1200 -p parallelHyperedge -t $E_SIM_STEPS -x $H_SIM_STEPS -s $SEED -k $MESSAGE_SIZE -o 2 -b $BM_FILE -W -c 0 -r 950
	#sleep 1
	
	aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"_zoltanVertex" -h $HYPERGRAPH_FILE -i 100 -m 1200 -p zoltanVertex -t $E_SIM_STEPS -x $H_SIM_STEPS -s $SEED -k $MESSAGE_SIZE -b $BM_FILE
	sleep 1
	#aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"_zoltanHyperedge" -h $HYPERGRAPH_FILE -i 100 -m 1200 -p zoltanHyperedge -t $E_SIM_STEPS -x $H_SIM_STEPS -s $SEED -k $MESSAGE_SIZE -b $BM_FILE
	#sleep 1

	#aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"_sequential_default" -h $HYPERGRAPH_FILE -i 100 -m 1200 -p sequentialVertex -t $E_SIM_STEPS -x $H_SIM_STEPS -s $SEED -k $MESSAGE_SIZE -H -b $BM_FILE
	#sleep 1
	#aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"_sequential_bandwidth" -h $HYPERGRAPH_FILE -i 100 -m 1200 -p sequentialVertex -t $E_SIM_STEPS -x $H_SIM_STEPS -s $SEED -k $MESSAGE_SIZE -H -W -b $BM_FILE
	#sleep 1

	aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"_simpleParallelVertex" -h $HYPERGRAPH_FILE -i 100 -m 1200 -p simpleParallelVertex -t $E_SIM_STEPS -x $H_SIM_STEPS -s $SEED -k $MESSAGE_SIZE -W -H -b $BM_FILE
	sleep 1
	aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"_rHDRF_bandwidth" -h $HYPERGRAPH_FILE -i 100 -m 1200 -p rHDRF -t $E_SIM_STEPS -x $H_SIM_STEPS -s $SEED -k $MESSAGE_SIZE -W -H -b $BM_FILE -e $GRAPH_STREAM
	sleep 1
	aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"_rHDRF_default" -h $HYPERGRAPH_FILE -i 100 -m 1200 -p rHDRF -t $E_SIM_STEPS -x $H_SIM_STEPS -s $SEED -k $MESSAGE_SIZE -H -b $BM_FILE -e $GRAPH_STREAM
	sleep 1
}

for i in $(seq 1 $TEST_REPETITIONS)
do
	SEED=$RANDOM

	#large graphs (> 3M)
	#run_experiment "sat14_11pipe_k.cnf.dual.hgr" $SEED 0 0
	#run_experiment "sat14_atco_enc3_opt1_04_50.cnf.hgr" $SEED 0 0
	#run_experiment "sat14_blocks-blocks-37-1.130-NOTKNOWN.cnf.dual.hgr" $SEED 0 0
	#run_experiment "sat14_SAT_dat.k100-24_1_rule_1.cnf.dual.hgr" $SEED 0 0
	#run_experiment "sat14_q_query_3_L200_coli.sat.cnf.dual.hgr" $SEED 0 0

	#medium graphs (> 800K)
	run_experiment "atmosmodj.mtx.hgr" $SEED 5 5
	run_experiment "kkt_power.mtx.hgr" $SEED 5 5
	run_experiment "sat14_velev-vliw-uns-2.0-uq5.cnf.dual.hgr" $SEED 5 5

	#run_experiment "sat14_itox_vc1130.cnf.dual.hgr" $SEED 1 10
	#run_experiment "2cubes_sphere.mtx.hgr" $SEED 4 3
	#run_experiment "ABACUS_shell_hd.mtx.hgr" $SEED 50 40
	#run_experiment "sparsine.mtx.hgr" $SEED 2 2
	
	#large graphs
	#run_experiment "pdb1HYS.mtx.hgr" $SEED 2 2 #
	#run_experiment "sat14_10pipe_q0_k.cnf.primal.hgr" $SEED 2 2
	#run_experiment "sat14_E02F22.cnf.hgr" $SEED 2 2
	#run_experiment "webbase-1M.mtx.hgr" $SEED 2 2
	#run_experiment "ship_001.mtx.hgr" $SEED 2 2
	#run_experiment "sat14_atco_enc1_opt1_05_21.cnf.dual.hgr" $SEED 2 2
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




