# Create ARCHER job files based on parameters passed

# parallel: This experiment characterises hyper PRAW streaming partitioning
# Factors to characterise:
#	impact of window-based streaming (updating only in batches)
#	impact of local updates only (only send remote updates when a pin has never been seen before, otherwise only update local datastructure)

# fixed prameters
# 	lambda = 0.5
# 	impalance tolerance = 1.2
#	staggered start



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
#PBS -l walltime=24:00:0
# budget code
#PBS -A e582

REPETITIONS=1
SIMS_PER_TRIAL=3
PROCESSES='''
template_5='''
EXPERIMENT_NAME='''
template_6='''
MESSAGE_SIZE='''
template_7='''
# This shifts to the directory that you submitted the job from
cd $PBS_O_WORKDIR

# bandwidth matrix creation
# bandwidth probing parameters
SIZE=512
ITERATIONS=20
WINDOW=10
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

	# run parallel versions
	NUM_PARALLEL_EXPERIMENTS=3
	MAX_PROCESSES="1"
	FACTOR="8"
	for p in $(seq 1 $NUM_PARALLEL_EXPERIMENTS)
	do
		# window based streaming tests
		aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"_hyperPraw_bandwidth_w1_"$MAX_PROCESSES -h $HYPERGRAPH_FILE -i 100 -m 1200 -p hyperPrawVertex -t $E_SIM_STEPS -x $H_SIM_STEPS -s $SEED -k $MESSAGE_SIZE -e $GRAPH_STREAM -P -K $MAX_PROCESSES -g 1 -b $BM_FILE -W -r 500 -q $SIMS_PER_TRIAL -H
		sleep 1
		aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"_hyperPraw_bandwidth_w3_"$MAX_PROCESSES -h $HYPERGRAPH_FILE -i 100 -m 1200 -p hyperPrawVertex -t $E_SIM_STEPS -x $H_SIM_STEPS -s $SEED -k $MESSAGE_SIZE -e $GRAPH_STREAM -P -K $MAX_PROCESSES -g 3 -b $BM_FILE -W -r 500 -q $SIMS_PER_TRIAL -H
		sleep 1
		aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"_hyperPraw_bandwidth_w10_"$MAX_PROCESSES -h $HYPERGRAPH_FILE -i 100 -m 1200 -p hyperPrawVertex -t $E_SIM_STEPS -x $H_SIM_STEPS -s $SEED -k $MESSAGE_SIZE -e $GRAPH_STREAM -P -K $MAX_PROCESSES -g 10 -b $BM_FILE -W -r 500 -q $SIMS_PER_TRIAL -H
		sleep 1
		#aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"_hyperPraw_bandwidth_w50_"$MAX_PROCESSES -h $HYPERGRAPH_FILE -i 100 -m 1200 -p hyperPrawVertex -t $E_SIM_STEPS -x $H_SIM_STEPS -s $SEED -k $MESSAGE_SIZE -e $GRAPH_STREAM -P -K $MAX_PROCESSES -g 50 -b $BM_FILE -W -r 500 -q $SIMS_PER_TRIAL -H
		#sleep 1

		# only local updates
		aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"_hyperPraw_bandwidth_w1_local_"$MAX_PROCESSES -h $HYPERGRAPH_FILE -i 100 -m 1200 -p hyperPrawVertex -t $E_SIM_STEPS -x $H_SIM_STEPS -s $SEED -k $MESSAGE_SIZE -e $GRAPH_STREAM -P -K $MAX_PROCESSES -g 1 -b $BM_FILE -W -r 500 -q $SIMS_PER_TRIAL -H -L
		sleep 1
		aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"_hyperPraw_bandwidth_w3_local_"$MAX_PROCESSES -h $HYPERGRAPH_FILE -i 100 -m 1200 -p hyperPrawVertex -t $E_SIM_STEPS -x $H_SIM_STEPS -s $SEED -k $MESSAGE_SIZE -e $GRAPH_STREAM -P -K $MAX_PROCESSES -g 3 -b $BM_FILE -W -r 500 -q $SIMS_PER_TRIAL -H -L
		sleep 1
		aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"_hyperPraw_bandwidth_w10_local_"$MAX_PROCESSES -h $HYPERGRAPH_FILE -i 100 -m 1200 -p hyperPrawVertex -t $E_SIM_STEPS -x $H_SIM_STEPS -s $SEED -k $MESSAGE_SIZE -e $GRAPH_STREAM -P -K $MAX_PROCESSES -g 10 -b $BM_FILE -W -r 500 -q $SIMS_PER_TRIAL -H -L
		sleep 1
		#aprun -n $PROCESSES hyperPraw -n $EXPERIMENT_NAME"_hyperPraw_bandwidth_w50_local_"$MAX_PROCESSES -h $HYPERGRAPH_FILE -i 100 -m 1200 -p hyperPrawVertex -t $E_SIM_STEPS -x $H_SIM_STEPS -s $SEED -k $MESSAGE_SIZE -e $GRAPH_STREAM -P -K $MAX_PROCESSES -g 50 -b $BM_FILE -W -r 500 -q $SIMS_PER_TRIAL -H -L
		#sleep 1

		MAX_PROCESSES=$(($MAX_PROCESSES * $FACTOR))
	done
	
}


for p in $(seq 1 $REPETITIONS)
do
	SEED=$RANDOM
	#synthetic graphs
	run_experiment "small_uniform_dense_c96.hgr" $SEED 2 30
	run_experiment "small_uniform_sparse_c96.hgr" $SEED 18 30
	run_experiment "large_uniform_sparse_c96.hgr" $SEED 11 30
	run_experiment "large_powerlaw_sparse_c96.hgr" $SEED 15 30
	run_experiment "small_powerlaw_dense_c96.hgr" $SEED 2 30
	run_experiment "small_uniform_dense_c192.hgr" $SEED 2 30
	run_experiment "small_uniform_sparse_c48.hgr" $SEED 19 30
	#run_experiment "huge_uniform_dense_c96.hgr" $SEED 1 5
	
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




