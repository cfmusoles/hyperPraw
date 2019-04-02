#!/bin/bash

APP_NAME="hyperPraw"
EXPERIMENT_NAME="$1"
MIN_PROCESSORS="$2"
NUM_EXPERIMENTS="$3"
BIG_MEM_NODE="$4"
SIM_STEPS="$5"
MESSAGE_SIZE="$6"
EXPERIMENT_TYPE="$7"
#iterations of each experiment (number of jobs)
REPETITIONS="$8"
HGRAPHS_DIR="$9"

GEOMETRIC_STEP=4

#create working directory
mkdir $WORK_DIR/$EXPERIMENT_NAME

#generate job files
python $EXPERIMENT_TYPE $EXPERIMENT_NAME $MIN_PROCESSORS $NUM_EXPERIMENTS $GEOMETRIC_STEP $BIG_MEM_NODE $SIM_STEPS $MESSAGE_SIZE
sleep 1

#geometric progression
for p in $(seq 1 $NUM_EXPERIMENTS)
do
	PROCESS_COUNT=$(($MIN_PROCESSORS * $GEOMETRIC_STEP ** ($p-1)))
	FILENAME="archer_job_"$EXPERIMENT_NAME"_"$PROCESS_COUNT".sh"
	cp $FILENAME $WORK_DIR/$EXPERIMENT_NAME/
	rm $FILENAME
done

#copy necessary files
cp $APP_NAME $WORK_DIR/$EXPERIMENT_NAME/
cp ../mpi_perf/mpi_perf $WORK_DIR/$EXPERIMENT_NAME/
cp resources/$HGRAPHS_DIR/* $WORK_DIR/$EXPERIMENT_NAME/
chmod +x archer_retrieve_results.sh
cp archer_retrieve_results.sh $WORK_DIR/$EXPERIMENT_NAME/

#submit job files
cd $WORK_DIR/$EXPERIMENT_NAME

# untar hgraph files
for f in *.tar.gz
do
  tar zxvf "$f" -C ./
done
rm *.tar.gz

#geometric progression
for p in $(seq 1 $NUM_EXPERIMENTS)
do
	PROCESS_COUNT=$(($MIN_PROCESSORS * $GEOMETRIC_STEP ** ($p-1)))
	FILENAME="archer_job_"$EXPERIMENT_NAME"_"$PROCESS_COUNT".sh"
	echo "Launching job: "$FILENAME
	for a in $(seq 1 $REPETITIONS)
	do
		qsub $FILENAME
		#qsub -q short $FILENAME
		rm $FILENAME
	done
done

