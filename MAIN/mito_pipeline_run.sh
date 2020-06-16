#! /bin/bash

#$ -S /bin/bash
#$ -N mito_pipeline_run
#$ -cwd
#$ -l h_rt=1:00:00
#$ -l h_vmem=20G
#$ -l tmpspace=20G
#$ -m beas

module load autossh/1.4e
autossh -M 13306 horus -L 3306:localhost:3306 -N -f &
autossh -M 18009 horus -R 8007:localhost:8007 -N -f &
pid=$!
sleep 5

# Take previous output
MAIN_LOC=$1
CROMWELL_LOC=$2
WDL_LOC=$3
INPUT_LOC=$4
NAME=$5

# Check whether sample has been successfully run already 
java -Xms4G -Xmx16G -Dconfig.file=${MAIN_LOC}application.conf.31 -jar ${CROMWELL_LOC}cromwell-29.jar \
  run ${MAIN_LOC}${WDL_LOC}MitochondriaPipeline.wdl \
  -i ${MAIN_LOC}${INPUT_LOC}mitochondrial_workflow_inputs_${NAME}.json \
  -o ${MAIN_LOC}${INPUT_LOC}mitochondrial_workflow_options_${NAME}.json \
  -p ${MAIN_LOC}${WDL_LOC}modules.zip
fi
