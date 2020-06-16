#! /bin/bash

#$ -N mitopipeline_startup
#$ -l h_rt=12:00:00
#$ -l h_vmem=8G
#$ -cwd
#$ -m e

set -o errexit

# Check existence of CONFIG file in root dir
source mitopipeline.config

##Load required modules
module load Java/1.8.0_60
module load python/3.6.1

cd $HG38_MAP_LOC

# Make sure we can look for files recursively
shopt -s globstar

FLOW=(q)
for BAMLOC in **/*.bam; do

# Remove BAM files created by structural variant calling
    if ! [[ "${BAMLOC}" =~ "gridss" ]]; then

# Only pick BAM files located in /mapping/
        if [[ "${BAMLOC}" =~ "mapping" ]]; then
            LOC=${BAMLOC%.*}/
            ARRAY=(` echo $LOC | sed "s:/: :g"`)
            BAM=${ARRAY[@]: -1:1}.bam

# Grab Sample ID from BAM file
            if ! [[ "$BAM" =~ "sorted" ]]; then
                if [[ "$BAM" =~ "_dedup" ]]; then
                    SAMPLE=`echo ${ARRAY[@]: -1:1} | sed -e "s/_dedup//"`
                    if [[ "$SAMPLE" =~ ".realigned" ]]; then
                        SAMPLE=`echo $SAMPLE | sed -e "s/.realigned//"`
                        if ! [[ "${FLOW[*]}" =~ "${SAMPLE}" ]]; then 
                            FLOW+=($SAMPLE)
                            
                        fi

                    else
                        SAMPLE=$SAMPLE
                        if ! [[ "${FLOW[*]}" =~ "${SAMPLE}" ]]; then 
                            FLOW+=($SAMPLE)
                        fi
                    fi

                elif [[ "$BAM" =~ ".realigned" ]]; then
                    SAMPLE=`echo ${ARRAY[@]: -1:1} | sed -e "s/.realigned//"`
                    if ! [[ "${FLOW[*]}" =~ "${SAMPLE}" ]]; then 
                        FLOW+=($SAMPLE)
                    fi

                else
                    SAMPLE=`echo ${ARRAY[@]: -1:1}`
                    if ! [[ "${FLOW[*]}" =~ "${SAMPLE}" ]]; then 
                        FLOW+=($SAMPLE)
                    fi
                fi
            for QC in **/*_dedup_WGSMetrics.txt; do
                if [[ "$QC" =~ "$SAMPLE" ]]; then
                    QCLOC=$MAINLOC$QC
                    cd $MAIN_LOC$INPUT_LOC
# Finish input.json with specified files
                    python3 finish_json_slurm.py -b $BAMLOC -n $SAMPLE -q $QCLOC -l $HG38_MAP_LOC -wl $WORK_LOC
                    cd $HG38_MAP_LOC
                fi
            done


            fi
        fi
    fi
done

cd $MAIN_LOC$WDL_LOC

zip modules.zip MitochondriaPipeline.wdl AlignmentPipeline.wdl AlignAndCall.wdl

cd $MAIN_LOC

## RUN cromwell for each specified sample ; skip first variable ( q )
for NAME in ${FLOW[@]:1}; do
    if ! [[ -d "$MAIN_LOC$WORK_LOC/final_output/$NAME" ]]
        then
            qsub mito_pipeline_run.sh $MAIN_LOC $CROMWELL_LOC $WDL_LOC $INPUT_LOC $NAME
            sleep 10
done
