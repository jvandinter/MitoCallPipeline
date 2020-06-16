#!/bin/bash

#$ -N annotate_chrM_to_MT
#$ -o vcf_convert_log.out
#$ -e vcf_convert_errlog.out
#$ -cwd

source mitopipeline.config

cd $MAIN_LOC

# Make sure we can look for files recursively
shopt -s globstar

for VCF in */*/*/*.vcf; do
    LOC=${VCF%.*}/
    ARRAY=(` echo $LOC | sed "s:/: :g"`)
    NEWVCF=${ARRAY[@]: -1:1}.vcf
    SAMPLE=${ARRAY[@]: -1:1}

# Convert chrM to MT in all VCFs
    if ! [[ -d "$MAIN_LOC/annotated_vcfs/$SAMPLE" ]]
    then
        awk '{gsub(/^chrM/,"MT"); print}' $MAIN_LOC$VCF > ${MAIN_LOC}annotated_input/$NEWVCF

# Create settings.config
        mkdir $MAIN_LOC/annotated_vcfs/$SAMPLE
        cd $MAIN_LOC/annotated_vcfs/$SAMPLE
        python3 /hpc/pmc_vanboxtel/projects/jip_mito/create_annotate_vcf_config.py -v ${MAIN_LOC}annotated_input/$NEWVCF -e $EMAIL -o $MAIN_LOC/annotated_vcfs/$SAMPLE

# Run IAP annotation on a file-to-file basis

        /hpc/pmc_vanboxtel/tools/IAP_v2.7.1_GRIDSS/illumina_pipeline.pl settings.config
    fi
done
