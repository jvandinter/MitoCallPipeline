#!/bin/bash

source mitopipeline.config

HOME_LOC=$MAIN_LOC$WORK_LOC

cd $HOME_LOC
for FOLDER in final_output/*; do
    echo $HOME_LOC$FOLDER
    cd $HOME_LOC$FOLDER
    find . -mindepth 2 -type f -print -exec mv {} . \;
    cd $HOME_LOC
done
