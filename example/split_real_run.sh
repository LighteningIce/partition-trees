#!/bin/bash

function ergodic(){

for file in ` ls $1`
do
    echo $1"/"$file
    ./split_real_tree  -d $1"/"$file>>test3.txt

done;

 }

INIT_PATH="../real_tree_every"
#RES_PATH="test_real_tree_every"

# INIT_PATH="../randTrees/randTrees_20001_60000"
# RES_PATH="randTrees/"
NPR=1
#INIT_PATH="../real_trees"
ergodic $INIT_PATH


