#!/bin/bash

function ergodic(){
for file in ` ls $1`
do
    echo $1"/"$file
    ./hetro-heuristics -c 0.01 -n 1 -m 1 -d $1"/"$file
done }

# INIT_PATH="../randTrees_15_100"
INIT_PATH="../realTrees"
ergodic $INIT_PATH



#!/bin/bash

function ergodic(){

 org_path="../result_hetro_thin_noconstaint/result_p$3"
 mkdir $org_path
 dir_path="../result_hetro_thin_noconstaint/result_p$3/$2";
 mkdir  $dir_path;
#  0.10,1.00,5.00,2.00,3.00,8.00,
for i in {0.01,10.00}
do
    #i=10.00
    echo $i;
    
    for file in ` ls $1`
    do
        echo $1"/"$file
       ./hetro-heuristics -c 0.01 -n 1 -m 0 -d $1"/"$file   >>$dir_path/$2_c$i.txt;
    done;
done
 }

#,  ,20,30,40    1,10,50,100,
#for  p in {1000,10000}
#do
    p=1
    #        "20","40","60","80","120","160","240","360","512",  ,"4096","8192"
     for  x in {"20","40","60","80","120","160","240","360","512","1024","2048","4096","8192"}
     do 
        #x=2048
        INIT_PATH="../randTrees_thin/randTrees_"$x
        RES_PATH="test_thin_"$x
        NPR=$p
        ergodic $INIT_PATH $RES_PATH $NPR

   done
#done
