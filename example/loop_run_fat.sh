#!/bin/bash

function ergodic(){
#dir_path="../result_fat_noconstraint/result_p1/$2";
 #test_lessprocessor
org_path="../result_fat_noconstraint/test_lessprocessor_p$3"
mkdir $org_path
dir_path="../result_fat_noconstraint/test_lessprocessor_p$3/$2";
mkdir  $dir_path;
#
for i in {0.01,0.10,1.00,2.00,3.00,5.00,8.00,10.00}
do
     #i=2.00
    echo $i;
    
    for file in ` ls $1`
    do
        echo $1"/"$file

        ./call-heuristics -s -l  -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt;
        ./select-heuristics -s -l -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt;
        ./call-heuristics-main -s -l  -q  -c $i -n $3 -m 0 -d   $1"/"$file >>$dir_path/$2_c$i.txt
        # ./memory-heuristics-test -s -l -q -c $i -n $3 -m 2 -d   $1"/"$file >>$dir_path/$2_c$i.txt
        #  ./call-heuristics -s -l  -q  -c $i -n 50 -d   $1"/"$file ;
        # ./select-heuristics -s -l -q  -c $i -n 50 -d   $1"/"$file ;
        # ./call-heuristics-main -s -l  -q  -c $i -n 50 -d   $1"/"$file
    done;
done
 }

#,
for  p in {10,20,30,40,50}
do
    #p=1
    # for x in {"15_100","101_500","501_1000","1001_2000","2001_10000","10001_20000","20001_60000","20","40","80","256"}
    # ,
   for  x in {"20","40","60","80","120","160","240","360","512","1024","2048","4096","8192"}
   do 
        #x=8192
        INIT_PATH="../randTrees_fat/randTrees_"$x
        RES_PATH="test_fat_"$x
        NPR=$p
        ergodic $INIT_PATH $RES_PATH $NPR
        # echo $RES_PATH
        # echo $INIT_PATH
    done
done
# INIT_PATH="../randTrees_fat/randTrees_120"
# RES_PATH="test_fat_120"
# NPR=1
# ergodic $INIT_PATH $RES_PATH $NPR


