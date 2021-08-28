#!/bin/bash

function ergodic(){
dir_path="../result/$2";
# test_lessprocessor
# org_path="../result/test_lessprocessor_p$3"
 #mkdir $org_path
 #dir_path="../result/test_lessprocessor_p$3/$2";
mkdir  $dir_path;

 for i in {0.01,0.10,1.00,2.00,3.00}
 do
    echo $i;
    # i=2.00
    for file in ` ls $1`
    do
        echo $1"/"$file

        ./call-heuristics -s -l  -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt;
        ./select-heuristics -s -l -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt;
        ./call-heuristics-main -s -l  -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt

        #  ./call-heuristics -s -l  -q  -c $i -n 50 -d   $1"/"$file ;
        # ./select-heuristics -s -l -q  -c $i -n 50 -d   $1"/"$file ;
        # ./call-heuristics-main -s -l  -q  -c $i -n 50 -d   $1"/"$file
    done;
 done
 }

#10,20,30,
 #for  p in {40,50}
 #do
      p=1
     #
     #x=20001_60000
     for x in {"20","40","80","256","15_100","101_500","501_1000","1001_2000","2001_10000","10001_20000","20001_60000"}
     do 
         INIT_PATH="../randTrees/randTrees_"$x
         RES_PATH="test_"$x
         NPR=$p
         ergodic $INIT_PATH $RES_PATH $NPR
#         # echo $RES_PATH
#         # echo $INIT_PATH
     done
 #done
#INIT_PATH="../randTrees_15_100"
#RES_PATH="testtest_15_100"
#NPR=1
#ergodic $INIT_PATH $RES_PATH $NPR


