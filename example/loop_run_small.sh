#!/bin/bash

function ergodic(){
#dir_path="../result_noconstraint/$2";
# test_lessprocessor
 org_path="../result_outMdeg/result_p$3"
 mkdir $org_path
 dir_path="../result_outMdeg/result_p$3/$2";
 mkdir  $dir_path;
#,2.00,3.00,8.00
 for i in {0.01,0.10,1.00,5.00,10.00}
 do

    #i=1.00
     echo $i;
    for file in ` ls $1`
    do
        echo $1"/"$file


        ./call-heuristics -s -l  -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt;
       ./call-heuristics -i -l  -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt;
        ./call-heuristics -a -l   -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt;
        ./select-heuristics -l  -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt;
        ./call-heuristics-main -l   -q  -c $i -n $3 -m 0 -d   $1"/"$file >>$dir_path/$2_c$i.txt   

        #./call-heuristics -s -l  -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt;
        #./select-heuristics -s -l -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt;
        #./call-heuristics-main -s -l  -q  -c $i -n $3 -m 2 -d   $1"/"$file >>$dir_path/$2_c$i.txt

        #  ./call-heuristics -s -l  -q  -c $i -n 50 -d   $1"/"$file ;
        # ./select-heuristics -s -l -q  -c $i -n 50 -d   $1"/"$file ;
        # ./call-heuristics-main -s -l  -q  -c $i -n 50 -d   $1"/"$file
    done;
 done
 }

#

for  p in {1,10,20,30,40,50}
do
     # p=50
      #        "20","40","80","256","15_100","101_500","501_1000","1001_2000",
      for x in {"2001_10000","10001_20000","20001_60000"}
      do 
         #x=20001_60000
         INIT_PATH="../randTrees/randTrees_"$x
         RES_PATH="test_"$x
         NPR=$p
         ergodic $INIT_PATH $RES_PATH $NPR
#         # echo $RES_PATH
#         # echo $INIT_PATH
     done
done
#INIT_PATH="../randTrees/randTrees_15_100"
#RES_PATH="testtest_15_100"
#NPR=1
#ergodic $INIT_PATH $RES_PATH $NPR


