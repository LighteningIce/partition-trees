#!/bin/bash

function ergodic(){
# dir_path="../result_thin_outMdeg/result_p$3/$2";
# test_lessprocessor
 org_path="../result_thin_outMdeg/result_p$3"
 mkdir $org_path
 dir_path="../result_thin_outMdeg/result_p$3/$2";
 mkdir  $dir_path;
#  0.10,1.00,5.00,2.00,3.00,8.00,
for i in {0.01,10.00}
do
    #i=10.00
    echo $i;
    
    for file in ` ls $1`
    do
        echo $1"/"$file

        #./call-heuristics -s  -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt;
        ./call-heuristics -s  -l  -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt;
       ./call-heuristics -i -l  -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt;
        ./call-heuristics -a  -l  -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt;
        #./select-heuristics -l  -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt;
        ./call-heuristics-main  -l  -q  -c $i -n $3 -m 0 -d   $1"/"$file >>$dir_path/$2_c$i.txt        

#        ./call-heuristics -s -l  -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt;
#       ./call-heuristics -s -f  -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt;
#        ./call-heuristics -a -l  -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt;
#        ./call-heuristics -a -l  -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt;
#        ./call-heuristics-main -s -l  -q  -c $i -n $3 -m 0 -d   $1"/"$file >>$dir_path/$2_c$i.txt        
#        ./call-heuristics-main -s -f  -q  -c $i -n $3 -m 0 -d   $1"/"$file >>$dir_path/$2_c$i.txt
#        ./select-heuristics -s -l -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt;
        # ./memory-heuristics-test -s -l -q -c $i -n $3 -m 2 -d   $1"/"$file >>$dir_path/$2_c$i.txt
        #  ./call-heuristics -s -l  -q  -c $i -n 50 -d   $1"/"$file ;
        # ./select-heuristics -s -l -q  -c $i -n 50 -d   $1"/"$file ;
        # ./call-heuristics-main -s -l  -q  -c $i -n 50 -d   $1"/"$file
    done;
done
 }

#,  ,20,30,40    1,10,50,100,
for  p in {1000,10000}
do
     #p=100
     #for x in {"15_100","101_500","501_1000","1001_2000","2001_10000","10001_20000","20001_60000","20","40","80","256"}
    #
    #        "20","40","60","80","120","160","240","360","512",  ,"4096","8192"
     for  x in {"1024","2048"}
     do 
        #x=2048
        INIT_PATH="../randTrees_thin/randTrees_"$x
        RES_PATH="test_thin_"$x
        NPR=$p
        ergodic $INIT_PATH $RES_PATH $NPR
        # echo $RES_PATH
        # echo $INIT_PATH
   done
done
# INIT_PATH="../randTrees_thin/randTrees_120"
# RES_PATH="test_thin_120"
# NPR=1
# ergodic $INIT_PATH $RES_PATH $NPR
