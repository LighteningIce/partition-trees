#!/bin/bash

function ergodic(){
for file in ` ls $1`
do
    echo $1"/"$file
#     if [ -d $1"/"$file] #如果 file存在且是一个目录则为真
#     then                       ergodic $1"/"$file
# 　　else                       local path=$1"/"$file #得到文件的完整的目录
# 　　　　　　　　　　　　　　　　　local name=$file       #得到文件的名字
# 　　 fi
    # ./call-heuristics -s -q  -c 0.001 -n 0.7 -d   $1"/"$file  >>  test.txt
    # ./call-heuristics -s -f  -q  -c 0.001 -n 100 -d   $1"/"$file >>test2.txt;
    # ./call-heuristics -s -l  -q   -c  0.01  -n 50 -d   $1"/"$file >>test3.txt;
    # ./select-heuristics -s -l   -q -c 0.01  -n 50 -d   $1"/"$file >>test3.txt;
    
    # ./memory-heuristics-test -s -l -q -c 1-n 1 -m 0 -d $1"/"$file >>test3.txt
    # ./call-heuristics -s   -c 0.01 -n 1 -d   $1"/"$file >>test3.txt;
    ./call-heuristics -i  -f  -c 1 -n 1 -d   $1"/"$file >>test3.txt;
    #./call-heuristics -a  -l -c 1 -n 1 -d   $1"/"$file >>test3.txt;
    #./select-heuristics -s    -q -c 1 -n 1 -d   $1"/"$file >>test3.txt;
    #./call-heuristics-main -s -l     -c 1  -n 1 -m 0 -d   $1"/"$file >>test3.txt
    ./hetro-heuristics -c 1 -n 1 -q  -m 0 -d $1"/"$file   >>test3.txt;

done }

# INIT_PATH="../randTrees_15_100"
INIT_PATH="../re1tree"
ergodic $INIT_PATH



