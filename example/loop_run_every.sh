#!/bin/bash

function ergodic(){
dir_path="../result/$2";
# test_lessprocessor
# org_path="../result/test_real_tree_every/test_realtree_lessprocessor_p$3"
# mkdir $org_path
# dir_path="../result/test_real_tree_every/test_realtree_lessprocessor_p$3/$2";
mkdir  $dir_path;
for i in {0.01,0.10,1.00,2.00}
do
    echo $i;
    for file in ` ls $1`
    do
        echo $1"/"$file
    #     if [ -d $1"/"$file] #如果 file存在且是一个目录则为真
    #     then                       ergodic $1"/"$file
    # 　　else                       local path=$1"/"$file #得到文件的完整的目录
    # 　　　　　　　　　　　　　　　　　local name=$file       #得到文件的名字
    # 　　 fi
        # ./call-heuristics -s -q  -c 0.001 -n 0.7 -d   $1"/"$file  >>  test.txt
        
        ./call-heuristics -s -l  -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt;
        ./select-heuristics -s -l -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt;
        ./call-heuristics-main -s -l  -q  -c $i -n $3 -d   $1"/"$file >>$dir_path/$2_c$i.txt

        #  ./call-heuristics -s -l  -q  -c $i -n 50 -d   $1"/"$file ;
        # ./select-heuristics -s -l -q  -c $i -n 50 -d   $1"/"$file ;
        # ./call-heuristics-main -s -l  -q  -c $i -n 50 -d   $1"/"$file
    done;
done
 }


# for  p in {10,20,30,40,50}
# do
#     for x in {"15_100","101_500","501_1000","1001_2000","2001_10000","10001_20000","20001_60000","20","40","80","256"}
#     do 
#         INIT_PATH="../real_tree_every"
#         RES_PATH="test_real_tree_every"
#         NPR=$p
#         ergodic $INIT_PATH $RES_PATH $NPR
#         # echo $RES_PATH
#         # echo $INIT_PATH
#     done
# done
INIT_PATH="../real_tree_every"
RES_PATH="test_real_tree_every"
NPR=1
ergodic $INIT_PATH $RES_PATH $NPR


