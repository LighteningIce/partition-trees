#!/bin/bash
function rand(){
    min=$1
    max=$(($2-$min+1))
    num=$(date +%s%N)
    echo $(($num%$max+$min))
}
function ergodic(){
for ((i=1;i<=100;i++));
do
    # tree_num=$(rand  20001 60000);#结点个数
    tree_num=8192
    # echo $tree_num;
    if [ $i -eq 1 ]
    # then mkdir ../randTrees_20001_60000
    then mkdir ../randTrees_thin/randTrees_$tree_num;
    fi
    
    rate=0.8
    #结点最大度
    degree=$(rand 2 5);
    #树最大高度
    height=`echo "$tree_num * $rate"|bc`;
    max_node=`echo "$height * $degree"|bc`
    f_num1=${tree_num%.*}
     f_num2=${max_node%.*}
    echo "height:"$f_num2;
    while [  $f_num2 -lt $f_num1 ]
    do
        let degree=$(rand 2 4);
        height=`echo "$tree_num * $rate"|bc`;
        height=${height%.*}
        max_node=`echo "$height * $degree"|bc`
        f_num2=${max_node%.*}
        
        echo "f_num2:"$f_num2;
        # echo "f_num1:"$f_num1;
    done;  
    # echo "degree=$degree"
    # echo "height=$height"
    # echo `expr $height \* $degree`;
    # echo $i
    if [ $((i % 5)) -eq 0 ]
    then sleep 3;
    fi
    ./main_gen $tree_num  $degree   $height  >> ../randTrees_thin/randTrees_$tree_num/gen_data_$tree_num'_'$i.txt;
    # ./main_gen $tree_num  $degree   $height>> ../randTrees_20001_60000/gen_data_20001_60000_$i.txt;
    # echo -n "随机数:"   ;expr $(date +%s%N)%$[$max - $min  + 1] + $min
    sleep 1;
done 
}


ergodic
