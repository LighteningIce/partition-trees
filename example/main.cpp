// Created by changjiang GOU on 10/05/2018.
// Copyright © 2018 Changjiang GOU. All rights reserved.

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "lib-io-tree.h"
#include "heuristics.h"

/**
 * optarg：表示当前选项对应的参数值
 * optind：表示的是下一个将被处理到的参数在argv中的下标值。
 * optopt：表示没有被未标识的选项
 * */

int main(int argc, char** argv) {
    clock_t beginTime = clock();
   int c;
   string stage1, stage2="NA", stage3;
   bool broken_already=false;

   int memory_constraint;
   double CCR, NPR;

   string dir;
   string brokenEdgesFile_dir;
   bool quiet=false;

    // printf("s=%d\n",'s');
    // printf("c=%d\n",'c');
    // printf("n=%d\n",'n');
    // printf("d=%d\n",'d');

   while (1) {
       //int this_option_optind = optind ? optind : 1;
       int option_index = 0;
       static struct option long_options[]={
           {"split",0,0,0},
           {"im",0,0,0},
           {"avoidchain",0,0,0},
           {"firstfit",0,0,0},
           {"immediately",0,0,0},
           {"largestfirst",0,0,0},
           {"ccr",1,0,0},
           {"npr",1,0,0},
           {"m",1,0,0},
           {"b",1,0,0},
           {"d",1,0,0},
           {0,0,0,0}//!!!
       };
    //getopt_long
    /**
     * int getopt_long(int argc, char * const argv[], const char *optstring, const struct option *longopts, int *longindex);  
     * argc和argv和main函数的两个参数一致
     * opstring ：表示【短选项字】符串
     *      1. 只有一个字符，不带冒号--只表示选项    -a
     *      2. 一个字符，后接一个冒号--表示选项后面带一个参数  -m 100
     *      3. 一个字符，两个冒号 -- 表示选项后面带一个可选参数，参数可有可无， 
     *      如果带参数，则选项与参数直接不能有空格，形式应该如-b200 
     * longopts:表示【长选项】结构体
     * longindex：longindex非空，它指向的变量将记录当前找到参数符合longopts里的第几个元素的描述，即是longopts的下标值。
     * */
       c = getopt_long(argc, argv, ":siafelqc:n:m:b:d:", long_options, &option_index);
        // printf("c=%d\n",c);


       if (c==-1) {
           break;
       }

       switch (c) {
           case 's':
               stage1 = "SplitSubtrees";
               break;

           case 'i':
               stage1 = "ImprovedSplit";
               break;

           case 'a':
               stage1 = "AvoidChain";
               break;

           case 'f':
               stage2 = "FirstFit";
               break;

           case 'e':
               stage2 = "Immediately";
               break;

           case 'l':
               stage2 = "LargestFirst";
               break;

           case 'q':
               quiet = true;
               break;

           case 'c':
           {//double CCRs[] = {1, 0.1, 0.001};//communication to computation
           //optarg：表示当前选项对应的参数值
               CCR=atof(optarg);
               break;
           }

           case 'n':
           {//double NPR[] = {1000, 100, 10};//amonut of nodes to processors
               NPR=atof(optarg);
               break;
           }

           case 'm':
           {
               memory_constraint = atoi(optarg);
           }
               break;

           case 'b':
           {
               broken_already=true;
               brokenEdgesFile_dir=optarg;
           }
               break;

           case 'd':
           {
               dir = optarg;
               std::cout<<dir<<endl;
           }
               break;

           case '?':
               printf("Usage call-heuristics -s -i -a -f -e -l -c [ccr] -n [npr] -m [0|1|2] -b [brokenfile] Directory Treelist\n");

           default:
               break;
       }
   }
//optind：表示的是下一个将被处理到的参数在argv中的下标值。
    // printf("optind=%d\n",optind);
    // printf("argc=%d\n",argc);
    //???有问题
//    if (optind >= argc) {
//        printf("Please provide a tree.\n");
//        exit(1);
//    }
    if (optind < argc) {
       printf("Please provide a tree.\n");
       exit(1);
   }

   int tree_size=0;
   int *prnts;
   //spacewghts=nwghts；timewghts=mswghts
   double *ewghts, *spacewghts, *timewghts;
   list<Cnode*> parallelSubtrees;
   int *chstart,*chend,*children,root=1;
   string treename,buffer;
   char cur_char;
   unsigned int num_processors;
   double makespan,maxoutd,minMem,memorySize;
   uint64_t count;
   clock_t time;

   unsigned int number_subtrees;

   cout.precision(20);

   std::cout<<"TreeName "<<"NPR "<<"CCR "<<"MemoryConstraint "<<"AmountSubtrees "<<"AmountProcessors "<<"Makespan "<<"Stage1 "<<"Stage2 "<<"Stage3 "<<"TimeConsuming"<<std::endl;
//    ifstream OpenFile(dir+argv[optind]);
    ifstream OpenFile(dir);
   ifstream BrokenEdgesFile(brokenEdgesFile_dir);
   int memory_constraint_options[3]={1,2,3};
   std::vector<int> brokenEdges;
//    do{
       treename="Tree1";
    //    OpenFile>>treename;

       parse_tree((dir).c_str(), &tree_size, &prnts, &spacewghts, &ewghts,&timewghts);
    // parse_tree_Seq((dir).c_str(), &tree_size, &prnts, &spacewghts, &ewghts,&timewghts);
       Ctree *treeobj = new Ctree(tree_size,prnts,spacewghts,ewghts,timewghts);
    //    treeobj->Print(cout);

       num_processors=ceil(tree_size/NPR);
       if (num_processors<3) {
           num_processors=3;
       }
       SetBandwidth(CCR, tree_size, ewghts, timewghts);

       if (broken_already) {
           time = 0;
           brokenEdges.clear();
           do{
               BrokenEdgesFile>>buffer;
               brokenEdges.push_back(stoi(buffer));
               treeobj->GetNode(stoi(buffer))->BreakEdge();//Edges already broken
               cur_char = BrokenEdgesFile.get();
           }while(cur_char != '\n' && BrokenEdgesFile.good());
           makespan=treeobj->GetRoot()->GetMSCost(true, true);
       }else{
           time = clock();
           if (stage1=="SplitSubtrees") {
               unsigned long sequentialLen;
               //？twolevel?=true
               makespan = SplitSubtrees(treeobj->GetRoot(), num_processors, false, parallelSubtrees,sequentialLen);// for counting how many subtrees produced, twolevel is set as false
           
           }else if(stage1=="ImprovedSplit"){
               po_construct(tree_size, prnts, &chstart, &chend, &children, &root);
               makespan = ImprovedSplit(treeobj, num_processors, chstart, children);
               //makespan = ImprovedSplit(treeobj, num_processors);
           }else if(stage1=="AvoidChain"){
               makespan = ASAP(treeobj, num_processors);
               number_subtrees = HowmanySubtrees(treeobj,true);
               time = clock()-time;
               std::cout<<treename<<" "<<NPR<<" "<<CCR<<" NA "<<number_subtrees<<" "<<num_processors<<" "<<makespan<<" "<<"ASAP NA NA "<<time<<std::endl;

               time = clock();
               number_subtrees = AvoidChain(treeobj);
               makespan=treeobj->GetRoot()->GetMSCost(true, true);
           }
           time = clock()-time;
       }
       number_subtrees = HowmanySubtrees(treeobj, quiet);
       std::cout<<treename<<" "<<NPR<<" "<<CCR<<" NA "<<number_subtrees<<" "<<num_processors<<" "<<makespan<<" "<<stage1<<" NA NA "<<time<<std::endl;

       if (stage2 != "NA") {
           maxoutd = MaxOutDegree(treeobj, true);
           schedule_t * schedule_f = new schedule_t();
           count = 0;
           MinMem(treeobj, maxoutd, minMem, *schedule_f, true, count);
           delete schedule_f;
           po_construct(tree_size, prnts, &chstart, &chend, &children, &root);

           vector<Ctree*> trees;
           Ctree *treeobj2 = new Ctree(tree_size,prnts,spacewghts,ewghts,timewghts);
           Ctree *treeobj3 = new Ctree(tree_size,prnts,spacewghts,ewghts,timewghts);
           trees.push_back(treeobj3);
           trees.push_back(treeobj2);
           trees.push_back(treeobj);

           for (int i=0; i<3; ++i) {
               treeobj=trees.back();
               trees.pop_back();
               for (std::vector<int>::iterator it = brokenEdges.begin(); it!=brokenEdges.end(); ++it) {
                   treeobj->GetNode(*(it))->BreakEdge();//Edges already broken
               }

               memory_constraint = memory_constraint_options[i];
               if (memory_constraint==1) {
                   memorySize = maxoutd;
               }else if (memory_constraint==2){
                   memorySize = (maxoutd + minMem)/2;
               }else{
                   memorySize = minMem;
               }

               time=clock();
               if (stage2=="LargestFirst") {
                   MemoryCheck(treeobj, chstart, children, memorySize, LARGEST_FIT);
               }else if(stage2=="FirstFit"){
                   MemoryCheck(treeobj, chstart, children, memorySize, FIRST_FIT);
               }else if(stage2=="Immediately"){
                   MemoryCheck(treeobj, chstart, children, memorySize, IMMEDIATELY);
               }

               time=clock()-time;
               number_subtrees = HowmanySubtrees(treeobj, true);
               makespan = treeobj->GetRoot()->GetMSCost(true,true);
               std::cout<<treename<<" "<<NPR<<" "<<CCR<<" "<<memory_constraint<<" "<<number_subtrees<<" "<<num_processors<<" "<<makespan<<" "<<stage1<<" "<<stage2<<" NA "<<time<<std::endl;

               time=clock();
               if (number_subtrees>num_processors) {
                   stage3 = "Merge";
                   makespan = Merge(treeobj, number_subtrees, num_processors, memorySize, chstart, children, true);
               } else if (number_subtrees==num_processors) {
                   stage3 = "Nothing";
               } else {
                   stage3 = "SplitAgain";
                   makespan = SplitAgain(treeobj, num_processors, number_subtrees);
               }
               time=clock()-time;

               number_subtrees = HowmanySubtrees(treeobj, true);
               std::cout<<treename<<" "<<NPR<<" "<<CCR<<" "<<memory_constraint<<" "<<number_subtrees<<" "<<num_processors<<" "<<makespan<<" "<<stage1<<" "<<stage2<<" "<<stage3<<" "<<time<<std::endl;
               
               delete treeobj;
           }
        cout << "total time: " << clock() - beginTime << endl;
           delete[] chstart;
           delete[] chend;
           delete[] children;
       }

       //delete treeobj;
       delete[] prnts;
       delete[] ewghts;
       delete[] spacewghts;
       delete[] timewghts;
//    }while (OpenFile.good());
   OpenFile.close();
   BrokenEdgesFile.close();

   exit(EXIT_SUCCESS);
}