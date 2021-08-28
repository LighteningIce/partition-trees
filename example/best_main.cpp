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
class Solution
{
public:
    vector<vector<int>> ans;
    vector<int> tmp;
    void find(int dep, vector<int> &nums)
    {
        if (dep <= 0)
        {
            ans.push_back(tmp);
            return;
        }
        tmp.push_back(nums[dep - 1]);
        find(dep - 1, nums);
        tmp.pop_back();
        find(dep - 1, nums);
    }

    vector<vector<int>> subsets(vector<int> &nums)
    {
        find(nums.size(), nums);
        return ans;
    }
};

int main(int argc, char **argv)
{
    clock_t beginTime = clock();
    int c;
    string stage1, stage2 = "NA", stage3;
    bool broken_already = false;

    int memory_constraint;
    double CCR, NPR;

    string dir;
    string brokenEdgesFile_dir;
    bool quiet = false;

    // printf("s=%d\n",'s');
    // printf("c=%d\n",'c');
    // printf("n=%d\n",'n');
    // printf("d=%d\n",'d');

    while (1)
    {
        //int this_option_optind = optind ? optind : 1;
        int option_index = 0;
        static struct option long_options[] = {
            {"split", 0, 0, 0},
            {"im", 0, 0, 0},
            {"avoidchain", 0, 0, 0},
            {"firstfit", 0, 0, 0},
            {"immediately", 0, 0, 0},
            {"largestfirst", 0, 0, 0},
            {"ccr", 1, 0, 0},
            {"npr", 1, 0, 0},
            {"m", 1, 0, 0},
            {"b", 1, 0, 0},
            {"d", 1, 0, 0},
            {0, 0, 0, 0} //!!!
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

        if (c == -1)
        {
            break;
        }

        switch (c)
        {
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
        { //double CCRs[] = {1, 0.1, 0.001};//communication to computation
            //optarg：表示当前选项对应的参数值
            CCR = atof(optarg);
            break;
        }

        case 'n':
        { //double NPR[] = {1000, 100, 10};//amonut of nodes to processors
            NPR = atof(optarg);
            break;
        }

        case 'm':
        {
            memory_constraint = atoi(optarg);
        }
        break;

        case 'b':
        {
            broken_already = true;
            brokenEdgesFile_dir = optarg;
        }
        break;

        case 'd':
        {
            dir = optarg;
            // std::cout << dir << endl;
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
    if (optind < argc)
    {
        printf("Please provide a tree.\n");
        exit(1);
    }

    int tree_size = 0;
    int *prnts;
    //spacewghts=nwghts；timewghts=mswghts
    double *ewghts, *spacewghts, *timewghts;
    list<Cnode *> parallelSubtrees;
    int *chstart, *chend, *children, root = 1;
    string treename, buffer;
    char cur_char;
    unsigned int num_processors;
    double makespan, maxoutd, minMem, memorySize;
    uint64_t count;
    clock_t time, time_2, time_3;

    unsigned int number_subtrees;

    cout.precision(20);

    // std::cout << "TreeName "
    //           << "NPR "
    //           << "CCR "
    //           << "MemoryConstraint "
    //           << "AmountSubtrees "
    //           << "AmountProcessors "
    //           << "Makespan "
    //           << "Stage1 "
    //           << "Stage2 "
    //           << "Stage3 "
    //           << "TimeConsuming" << std::endl;
    //    ifstream OpenFile(dir+argv[optind]);
    ifstream OpenFile(dir);
    ifstream BrokenEdgesFile(brokenEdgesFile_dir);
    int memory_constraint_options[3] = {1, 2, 3};
    std::vector<int> brokenEdges;
    //    do{
    treename = "Tree1";
    //    OpenFile>>treename;

    //    parse_tree((dir).c_str(), &tree_size, &prnts, &spacewghts, &ewghts,&timewghts);
    parse_tree_Seq((dir).c_str(), &tree_size, &prnts, &spacewghts, &ewghts, &timewghts);
    Ctree *treeobj = new Ctree(tree_size, prnts, spacewghts, ewghts, timewghts);
    treeobj->Print(cout);
    const vector<Cnode *> *NChildren = treeobj->GetNodes();

    bool *broken = new bool[tree_size];
    broken[0] = true;
    double ms_min = __DBL_MAX__;

    // for (int i = 1; i < tree_size; i++)
    // {
    //     // cout << "i=" << i << endl;
    //     treeobj->GetNode(1)->BreakEdge();
    //     for (int j = i; j < tree_size; j++)
    //     {
    //         for (int l = i; l <= j; l++)
    //         {
    //             cout << l;
    //             broken[l] = true;
    //             treeobj->GetNode(l + 1)->BreakEdge();
    //         }
    //         cout << endl;
    //         double MS=treeobj->GetRoot()->GetMSCost(true, true) ;
    //         cout << "MS=" << MS<< endl;
    //         if(MS<ms_min)
    //         {
    //             ms_min=MS;
    //         }
    //         for (int k = 0; k < tree_size; k++)
    //         {
    //             // cout << broken[k];
    //             if (broken[k] == true)
    //             {
    //                 treeobj->GetNode(k + 1)->RestoreEdge();
    //             }
    //             broken[k] = 0;
    //         }
    //         cout << endl;
    //     }
    //     //
    // }
    // cout<<"min MS:"<<ms_min<<endl;
    Solution S;
    vector<int> nums;
    nums.push_back(1);
    nums.push_back(2);
    nums.push_back(3);
    nums.push_back(4);
    nums.push_back(5);
    nums.push_back(6);
    nums.push_back(7);
    nums.push_back(8);
    nums.push_back(9);
    nums.push_back(10);
    nums.push_back(11);
    nums.push_back(12);
    nums.push_back(13);
    nums.push_back(14);

    vector<vector<int>> subset = S.subsets(nums);
    for (vector<vector<int>>::iterator iter = subset.begin(); iter != subset.end(); iter++)
    {
        treeobj->GetNode(1)->BreakEdge();
        for (vector<int>::iterator iiter = (*iter).begin(); iiter != (*iter).end(); iiter++)
        {
            cout << (*iiter);
            broken[(*iiter)-1] = true;
            treeobj->GetNode((*iiter))->BreakEdge();
        }
        cout << endl;
         for (int k = 0; k < tree_size; k++)
        {
            cout << broken[k];
        }
        cout << endl;
        double MS = treeobj->GetRoot()->GetMSCost(true, true);
        cout << "MS=" << MS << endl;
        if (MS < ms_min)
        {
            ms_min = MS;
        }

        //             if (broken[k] == true)
        //             {
        //                 treeobj->GetNode(k + 1)->RestoreEdge();
        //             }
        //             broken[k] = 0;

        for (vector<int>::iterator iiter = (*iter).begin(); iiter != (*iter).end(); iiter++)
        {
            if (broken[(*iiter)-1] == true)
            {
                treeobj->GetNode((*iiter) )->RestoreEdge();
            }
            broken[(*iiter)-1] = 0;
        }
        for (vector<int>::iterator iiter = (*iter).begin(); iiter != (*iter).end(); iiter++)
        {
            cout<<treeobj->GetNode((*iiter) )->IsBorken();
        }
        // for (int k = 0; k < tree_size; k++)
        // {
        //     cout << broken[k];
        // }
        cout << endl;
        cout << "=============================================" << endl;
    }
    cout << "min MS:" << ms_min << endl;
}
