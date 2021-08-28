#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "heuristics.h"
#include "lib-io-tree.h"
#include <sys/stat.h>

/**
 * optarg：表示当前选项对应的参数值
 * optind：表示的是下一个将被处理到的参数在argv中的下标值。
 * optopt：表示没有被未标识的选项
 * */
int j=0;
void get_small_tree(Ctree*treeobj,int *chstart,int *children,string filename);
int main(int argc, char **argv)
{
    clock_t beginTime = clock();
    int c;


    string dir;
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
        clock_t time;

        unsigned int number_subtrees;

        cout.precision(20);

        ifstream OpenFile(dir);

        int memory_constraint_options[3] = {1, 2, 3};
        std::vector<int> brokenEdges;
        // cout<<"test_1"<<endl;
        // dir ="../randTrees/randTrees_20001_60000/gen_data_20001_60000_23.txt";
        // cout<<dir<<endl;
        //  cout<<"test_2"<<endl;
        string fileName=dir.substr(dir.find_last_of('/')+1);
        cout<<fileName<<endl;
         string fileName_2=fileName.substr(0,fileName.length()-8);
         cout<<fileName_2<<endl;


        // dir="../real_tree_every/dump.0.1.amd.CEMW.t2em-1897.tree.nf";
        
        parse_tree((dir).c_str(), &tree_size, &prnts, &spacewghts, &ewghts, &timewghts);
        // parse_tree_Seq((dir).c_str(), &tree_size, &prnts, &spacewghts, &ewghts,&timewghts);
        Ctree *treeobj = new Ctree(tree_size, prnts, spacewghts, ewghts, timewghts);
        // tree_size=treeobj->GetNodes()->size();
        po_construct(tree_size, prnts, &chstart, &chend, &children, &root);
        vector<Cnode*>*Children=treeobj->GetRoot()->GetChildren();
        cout<<"tree size:"<<tree_size<<endl;
        cout<<"children size:"<<Children->size()<<endl;
        get_small_tree(treeobj,chstart,children,fileName_2);
        cout<<endl;
        cout<<"====================================================="<<endl;
}

void get_small_tree(Ctree*treeobj,int *chstart,int *children,string filename)
{
    int tree_size = treeobj->GetNodes()->size();
    double *ewght, *timewght, *spacewght;
     int *prnt;
     vector<Cnode*>*Children=treeobj->GetRoot()->GetChildren();
    
     for(vector<Cnode*>::iterator i=Children->begin();i!=Children->end();i++)
    {
       
        Ctree* subtree = BuildSubtree(treeobj, (*i), tree_size, &prnt, &ewght, &timewght, &spacewght, chstart, children);
        int *chstart_1, *chend_1, *children_1, root = 1;
        po_construct(tree_size, prnt, &chstart_1, &chend_1, &children_1, &root);
        int subTree_size=subtree->GetNodes()->size();
        // cout<<subTree_size<<endl;
        if(subTree_size>1000&&subTree_size<15000){
            cout<<subTree_size<<" ";
            string file_name_2=filename+"-"+to_string(j)+".tree.nf";
            j++;
            cout<<"file_name_2:"<<file_name_2<<endl;
            subtree->OutputToFile("/home/jgwu2_jsjxy/hesn/MemComJournal/realTrees_1",file_name_2);
        }else if(subTree_size<1000)
        {
            continue;
        }else
        {
            get_small_tree(subtree,chstart_1,children_1,filename);
        }
        
    }
}