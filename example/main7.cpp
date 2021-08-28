#include <iostream>
#include <math.h>
#include <fstream>
#include <list>
#include <cmath>
#include <map>
#include <getopt.h>
#include "lib-io-tree-utils.h"
#include "heuristics.h"
#include <sys/stat.h>
#include <string.h>

using namespace std;

int main(int argc, char **argv)
{
    double CCR = 0.01, NPR = 100;
    bool quiet = false;
    int memory_constraint = 0;
    int memory_constraint_options[3] = {0, 1};
    string dir;
    int c;
    // string stage =
    while (1)
    {
        int option_index = 0;
        static struct option long_options[] = {
            {"ccr", 1, 0, 0},
            {"npr", 1, 0, 0},
            {"b", 1, 0, 0},
            {"d", 1, 0, 0},
            {"m", 1, 0, 0},
            {0, 0, 0, 0}};
        c = getopt_long(argc, argv, "qc:n:m:d:", long_options, &option_index);

        if (c == -1)
        {
            break;
        }

        switch (c)
        {
        case 'c':
            CCR = atof(optarg);
            break;
        case 'n':
            NPR = atof(optarg);
            break;
        case 'm':
            memory_constraint = atoi(optarg);
            break;
        case 'q':
            quiet = true;
            break;
        case 'd':
            dir = optarg;
            if (!quiet)
            {
                std::cout << dir << endl;
            }
            break;
        case '?':
            printf("Usage call-heuristics -s -i -a -f -e -l -c [ccr] -n [npr] -m [0|1|2] -b [brokenfile] Directory Treelist\n");

        default:
            break;
        }
    }
    if (optind < argc)
    {
        cout << "Please provide a tree.." << endl;
        exit(0);
    }
    cout.precision(15);
    ifstream OpenFile(dir);
    int tree_size;
    int *prnts;
    double *ewghts, *spacewghts, *timewghts;
    parse_tree_Seq((dir).c_str(), &tree_size, &prnts, &spacewghts, &ewghts, &timewghts);
    Ctree *treeobj = new Ctree(tree_size, prnts, spacewghts, ewghts, timewghts);
    // treeobj->Print(cout);
    unsigned int num_processors;
    unsigned int number_subtrees;

    num_processors = ceil(tree_size / NPR);
    if (num_processors < 3)
    {
        num_processors = 3;
    }
    SetBandwidth(CCR, tree_size, ewghts, timewghts);

    /**
     * 输入处理器的信息
     * */
    // memory_constraint = memory_constraint_options[i];
    double maxoutd = MaxOutDegree(tree_size, prnts, spacewghts, ewghts);
    schedule_t *schedule_f = new schedule_t();
    int count = 0;
    double minMem;
    MinMem(treeobj, maxoutd, minMem, *schedule_f, true, count);
    // cout << "ccr:" << CCR << endl;
    // cout << "num_processors:" << num_processors << endl;
    Procs *procs=new Procs(num_processors,minMem,memory_constraint);
    
    // procs->PrintProcs();
    
    // cout << "==========================partitioning=================================" << endl;
    io_method_t method = LARGEST_FIT;
    PartitionTree(treeobj, num_processors, quiet, method, memory_constraint);
    cout << "the subtree size:" << HowmanySubtrees(treeobj, false) << endl;
    cout  << treeobj->GetRoot()->GetMSCost(true, true)<< "   ";
    // cout << "====================================================================" << endl;
    // cout << endl;
    // cout << "==========================allocation=================================" << endl;
    double lastms = AllocPower(treeobj,prnts, procs, num_processors, quiet);
    // cout << "the subtree size:" << HowmanySubtrees(treeobj, true) << endl;
    cout << lastms << endl;
    // cout << "====================================================================" << endl;
    // cout << endl;
}
