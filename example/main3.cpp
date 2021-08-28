#include <iostream>
 #include <fstream>
 #include <string>
 #include <math.h>
 #include <stdlib.h>
 #include "lib-io-tree.h"
 #include "heuristics.h"

//./call-heuristics -s -c 0.001 -n 100 -d   .././realTrees/example10.tree.nf  >>  ../result/result_split_1
 int main(int argc, const char * argv[]) {
    //  for(int i=0;i<argc;i++){
    //      cout<<i<<" : "<<argv[i]<<endl;
    //  }
    int tree_size=0;
    int *prnts;
    double *ewghts, *spacewghts, *timewghts;
    int * chstart,*chend,*children;
    int root;
    string dir=argv[7]; //".././realTrees/example10.tree.nf";
    // cout<<"dir:"<<dir<<endl;
    string treename="Tree1";
    double makespan;
    double maxoutd;
    //double CCRs[] = {1, 0.1, 0.001};//commumakenication to computation
    double CCR=atof(argv[3]);//0.001;//
    //double NPR[] = {100000, 10000, 1000};//ratio of nodes' amount to processors
    double NPR=atof(argv[5]);//100;//;
    clock_t time;

    unsigned int number_subtrees;
    unsigned int num_processors;
    int memory_constraint_options[3]={1,2,3};
    int memory_constraint;
    double memorySize, minMem;
    uint64_t count;
    string stage2heuristic;

    cout.precision(20 );

    std::cout<<"TreeName NPR CCR MemoryConstraint AmountSubtrees AmountProcessors Makespan Heuristic TimeConsuming"<<std::endl;

    ifstream OpenFile(dir);//+argv[2]);
    // do{
	// std::cout<<"this a test1！"<<std::endl;
    // //     OpenFile>>treename;
	// std::cout<<"this a test2！"+treename+"uu"<<std::endl;
    
        parse_tree((dir).c_str(), &tree_size, &prnts, &spacewghts, &ewghts,&timewghts);
	 // std::cout<<"this a test3！"<<std::endl;
        num_processors=ceil(tree_size/NPR);
        if(num_processors<3){
            num_processors=3;
        }
        SetBandwidth(CCR, tree_size, ewghts, timewghts);

        Ctree *treeobj = new Ctree(tree_size,prnts,spacewghts,ewghts,timewghts);
        maxoutd = MaxOutDegree(treeobj, true);
        po_construct(tree_size, prnts, &chstart,&chend,&children, &root);

        time = clock();
        makespan = treeobj->GetRoot()->GetMSCost();
        number_subtrees = 1;
        time = clock()-time;
        cout<<treename<<" "<<NPR<<" "<<CCR<<" NA "<<number_subtrees<<" "<<num_processors<<" "<<makespan<<" Sequence "<<time<<endl;

        schedule_t * schedule_f = new schedule_t();
        count = 0;
        MinMem(treeobj, maxoutd, minMem, *schedule_f, true, count);
        delete schedule_f;
        delete treeobj;

	  //step2:fiting into memory
        for (int stage2Method=0; stage2Method<3; ++stage2Method) {
            for (int i=0; i<3; ++i){
                Ctree *treeobj = new Ctree(tree_size,prnts,spacewghts,ewghts,timewghts);

                memory_constraint = memory_constraint_options[i];
                if (memory_constraint==1) {
                    memorySize = maxoutd;
                }else if (memory_constraint==2){
                    memorySize = (maxoutd + minMem)/2;
                }else{
                    memorySize = minMem;
                }

                time=clock();
                switch (stage2Method) {
                    case 0:
                        stage2heuristic = "FIRST_FIT";
                        MemoryCheck(treeobj, chstart, children, memorySize, FIRST_FIT);
                        break;
                    case 1:
                        stage2heuristic = "LARGEST_FIT";
                        MemoryCheck(treeobj, chstart, children, memorySize, LARGEST_FIT);
                        break;
                    case 2:
                        stage2heuristic = "IMMEDIATELY";
                        MemoryCheck(treeobj, chstart, children, memorySize, IMMEDIATELY);
                        break;

                    default:
                        stage2heuristic = "FIRST_FIT";
                        MemoryCheck(treeobj, chstart, children, memorySize, IMMEDIATELY);
                        break;
                }
                time=clock()-time;

                makespan = treeobj->GetRoot()->GetMSCost(true,true);
                number_subtrees = HowmanySubtrees(treeobj, true);

               std::cout<<treename<<" "<<NPR<<" "<<CCR<<" "<<memory_constraint<<" "<<number_subtrees<<" "<<num_processors<<" "<<makespan<<" "<<stage2heuristic<<" "<<time<<endl;

		//step 3:
                if (number_subtrees > num_processors) {
                    time = clock();
                    makespan = MergeV2(treeobj, number_subtrees, num_processors, memorySize, chstart, children, true);
                    time = clock() - time;
                    number_subtrees = HowmanySubtrees(treeobj, true);
                    std::cout<<treename<<" "<<NPR<<" "<<CCR<<" "<<memory_constraint<<" "<<number_subtrees<<" "<<num_processors<<" "<<makespan<<" "<<stage2heuristic<<"+Merge "<<time<<endl;
                } else if (number_subtrees == num_processors){
                    std::cout<<treename<<" "<<NPR<<" "<<CCR<<" "<<memory_constraint<<" "<<number_subtrees<<" "<<num_processors<<" "<<makespan<<" "<<stage2heuristic<<"+Nothing "<<0<<endl;
                } else {
                    time = clock();
                    makespan = SplitAgain(treeobj, num_processors, number_subtrees);
                    time = clock() - time;
                    number_subtrees = HowmanySubtrees(treeobj, true);
                    std::cout<<treename<<" "<<NPR<<" "<<CCR<<" "<<memory_constraint<<" "<<number_subtrees<<" "<<num_processors<<" "<<makespan<<" "<<stage2heuristic<<"+SplitAgain "<<time<<endl;
                }

                delete treeobj;
            }
        }

        delete[] prnts;
        delete[] ewghts;
        delete[] spacewghts;
        delete[] timewghts;
        delete [] chstart;
        delete [] chend;
        delete [] children;
    // }while (OpenFile.good());
    OpenFile.close();

    return 0;
 }
