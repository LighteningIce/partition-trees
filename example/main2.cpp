#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <stdlib.h>
#include "lib-io-tree.h"
#include "heuristics.h"

void PrintTree(Ctree* tree){
    const vector<Cnode*>* children = tree->GetNodes();
    //NW->me_weight spaceweight  
    //MSW->ms_weight   timeweigt
    cout<<"nodeId   parentId   ms_weight   me_weight   edge_weight"<<endl;
    for (vector<Cnode*>::const_iterator it=children->begin(); it!=children->end(); ++it) {
        cout<<(*it)->GetId()<<" "<<(*it)->GetParentId()<<" "<<(*it)->GetMSW()<<" "<<(*it)->GetNW()<<" "<<(*it)->GetEW()<<endl;
    }
    cout<<endl;
}


int main(int argc, const char * argv[]) {
    int tree_size=0;
    int *prnts;
    double *ewghts, *spacewghts, *timewghts;
    int *chstart,*chend,*children,root=1;
    string treename, buffer;;

    cout.precision(10);

    string tree = ".././realTrees/example12.tree.nf";

    parse_tree(tree.c_str(), &tree_size, &prnts, &spacewghts, &ewghts,&timewghts);
    Ctree *treeobj = new Ctree(tree_size,prnts,spacewghts,ewghts,timewghts);
    
    PrintTree(treeobj);

    double NPR = 100;
    double CCR = 1;
    SetBandwidth(CCR, tree_size, ewghts, timewghts);

    unsigned int numberProcessors=ceil(tree_size/NPR);
    if (numberProcessors<3) {
        numberProcessors=3;
    }

//    po_construct(tree_size, prnts, &chstart, &chend, &children, &root);
//    double makespan = ImprovedSplit(treeobj, numberProcessors, chstart, children);
//    cout<<"makespan of ImprovedSplit: "<<makespan<<endl;

    double makespan = treeobj->GetRoot()->GetMSCost();
    cout<<"number of processor:"<<numberProcessors<<endl;
    cout<<"number of tree: "<<tree_size<<endl;
    cout<<"makespan before "<<makespan<<endl;
    
    double maxoutd = MaxOutDegree(treeobj, true);
    po_construct(tree_size, prnts, &chstart, &chend, &children, &root);
    MemoryCheck(treeobj, chstart, children, maxoutd, LARGEST_FIT);
    makespan = treeobj->GetRoot()->GetMSCost(true,true);
    cout<<"makespan now "<<makespan<<endl;
    
    delete[] chstart;
    delete[] chend;
    delete[] children;
    delete treeobj;
    delete[] prnts;
    delete[] ewghts;
    delete[] spacewghts;
    delete[] timewghts;
    return 0;
}
