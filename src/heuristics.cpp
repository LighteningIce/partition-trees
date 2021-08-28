//
//  heuristics.cpp
//  memCom2
//
//  Created by changjiang GOU on 11/05/2018.
//  Copyright © 2018 ROMA. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <list>
#include <algorithm>
#include "heuristics.h"
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>
#include <random>
#include <chrono>
//#include <omp.h>

extern double BANDWIDTH;

bool cmp_noincreasing(Cnode *a, Cnode *b) { return (a->GetMSCost(true, false) >= b->GetMSCost(true, false)); };
bool cmp_nodecreasing(Cnode *a, Cnode *b) { return (a->GetMSCost(true, false) < b->GetMSCost(true, false)); };
//bool cmp_noIn_minusCommu (Cnode* a, Cnode* b){return (a->GetMSminusComu()>=b->GetMSminusComu());};
bool cmp_noIn_noCommu(Cnode *a, Cnode *b) { return (a->GetMSCost(false, false) >= b->GetMSCost(false, false)); };
bool cmp_c_nodecrea(Cnode *a, Cnode *b) { return (a->GetMSW() < b->GetMSW()); };
bool cmp_nodecreasing_double(double a, double b) { return (a < b); }; //double 型非增数据
bool cmp_noincreasing_msw(Cnode *a, Cnode *b)
{
    return (a->GetMSW() <= b->GetMSW());
};
//void GetTwoLargestElementTypeone(vector<Cnode*>* container, Cnode* & Largest, Cnode* & secondLargest) {
//    if (container->front()->GetMSminusComu()>=container->at(1)->GetMSminusComu()) {
//        Largest=container->front();
//        secondLargest=container->at(1);
//    }else{
//        Largest=container->at(1);
//        secondLargest=container->front();
//    }
//
//    if (container->size()>2) {
//        vector<Cnode*>::iterator iter=container->begin();
//        iter=iter+2;
//        for (;iter!=container->end(); ++iter) {
//            if ((*iter)->GetMSminusComu()>Largest->GetMSminusComu()) {
//                secondLargest=Largest;
//                Largest=*iter;
//            } else if((*iter)->GetMSminusComu()>secondLargest->GetMSminusComu()){
//                secondLargest=*iter;
//            }
//        }
//    }
//}

template <class T, class U>
void GetTwoLargestElementTypetwo(T container, U &Largest, U &secondLargest)
{
    if (container->front()->GetMSCost(true, false) > container->at(1)->GetMSCost(true, false))
    {
        Largest = container->front();
        secondLargest = container->at(1);
    }
    else
    {
        Largest = container->at(1);
        secondLargest = container->front();
    }

    if (container->size() > 2)
    {
        vector<Cnode *>::iterator iter = container->begin();
        iter = iter + 2;
        for (; iter != container->end(); ++iter)
        {
            if ((*iter)->GetMSCost(true, false) > Largest->GetMSCost(true, false))
            {
                secondLargest = Largest;
                Largest = *iter;
            }
            else if ((*iter)->GetMSCost(true, false) > secondLargest->GetMSCost(true, false))
            {
                secondLargest = *iter;
            }
        }
    }
}

void GetTwoLargestElementTypethree(vector<Cnode *> *container, vector<Cnode *>::iterator &Largest, vector<Cnode *>::iterator &secondLargest)
{
    if (container->front()->GetMSCost(false, false) >= container->back()->GetMSCost(false, false))
    {
        Largest = container->begin();
        secondLargest = Largest;
        advance(secondLargest, 1);
    }
    else
    {
        secondLargest = container->begin();
        Largest = secondLargest;
        advance(Largest, 1);
    }

    if (container->size() > 2)
    {
        vector<Cnode *>::iterator iter = container->begin();
        advance(iter, 2);
        for (; iter != container->end(); ++iter)
        {
            if ((*iter)->GetMSCost(false, false) > (*Largest)->GetMSCost(false, false))
            {
                secondLargest = Largest;
                Largest = iter;
            }
            else if ((*iter)->GetMSCost(false, false) > (*secondLargest)->GetMSCost(false, false))
            {
                secondLargest = iter;
            }
        }
    }
}

void GetTwoSmallestElement(list<Cnode *> *container, list<Cnode *>::iterator &Smallest, list<Cnode *>::iterator &secondSmallest)
{
    if (container->front()->GetMSW() <= container->back()->GetMSW())
    {
        Smallest = container->begin();
        secondSmallest = Smallest;
        advance(secondSmallest, 1);
    }
    else
    {
        secondSmallest = container->begin();
        Smallest = secondSmallest;
        advance(Smallest, 1);
    }

    if (container->size() > 2)
    {
        list<Cnode *>::iterator iter = container->begin();
        advance(iter, 2);
        for (; iter != container->end(); ++iter)
        {
            if ((*iter)->GetMSW() < (*Smallest)->GetMSW())
            {
                secondSmallest = Smallest;
                Smallest = iter;
            }
            else if ((*iter)->GetMSW() < (*secondSmallest)->GetMSW())
            {
                secondSmallest = iter;
            }
        }
    }
}

double SplitSubtrees(Cnode *root, unsigned long num_processor, double twolevel, list<Cnode *> &parallelRoots, unsigned long &sequentialLength)
{
    // cout<<"this is SplitSubtrees "<<endl;
    parallelRoots.clear();
    parallelRoots.emplace_front(root); //在开头添加root
    //cout<<"   insert root"<<endl;
    vector<double> MS(1, root->GetMSCost(true, true)); // take communication cost into account
    double MS_sequential = root->GetEW() / BANDWIDTH, Weight_more, Weight_PQ;
    unsigned long amountSubtrees;
    vector<Cnode *> *children;

    Cnode *currentNode = root;
    double temp;
    unsigned int mergetime;
    //list<Cnode*>::iterator iter;
    //unsigned int target;
    //unsigned int round=0;
    while (!currentNode->IsLeaf())
    {
        MS_sequential = MS_sequential + currentNode->GetMSW();

        Weight_PQ = 0;
        parallelRoots.remove(currentNode);
        // cout<<"pop up "<<currentNode->GetId()<<endl;

        children = currentNode->GetChildren();
        for (vector<Cnode *>::iterator iter = children->begin(); iter != children->end(); iter++)
        {
            if ((*iter)->IsBorken())
            {
                temp = (*iter)->GetMSCost(true, false);
                if (temp > Weight_PQ)
                {
                    Weight_PQ = temp;
                }
            }
            else
            {
                parallelRoots.push_back(*iter);
                // cout<<"   insert "<<(*iter)->GetId()<<endl;
            }
        }

        if (parallelRoots.empty())
        {
            break;
        }
        else
        {
            //找parallelRoots中MS最大的元素
            currentNode = *max_element(parallelRoots.begin(), parallelRoots.end(), cmp_nodecreasing); //non-decreasing
        }

        temp = currentNode->GetMSCost(true, false);
        if (temp > Weight_PQ)
        {
            Weight_PQ = temp;
        }

        Weight_more = 0;
        amountSubtrees = parallelRoots.size() + 1;
        //计算超出处理器个数的那部分子树的MS之和Weight_more
        if (amountSubtrees > num_processor)
        {
            parallelRoots.sort(cmp_noIn_noCommu); //non-increasing sort, computation weight, no communication
            list<Cnode *>::reverse_iterator iter = parallelRoots.rbegin();
            mergetime = amountSubtrees - num_processor;
            //计算最小的mergetime个子树的权重之和
            for (unsigned int i = 0; i < mergetime; ++i, ++iter)
            {
                Weight_more += (*iter)->GetMSCost(false, false); // no comunication cost, ImprovedSplit never goes to here.
            }
        }

        //round++;
        //cout<<"round "<<round<<"---MS Now "<<MS_sequential+Weight_more+Weight_PQ<<", edges broken{ ";
        //        iter=parallelRoots.begin();
        //        if (parallelRoots.size()>(num_processor-1)) {
        //            target = num_processor-1;
        //        }else{
        //            target = parallelRoots.size();
        //        }
        //        if (!parallelRoots.empty()) {
        //            for (unsigned int count=0; count<target; ++iter, ++count) {
        //                cout<<(*iter)->GetId()<<" ";
        //            }
        //        }
        //        cout<<"}"<<endl;

        //cout<<"makespan "<<MS_sequential+Weight_more+Weight_PQ<<endl;
        //本次切割的MS
        MS.push_back(MS_sequential + Weight_more + Weight_PQ);
    }

    //只是计算两层划分的MS，从这返回最小的MS
    if (twolevel == true)
    {

        double makespan;
        makespan = *std::min_element(MS.begin(), MS.end());
        return makespan;
    }

    //return broken edges, i.e., root of subtrees
    vector<double>::iterator smallestMS_iter = min_element(MS.begin(), MS.end());
    //需要划分minMS_step才能得到最小的MS
    unsigned long minMS_step = smallestMS_iter - MS.begin();
    //cout<<"minMS_step "<<minMS_step<<endl;
    sequentialLength = minMS_step;
    unsigned int i = 0;
    parallelRoots.clear();
    parallelRoots.push_back(root);
    currentNode = root;
    while (i < minMS_step)
    {
        parallelRoots.remove(currentNode);

        children = currentNode->GetChildren();
        //parallelRoots==PQ
        for (vector<Cnode *>::iterator iter = children->begin(); iter != children->end(); iter++)
        {
            if (!(*iter)->IsBorken())
            {
                parallelRoots.push_back(*iter);
            }
        }

        currentNode = *max_element(parallelRoots.begin(), parallelRoots.end(), cmp_nodecreasing); //non-decreasing
        i++;
    }

    if (parallelRoots.size() > 1)
    {
        amountSubtrees = parallelRoots.size() + 1;
    }
    else
    {
        amountSubtrees = 1;
    }

    if (amountSubtrees > num_processor)
    {
        parallelRoots.sort(cmp_noIn_noCommu); //non-increasing sort, computation weight, no communication cost
        mergetime = amountSubtrees - num_processor;
        for (unsigned int i = 0; i < mergetime; ++i)
        {
            parallelRoots.pop_back();
        }
    }

    root->BreakEdge(); //root should always be broken
    for (list<Cnode *>::iterator iter = parallelRoots.begin(); iter != parallelRoots.end(); ++iter)
    {
        (*iter)->BreakEdge();
    }

    //    cout<<"   broken edges: ";
    //    for (list<Cnode*>::iterator iter=parallelRoots.begin(); iter!=parallelRoots.end(); ++iter) {
    //        cout<<(*iter)->GetId()<<" ";
    //        if (!(*iter)->IsBorken()) {
    //            cout<<"(error) ";
    //        }
    //    }
    //    cout<<endl;

    //    cout<<"makespan from the tree root "<<root->GetMSCost(true,true)<<endl;

    return *smallestMS_iter;
}

void ISCore(Cnode *root, unsigned long num_processors, bool sequentialPart)
{ //number of processors here assumed to the same as tree'size
    // cout<<"this is ISCore "<<endl;
    list<Cnode *> parallelRoots;
    double MS_before;
    double MS_now;
    unsigned long SF_now; //avoid dead lock
    double makespan;

    if (root->IsLeaf())
    {
        //cout<<"root is leaf, return."<<endl;
        return;
    }

    makespan = SplitSubtrees(root, num_processors, false, parallelRoots, SF_now); //SF_now will be modified in SplitSubtrees, it represents the length of sequential part, 0 means the subtree no need to partition

    if (sequentialPart == true)
    {
        if (SF_now == 0)
        {
            //cout<<"this subtree has already been fully checked, return."<<endl;
            return;
        }
    }

    parallelRoots.sort(cmp_noincreasing); //non-increasing sort, communication counted

    Cnode *frontNode;
    if (parallelRoots.size() > 1)
    { //==1 means there is no parallel part
        while (true)
        {
            frontNode = parallelRoots.front();
            parallelRoots.pop_front();
            MS_before = frontNode->GetMSCost(true, false);

            //cout<<"---ISCore works on Parallel root "<<frontNode->GetId()<<endl;
            ISCore(frontNode, num_processors, false);

            MS_now = frontNode->GetMSCost(true, true); //Makespan updated enforced

            if (MS_now >= MS_before)
            {
                break;
            }

            if (parallelRoots.empty())
            {
                break;
            }

            if (parallelRoots.front()->GetMSCost(true, false) <= MS_now)
            {
                break;
            }
        }

        //cout<<"---ISCore works on Sequential root "<<root->GetId()<<endl;
        ISCore(root, num_processors, true);
    }

    return;
}

double ImprovedSplit(Ctree *tree, unsigned int number_processor, int *chstart, int *childrenID)
{
    //double ImprovedSplit(Ctree* tree, unsigned int number_processor){
    // cout<<"this is ImprovedSplit "<<endl;
    unsigned long tree_size = tree->GetNodes()->size();
    Cnode *root = tree->GetRoot();
    //cout<<"---ISCore works on the root"<<endl;
    clock_t c = clock();
    ISCore(root, tree_size, false);
    c = clock() - c;
    // cout<<"IScore time:"<<c<<endl;
    //    Ctree* Qtreeobj = BuildQtree(tree);
    //    long index=Qtreeobj->GetNodes()->size()-number_processor;
    //    if (index>0) {
    //        list<Cnode*> C;
    //        C.assign(Qtreeobj->GetNodes()->begin(),Qtreeobj->GetNodes()->end());
    //        C.pop_front();//pop up the root
    //        C.sort(cmp_c_nodecrea);
    //        list<Cnode*>::iterator iter=C.begin();
    //        while (index>0) {
    //            tree->GetNode((*iter)->GetothersideID())->RestoreEdge();
    //            advance(iter, 1);
    //            index--;
    //        }
    //    }

    unsigned int numberSubtrees = HowmanySubtrees(tree, true);
    double makespan = Merge(tree, numberSubtrees, number_processor, 0, chstart, childrenID, false);

    //double makespan = root->GetMSCost(true,true);
    //delete Qtreeobj;
    return makespan;
}

//double ImprovedSplit(Ctree* tree, unsigned int processor_number){//first implementation
//    unsigned long tree_size=tree->GetNodes()->size();
//    Cnode* root=tree->GetRoot();
//    ISCore(root, tree_size, false, 0);
//
////    cout<<"Broken Edges: ";
////    for (unsigned int i=1; i<=tree_size; ++i) {
////        if (tree->GetNode(i)->IsBorken()) {
////            cout<<i<<" ";
////        }
////    }
////    cout<<endl;
//
//    unsigned int num_subtrees=0;
//    num_subtrees=HowmanySubtrees(tree,true);
//
//    if (processor_number>=num_subtrees) {
//        return root->GetMSCost(true, true);
//    }
//
//    Ctree* Qtreeobj = BuildQtree(tree);//makespan will also be updated in BuildQtree //root->GetMSCost(true, true);
//    Cnode* currentNode;
//    double temp;
//    unsigned int shortage=num_subtrees-processor_number;
//    unsigned int subtree_root_id;
//    double smallestParaPart;
//    vector<Cnode*>* Children;
//    while (shortage>0) {
//        currentNode=Qtreeobj->GetRoot();
//        smallestParaPart=currentNode->GetMSCost(true, true);
//        while (!currentNode->IsLeaf()) {
//            Children=currentNode->GetChildren();
//            for (vector<Cnode*>::iterator iter=Children->begin(); iter!=Children->end(); iter++) {
//                temp=(*iter)->GetMSCost(true, false);
//                if (temp<smallestParaPart) {
//                    smallestParaPart=temp;
//                    currentNode=(*iter);
//                }
//            }
//        }//now current node is the leaf node and is the smallest one from the smallest one iteratively
//
//        subtree_root_id = currentNode->GetothersideID();
//        tree->GetNode(subtree_root_id)->RestoreEdge();
//        currentNode->MergetoParent();
//        shortage--;
//
////        for (unsigned int i=1; i<=Qtreeobj->GetNodes()->size(); ++i) {
////            cout<<i<<" "<<Qtreeobj->GetNode(i)->GetParentId()<<" "<<Qtreeobj->GetNode(i)->GetMSW()<<" "<<Qtreeobj->GetNode(i)->GetEW()<<endl;
////        }
////        cout<<"-----------------------------"<<endl;
////        cout<<"Restore edge "<<subtree_root_id<<", shortage "<<shortage<<endl;
//    }
//
//    double Makespan=root->GetMSCost(true, true);
//
//    delete Qtreeobj;
//
//    return Makespan;
//}

//MemoryEnough(tree, (*smallest)->GetParent(), (*smallest), leaf, memory_size, chstart, childrenID);
// MemoryEnough(tree, currentQNode->GetParent(), currentQNode, leaf, memory_size, chstart, childrenID);
bool MemoryEnough(Ctree *tree, Cnode *Qrootone, Cnode *Qroottwo, bool leaf, double memory_size, int *chstart, int *children)
{
    bool enough = false;
    unsigned long new_tree_size = tree->GetNodes()->size();

    Cnode *SubtreeRoot = tree->GetNode(Qrootone->GetothersideID());

    vector<Cnode *> *childrenvector = Qrootone->GetChildren();
    //如果该子树（MS最小）是叶子结点且只有一个兄弟，则将该子树、兄弟子树与它的父结点合并
    if ((leaf == true) & (childrenvector->size() == 2))
    {
        tree->GetNode(childrenvector->front()->GetothersideID())->RestoreEdge();
        tree->GetNode(childrenvector->back()->GetothersideID())->RestoreEdge();
    }
    else
    {
        //否则，将该子树与父结点合并
        tree->GetNode(Qroottwo->GetothersideID())->RestoreEdge(); //restore edge temporarilly
    }

    //clock_t time;
    //time=clock();

    double *ewghts, *timewghts, *spacewghts;
    int *prnts;
    Ctree *subtree = BuildSubtree(tree, SubtreeRoot, new_tree_size, &prnts, &ewghts, &timewghts, &spacewghts, chstart, children);
    delete[] ewghts;
    delete[] timewghts;
    delete[] spacewghts;
    delete[] prnts;

    //time=clock()-time;
    //time=time/10000;
    //printf("  building subtree took me %d *10^4 clicks. \n",time);

    //time=clock();
    double maxout, requiredMemory;
    uint64_t count = 0;
    schedule_t *schedule_f = new schedule_t(); //list<int>
    maxout = MaxOutDegree(subtree, true);
    //重新遍历该子树
    MinMem(subtree, maxout, requiredMemory, *schedule_f, true, count);

    //time=clock()-time;
    //time=time/10000;
    //printf("  memory checking took me %d *10^4 clicks. \n",time);

    if (requiredMemory <= memory_size)
    {
        enough = true;
    }

    //复原
    if ((leaf == true) & (childrenvector->size() == 2))
    {
        tree->GetNode(childrenvector->front()->GetothersideID())->BreakEdge();
        tree->GetNode(childrenvector->back()->GetothersideID())->BreakEdge();
    }
    else
    {
        tree->GetNode(Qroottwo->GetothersideID())->BreakEdge();
    }

    delete subtree;
    delete schedule_f;

    return enough;
}

///Qtree corresponds to a whole original tree  创建商树
Ctree *BuildQtree(Ctree *tree)
{ //Qtree is for makespan side, so do not use it for space side
    Cnode *root = tree->GetRoot();
    root->BreakEdge();
    tree->GetRoot()->GetMSCost(true, true); //update
    size_t tree_size = tree->GetNodes()->size();
    unsigned long num_subtrees = HowmanySubtrees(tree, true);

    int *prnts = new int[num_subtrees + 1];
    double *ewghts = new double[num_subtrees + 1];
    double *timewghts = new double[num_subtrees + 1];
    int *brokenEdges = new int[num_subtrees + 1];

    //creat Quotient tree
    brokenEdges[1] = 1; //root node
    prnts[1] = 0;
    ewghts[1] = 0;
    timewghts[1] = root->GetSequentialPart();
    unsigned int j = 2;
    root->SetothersideID(1);

    Cnode *currentNode;
    //ℹ=1为根结点  所以根结点的孩子结点从2开始
    for (unsigned int i = 2; i <= tree_size; ++i)
    {
        currentNode = tree->GetNode(i);
        if (currentNode->IsBorken())
        {
            currentNode->SetothersideID(j); //corresponding node's ID on Qtree
            brokenEdges[j] = i;
            timewghts[j] = currentNode->GetSequentialPart();
            ewghts[j] = currentNode->GetEW();
            ++j;
        }
    }
    //为什么i=2开始？？  ℹ=1为根结点  所以根结点的孩子结点从2开始，包含根结点的子树GetothersideID=1
    //
    for (unsigned int i = 2; i <= num_subtrees; ++i)
    {
        currentNode = tree->GetNode(brokenEdges[i])->GetParent();
        while (!currentNode->IsBorken())
        {
            currentNode = currentNode->GetParent();
        }
        prnts[i] = currentNode->GetothersideID();
    }

    Ctree *Qtreeobj = new Ctree(num_subtrees, prnts, timewghts, ewghts, timewghts); //Qtree only reprents makespan, not memory consumption

    for (unsigned int i = 1; i <= num_subtrees; i++)
    {
        Qtreeobj->GetNode(i)->BreakEdge();                    //break edge
        Qtreeobj->GetNode(i)->SetothersideID(brokenEdges[i]); //corresponding node's ID on tree
    }

    delete[] prnts;
    delete[] ewghts;
    delete[] timewghts;
    delete[] brokenEdges;

    return Qtreeobj;
}

bool increaseMS(Ctree *tree, Ctree *Qtree, Cnode *&smallestNode, int *chstart, int *childrenID, double memory_size, bool CheckMemory)
{
    //cout<<"   ---start compute the minimum combination"<<endl;
    //vector<Cnode*> que;
    //vector<Cnode*>* children=Qrooot->GetChildren();
    //que.insert(que.end(), children->begin(),children->end());

    Cnode *currentNode;
    double diff, increase, temp;
    bool memoryEnough;
    bool feasible = false;
    Cnode *LargestNode;
    Cnode *secondLargest;
    Cnode *parent;
    double smallestIncrease = tree->GetRoot()->GetMSCost(true, false);
    bool leaf = false;
    const vector<Cnode *> *subtrees = Qtree->GetNodes();
    vector<Cnode *> *children;

    if (subtrees->front()->GetId() != 1)
    {
        cout << "error in function increaseMs" << endl;
        return false;
    }

    vector<Cnode *>::const_iterator iter = subtrees->begin();
    ++iter;
    for (; iter != subtrees->end(); ++iter)
    {
        currentNode = (*iter);

        if (tree->GetNode(currentNode->GetothersideID())->IsBorken() == true)
        { //this subtree has not been merged yet
            children = currentNode->GetChildren();

            if (children->empty())
            {
                leaf = true;
            }

            //check the memory cost, if merge itself to its parent subtree
            if (CheckMemory == true)
            {
                memoryEnough = MemoryEnough(tree, currentNode->GetParent(), currentNode, leaf, memory_size, chstart, childrenID);
            }
            else
            {
                memoryEnough = true;
            }

            //cout<<"   subtree "<<currentNode->GetothersideID()<<" ";//print subtree's root id
            if (memoryEnough == true)
            {
                feasible = true;
                //cout<<"memory fit."<<endl;
                if (children->empty())
                {
                    children = currentNode->GetParent()->GetChildren();
                    if (children->size() == 2)
                    {
                        increase = children->front()->GetMSCost(false, false) + children->back()->GetMSCost(false, false) - currentNode->GetParent()->GetParallelPart();
                    }
                    else if (children->size() == 1)
                    {
                        increase = -currentNode->GetEW() / BANDWIDTH;
                    }
                    else
                    {
                        //cout<<"   current subtree "<<currentNode->GetothersideID()<<", parent id "<<currentNode->GetParent()->GetothersideID()<<", number of siblings "<<children->size()-1<<endl;

                        GetTwoLargestElementTypetwo(children, LargestNode, secondLargest); //no-increasing, communication counted
                        if (currentNode->GetMSCost(true, false) == LargestNode->GetMSCost(true, false))
                        {
                            if (currentNode->GetMSCost(true, false) == secondLargest->GetMSCost(true, false))
                            {
                                increase = currentNode->GetMSW();
                            }
                            else
                            {
                                increase = -currentNode->GetEW() / BANDWIDTH + secondLargest->GetMSCost(true, false);
                            }
                        }
                        else
                        {
                            increase = currentNode->GetMSW();
                        }
                    }
                }
                else
                {
                    children = currentNode->GetParent()->GetChildren();
                    //cout<<"   current subtree "<<currentNode->GetothersideID()<<", parent id "<<currentNode->GetParent()->GetothersideID()<<", number of siblings "<<children->size()-1<<endl;
                    diff = currentNode->GetMSCost(true, false) - currentNode->GetParent()->GetParallelPart();
                    if (diff < 0)
                    {
                        increase = currentNode->GetMSW();
                    }
                    else
                    { //diff=0
                        if (children->size() == 1)
                        {
                            increase = -currentNode->GetEW() / BANDWIDTH;
                        }
                        else
                        {                                                                      //children's size larger than 1
                            GetTwoLargestElementTypetwo(children, LargestNode, secondLargest); //no-increasing, communication counted
                            temp = currentNode->GetParallelPart() - secondLargest->GetMSCost(true, false);
                            if (temp >= 0)
                            {
                                increase = -currentNode->GetEW() / BANDWIDTH;
                            }
                            else
                            {
                                increase = -temp - currentNode->GetEW() / BANDWIDTH;
                            }
                        }
                    }
                }

                parent = currentNode->GetParent();
                while (parent->GetId() != 1)
                { //not the root node
                    temp = parent->GetParent()->GetParallelPart() - (parent->GetMSCost(true, false) + increase);
                    if (temp >= 0)
                    {
                        increase = 0;
                        break;
                    }
                    else
                    {
                        increase = -temp;
                        parent = parent->GetParent();
                    }
                }

                //cout<<"   merge, increase in MS(r) "<<increase<<endl;
                if (increase < smallestIncrease)
                {
                    smallestIncrease = increase;
                    smallestNode = currentNode;
                }
            }
            else
            {
                //cout<<"memory does not fit!!!"<<endl;
            }
        }
    }

    //cout<<"   ---end compute the minimum combination"<<endl;
    return feasible;
}

bool cmp_merge_smallest(const pair<double, Cnode *> &a, const pair<double, Cnode *> &b) { return a.first < b.first; };

bool estimateMS(Ctree *tree, Ctree *Qtree, Cnode *&smallestNode, int *chstart, int *childrenID, double memory_size, bool CheckMemory)
{
    //cout<<"   ---start compute the minimum combination"<<endl;

    Cnode *currentQNode;
    double increase;
    bool memoryEnough;
    Cnode *LargestNode;
    Cnode *secondLargest;
    bool leaf = false;
    const vector<Cnode *> *subtrees = Qtree->GetNodes();
    vector<Cnode *> *children;

    if (subtrees->front()->GetId() != 1)
    { //the root is supposed to be the first element in vector nodes
        cout << "error in function estimateMS" << endl;
        return false;
    }

    vector<Cnode *> tempQue;
    currentQNode = Qtree->GetRoot();
    //SetMSDiff   ==slack (di)the threshold such that MS(r) is not impacted by the increase of MS(i) up to MS(i) + di .
    currentQNode->SetMSDiff(0);
    children = currentQNode->GetChildren();
    tempQue.insert(tempQue.end(), children->begin(), children->end());
    //cout<<"   ---compute makespan difference---"<<endl;
    while (!tempQue.empty())
    {
        currentQNode = tempQue.back();
        tempQue.pop_back();

        currentQNode->SetMSDiff(currentQNode->GetParent()->GetMSDiff() + currentQNode->GetParent()->GetParallelPart() - currentQNode->GetMSCost(true, false));
        //cout<<"   subtree "<<currentQNode->GetothersideID()<<", makespan difference: "<<currentQNode->GetMSDiff()<<endl;
        children = currentQNode->GetChildren();
        tempQue.insert(tempQue.end(), children->begin(), children->end());
    }
    //cout<<"   ----------------------------------"<<endl;

    list<pair<double, Cnode *>> list_increase_id;
    vector<Cnode *>::const_iterator iter = subtrees->begin();
    ++iter;
    unsigned long size = subtrees->size() - 1;
    //  #pragma omp parallel for
    for (unsigned int step = 0; step < size; ++step)
    {
        currentQNode = *(iter + step);

        if (tree->GetNode(currentQNode->GetothersideID())->IsBorken() == true)
        { //this subtree has not been merged yet
            children = currentQNode->GetChildren();
            if (children->empty())
            { //this is a leaf node
                children = currentQNode->GetParent()->GetChildren();
                if (children->size() == 2)
                { //叶子结点且只有一个兄弟，则将该结点和兄弟结点与父结点合并，增量为
                    increase = children->front()->GetMSW() + children->back()->GetMSW() - currentQNode->GetParent()->GetParallelPart();
                }
                else if (children->size() == 1)
                {
                    increase = -currentQNode->GetEW() / BANDWIDTH;
                }
                else
                {

                    GetTwoLargestElementTypetwo(children, LargestNode, secondLargest); //no-increasing, communication counted

                    if (currentQNode->GetId() == LargestNode->GetId())
                    {
                        increase = currentQNode->GetMSW() + secondLargest->GetMSCost(true, false) - currentQNode->GetMSCost(true, false);
                    }
                    else
                    {
                        increase = currentQNode->GetMSW();
                    }
                }
            }
            else
            { //not a leaf node
                children = currentQNode->GetParent()->GetChildren();

                if (children->size() == 1)
                {
                    increase = -currentQNode->GetEW() / BANDWIDTH;
                }
                else
                {
                    GetTwoLargestElementTypetwo(children, LargestNode, secondLargest); //no-increasing, communication counted

                    if (currentQNode->GetId() == LargestNode->GetId())
                    {
                        increase = currentQNode->GetMSW() + max(secondLargest->GetMSCost(true, false), currentQNode->GetParallelPart()) - currentQNode->GetMSCost(true, false);
                    }
                    else
                    {
                        increase = currentQNode->GetMSW();
                    }
                }
            }

            //now consider the "slack" between the parent of current node with siblings of its parent
            increase = increase - currentQNode->GetParent()->GetMSDiff();

            //cout<<"merge, increase in MS(r) "<<increase<<endl;
            list_increase_id.push_back(pair<double, Cnode *>(increase, currentQNode));
        }
    }

    bool feasible = false;
    list<pair<double, Cnode *>>::iterator smallest_iter;
    //选择最小的增量的合并，并检查内存可否匹配，如果内存不够，则放弃该合并，选择增量次小的
    while (feasible == false && list_increase_id.empty() == false)
    {
        smallest_iter = min_element(list_increase_id.begin(), list_increase_id.end(), cmp_merge_smallest);
        currentQNode = (*smallest_iter).second;
        //cout<<"   increase in MS(r) estimated: "<<(*smallest_iter).first<<endl;

        children = currentQNode->GetChildren();
        if (children->empty())
        {
            leaf = true;
        }

        if (CheckMemory == true)
        {
            //检查currentQNode与父结点合并时，内存是否足够
            memoryEnough = MemoryEnough(tree, currentQNode->GetParent(), currentQNode, leaf, memory_size, chstart, childrenID);
        }
        else
        {
            memoryEnough = true;
        }

        if (memoryEnough == true)
        {
            feasible = true;
            smallestNode = currentQNode;
        }
        else
        {
            list_increase_id.erase(smallest_iter);
        }
    }

    //cout<<"   ---end compute the minimum combination"<<endl;
    return feasible;
}

double Merge(Ctree *tree, unsigned int num_subtrees, unsigned int processor_number, double const memory_size, int *chstart, int *childrenID, bool CheckMemory)
{
    Cnode *root = tree->GetRoot();

    if (processor_number >= num_subtrees)
    {
        return root->GetMSCost(true, true);
    }

    Ctree *Qtreeobj = BuildQtree(tree);

    Cnode *node_smallest_increase;
    Cnode *parent;
    int shortage = num_subtrees - processor_number;
    double temp;
    Cnode *nodeone;
    Cnode *nodetwo;
    bool memoryEnough;

    while (shortage > 0)
    { //merge subtree
        //cout<<"shortage "<<shortage<<endl;
        temp = Qtreeobj->GetRoot()->GetMSCost(true, true); //initilize ms
        temp = tree->GetRoot()->GetMSCost(true, true);     //update ms

        //memoryEnough=increaseMS(tree, Qtreeobj, node_smallest_increase, chstart, childrenID, memory_size, CheckMemory);
        memoryEnough = estimateMS(tree, Qtreeobj, node_smallest_increase, chstart, childrenID, memory_size, CheckMemory);

        //when parameter checkMemory is false, memoryEnough will always be true;
        if (memoryEnough == true)
        {
            //merge currentNode (or and its sibling) to its parent
            if (node_smallest_increase->IsLeaf())
            {
                parent = node_smallest_increase->GetParent();
                if (parent->GetChildren()->size() == 2)
                {
                    nodeone = parent->GetChildren()->front();
                    nodetwo = parent->GetChildren()->back();
                    //cout<<"Merge node "<<nodeone->GetothersideID()<<" and its sibling "<<nodetwo->GetothersideID()<<endl;
                    nodeone->MergetoParent();
                    nodetwo->MergetoParent();
                    shortage = shortage - 2;
                    tree->GetNode(nodeone->GetothersideID())->RestoreEdge();
                    tree->GetNode(nodetwo->GetothersideID())->RestoreEdge();
                }
                else
                {
                    //cout<<"Merge node "<<node_smallest_increase->GetothersideID()<<endl;
                    node_smallest_increase->MergetoParent();
                    shortage--;
                    tree->GetNode(node_smallest_increase->GetothersideID())->RestoreEdge();
                }
            }
            else
            {
                //cout<<"Merge node "<<node_smallest_increase->GetothersideID()<<endl;
                node_smallest_increase->MergetoParent();
                shortage--;
                tree->GetNode(node_smallest_increase->GetothersideID())->RestoreEdge();
            }
            //cout<<"------------------------"<<endl;
        }
        else
        {
            break;
        }
    }

    if (shortage > 0)
    { //failure
        temp = -1;
    }
    else
    {
        temp = root->GetMSCost(true, true);
    }
    delete Qtreeobj;

    return temp;
}
//当处理器数目大于子树数目时，合并子树
//与Merge的区别：MergeV2先选择非关键路径上的第二层以上结点进行合并
//merge 对所有商树结点进行试探合并，选择增量最小且内存匹配的合并
double MergeV2(Ctree *tree, unsigned int num_subtrees, unsigned int processor_number, double const memory_size, int *chstart, int *childrenID, bool CheckMemory)
{
    if (processor_number >= num_subtrees)
    {
        return tree->GetRoot()->GetMSCost(true, true);
    }

    tree->GetRoot()->GetMSCost(true, true); //update makespan

    Ctree *Qtreeobj = BuildQtree(tree); //建立商树

    Cnode *currentNode;
    Cnode *Qroot = Qtreeobj->GetRoot();
    //shortage缺少处理器的个数
    int shortage = num_subtrees - processor_number;
    list<Cnode *> Llist;
    vector<unsigned int> CriticalPath;
    double temp;
    Cnode *largestNode;
    //list<Cnode*>::iterator largest;
    //list<Cnode*>::iterator secondLargest;
    list<Cnode *>::iterator smallest;
    list<Cnode *>::iterator secondSmallest;
    vector<Cnode *> *Children;
    long pathlength;
    vector<Cnode *> queue;
    bool memoryCheckPass = false, leaf = false;
    Cnode *nodeone;
    Cnode *nodetwo;
    bool DeadBreak, firstTime;

    //clock_t time;
    //time = clock();
    while (shortage > 0)
    {
        DeadBreak = true;
        firstTime = true;
        Llist.clear();
        CriticalPath.clear();
        temp = Qroot->GetMSCost(true, true);           //update ms
        temp = tree->GetRoot()->GetMSCost(true, true); //update ms  ？？？为什么重复

        CriticalPath.push_back(1); //关键路径
        largestNode = Qroot;
        Children = largestNode->GetChildren();

        //time=clock();
        //  关键路径： 找到当前结点的孩子结点并行处理的MS最大的结点，即为关键路径上的元素，再从该结点出发，寻找下一个关键结点
        while (!Children->empty())
        { //initialize critical path
            temp = largestNode->GetParallelPart();
            for (vector<Cnode *>::iterator iter = Children->begin(); iter != Children->end(); ++iter)
            {
                if ((*iter)->GetMSCost(true, false) == temp)
                {
                    largestNode = (*iter);
                    break;
                }
            }
            CriticalPath.push_back(largestNode->GetId());
            Children = largestNode->GetChildren();
        }
        //time=clock()-time;
        //time=time/10000;
        //printf("  initializing critical path took me %d *10^4 clicks. \n",time);

        //        cout<<"Critical Path: ";
        //        for (vector<unsigned int>::iterator iter=CriticalPath.begin(); iter!=CriticalPath.end(); ++iter) {
        //            cout<<(*iter)<<" ";
        //        }
        //        cout<<endl;
        //        cout<<endl;

        Children = Qroot->GetChildren();
        pathlength = CriticalPath.size();

        //time=clock();
        //queue为商树关键路径上的结点的    【非关键路径】上的孩子结点，即第一层非关键路径上的所有结点
        for (unsigned int i = 1; i < pathlength; ++i)
        { //initialize vector L
            for (vector<Cnode *>::iterator iter = Children->begin(); iter != Children->end(); ++iter)
            {
                if ((*iter)->GetId() != CriticalPath[i])
                {
                    queue.push_back(*iter);
                }
            }
            Children = Qtreeobj->GetNode(CriticalPath[i])->GetChildren();
        }

        //Llist  商树关键路径上的结点的    【非关键路径】上的孩子结点 的所有孩子结点，即第二层及以上非关键路径上的所有结点
        //queue 非关键路径上的所有结点
        while (!queue.empty())
        { //initialize vector L
            currentNode = queue.back();
            queue.pop_back();
            Children = currentNode->GetChildren();
            for (vector<Cnode *>::iterator iter = Children->begin(); iter != Children->end(); ++iter)
            {
                Llist.push_back(*iter);
                queue.push_back(*iter);
            }
        }

    //当Llist为空时，list 根结点的所有孩子结点
    CheckOnCritical:
        if (Llist.empty())
        {
            queue.push_back(Qroot);
            while (!queue.empty())
            { //initialize vector L
                currentNode = queue.back();
                queue.pop_back();
                Children = currentNode->GetChildren();
                for (vector<Cnode *>::iterator iter = Children->begin(); iter != Children->end(); ++iter)
                {
                    Llist.push_back(*iter);
                    queue.push_back(*iter);
                }
            }
        }
        //time=clock()-time;
        //time=time/10000;
        //printf("  initializing L took me %d *10^4 clicks. \n",time);

        //        cout<<"List L: ";
        //        for (list<Cnode*>::iterator iter=Llist.begin(); iter!=Llist.end(); ++iter) {
        //            cout<<(*iter)->GetId()<<" ";
        //        }
        //        cout<<endl;
        //        cout<<endl;

        //cout<<"List L size: "<<Llist.size()<<", shortage: "<<shortage<<endl;
        do
        {
            /**
             * 每次选择对最小的子树进行合并操作，并检查是否内存是否足够。
             * 若最小子树内存足够，则进行合并操作； 若合并后内存不够，继续对第二小的子树进行合并操作。
             * 如果第二小的子树内存不够，则将该子树从Llist中删除；若第二小的子树内存足够，则进行合并操作
             * 
             * 直到内存足够或 Llist为空
            */
            //Llist  商树上的结点
            if (Llist.size() == 1)
            {
                secondSmallest = Llist.begin();
                smallest = Llist.begin();
            }
            else
            {
                //找到商树中最小的两个结点
                GetTwoSmallestElement(&Llist, smallest, secondSmallest);
            }

            if ((*smallest)->IsLeaf())
            {
                leaf = true;
            }
            else
            {
                leaf = false;
            }

            if (CheckMemory == true)
            {
                memoryCheckPass = MemoryEnough(tree, (*smallest)->GetParent(), (*smallest), leaf, memory_size, chstart, childrenID);
            }
            else
            {
                memoryCheckPass = true;
            }

            if (memoryCheckPass == false)
            {
                if ((*secondSmallest)->IsLeaf())
                {
                    leaf = true;
                }
                else
                {
                    leaf = false;
                }
                memoryCheckPass = MemoryEnough(tree, (*secondSmallest)->GetParent(), *secondSmallest, leaf, memory_size, chstart, childrenID);
                if (memoryCheckPass == true)
                {
                    currentNode = *secondSmallest;
                    DeadBreak = false;
                }
                else
                {
                    Llist.erase(secondSmallest);
                }
            }
            else
            {
                currentNode = *smallest;
                DeadBreak = false;
            }
        } while ((memoryCheckPass == false) && (!Llist.empty()));

        if (DeadBreak == true && firstTime == true)
        {
            Llist.clear();
            firstTime = false;
            goto CheckOnCritical;
        }

        if (DeadBreak == true)
        {
            delete Qtreeobj;
            return -1; //which means failure
        }

        //merge currentNode (or and its sibling) to its parent
        if (currentNode->IsLeaf())
        {
            /**
             * 如果当前结点是叶子结点且只有一个兄弟，则将它们与父结点合并
             * 否则该子树与父结点合并
            */
            if (currentNode->GetParent()->GetChildren()->size() == 2)
            {
                nodeone = currentNode->GetParent()->GetChildren()->front();
                nodetwo = currentNode->GetParent()->GetChildren()->back();
                //cout<<"Merge node "<<nodeone->GetId()<<"-"<<" and its sibling "<<nodetwo->GetId()<<"-"<<endl;
                nodeone->MergetoParent();
                nodetwo->MergetoParent();
                shortage = shortage - 2;
                tree->GetNode(nodeone->GetothersideID())->RestoreEdge();
                tree->GetNode(nodetwo->GetothersideID())->RestoreEdge();
            }
            else
            {
                //cout<<"Merge node "<<currentNode->GetId()<<"-"<<endl;
                //time=clock();
                currentNode->MergetoParent();
                //time=clock()-time;
                //time=time/10000;
                //printf("  merge node to its parent took me %d *10^4 clicks. \n",time);
                shortage--;
                tree->GetNode(currentNode->GetothersideID())->RestoreEdge();
            }
        }
        else
        {
            //cout<<"Merge node "<<currentNode->GetId()<<"-"<<endl;
            //当前结点是非叶结点则与父结点合并
            currentNode->MergetoParent();
            shortage--;
            tree->GetNode(currentNode->GetothersideID())->RestoreEdge();
        }
        //cout<<"------------------------"<<endl;
    }

    temp = tree->GetRoot()->GetMSCost(true, true);
    delete Qtreeobj;

    return temp;
}

bool cmp_asapc(Cnode *a, Cnode *b) { return (a->GetMSminusComu() < b->GetMSminusComu()); };

double ASAP(Ctree *tree, unsigned int num_processors, unsigned int depth)
{ //depth should be at least 1
    list<Cnode *> PriorityQue;
    vector<Cnode *> BrokenEdges;
    unsigned long step_minimumMS = 0;
    double minimumMS = tree->GetRoot()->GetMSCost(true, true);
    //cout<<"Excuting sequentially, makespan "<<minimumMS<<endl;
    double temp;
    list<Cnode *> children_buffer;
    Cnode *currentNode;
    unsigned int add_child_deepth;
    list<Cnode *>::iterator node_position;
    unsigned int i;
    Cnode *LargestNode;

    vector<Cnode *> *children = tree->GetRoot()->GetChildren();

    while (children->size() == 1)
    { //avoid the linear chain
        children = children->front()->GetChildren();
    }

    for (vector<Cnode *>::iterator iter = children->begin(); iter != children->end(); ++iter)
    {
        PriorityQue.push_back(*iter);
        //cout<<"   add node "<<(*iter)->GetId()<<" into PQ."<<endl;
        children_buffer.push_back(*iter);
    }

    unsigned long buffer_size;
    unsigned int generation = 1;
    while (generation < depth)
    {
        buffer_size = children_buffer.size();
        while (buffer_size > 0)
        {
            children = children_buffer.front()->GetChildren();
            children_buffer.pop_front();
            for (vector<Cnode *>::iterator child = children->begin(); child != children->end(); ++child)
            {
                children_buffer.push_back(*child);
                PriorityQue.push_back(*child);
                //cout<<"   add node "<<(*child)->GetId()<<" into PQ."<<endl;
            }
            --buffer_size;
        }
        ++generation;
    }

    while (num_processors > 1)
    { //Breaking an edge a time
        if (PriorityQue.empty())
        { //when having more processors than nodes
            break;
        }

        if (PriorityQue.size() == 1)
        {
            LargestNode = PriorityQue.front();
        }
        else
        {
            LargestNode = *max_element(PriorityQue.begin(), PriorityQue.end(), cmp_asapc); //computation weight minus communication cost
        }

        if (LargestNode->GetParent()->GetChildren()->size() > 1)
        {
            LargestNode->BreakEdge(); //break edge
            BrokenEdges.push_back(LargestNode);
            temp = tree->GetRoot()->GetMSCost(true, true);
            //cout<<"Break edge "<<LargestNode->GetId()<<", makespan now: "<<temp;
            num_processors--;
            if (temp < minimumMS)
            {
                minimumMS = temp;
                step_minimumMS = BrokenEdges.size();
                //cout<<", makespan decreased";
            }
            //cout<<endl;
        }

        currentNode = LargestNode;
        node_position = find(PriorityQue.begin(), PriorityQue.end(), currentNode);
        PriorityQue.erase(node_position);
        //cout<<"   pop up node "<<currentNode->GetId()<<endl;
        children_buffer.clear();
        children_buffer.push_back(currentNode);
        add_child_deepth = 1;

        node_position = find(PriorityQue.begin(), PriorityQue.end(), currentNode->GetParent());
        while (node_position != PriorityQue.end())
        {
            PriorityQue.erase(node_position);
            ++add_child_deepth;
            currentNode = currentNode->GetParent();
            node_position = find(PriorityQue.begin(), PriorityQue.end(), currentNode->GetParent());
        }

        i = 0;
        while (i < depth - add_child_deepth)
        { //avoid insert nodes that are already in PQ, when depth is larger than 1
            ++i;
            buffer_size = children_buffer.size();
            while (buffer_size > 0)
            {
                children = children_buffer.front()->GetChildren();
                children_buffer.pop_front();
                for (vector<Cnode *>::iterator child = children->begin(); child != children->end(); ++child)
                {
                    children_buffer.push_back(*child);
                }
                --buffer_size;
            }
        }

        while (add_child_deepth >= 1)
        {
            buffer_size = children_buffer.size();
            while (buffer_size > 0)
            {
                children = children_buffer.front()->GetChildren();
                children_buffer.pop_front();
                for (vector<Cnode *>::iterator child = children->begin(); child != children->end(); ++child)
                {
                    children_buffer.push_back(*child);
                    PriorityQue.push_back(*child);
                    //cout<<"   add node "<<(*child)->GetId()<<" into PQ."<<endl;
                }
                --buffer_size;
            }
            --add_child_deepth;
        }
    }

    //cout<<"resotre edge ";
    unsigned long restore_index = BrokenEdges.size();
    while (restore_index > step_minimumMS)
    {
        BrokenEdges[restore_index - 1]->RestoreEdge();
        //cout<<BrokenEdges[restore_index-1]->GetId()<<" ";
        restore_index--;
    }
    //cout<<endl;

    return minimumMS;
}

bool cmp_asap(Cnode *a, Cnode *b) { return (a->GetMSCost(false, false) < b->GetMSCost(false, false)); };

double ASAP(Ctree *tree, unsigned int num_processors)
{
    list<Cnode *> PriorityQue;
    vector<Cnode *> BrokenEdges;
    unsigned long step_minimumMS = 0;
    double minimumMS = tree->GetRoot()->GetMSCost(true, true);
    //cout<<"Excuting sequentially, makespan "<<minimumMS<<endl;
    double temp;
    Cnode *LargestNode;
    list<Cnode *>::iterator node_position;

    vector<Cnode *> *children = tree->GetRoot()->GetChildren();
    while (children->size() == 1)
    { //avoid the linear chain
        children = children->front()->GetChildren();
    }

    PriorityQue.insert(PriorityQue.end(), children->begin(), children->end());

    while (num_processors > 1)
    { //Breaking an edge a time
        if (PriorityQue.empty())
        { //when having more processors than nodes
            break;
        }

        if (PriorityQue.size() == 1)
        {
            LargestNode = PriorityQue.front();
        }
        else
        {
            LargestNode = *max_element(PriorityQue.begin(), PriorityQue.end(), cmp_asap); //computation weight, no communication
        }

        node_position = find(PriorityQue.begin(), PriorityQue.end(), LargestNode);
        PriorityQue.erase(node_position);

        if (LargestNode->GetParent()->GetChildren()->size() > 1)
        {
            LargestNode->BreakEdge(); //break edge
            BrokenEdges.push_back(LargestNode);
            temp = tree->GetRoot()->GetMSCost(true, true);
            //cout<<"Break edge "<<LargestNode->GetId()<<", makespan now: "<<temp;
            num_processors--;
            if (temp < minimumMS)
            {
                minimumMS = temp;
                step_minimumMS = BrokenEdges.size();
                //cout<<", makespan decreased";
            }
            //cout<<endl;
        }
        //cout<<"   pop up node "<<LargestNode->GetId()<<endl;
        children = LargestNode->GetChildren();
        PriorityQue.insert(PriorityQue.end(), children->begin(), children->end());
    }

    //cout<<"resotre edge ";
    unsigned long restore_index = BrokenEdges.size();
    while (restore_index > step_minimumMS)
    {
        BrokenEdges[restore_index - 1]->RestoreEdge();
        //cout<<BrokenEdges[restore_index-1]->GetId()<<" ";
        restore_index--;
    }
    //cout<<endl;

    return minimumMS;
}

//unsigned long AvoidChain(Ctree* tree){//it works, the first implementation
//    Cnode* root=tree->GetRoot();
//
//    root->BreakEdge();
//    unsigned long num_subtrees=0;
//    num_subtrees = HowmanySubtrees(tree, true);
//
//    Ctree* Qtreeobj = BuildQtree(tree);
//
//    //find the chain
//    Cnode* currentNode;
//    forward_list<Cnode*> Que;
//    Que.push_front(Qtreeobj->GetRoot());
//    vector<Cnode*>* children;
//    while (!Que.empty()) {
//        currentNode=Que.front();
//        Que.pop_front();
//
//        while (currentNode->GetChildren()->size()==1) {
//            currentNode=currentNode->GetChildren()->front();
//            tree->GetNode(currentNode->GetothersideID())->RestoreEdge();// Restore Edge
//            //cout<<"Restore edge "<<currentNode->GetothersideID()<<endl;
//            num_subtrees--;
//        }
//
//        children=currentNode->GetChildren();
//        for (vector<Cnode*>::iterator iter=children->begin();iter!=children->end();iter++){
//            Que.push_front((*iter));
//        }
//    }
//
//    delete Qtreeobj;
//
//    return num_subtrees;
//}

unsigned long AvoidChain(Ctree *tree)
{
    Cnode *root = tree->GetRoot();
    root->BreakEdge();
    unsigned long num_subtrees = 0;
    num_subtrees = HowmanySubtrees(tree, true);

    Ctree *Qtreeobj = BuildQtree(tree);
    const vector<Cnode *> *AllNodes = Qtreeobj->GetNodes();
    vector<Cnode *> *children;
    for (vector<Cnode *>::const_iterator iter = AllNodes->begin(); iter != AllNodes->end(); iter++)
    {
        children = (*iter)->GetChildren();
        if (children->size() == 1)
        {
            tree->GetNode(children->front()->GetothersideID())->RestoreEdge(); //restore edge
            //cout<<"resotre edge "<<children->front()->GetothersideID()<<endl;
            num_subtrees--;
        }
    }

    delete Qtreeobj;

    return num_subtrees;
}

bool cmp_larSav(Cnode *a, Cnode *b) { return (a->GetMSminusComu() > b->GetMSminusComu()); };

//double LarSav(Ctree* tree, unsigned int processor_number, unsigned int num_subtrees){
//    Cnode* root=tree->GetRoot();
//    if (processor_number<=num_subtrees) {
//        return root->GetMSCost(true, true);
//    }
//
//    root->GetMSCost(true, true);//update makespan
//
//    Ctree* Qtreeobj = BuildQtree(tree);
//
////    cout<<"id parentId ew msw nw"<<endl;
////    for (int i=1; i<=100; ++i) {
////        cout<<Qtreeobj->GetNode(i)->GetId()<<" "<<Qtreeobj->GetNode(i)->GetParentId()<<" "<<Qtreeobj->GetNode(i)->GetEW()<<" "<<Qtreeobj->GetNode(i)->GetMSW()<<" "<<Qtreeobj->GetNode(i)->GetNW()<<endl;
////        Qtreeobj->GetNode(i)->BreakEdge();
////    }
//
//    vector<Cnode*> CriticalPath;
//    list<Cnode*> listL;
//    list<Cnode*> que;
//    int idleProcessors=processor_number-num_subtrees;
//    Cnode* largestNode=root;
//    vector<Cnode*>* Children;
//    double temp;
//    bool childSubtreeEnd;
//    Cnode* Largest;
//    Cnode* secondLargest;
//    unsigned int tempid;
//    Cnode* currentNode;
//    Cnode* Qroot=Qtreeobj->GetRoot();
//    long QtreeSize;
////    clock_t time;
//    vector<Cnode*> tempStore;
//    vector<Cnode*> nodesLastSubtree;
//    bool DeadLock=true;
//    while (idleProcessors>0) {
//        DeadLock=true;
//        childSubtreeEnd=false;
//        CriticalPath.clear();
//        CriticalPath.push_back(root);
//        Children=Qroot->GetChildren();
//        Qroot->GetMSCost(true,true);//update makespan
//        largestNode=Qroot;
////        time = clock();
//        while (!Children->empty()) {//initialize critical path
//            temp=largestNode->GetParallelPart();
//            for (vector<Cnode*>::iterator iter=Children->begin(); iter!=Children->end(); ++iter) {
//                //cout<<(*iter)->GetId()<<" "<<endl;
//                if ((*iter)->GetMSCost(true, false)==temp) {
//                    largestNode=(*iter);
//                    break;
//                }
//            }
//            CriticalPath.push_back(tree->GetNode(largestNode->GetothersideID()));
//            Children=largestNode->GetChildren();
//        }
//
////        time=clock()-time;
////        time=time/10000;
////        printf("Build CriPat took me %d *10^4 clicks. \n",time);
//
////        cout<<endl;
////        cout<<"Critical path: ";
////        for (vector<Cnode*>::iterator iter=CriticalPath.begin(); iter!=CriticalPath.end(); ++iter) {
////            cout<<(*iter)->GetId()<<" ";
////        }
////        cout<<endl;
////        cout<<"Critical path size: "<<CriticalPath.size()<<endl;
////
////        cout<<"Idle processor number: "<<idleProcessors<<endl;
//
//        listL.clear();
////        time=clock();
//
//        QtreeSize=Qtreeobj->GetNodes()->size();
//        for (unsigned int i=2; i<=QtreeSize; ++i) {//remove ancestors of roots
//            currentNode=tree->GetNode(Qtreeobj->GetNode(i)->GetothersideID())->GetParent();
//            while (!currentNode->IsBorken()) {
//                currentNode->BreakEdge();
//                tempStore.push_back(currentNode);
//                que.push_back(currentNode);
//                currentNode=currentNode->GetParent();
//            }
//        }
//
//        for (unsigned int i=0; i<CriticalPath.size(); ++i) {//initialize set L
//            que.push_back(CriticalPath[i]);
//            while (!que.empty()) {
//                Children=que.front()->GetChildren();
//                que.pop_front();
//                for (vector<Cnode*>::iterator iter=Children->begin(); iter!=Children->end(); ++iter) {
//                    if (!(*iter)->IsBorken()) {
//                        listL.push_back((*iter));
//                        que.push_back((*iter));
//                    }
//                }
//            }
//        }
//
//        while (!tempStore.empty()) {
//            currentNode=tempStore.back();
//            currentNode->RestoreEdge();
//            tempStore.pop_back();
//        }
//
////        cout<<endl;
////        cout<<"List L: ";
////        for (list<Cnode*>::iterator iter=listL.begin(); iter!=listL.end(); ++iter) {
////            cout<<(*iter)->GetId()<<" ";
////        }
////        cout<<endl;
////        cout<<"List L size: "<<listL.size()<<endl;
//
//        if (listL.empty()) {
//            break;
//        }
//
////        time=clock()-time;
////        time=time/10000;
////        printf("Build listL took me %d *10^4 clicks. \n",time);
//
////        time=clock();
////        largestNode = *max_element(listL.begin(), listL.end(), cmp_larSav);//computation weight minus communication cost
//        listL.sort(cmp_larSav);//computation weight minus communication cost, non-increasing
////        cout<<endl;
////        for (list<Cnode*>::iterator iter=listL.begin(); iter!=listL.end(); advance(iter, 1)) {
////            cout<<(*iter)->GetId()<<"_"<<(*iter)->GetMSminusComu()<<" ";
////        }
////        cout<<endl;
////        time=clock()-time;
////        time=time/10000;
////        printf("Get max element took me %d *10^4 clicks. \n",time);
//
//        //cout<<"Nodes of the last subtree on Critical path: ";
//        Children=CriticalPath.back()->GetChildren();
//        que.clear();
//        nodesLastSubtree.clear();
//        que.assign(Children->begin(),Children->end());
//        while (!que.empty()) {
//            for (vector<Cnode*>::iterator iter=Children->begin(); iter!=Children->end(); ++iter) {
//                nodesLastSubtree.push_back((*iter));
//                que.push_back((*iter));
//                //cout<<(*iter)->GetId()<<" ";
//            }
//            Children=que.front()->GetChildren();
//            que.pop_front();
//        }
//        //cout<<endl;
//
//
//        while (!listL.empty()) {
//            tempid = listL.front()->GetId();
//            childSubtreeEnd=false;
//            for (vector<Cnode*>::iterator iter=nodesLastSubtree.begin(); iter!=nodesLastSubtree.end(); ++iter) {
//                if ((*iter)->GetId()==tempid) {
//                    childSubtreeEnd=true;
//                    break;
//                }
//            }
//
//            if (childSubtreeEnd==true) {//on the last subtree of critical path
//                currentNode=listL.front();
//                if ((idleProcessors>1)&&(currentNode->GetParent()->GetChildren()->size()>1)) {
//                    GetTwoLargestElementTypeone(currentNode->GetParent()->GetChildren(), Largest, secondLargest);
//
//                    temp=Qtreeobj->GetNode(CriticalPath.back()->GetothersideID())->GetMSW();
//                    Largest->BreakEdge();
//                    Largest->SetothersideID(Qtreeobj->GetNodes()->size()+1);
//                    secondLargest->BreakEdge();
//                    secondLargest->SetothersideID(Qtreeobj->GetNodes()->size()+2);
//
//                    //cout<<"Break edge "<<Largest->GetId()<<", "<<secondLargest->GetId()<<" of subtree "<<CriticalPath.back()->GetothersideID()<<", create Q node ";
//
//                    Cnode* newNodeone = new Cnode(CriticalPath.back()->GetothersideID(), 0, Largest->GetEW(), Largest->GetMSCost(false, false));//bug here
//                    newNodeone->SetId(Qtreeobj->GetNodes()->size()+1);
//                    //cout<<Qtreeobj->GetNodes()->size()+1<<", ";
//                    newNodeone->GetChildren()->clear();
//                    newNodeone->SetParent(Qtreeobj->GetNode(CriticalPath.back()->GetothersideID()));
//                    newNodeone->BreakEdge();
//                    newNodeone->SetothersideID(Largest->GetId());
//                    Qtreeobj->GetNode(CriticalPath.back()->GetothersideID())->AddChild(newNodeone);
//                    Qtreeobj->addNode(newNodeone);
//                    temp=temp-newNodeone->GetMSW();
//
//                    Cnode* newNodetwo = new Cnode(CriticalPath.back()->GetothersideID(), 0, secondLargest->GetEW(), secondLargest->GetMSCost(false, false));//bug here
//                    newNodetwo->SetId(Qtreeobj->GetNodes()->size()+1);
//                    //cout<<Qtreeobj->GetNodes()->size()+1<<endl;
//                    newNodetwo->GetChildren()->clear();
//                    newNodetwo->SetParent(Qtreeobj->GetNode(CriticalPath.back()->GetothersideID()));
//                    newNodetwo->BreakEdge();
//                    newNodetwo->SetothersideID(secondLargest->GetId());
//                    Qtreeobj->GetNode(CriticalPath.back()->GetothersideID())->AddChild(newNodetwo);
//                    Qtreeobj->addNode(newNodetwo);
//                    temp=temp-newNodetwo->GetMSW();
//                    Qtreeobj->GetNode(CriticalPath.back()->GetothersideID())->SetMSW(temp);
//
//                    //cout<<"Node "<<newNodeone->GetId()<<", msweight "<<newNodeone->GetMSW()<<", "<<"Node "<<newNodetwo->GetId()<<", msweight "<<newNodetwo->GetMSW()<<endl;
//
//                    idleProcessors=idleProcessors-2;
//                    DeadLock=false;
//                    break;
//                }else{
//                    //cout<<"list L pop node "<<listL.front()->GetId()<<endl;
//                    listL.pop_front();
//                }
//            }else{//not on the last subtree
//                largestNode=listL.front();
//                largestNode->BreakEdge();
//                largestNode->SetothersideID(Qtreeobj->GetNodes()->size()+1);
//                currentNode=largestNode->GetParent();
//                while (!currentNode->IsBorken()) {
//                    currentNode=currentNode->GetParent();
//                }
//                Cnode* newNode = new Cnode(currentNode->GetothersideID(), 0, largestNode->GetEW(), largestNode->GetMSCost(false, false));//bug here
//                //cout<<"Break edge "<<listL.front()->GetId()<<" of subtree "<<currentNode->GetothersideID()<<", ";
//                newNode->SetId(Qtreeobj->GetNodes()->size()+1);
//                //cout<<"create Q node "<<Qtreeobj->GetNodes()->size()+1<<endl;
//                newNode->GetChildren()->clear();
//                newNode->SetParent(Qtreeobj->GetNode(currentNode->GetothersideID()));
//                newNode->BreakEdge();
//                newNode->SetothersideID(largestNode->GetId());
//                Qtreeobj->GetNode(currentNode->GetothersideID())->AddChild(newNode);
//                Qtreeobj->addNode(newNode);
//                temp=Qtreeobj->GetNode(currentNode->GetothersideID())->GetMSW();
//                Qtreeobj->GetNode(currentNode->GetothersideID())->SetMSW(temp-newNode->GetMSW());
//
//                //cout<<"Node "<<newNode->GetId()<<", msweight "<<newNode->GetMSW()<<endl;
//
//                idleProcessors--;
//                DeadLock=false;
//                break;
//            }
//        }
//
//        if (DeadLock==true) {
//            break;
//        }
//
////        cout<<"---------------------------------------------"<<endl;
//
////        Qroot->GetMSCost(true,true);//update makespan
////        for (int i=1; i<=Qtreeobj->GetNodes()->size(); ++i) {
////            cout<<"node "<<i<<", MScost "<<Qtreeobj->GetNode(i)->GetMSCost(true, false)<<", SequentialPart "<<Qtreeobj->GetNode(i)->GetSequentialPart()<<", ParallelPart "<<Qtreeobj->GetNode(i)->GetParallelPart()<<endl;
////            for (vector<Cnode*>::iterator iter=Qtreeobj->GetNode(i)->GetChildren()->begin(); iter!=Qtreeobj->GetNode(i)->GetChildren()->end(); ++iter) {
////                cout<<"   child "<<(*iter)->GetId()<<", MScost "<<(*iter)->GetMSCost(true, false)<<endl;
////            }
////        }
//
//    }
//
//    temp=tree->GetRoot()->GetMSCost(true, true);
//    delete Qtreeobj;
//
//    return temp;
//}

bool EstimateDecrase(int idleP, Ctree *tree, vector<Cnode *> *criticalPath, bool *lastsubtree, Cnode **node_i, Cnode **node_j)
{
    //cout<<"   --------------estimate decrease in makespan-----------------"<<endl;
    *lastsubtree = false;
    bool MSdecreased = false;
    vector<Cnode *> *children;
    vector<double> decreaseSequence;
    double temp, decrease = -1;
    vector<Cnode *> tempQue;
    Cnode *lastSubtreeRoot = tree->GetNode(criticalPath->back()->GetothersideID());

    //cout<<"   Last subtree root "<<lastSubtreeRoot->GetId()<<endl;
    //nodes on the last subtree of critical path
    //
    if (idleP > 1)
    { //has at least 2 idle processor
        tempQue.push_back(lastSubtreeRoot);
        vector<Cnode *>::iterator largestNode, secondLargest;
        //cout<<"   work on the last subtree "<<criticalPath->back()->GetothersideID()<<endl;
        /**
         * 对关键路径上最后一棵子树进行切割：如果子树的孩子数目大于2，则选择最大子树
         */
        while (!tempQue.empty())
        {
            children = tempQue.back()->GetChildren();
            tempQue.pop_back();
            if (children->size() > 1)
            {                                                                        //has at least 2 children
                GetTwoLargestElementTypethree(children, largestNode, secondLargest); //node_i is the largest, in terms of W
                //temp: the decrease of MS(t)
                temp = min((*largestNode)->GetSequentialPart() - (*secondLargest)->GetEW() / BANDWIDTH, (*secondLargest)->GetSequentialPart() - (*largestNode)->GetEW() / BANDWIDTH);
                if (temp > decrease)
                {
                    decrease = temp;
                    *node_i = *largestNode;
                    *node_j = *secondLargest;
                }

                //for循环结束后  tempQue保存了最后一个子树的所有孩子结点
                for (vector<Cnode *>::iterator it = children->begin(); it != children->end(); ++it)
                {
                    tempQue.push_back(*it);
                    if (it != largestNode)
                    {
                        //？？？是不是应该计算的是largestNode和it的并行的时间，，而不是largestNode和secondLargest的并行时间？？？
                        //如果是，上面的if 语句是否多余？？？
                        temp = min((*largestNode)->GetSequentialPart() - (*secondLargest)->GetEW() / BANDWIDTH, (*secondLargest)->GetSequentialPart() - (*largestNode)->GetEW() / BANDWIDTH);
                        if (temp > decrease)
                        {
                            decrease = temp;
                            *node_i = *largestNode;
                            *node_j = *it;
                        }
                    }
                }
            }
            else
            {
                for (vector<Cnode *>::iterator it = children->begin(); it != children->end(); ++it)
                {
                    tempQue.push_back(*it);
                }
            }
        }
    }

    //    if (decrease>=0) {
    //        cout<<"   if cut edge "<<(*node_i)->GetId()<<" and "<<(*node_j)->GetId()<<" of last subtree, decrease is "<<decrease<<endl;
    //    }else{
    //        cout<<"   fail on the last subtree."<<endl;
    //    }
    //
    //
    if (criticalPath->size() == 1)
    { //only one subtree on the critical path
        if (decrease >= 0)
        {
            *lastsubtree = true;
            MSdecreased = true;
        }
        return MSdecreased;
    }

    //nodes on other subtrees
    double decrease_othersubtrees = -1;
    Cnode *output_node;
    Cnode *subtreeRoot;
    Cnode *currentNode; //current node is on the path composed of critial path nodes
    Cnode *nodeOnPath;
    Cnode *SubtreeT = criticalPath->back();
    double MS_t, W_t;

    //cout<<"   working on subtree ";
    do
    {
        //currentNode子树的根结点
        currentNode = tree->GetNode(SubtreeT->GetothersideID());
        //SubtreeT商树中子树的父结点
        SubtreeT = SubtreeT->GetParent();
        //cout<<"   "<<SubtreeT->GetothersideID()<<"{ "<<endl;
        //subtreeRoot  子树的父结点的根结点
        subtreeRoot = tree->GetNode(SubtreeT->GetothersideID());
        MS_t = SubtreeT->GetMSCost(true, false);
        W_t = SubtreeT->GetMSW();

        do
        {
            //nodeOnPath关键路径上子树的根结点
            nodeOnPath = currentNode;

            currentNode = currentNode->GetParent();
            tempQue.push_back(currentNode);
            while (!tempQue.empty())
            {
                children = tempQue.back()->GetChildren();
                tempQue.pop_back();
                for (vector<Cnode *>::iterator it = children->begin(); it != children->end(); ++it)
                {
                    if ((*it)->GetId() != nodeOnPath->GetId() && (!(*it)->IsBorken()))
                    {
                        //cout<<"    "<<(*it)->GetId()<<" W_i "<<(*it)->GetSequentialPart()<<", MS(t) "<<MS_t<<", W_t "<<W_t<<", MS_tj "<<(*it)->GetParallelPart()<<endl;
                        tempQue.push_back((*it));
                        temp = min((*it)->GetSequentialPart(), MS_t - W_t - (*it)->GetEW() / BANDWIDTH - (*it)->GetParallelPart());
                        if (temp > decrease_othersubtrees)
                        {
                            decrease_othersubtrees = temp;
                            output_node = (*it);
                        }
                    }
                }
            }
        } while (!currentNode->IsBorken());
        //cout<<"   }"<<endl;
    } while (subtreeRoot->GetId() != tree->GetRootId());
    //cout<<endl;

    //    if (decrease_othersubtrees>=0) {
    //        cout<<"   if cut edge "<<output_node->GetId()<<", decrease is "<<decrease_othersubtrees<<endl;
    //    }else{
    //        cout<<"   cut an edge, fail"<<endl;
    //    }

    if (decrease_othersubtrees >= 0)
    {
        MSdecreased = true;
        if (decrease_othersubtrees < decrease)
        {
            *lastsubtree = true;
        }
        else
        {
            *node_i = output_node;
        }
    }
    else
    {
        if (decrease >= 0)
        {
            *lastsubtree = true;
            MSdecreased = true;
        }
    }

    return MSdecreased;
}

//Cnode* GetLargestSibling(Cnode* node){
//    unsigned int id=node->GetId();
//    Cnode* node_return;
//    vector<Cnode*>* Children=node->GetParent()->GetChildren();
//    double temp = 0;
//    double cost_sibling;
//    for (vector<Cnode*>::iterator iter=Children->begin(); iter!=Children->end(); ++iter) {
//        cost_sibling = (*iter)->GetMSCost(false,false);
//        if ((cost_sibling>=temp)&&((*iter)->GetId()!=id)) {
//            node_return = (*iter);
//            temp = cost_sibling;
//        }
//    }
//
//    return node_return;
//}

//如果处理器的个数>子树的个数
double SplitAgain(Ctree *tree, unsigned int processor_number, unsigned int num_subtrees)
{
    double MS_now;
    Cnode *root = tree->GetRoot();
    Ctree *Qtreeobj = BuildQtree(tree);

    vector<Cnode *> CriticalPath; //Q nodes on Critical Path

    Cnode *Qroot = Qtreeobj->GetRoot();
    Cnode *largestNode;
    Cnode *node_i;
    Cnode *node_j;
    Cnode *parent;
    double temp;
    vector<Cnode *> *Children;
    bool MSReduced, onLastSubtree;
    vector<Cnode *> tempVector;

    //idleProcessors  多出的处理器的个数
    int idleProcessors = processor_number - num_subtrees;

    while (idleProcessors > 0)
    {
        //cout<<"******** root id "<<tree->GetRootId()<<" ********"<<endl;
        CriticalPath.clear();
        CriticalPath.push_back(Qroot);
        MS_now = root->GetMSCost(true, true); //update makespan
        Qroot->GetMSCost(true, true);         //update critical path
        largestNode = Qroot;
        Children = Qroot->GetChildren();
        //cout<<"critical path (subtres' roots){1 ";
        while (!Children->empty())
        { //initialize critical path
            temp = largestNode->GetParallelPart();
            for (vector<Cnode *>::iterator iter = Children->begin(); iter != Children->end(); ++iter)
            {
                if ((*iter)->GetMSCost(true, false) == temp)
                {
                    largestNode = (*iter);
                    break;
                }
            }
            //cout<<largestNode->GetothersideID()<<" ";
            CriticalPath.push_back(largestNode);
            Children = largestNode->GetChildren();
        }
        //cout<<"}"<<endl;

        //cout<<"Idle processor now: "<<idleProcessors<<endl;
        //检测切割是否会减少makespan
        MSReduced = EstimateDecrase(idleProcessors, tree, &CriticalPath, &onLastSubtree, &node_i, &node_j);

        if (MSReduced == true)
        {
            if (onLastSubtree == false)
            { //不是最后一个结点
                //cout<<"cut edge "<<node_i->GetId()<<endl;
                node_i->BreakEdge(); //C<-C\cup C_k
                idleProcessors--;

                //对新切割出来的结点初始化
                node_i->SetothersideID(Qtreeobj->GetNodes()->size() + 1); //设置商树上的编号
                parent = node_i->GetParent();
                while (!parent->IsBorken())
                {
                    parent = parent->GetParent();
                }
                Cnode *Qparent = Qtreeobj->GetNode(parent->GetothersideID());
                Cnode *Qchild;
                Cnode *newNode = new Cnode(parent->GetothersideID(), 0, node_i->GetEW(), node_i->GetSequentialPart());
                newNode->SetId(Qtreeobj->GetNodes()->size() + 1);
                newNode->SetParent(Qparent);
                newNode->BreakEdge();
                newNode->SetothersideID(node_i->GetId());
                Qparent->AddChild(newNode);
                Qtreeobj->addNode(newNode);
                temp = Qparent->GetMSW();
                Qparent->SetMSW(temp - newNode->GetMSW());
                //cout<<"create new Q node "<<newNode->GetId()<<", msw "<<newNode->GetMSW()<<", its parent "<<Qparent->GetId()<<", msw "<<Qparent->GetMSW()<<endl;

                //设置新商树结点的  孩子
                newNode->GetChildren()->clear();
                if (node_i->GetParallelPart() > 0)
                {
                    //cout<<"went to here1."<<endl;
                    tempVector.push_back(node_i);
                    while (!tempVector.empty())
                    {
                        Children = tempVector.back()->GetChildren();
                        tempVector.pop_back();
                        for (vector<Cnode *>::iterator iter = Children->begin(); iter != Children->end(); ++iter)
                        {
                            if ((*iter)->IsBorken())
                            {
                                //如果边被切割，则将该结点对应商树上的结点为新结点的孩子
                                //cout<<"went to here2."<<endl;
                                Qchild = Qtreeobj->GetNode((*iter)->GetothersideID());
                                newNode->AddChild(Qchild);
                                Qchild->SetParent(newNode);
                                Qchild->SetParentId(newNode->GetId());
                                Qparent->RemoveChild((*iter)->GetothersideID());
                            }
                            else
                            {
                                tempVector.push_back((*iter));
                            }
                        }
                    }
                }
            }
            else
            { //商树上最后 一个结点，切割成三部分，切割该结点和其最大兄弟结点
                //cout<<"cut edge "<<node_i->GetId()<<" and edge "<<node_j->GetId()<<endl;
                node_i->BreakEdge(); //C<-C\cup C_k
                node_j->BreakEdge(); //C<-C\cup C_k
                idleProcessors = idleProcessors - 2;

                //设置结点i,j在商树上的结点的编号
                node_i->SetothersideID(Qtreeobj->GetNodes()->size() + 1);
                node_j->SetothersideID(Qtreeobj->GetNodes()->size() + 2);

                //初始化商树上的第一个新结点
                Cnode *newNodeone = new Cnode(CriticalPath.back()->GetId(), 0, node_i->GetEW(), node_i->GetSequentialPart());
                newNodeone->SetId(Qtreeobj->GetNodes()->size() + 1); //设置ID
                newNodeone->GetChildren()->clear();
                newNodeone->SetParent(CriticalPath.back()); //设置父结点
                newNodeone->BreakEdge();
                newNodeone->SetothersideID(node_i->GetId()); //设置子树根结点
                CriticalPath.back()->AddChild(newNodeone);   //
                Qtreeobj->addNode(newNodeone);
                temp = CriticalPath.back()->GetMSW();
                temp = temp - newNodeone->GetMSW();

                Cnode *newNodetwo = new Cnode(CriticalPath.back()->GetId(), 0, node_j->GetEW(), node_j->GetSequentialPart());
                newNodetwo->SetId(Qtreeobj->GetNodes()->size() + 1);
                newNodetwo->GetChildren()->clear();
                newNodetwo->SetParent(CriticalPath.back());
                newNodetwo->BreakEdge();
                newNodetwo->SetothersideID(node_j->GetId());
                CriticalPath.back()->AddChild(newNodetwo);
                Qtreeobj->addNode(newNodetwo);
                temp = temp - newNodetwo->GetMSW();
                CriticalPath.back()->SetMSW(temp); //更新原商树中最后一个子树的时间MS
                //cout<<"create new Q node "<<newNodetwo->GetId()<<", msw "<<newNodetwo->GetMSW()<<" and new node "<<newNodeone->GetId()<<", msw "<<newNodeone->GetMSW()<<", their parent "<<CriticalPath.back()->GetId()<<", msw "<<CriticalPath.back()->GetMSW()<<endl;
            }
        }
        else
        {
            break;
        }
    }

    delete Qtreeobj;

    MS_now = root->GetMSCost(true, true); //更新MS
    return MS_now;
}

//num_para_subtrees 子树划分的次数
void Immediately(Ctree *tree, unsigned long N, double *nwghts, double *ewghts, int *chstart, int *children, int *schedule, double m_availble, unsigned int &num_para_subtrees, vector<unsigned int> *brokenEdges)
{
    double memory_occupied = ewghts[schedule[N - 1]];
    list<unsigned int> allNodes;
    list<unsigned int> queue;
    unsigned long subtree_size;
    list<unsigned int>::iterator iter;
    unsigned int com_freq;
    schedule_t *schedule_f = new schedule_t();
    list<int>::iterator ite_sche;
    double maxoutD, memory_required, node_cost, data_to_unload;
    uint64_t count;
    int rootid, cur_task_id;
    unsigned int child_start, child_end;
    vector<unsigned int> subtreeBrokenEdges;

    //cout<<"current task:";
    for (unsigned long rank = N - 1; rank >= 1; rank--)
    {
        cur_task_id = schedule[rank];
        if (cur_task_id != 0)
        { //=0 means this node has already been moved to another processor
            //cout<<" "<<cur_task_id;
            // node_cost 计算MemReq(cur_task_id)
            node_cost = ewghts[cur_task_id] + nwghts[cur_task_id];
            for (int j = chstart[cur_task_id]; j < chstart[cur_task_id + 1]; j++)
            {
                node_cost += ewghts[children[j]];
            }

            data_to_unload = memory_occupied + node_cost - ewghts[cur_task_id] - m_availble;
            //内存不足，切掉相应的边
            if (data_to_unload > 0)
            { // schedule the subtree that is rooted at this node onto another processor
                //cout<<"(break)"<<endl;
                tree->GetNode(cur_task_id)->BreakEdge(); // set it cut, used for building a quotient tree later
                brokenEdges->push_back(tree->GetNode(cur_task_id)->GetothersideID());
                ++num_para_subtrees;
                allNodes.clear();
                queue.clear();
                queue.push_back(cur_task_id);

                //将以cur_task_id为根的子树上的所有结点入队
                do
                {
                    child_start = *(chstart + queue.front());
                    child_end = *(chstart + queue.front() + 1);
                    allNodes.push_back(queue.front());
                    queue.pop_front();
                    //将queue.front()的所有孩子入队
                    for (unsigned int i = child_start; i < child_end; ++i)
                    {
                        queue.push_back(*(children + i));
                    }
                } while (!queue.empty());

                subtree_size = allNodes.size();
                for (long i = rank - 1; i >= 0; i--)
                {
                    //将属于子树的结点的遍历顺序置为0  == 将该子树移至其它处理器
                    iter = find(allNodes.begin(), allNodes.end(), schedule[i]);
                    if (iter != allNodes.end())
                    {
                        schedule[i] = 0; //IO counter will pass 0;
                        allNodes.erase(iter);
                    }
                    if (allNodes.size() == 1)
                    {
                        break;
                    }
                }

                double *ewghtssub, *timewghtssub, *spacewghtssub;
                int *prntssub, *chstartsub, *chendsub, *childrensub;
                Ctree *subtree = BuildSubtree(tree, tree->GetNode(cur_task_id), subtree_size, &prntssub, &ewghtssub, &timewghtssub, &spacewghtssub, chstart, children);

                subtree_size = subtree->GetNodes()->size();
                //                        for (unsigned int index=1; index<=subtree_size; ++index) {
                //                            cout<<index<<" "<<subtree->GetNode(index)->GetParentId()<<" "<<subtree->GetNode(index)->GetNW()<<" "<<subtree->GetNode(index)->GetMSW()<<" "<<subtree->GetNode(index)->GetEW()<<endl;
                //                        }

                int *schedule_copy = new int[subtree_size + 1];
                maxoutD = MaxOutDegree(subtree, true);
                schedule_f->clear();
                count = 0;
                //获得新的遍历顺序
                MinMem(subtree, maxoutD, memory_required, *schedule_f, true, count);
                ite_sche = schedule_f->begin();
                for (unsigned int i = subtree_size; i >= 1; --i)
                {
                    schedule_copy[i] = *ite_sche;
                    advance(ite_sche, 1);
                }
                schedule_copy[0] = subtree_size + 1;
                po_construct(subtree_size, prntssub, &chstartsub, &chendsub, &childrensub, &rootid);

                //继续对子树进行遍历
                if (memory_required > m_availble)
                {
                    Immediately(subtree, subtree_size, spacewghtssub, ewghtssub, chstartsub, childrensub, schedule_copy, m_availble, com_freq, &subtreeBrokenEdges);

                    for (vector<unsigned int>::iterator iter = subtreeBrokenEdges.begin(); iter != subtreeBrokenEdges.end(); ++iter)
                    {
                        brokenEdges->push_back(tree->GetNode(*iter)->GetothersideID());
                    }
                }

                delete[] ewghtssub;
                delete[] timewghtssub;
                delete[] spacewghtssub;
                delete[] prntssub;
                delete[] chstartsub;
                delete[] chendsub;
                delete[] childrensub;
                delete[] schedule_copy;
                delete subtree;

                memory_occupied -= ewghts[cur_task_id];
                memory_occupied = max(0.0, memory_occupied);
            }
            else
            { //memory is enough for executing this node
                //？？？-2*ewghts[cur_task_id]
                memory_occupied += node_cost - 2 * ewghts[cur_task_id] - nwghts[cur_task_id];
                memory_occupied = max(0.0, memory_occupied);
            }
        }
    }
    delete schedule_f;
    //cout<<endl;
}

void MemoryCheck(Ctree *tree, int *chstart, int *children, double const memory_size, io_method_t method)
{ //chstart, children are not modified
    vector<Cnode *> subtreeRoots;
    Cnode *currentnode;
    Cnode *subtreeRoot;
    int rootid;
    tree->GetRoot()->BreakEdge();

    //cout<<"Subtrees' roots: ";
    unsigned long treeSize = tree->GetNodes()->size();
    //得到各子树的根结点
    for (unsigned int i = treeSize; i >= 1; --i)
    {
        currentnode = tree->GetNode(i);
        if (currentnode->IsBorken())
        {
            //cout<<i<<" ";
            subtreeRoots.push_back(currentnode);
        }
    }
    //cout<<endl;

    double maxoutD, memory_required;
    schedule_t *schedule_f = new schedule_t(); // 遍历顺序
    uint64_t count;
    unsigned int com_freq; //切割子树的次数
    unsigned long subtreeSize;
    list<int>::iterator ite_sche;
    vector<unsigned int> BrokenEdgesID;
    double IO_volume;
    while (!subtreeRoots.empty())
    {
        //每颗子树分别处理
        subtreeRoot = subtreeRoots.back();
        subtreeRoots.pop_back();

        double *ewghts, *timewghts, *spacewghts;
        int *prnts;
        Ctree *subtree = BuildSubtree(tree, subtreeRoot, treeSize, &prnts, &ewghts, &timewghts, &spacewghts, chstart, children);

        subtreeSize = subtree->GetNodes()->size();
        int *schedule_copy = new int[subtreeSize + 1];
        maxoutD = MaxOutDegree(subtree, true);
        schedule_f->clear();
        count = 0;
        //获取子树遍历序列 schedule_f，内存需求memory_required
        MinMem(subtree, maxoutD, memory_required, *schedule_f, true, count);

        int *chstartsub, *chendsub, *childrensub;
        po_construct(subtreeSize, prnts, &chstartsub, &chendsub, &childrensub, &rootid);
        // cout<<"rootid :"<<rootid<<endl;
        // cout<<"Subtree "<<subtreeRoot->GetId()<<" needs memory "<<memory_required;

        if (memory_required > memory_size)
        { //内存不足时，需要切割子树
            // cout<<", larger than what is available: "<<memory_size<<endl;

            ite_sche = schedule_f->begin();
            //逆序
            for (unsigned int i = subtreeSize; i >= 1; --i)
            {
                schedule_copy[i] = *ite_sche;
                // cout<<" schedule :"  <<schedule_copy[i];
                advance(ite_sche, 1);
            }
            schedule_copy[0] = subtreeSize + 1;

            switch (method)
            {
            case FIRST_FIT:
                IO_volume = IOCounter(subtree, subtreeSize + 1, spacewghts, ewghts, chstartsub, childrensub, schedule_copy, memory_size, false, true, com_freq, &BrokenEdgesID, FIRST_FIT);
                break;
            case LARGEST_FIT:
                IO_volume = IOCounter(subtree, subtreeSize + 1, spacewghts, ewghts, chstartsub, childrensub, schedule_copy, memory_size, false, true, com_freq, &BrokenEdgesID, LARGEST_FIT);
                break;
            case IMMEDIATELY:
                Immediately(subtree, subtreeSize + 1, spacewghts, ewghts, chstartsub, childrensub, schedule_copy, memory_size, com_freq, &BrokenEdgesID);
                break;

            default:
                break;
            }
        }
        // else
        // {
        //     cout<<"memory  available "<<endl;
        // }

        //cout<<endl;

        delete[] ewghts;
        delete[] timewghts;
        delete[] spacewghts;
        delete[] prnts;
        delete[] schedule_copy;
        delete[] chstartsub;
        delete[] chendsub;
        delete[] childrensub;
        delete subtree;
    }
    delete schedule_f;

    for (vector<unsigned int>::iterator iter = BrokenEdgesID.begin(); iter != BrokenEdgesID.end(); ++iter)
    {
        tree->GetNode(*iter)->BreakEdge();
    }
}

bool MemoryCheck_hetro(Ctree *tree, Cnode *subtreeRoot, int *chstart, int *children, double const memory_size, io_method_t method, double memory_required, vector<unsigned int> *BrokenEdgesID)
{                          //chstart, children are not modified
    unsigned int com_freq; //切割子树的次数

    vector<Cnode *> subtreeRoots;
    Cnode *currentnode;
    int rootid;
    double maxoutD;
    double IO_volume;
    unsigned long subtreeSize;
    list<int>::iterator ite_sche;
    double *ewghts, *timewghts, *spacewghts;
    int *prnts;
    uint64_t count;
    schedule_t *schedule_f = new schedule_t(); // 遍历顺序
    unsigned long treeSize = tree->GetNodes()->size();
    tree->GetRoot()->BreakEdge();

    Ctree *subtree = BuildSubtree(tree, subtreeRoot, treeSize, &prnts, &ewghts, &timewghts, &spacewghts, chstart, children);

    subtreeSize = subtree->GetNodes()->size();
    int *schedule_copy = new int[subtreeSize + 1];
    maxoutD = MaxOutDegree(subtree, true);
    schedule_f->clear();
    count = 0;
    //获取子树遍历序列 schedule_f，内存需求memory_required
    MinMem(subtree, maxoutD, memory_required, *schedule_f, true, count);

    int *chstartsub, *chendsub, *childrensub;
    po_construct(subtreeSize, prnts, &chstartsub, &chendsub, &childrensub, &rootid);
    // cout<<"rootid :"<<rootid<<endl;
    // cout<<"Subtree "<<subtreeRoot->GetId()<<" needs memory "<<memory_required;

    if (memory_required > memory_size)
    { //内存不足时，需要切割子树
        // cout<<", larger than what is available: "<<memory_size<<endl;

        ite_sche = schedule_f->begin();
        //逆序
        for (unsigned int i = subtreeSize; i >= 1; --i)
        {
            schedule_copy[i] = *ite_sche;
            // cout<<" schedule :"  <<schedule_copy[i];
            advance(ite_sche, 1);
        }
        schedule_copy[0] = subtreeSize + 1;

        switch (method)
        {
        case FIRST_FIT:
            IO_volume = IOCounter(subtree, subtreeSize + 1, spacewghts, ewghts, chstartsub, childrensub, schedule_copy, memory_size, false, true, com_freq, BrokenEdgesID, FIRST_FIT);
            break;
        case LARGEST_FIT:
            IO_volume = IOCounter(subtree, subtreeSize + 1, spacewghts, ewghts, chstartsub, childrensub, schedule_copy, memory_size, false, true, com_freq, BrokenEdgesID, LARGEST_FIT);
            break;
        case IMMEDIATELY:
            Immediately(subtree, subtreeSize + 1, spacewghts, ewghts, chstartsub, childrensub, schedule_copy, memory_size, com_freq, BrokenEdgesID);
            break;

        default:
            break;
        }
    }
    else
    {
        return true;
    }

    return false;
}

//计算子树的数目
unsigned int HowmanySubtrees(const Ctree *tree, bool quiet)
{
    unsigned int number_subtrees = 0;
    tree->GetRoot()->BreakEdge();
    const vector<Cnode *> *Nodes = tree->GetNodes();
    if (quiet == false)
    {
        cout << "Broken Edges { ";
    }
    for (auto it = Nodes->begin(); it != Nodes->end(); ++it)
    {
        if ((*it)->IsBorken())
        {
            number_subtrees++;
            if (quiet == false)
            {
                cout << (*it)->GetId() << " ";
            }
        }
    }
    if (quiet == false)
    {
        cout << "}" << endl;
    }
    return number_subtrees;
}

void SetBandwidth(double CCR, unsigned long tree_size, double *ewghts, double *timewghts)
{
    double sum_edges = 0;
    double sum_weights = 0;
    for (unsigned int i = 1; i <= tree_size; ++i)
    {
        sum_edges = sum_edges + ewghts[i];
        sum_weights = sum_weights + timewghts[i];
    }
    BANDWIDTH = sum_edges / (sum_weights * CCR); //  CCR = (sum_edges/BANDWIDTH) / sum_weights       ------the average communication to computation ratio
}

double Sequence(Cnode *root)
{
    return root->GetMSCost();
}

//初始化时破除所有边，再合并
void brokenAllEdge(Ctree *treeobj, int *prnts, unsigned int processor_number)
{
    if (processor_number <= 0)
    {
        return;
    }
    double maxoutd, minMem, memorySize;
    int N = treeobj->GetNodes()->size();
    int *chstart, *chend, *children, root = 1;
    // int *prnts = new int[N + 1];
    // for (int i = 1; i <= treeobj->getTreeSize(); i++)
    // {
    //     prnts[i] = treeobj->GetNode(i)->GetParentId();
    // }

    //破边
    const vector<Cnode *> *allnodes = treeobj->GetNodes();
    for (vector<Cnode *>::const_iterator i = allnodes->begin(); i != allnodes->end(); i++)
    {
        (*i)->BreakEdge();
    }
    int num_subtrees = treeobj->GetNodes()->size();

    //合并
    maxoutd = MaxOutDegree(treeobj, true); //树中 所有结点中所需的最大内存
    schedule_t *schedule_f = new schedule_t();
    int count = 0;
    MinMem(treeobj, maxoutd, minMem, *schedule_f, true, count);
    po_construct(N, prnts, &chstart, &chend, &children, &root);
    memorySize = minMem; //minMem为执行树所需最小内存,当内存大小等于minMem时,相当于没有内存限制

    MergeV2(treeobj, num_subtrees, processor_number, memorySize, chstart, children, false);
}

void GetFileNames(string path, vector<string> &filenames, int quite)
{
    DIR *pDir;
    struct dirent *ptr;
    if (!(pDir = opendir(path.c_str())))
    {
        cout << "Folder doesn't Exist!" << endl;
        return;
    }
    while ((ptr = readdir(pDir)) != 0)
    {
        if (strcmp(ptr->d_name, ".") != 0 && strcmp(ptr->d_name, "..") != 0)
        {
            filenames.push_back(path + "/" + ptr->d_name);
        }
    }
    if ((!quite))
    {
        cout << "all file:" << endl;
        for (vector<string>::iterator i = filenames.begin(); i != filenames.end(); i++)
        {
            cout << (*i) << endl;
        }
        cout << "============================================================================" << endl;
    }

    closedir(pDir);
}

Ctree *PartitionTree(Ctree *treeobj, unsigned int num_processors, bool quiet, io_method_t method, int memorySize_flag)
{
    if (num_processors == 0)
    {
        return treeobj;
    }
    int N = treeobj->GetNodes()->size();

    list<Cnode *> parallelSubtrees; //the subtree rooted in i which i is the root of tree.
    unsigned long sequentialLen;    //the length of parallelSubtrees
    clock_t time;
    unsigned int number_subtrees; //the number of subtrees
    double makespan = treeobj->GetRoot()->GetMSCost(true, true);
    double msbefore = makespan;
    double maxoutd, minMem, memorySize;
    int *chstart, *chend, *children, root = 1;
    Ctree *Qtree;
    uint64_t count;
    double memo[2] = {0, 0};

    int *prnts = new int[N + 1];

    // cout<<"test 1  in PartitionTree!"<<endl;
    for (int i = 1; i <= N; i++)
    {
        prnts[i] = treeobj->GetNode(i)->GetParentId();
    }

    time = clock();
    // cout<<"test 3  in PartitionTree!"<<endl;
    SplitSubtrees(treeobj->GetRoot(), num_processors, false, parallelSubtrees, sequentialLen);
    // cout<<"test 4  in PartitionTree!"<<endl;

    // time = clock() - time;
    number_subtrees = HowmanySubtrees(treeobj, quiet);

    vector<Cnode *> cut_Node, restoreNode;
    // cout<<"test 2  in PartitionTree!"<<endl;

    while (1)
    {
        msbefore = makespan;
        treeobj->GetRoot()->GetMSCost(true, true);
        Qtree = BuildQtree(treeobj);
        // vector<Cnode*>criQNode=getCriticalQNode(treeobj);
        vector<Cnode *> criQNode = getCriticalQNodeByQtree(Qtree);
        cut_Node.clear();
        // cout << "THE CRITICAL NODE:";
        // for (vector<Cnode *>::iterator iter = criQNode.begin(); iter != criQNode.end(); iter++)
        // {
        //     cout << (*iter)->GetothersideID() << " " << (*iter)->GetMSW() << endl;
        // }
        // cout << endl;
        // cout << "the number of critical path: " << CountCriPathByQtree(treeobj) << endl;

        Cnode *Qtemp, *temp;
        int pathNumbefore = CountCriPathByQtree(Qtree); //切割前关键路径的数目
        //  int pathNumbefore = CountCriPath(treeobj);
        int i = pathNumbefore;
        while (i > 0)
        {
            while (!criQNode.empty())
            {
                Qtemp = criQNode.back();
                criQNode.pop_back();

                temp = treeobj->GetNode(Qtemp->GetothersideID());
                SplitSubtrees(temp, N, false, parallelSubtrees, sequentialLen);
                cut_Node.insert(cut_Node.end(), parallelSubtrees.begin(), parallelSubtrees.end());
            }
            i--;
            // i = CountCriPathByQtree(Qtree);
            if (i > 0)
            {
                treeobj->GetRoot()->GetMSCost(true, true);
                Qtree = BuildQtree(treeobj);
                criQNode = getCriticalQNodeByQtree(Qtree);
            }
        }

        number_subtrees = HowmanySubtrees(treeobj, quiet);
        maxoutd = MaxOutDegree(treeobj, true); //树中 所有结点中所需的最大内存
        schedule_t *schedule_f = new schedule_t();
        count = 0;
        MinMem(treeobj, maxoutd, minMem, *schedule_f, true, count);
        po_construct(N, prnts, &chstart, &chend, &children, &root);
        // memorySize = minMem; //minMem为执行树所需最小内存,当内存大小等于minMem时,相当于没有内存限制
        // if (memorySize_flag == 1)
        // {
        //     memorySize = maxoutd;
        // }
        // else if (memorySize_flag == 2)
        // {
        //     memorySize = (maxoutd + minMem) / 2;
        // }
        // else
        // {
        //     memorySize = minMem;
        // }
        // // memorySize = maxoutd;
        // MemoryCheck(treeobj, chstart, children, memorySize, method);
        //
        //子树的个数大于处理器的个数，合并
        if (number_subtrees > num_processors)
        {

            Merge(treeobj, number_subtrees, num_processors, memorySize, chstart, children, false);
            // MergeV2(treeobj, number_subtrees, num_processors, memorySize, chstart, children, true);
        }
        number_subtrees = HowmanySubtrees(treeobj, quiet);
        double msafter = treeobj->GetRoot()->GetMSCost(true, true);
        treeobj->GetRoot()->GetMSCost(true, true);
        // cout << "msbefore :" << msbefore << endl;
        // cout << "msafter:" << msafter << endl;
        if (memo[0] > msafter)
        {
            memo[0] = msafter;
            memo[1] = 1;
        }
        else if (memo[0] == msafter)
        {
            memo[1]++;
        }
        else if (memo[0] == 0)
        {
            memo[0] = msafter;
        }

        if (msbefore == msafter)
        {
            int pathNumAfter = CountCriPath(treeobj);
            // cout << "pathNumAfter:" << pathNumAfter << endl;
            // cout << "pathNumbefore:" << pathNumbefore << endl;

            if (pathNumAfter == pathNumbefore)
            {

                break;
            }
        }
        if (msbefore > msafter)
        {
            if (memo[0] == msafter && memo[1] >= 2)
            {
                break;
            }
        }
        // time = clock() - time;
        // number_subtrees = HowmanySubtrees(treeobj, quiet);
        makespan = treeobj->GetRoot()->GetMSCost(true, true);
        // cout << PNR << " " << CCR << " NA " << number_subtrees << " " << num_processors << " " << makespan << " "
        //      << stage1 << "  splitQtree NA  " << time << endl;
    }
    time = clock() - time;
    number_subtrees = HowmanySubtrees(treeobj, quiet);
    makespan = treeobj->GetRoot()->GetMSCost(true, true);

    delete[] prnts;

    return treeobj;
}

Ctree *PartitionTree_memory(Ctree *treeobj, unsigned int num_processors, bool quiet, int memory_flag, io_method_t method)
{
    if (num_processors == 0)
    {
        return treeobj;
    }
    int N = treeobj->GetNodes()->size();

    list<Cnode *> parallelSubtrees; //the subtree rooted in i which i is the root of tree.
    unsigned long sequentialLen;    //the length of parallelSubtrees
    clock_t time;
    unsigned int number_subtrees; //the number of subtrees
    double makespan = treeobj->GetRoot()->GetMSCost(true, true);
    double msbefore = makespan;
    double maxoutd, minMem, memorySize;
    int *chstart, *chend, *children, root = 1;
    Ctree *Qtree;
    uint64_t count;
    double memo[2] = {0, 0};

    int *prnts = new int[N + 1];
    //     double *ewghts=new double[N+1], *nwghts=new double[N+1], *mswghts=new double[N+1];
    // cout<<"test 1  in PartitionTree!"<<endl;
    for (int i = 1; i <= N; i++)
    {
        prnts[i] = treeobj->GetNode(i)->GetParentId();
        //     ewghts[i]=treeobj->GetNode(i)->GetEW();
        //     nwghts[i]=treeobj->GetNode(i)->GetNW();
        //     mswghts[i]=treeobj->GetNode(i)->GetMSW();
    }

    time = clock();
    // cout<<"test 3  in PartitionTree!"<<endl;
    SplitSubtrees(treeobj->GetRoot(), num_processors, false, parallelSubtrees, sequentialLen);
    // cout<<"test 4  in PartitionTree!"<<endl;
    // cout<<"ks="<<ks<<endl;
    // time = clock() - time;
    number_subtrees = HowmanySubtrees(treeobj, quiet);
    // makespan = treeobj->GetRoot()->GetMSCost(true, true);
    // cout << PNR << " " << CCR << " NA " << number_subtrees << " " << num_processors << " " << makespan << " "
    //      << "splitSubtree"
    //      << "  NA  NA  " << time << endl;
    // time = clock();
    vector<Cnode *> cut_Node, restoreNode;
    cout << "test 2  in PartitionTree!" << endl;

    while (1)
    {
        msbefore = makespan;
        treeobj->GetRoot()->GetMSCost(true, true);
        Qtree = BuildQtree(treeobj);
        // vector<Cnode*>criQNode=getCriticalQNode(treeobj);
        vector<Cnode *> criQNode = getCriticalQNodeByQtree(Qtree);
        cut_Node.clear();
        // cout << "THE CRITICAL NODE:";
        // for (vector<Cnode *>::iterator iter = criQNode.begin(); iter != criQNode.end(); iter++)
        // {
        //     cout << (*iter)->GetothersideID() << " " << (*iter)->GetMSW() << endl;
        // }
        // cout << endl;
        // cout << "the number of critical path: " << CountCriPathByQtree(treeobj) << endl;

        Cnode *Qtemp, *temp;
        int pathNumbefore = CountCriPathByQtree(Qtree); //切割前关键路径的数目
        //  int pathNumbefore = CountCriPath(treeobj);
        int i = pathNumbefore;
        while (i > 0)
        {
            while (!criQNode.empty())
            {
                Qtemp = criQNode.back();
                criQNode.pop_back();

                temp = treeobj->GetNode(Qtemp->GetothersideID());
                SplitSubtrees(temp, N, false, parallelSubtrees, sequentialLen);
                cut_Node.insert(cut_Node.end(), parallelSubtrees.begin(), parallelSubtrees.end());
            }
            Qtree = BuildQtree(treeobj);
            i = CountCriPathByQtree(Qtree); //切割后关键路径的数目
            if (i > 0)
            {
                treeobj->GetRoot()->GetMSCost(true, true);
                //Qtree = BuildQtree(treeobj);
                criQNode = getCriticalQNodeByQtree(Qtree);
            }
        }
        cout << "test 3  in PartitionTree!" << endl;
        number_subtrees = HowmanySubtrees(treeobj, quiet);
        // the second step: 内存检测
        maxoutd = MaxOutDegree(treeobj, true); //树中 所有结点中所需的最大内存
        schedule_t *schedule_f = new schedule_t();
        count = 0;
        MinMem(treeobj, maxoutd, minMem, *schedule_f, true, count);
        po_construct(N, prnts, &chstart, &chend, &children, &root);
        if (memory_flag == 1)
        {
            memorySize = maxoutd;
        }
        else if (memory_flag == 2)
        {
            memorySize = (maxoutd + minMem) / 2;
        }
        else
        {
            memorySize = minMem;
        }
        cout << "test 4  in PartitionTree!" << endl;
        MemoryCheck(treeobj, chstart, children, memorySize, method);
        cout << "test 5  in PartitionTree!" << endl;
        //子树的个数大于处理器的个数，合并
        if (number_subtrees > num_processors)
        {
            //memorySize = minMem; //minMem为执行树所需最小内存,当内存大小等于minMem时,相当于没有内存限制

            // Merge(treeobj, restoreNode, number_subtrees, num_processors, memorySize, chstart, children, false);
            MergeV2(treeobj, number_subtrees, num_processors, memorySize, chstart, children, true);
        }
        cout << "test 6  in PartitionTree!" << endl;
        number_subtrees = HowmanySubtrees(treeobj, quiet);
        double msafter = treeobj->GetRoot()->GetMSCost(true, true);
        treeobj->GetRoot()->GetMSCost(true, true);
        // cout << "msbefore :" << msbefore << endl;
        // cout << "msafter:" << msafter << endl;
        if (memo[0] > msafter)
        {
            memo[0] = msafter;
        }
        else if (memo[0] == msafter)
        {
            memo[1]++;
        }
        else if (memo[0] == 0)
        {
            memo[0] = msafter;
        }

        if (msbefore == msafter)
        {
            int pathNumAfter = CountCriPath(treeobj);
            cout << "pathNumAfter:" << pathNumAfter << endl;
            cout << "pathNumbefore:" << pathNumbefore << endl;
            cout << endl;

            if (pathNumAfter == pathNumbefore)
            {

                break;
            }
        }
        if (msbefore > msafter)
        {
            if (memo[1] == msafter && memo[1] >= 2)
            {
                break;
            }
        }
        // time = clock() - time;
        // number_subtrees = HowmanySubtrees(treeobj, quiet);
        makespan = treeobj->GetRoot()->GetMSCost(true, true);
        // cout << PNR << " " << CCR << " NA " << number_subtrees << " " << num_processors << " " << makespan << " "
        //      << stage1 << "  splitQtree NA  " << time << endl;
    }
    time = clock() - time;
    number_subtrees = HowmanySubtrees(treeobj, quiet);
    makespan = treeobj->GetRoot()->GetMSCost(true, true);

    // cout << "MS:" << treeobj->GetRoot()->GetMSCost(true, true) << endl;
    // cout << PNR << " " << CCR << " NA " << number_subtrees << " " << num_processors << " " << makespan << " "
    //      << "splitSubtree"
    //      << "  splitQtree NA  " << time << endl;
    // cout << "hahah" << sizeof(Ctree) << endl;
    // cout << &treeobj << endl;
    // cout << "测试------------------" << endl;
    // cout << "测试------------------" << endl;
    delete[] prnts;

    return treeobj;
}

//直接通过商树获得关键路径上的结点
vector<Cnode *> getCriticalQNodeByQtree(Ctree *Qtree)
{

    vector<Cnode *> CriticalPath;
    Cnode *Qroot = Qtree->GetRoot();
    Cnode *largestNode, *tempNode;
    double temp;
    vector<Cnode *> *Children;

    list<Cnode *> parallelSubtrees;

    CriticalPath.clear();
    CriticalPath.push_back(Qroot);
    // MS_now=root->GetMSCost(true, true);//update makespan
    Qroot->GetMSCost(true, true); //update critical path
    largestNode = Qroot;
    Children = Qroot->GetChildren();
    //cout<<"critical path (subtres' roots){1 ";

    //get the critical [path]
    while (!Children->empty())
    { //initialize critical path
        temp = largestNode->GetParallelPart();
        for (vector<Cnode *>::iterator iter = Children->begin(); iter != Children->end(); ++iter)
        {
            if ((*iter)->GetMSCost(true, false) == temp)
            {
                largestNode = (*iter);
                break;
            }
        }
        //cout<<largestNode->GetothersideID()<<" ";
        CriticalPath.push_back(largestNode);
        Children = largestNode->GetChildren();
    }
    return CriticalPath;
}

//直接使用商树计算关键路径的数目
int CountCriPathByQtree(Ctree *Qtree)
{

    Cnode *Qroot = Qtree->GetRoot();
    double temp;
    vector<Cnode *> *Children = Qroot->GetChildren();
    vector<Cnode *> criPath;
    if (Children->size() < 1)
    {
        return 1;
    }

    Qroot->GetMSCost(true, true); //update critical path
    temp = Qroot->GetParallelPart();
    // cout<<"ttemp:"<<temp<<endl;
    for (vector<Cnode *>::iterator iter = Children->begin(); iter != Children->end(); ++iter)
    {
        if ((*iter)->GetMSCost(true, false) == temp)
        {
            criPath.push_back((*iter));
        }
    }

    return criPath.size();
}

//计算关键路径的数目
int CountCriPath(Ctree *tree)
{

    Cnode *root = tree->GetRoot();
    Ctree *Qtreeobj = BuildQtree(tree);
    Cnode *Qroot = Qtreeobj->GetRoot();
    double temp;
    vector<Cnode *> *Children = Qroot->GetChildren();
    vector<Cnode *> criPath;
    if (Children->size() <= 0)
    {
        return 1;
    }

    Qroot->GetMSCost(true, true); //update critical path
    temp = Qroot->GetParallelPart();
    for (vector<Cnode *>::iterator iter = Children->begin(); iter != Children->end(); ++iter)
    {
        // cout<<(*iter)->GetothersideID()<<" ";
        if ((*iter)->GetMSCost(true, false) == temp)
        {
            criPath.push_back((*iter));
        }
    }
    return criPath.size();
}

// double *getRandom(double E, double S, int num)
// {
//     unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); //每次程序运行产生不同的随机数
//     std::default_random_engine gen(seed);
//     std::normal_distribution<double> n(E, S); //均值为0，标准差为1
//     double *number = new double[num];
//     double rnd;
//     for (int i = 0; i < num; i++)
//     {
//         while ((rnd = n(gen)) < 0)
//             ;
//         number[i] = rnd;
//         sleep(0.1);
//     }
//     return number;
// }

// void swap(double *prower, double *mems, int low, int high)
// {
//     // double temp1,temp2;
//     swap(prower[low], prower[high]);
//     swap(mems[low], mems[high]);
// }
// //快速排序（从大到小）
// void quickSort(double *prower, double *mems, int left, int right)
// {

//     int low, high, m;
//     double pivot, temp;
//     low = left;
//     high = right;

//     m = (low + high) / 2;

//     // if ((prower[low] / mems[low]) > prower[high] / mems[high])
//     // {
//     //     swap(prower, mems, low, high);
//     // }
//     // if ((prower[m] / mems[m]) > (prower[high] / mems[high]))
//     // {
//     //     swap(prower, mems, m, high);
//     // }
//     // if ((prower[m] / mems[m]) > (prower[low] / mems[low]))
//     // {
//     //     swap(prower, mems, m, low);
//     // }
//     //<>
//     pivot = prower[low] / mems[low];
//     while (low < high)
//     {
//         while (low < high && (prower[high] / mems[high]) <= pivot)
//             high--;
//         while (low < high && (prower[low] / mems[low]) >= pivot)
//             low++;
//         if (low < high)
//             swap(prower, mems, low, high);
//     }

//     if (left < low)
//     {
//         swap(prower, mems, left, low);
//         quickSort(prower, mems, left, low - 1);
//         quickSort(prower, mems, low + 1, right);
//     }
// }

Cnode *getVectorMax(vector<Cnode *> nodes)
{
    Cnode *max = *min_element(nodes.begin(), nodes.end(), cmp_noincreasing_msw);
    // Cnode *max = nodes.front();
    // cout<<"min id :"<<max->GetId()<<endl;
    for (vector<Cnode *>::iterator i = nodes.begin(); i != nodes.end(); i++)
    {
        //返回关键路径上没有分配超算的结点的MSW最大的结点
        if (max->GetMSW() <= (*i)->GetMSW() && !(*i)->IsAlloc())
        {
            max = (*i);
        }
    }
    return max;
}
//获得nodes中的最大的节点，不考虑是否被分配
Cnode *getVectorMax_noAlloc(vector<Cnode *> *nodes)
{
    Cnode *max = *min_element(nodes->begin(), nodes->end(), cmp_noincreasing_msw);
    for (vector<Cnode *>::iterator i = nodes->begin(); i != nodes->end(); i++)
    {
        //返回关键路径上没有分配超算的结点的MSW最大的结点
        if (max->GetMSW() <= (*i)->GetMSW() )
        {
            max = (*i);
        }
    }
    return max;
}
//判断是否关键上的结点都分配完成
bool IsAllCriNodeAlloc(vector<Cnode *> nodes)
{
    for (vector<Cnode *>::iterator i = nodes.begin(); i != nodes.end(); i++)
    {
        if (!(*i)->IsAlloc())
            return false;
    }
    return true;
}
//获取vector中的最大值，主要用于求商树点中关键路径的最大值
Cnode *getVectorMax_crital(vector<Cnode *> nodes)
{

    if (IsAllCriNodeAlloc(nodes))
        return NULL; //如果关键路径上的所有结点都分配完，则返回NULL
    Cnode *max = *min_element(nodes.begin(), nodes.end(), cmp_noincreasing_msw);
    for (vector<Cnode *>::iterator i = nodes.begin(); i != nodes.end(); i++)
    {
        //返回关键路径上没有分配超算的结点的MSW最大的结点
        if (max->GetMSW() <= (*i)->GetMSW() && !(*i)->IsAlloc())
        {
            max = (*i);
        }
    }
    return max;
}

double AllocPower(Ctree *treeobj, int *prnts, Procs *procs, int num_processors, bool quite)
{

    if (HowmanySubtrees(treeobj, true) > num_processors)
    {
        return __DBL_MAX__;
    }
    Ctree *Q = BuildQtree(treeobj);
    vector<Cnode *> Children = *(Q->GetNodes());

    Cnode *node;
    vector<Cnode *> criQnodes;
    Cnode *maxNode;

    int N = Q->GetNodes()->size();
    // int n_p = num_processors;
    int ori_N = treeobj->GetNodes()->size();
    double minMem, maxoutd = MaxOutDegree(treeobj, true); //树中 所有结点中所需的最大内存

    schedule_t *schedule_f = new schedule_t();
    uint64_t count;
    MinMem(treeobj, maxoutd, minMem, *schedule_f, true, count);

    int *chstart, *chend, *children, root = 1;
    po_construct(ori_N, prnts, &chstart, &chend, &children, &root);
    Cnode *currentNode;
    int p_index = 0;
    // for (vector<Cnode *>::const_iterator iter = Q->GetNodes()->begin(); iter != Q->GetNodes()->end(); iter++)
    // {
    //     cout << (*iter)->GetMSW() << endl;
    // }
    while (procs->IsRemain() && N > 0)
    {
        // cout<<"remain processors:"<<procs->getRemainProcs()<<endl;
        p_index = 0;
        while (procs->IsUsedForPro(p_index)) //查找空闲处理器
            p_index++;
        node = Q->GetRoot();
        double ms = Q->GetRoot()->GetMSCost(true, true);

        criQnodes = getCriticalQNodeByQtree(Q);
        maxNode = getVectorMax_crital(criQnodes);
        if (maxNode == NULL)
        {
            maxNode = getVectorMax(Children);
        }
        if (maxNode != NULL && !maxNode->IsAlloc())
        {

            Cnode *subtreeRoot = treeobj->GetNode(maxNode->GetothersideID());
            vector<unsigned int> brokenEdgeID;
            double memory_required;
            bool flag = MemoryCheck_hetro(treeobj, subtreeRoot, chstart, children, procs->getMem(p_index), LARGEST_FIT, memory_required, &brokenEdgeID);
            //先查找符合内存的处理器
            //。。。。。
            //找不到再切割
            if (flag == false)
            {
                p_index = 0;
                while (p_index < procs->getProcNum() && !procs->EnoughMem(p_index, memory_required))
                    p_index++;
                // cout << "memory constaint:" << p_index << endl;
                if (p_index >= procs->getProcNum()) //找不到再切割
                {
                    if (N + brokenEdgeID.size() <= procs->getProcNum()) //若切割后的子树个数大于处理器个数，则不切割
                    {
                        for (vector<unsigned int>::iterator iter = brokenEdgeID.begin(); iter != brokenEdgeID.end(); ++iter)
                        {
                            treeobj->GetNode(*iter)->BreakEdge();
                            currentNode = new Cnode();
                            Cnode *oriNode = treeobj->GetNode((*iter));
                            oriNode->SetothersideID(++N);
                            currentNode->SetothersideID((*iter));
                            currentNode->SetEW(oriNode->GetEW());
                            // currentNode.SetNW(nwghts[i]);
                            currentNode->SetMSW(oriNode->GetSequentialPart());
                            currentNode->SetId(N);
                            currentNode->SetLabel(N);
                            currentNode->SetParentId(maxNode->GetId());
                            vector<int> new_add;
                            for (vector<Cnode *>::const_iterator i = treeobj->GetNodes()->begin(); i != treeobj->GetNodes()->end(); i++) //修改新子树的节点的所属的商树ID
                            {
                                if ((*i)->GetothersideID() == maxNode->GetId() && (*i)->GetParent()->GetothersideID() == N)
                                {
                                    (*i)->SetothersideID(N);
                                    new_add.push_back((*i)->GetId());
                                }
                            }
                            for (vector<Cnode *>::const_iterator i = maxNode->GetChildren()->begin(); i != maxNode->GetChildren()->end(); i++) //修改新子树的子节点
                            {
                                int childID = (*i)->GetothersideID();
                                vector<int>::iterator k = find(new_add.begin(), new_add.end(), childID);
                                if (k != new_add.end())
                                {
                                    Q->GetNode(childID)->SetParentId(N);
                                    Q->GetNode(childID)->SetParent(currentNode);
                                    currentNode->AddChild(Q->GetNode(childID));
                                }
                            }

                            Q->addNode(currentNode);
                            Children.push_back(currentNode);
                        }
                    }
                    else
                    {
                        return __DBL_MAX__;
                    }
                    treeobj->GetRoot()->GetMSCost(true, true);
                }
                else
                {
                    maxNode->setPower(procs->getPower(p_index));
                    procs->setUsed(p_index);
                    N--;

                    // cout<<"original:"<<maxNode->GetMSW()<<"     "<<powers[p_index]<<endl;
                    maxNode->SetMSW(maxNode->GetMSW() / procs->getPower(p_index));
                }
            }
            else //内存满足要求
            {
                if (!quite)
                {
                    cout << "the Qnode is :" << maxNode->GetId() << endl;
                    cout << "the MSW before: " << maxNode->GetMSW() << endl;
                    cout << "Is Have sc: " << maxNode->IsAlloc() << endl;
                }

                maxNode->setPower(procs->getPower(p_index));
                procs->setUsed(p_index);
                N--;

                // cout<<"original:"<<maxNode->GetMSW()<<"     "<<powers[p_index]<<endl;
                maxNode->SetMSW(maxNode->GetMSW() / procs->getPower(p_index));
                // cout<<"after:"<<maxNode->GetMSW()<<endl;

                if (!quite)
                {
                    cout << "the power is :" << procs->getPower(p_index) << endl;
                    cout << "Is Have sc: " << maxNode->IsAlloc() << endl;
                    cout << "the MSW after:" << maxNode->GetMSW() << endl;
                    cout << "=====================================================" << endl;
                }
                //
                // Children.insert(Children.end(), maxNode->GetChildren()->begin(), maxNode->GetChildren()->end());
                // Q->GetRoot()->GetMSCost(true, true);
            }
        }
    }
    // Q->GetNodes()
    // for (vector<Cnode *>::const_iterator iter = Q->GetNodes()->begin(); iter != Q->GetNodes()->end(); iter++)
    // {
    //     cout << (*iter)->GetMSW() << endl;
    // }
    // cout<<"n_p:"<<p_index<<endl;
    HowmanySubtrees(treeobj, true);
    double lastms = Q->GetRoot()->GetMSCost(true, true);
    // cout << "the end: " << lastms << endl;
    return lastms;
}

double AllocPowerV2(Ctree *treeobj, int *prnts, Procs *procs, int num_processors, bool quite)
{
    if (HowmanySubtrees(treeobj, true) > num_processors)
    {
        return __DBL_MAX__;
    }
    // Cnode *maxNode;
    Cnode *node;
    Ctree *Q = BuildQtree(treeobj);
    vector<Cnode *> allNodes = *(Q->GetNodes());
    int N = Q->GetNodes()->size();
    int ori_N = treeobj->GetNodes()->size();
    double minMem, maxoutd = MaxOutDegree(treeobj, true); //树中 所有结点中所需的最大内存

    schedule_t *schedule_f = new schedule_t();
    uint64_t count;
    MinMem(treeobj, maxoutd, minMem, *schedule_f, true, count);

    int *chstart, *chend, *children, root = 1;
    po_construct(ori_N, prnts, &chstart, &chend, &children, &root);
    int p_index = 0;
    Cnode *currentNode;
    while (procs->IsRemain() && N > 0)
    {
        p_index = 0;
        while (procs->IsUsedForPro(p_index))
            p_index++;
        node = Q->GetRoot();
        while (node->IsAlloc() && (!node->IsLeaf())) //查找
        {
            // double maxMS = node->GetParallelPart();
            vector<Cnode *> *ItsChildren = node->GetChildren();
           node=getVectorMax_noAlloc(ItsChildren);
        }
        if (!node->IsAlloc())
        {
            Cnode *subtreeRoot = treeobj->GetNode(node->GetothersideID());
            vector<unsigned int> brokenEdgeID;
            double memory_required;
            bool flag = MemoryCheck_hetro(treeobj, subtreeRoot, chstart, children, procs->getMem(p_index), LARGEST_FIT, memory_required, &brokenEdgeID);
            if (flag == false)
            {
                p_index = 0;
                while (p_index < procs->getProcNum() && !procs->EnoughMem(p_index, memory_required))
                    p_index++;
                // cout << "memory constaint:" << p_index << endl;
                if (p_index >= procs->getProcNum()) //找不到再切割
                {
                    if (N + brokenEdgeID.size() <= procs->getProcNum()) //若切割后的子树个数大于处理器个数，则不切割
                    {
                        for (vector<unsigned int>::iterator iter = brokenEdgeID.begin(); iter != brokenEdgeID.end(); ++iter)
                        {
                            treeobj->GetNode(*iter)->BreakEdge();
                            currentNode = new Cnode();
                            Cnode *oriNode = treeobj->GetNode((*iter));
                            oriNode->SetothersideID(++N);
                            currentNode->SetothersideID((*iter));
                            currentNode->SetEW(oriNode->GetEW());
                            // currentNode.SetNW(nwghts[i]);
                            currentNode->SetMSW(oriNode->GetSequentialPart());
                            currentNode->SetId(N);
                            currentNode->SetLabel(N);
                            currentNode->SetParentId(node->GetId());
                            vector<int> new_add;
                            for (vector<Cnode *>::const_iterator i = treeobj->GetNodes()->begin(); i != treeobj->GetNodes()->end(); i++) //修改新子树的节点的所属的商树ID
                            {
                                if ((*i)->GetothersideID() == node->GetId() && (*i)->GetParent()->GetothersideID() == N)
                                {
                                    (*i)->SetothersideID(N);
                                    new_add.push_back((*i)->GetId());
                                }
                            }
                            for (vector<Cnode *>::const_iterator i = node->GetChildren()->begin(); i != node->GetChildren()->end(); i++) //修改新子树的子节点
                            {
                                int childID = (*i)->GetothersideID();
                                vector<int>::iterator k = find(new_add.begin(), new_add.end(), childID);
                                if (k != new_add.end())
                                {
                                    Q->GetNode(childID)->SetParentId(N);
                                    Q->GetNode(childID)->SetParent(currentNode);
                                    currentNode->AddChild(Q->GetNode(childID));
                                }
                            }

                            Q->addNode(currentNode);
                            allNodes.push_back(currentNode);
                        }
                    }
                    else
                    {
                        return __DBL_MAX__;
                    }
                    treeobj->GetRoot()->GetMSCost(true, true);
                }
                else
                {
                    node->setPower(procs->getPower(p_index));
                    procs->setUsed(p_index);
                    N--;

                    // cout<<"original:"<<maxNode->GetMSW()<<"     "<<powers[p_index]<<endl;
                    node->SetMSW(node->GetMSW() / procs->getPower(p_index));
                }
            }
            else //内存满足要求
            {
                if (!quite)
                {
                    cout << "the Qnode is :" << node->GetId() << endl;
                    cout << "the MSW before: " << node->GetMSW() << endl;
                    cout << "Is Have sc: " << node->IsAlloc() << endl;
                }

                node->setPower(procs->getPower(p_index));
                procs->setUsed(p_index);
                N--;

                // cout<<"original:"<<maxNode->GetMSW()<<"     "<<powers[p_index]<<endl;
                node->SetMSW(node->GetMSW() / procs->getPower(p_index));
                // cout<<"after:"<<maxNode->GetMSW()<<endl;

                if (!quite)
                {
                    cout << "the power is :" << procs->getPower(p_index) << endl;
                    cout << "Is Have sc: " << node->IsAlloc() << endl;
                    cout << "the MSW after:" << node->GetMSW() << endl;
                    cout << "=====================================================" << endl;
                }
                //
                // Children.insert(Children.end(), maxNode->GetChildren()->begin(), maxNode->GetChildren()->end());
                // Q->GetRoot()->GetMSCost(true, true);
            }
        }
    
    }
