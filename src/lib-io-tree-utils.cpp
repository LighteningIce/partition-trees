/*
 *  lib-io-tree-utils.cpp
 *  lib-io-tree
 *
 *  Created by defbond on 8/10/10.
 *  Copyright 2010 LIP/ENS-Lyon. All rights reserved.
 *
 */
#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <cmath>

#include "lib-io-tree-utils.h"
#include "lib-io-tree.h"
#include "lib-io-tree-minmem.h"
#include <sys/time.h>
#include <algorithm>
#include <random>
#include <chrono>
#include <unistd.h>
using namespace std;

double BANDWIDTH=1;

bool sort_sche(node_sche a, node_sche b){
    return (a.second>b.second);
}

bool sort_ew(node_ew a, node_ew b){
    return (a.second>b.second);
}

double u_wseconds(void) {
    struct timeval tp;
    
    gettimeofday(&tp, NULL);
    
    return (double) tp.tv_sec + (double) tp.tv_usec / 1000000.0;
};

void parse_tree(const char *filename,Ctree * tree) {
    ifstream OpenFile(filename);
    char begin;
    char cur_char;
    unsigned int line_index = 1;
    bool nodes_cnt_read = false;
    
    do
    {
        /*skip commentary lines*/
        begin = OpenFile.peek();
        if(OpenFile.good() ){
            
            if(begin=='%'){
                do{
                    cur_char = OpenFile.get();
                }while(cur_char != '\n' && OpenFile.good());
            }
            else{
                if(!nodes_cnt_read){
                    /* get the number of nodes and skip last trailing character*/
                    int nb_of_nodes;
                    OpenFile >> nb_of_nodes;
                    do{cur_char = OpenFile.get();}while(cur_char != '\n' && OpenFile.good());
                    
                    
                    //OpenFile.get(cur_char);
                    nodes_cnt_read = true;
                    /*allocate space for nodes*/
                    tree->AllocateNodes(nb_of_nodes);
                    
                }
                else{
                    /*parse actual nodes*/
                    unsigned int parent;
                    double ew, nw;
                    
                    OpenFile >> parent >> nw >> ew;
                    do{cur_char = OpenFile.get();}while(cur_char != '\n' && OpenFile.good());
                    
                    
                    
                    
                    tree->GetNode(line_index-1)->SetParentId(parent-1);
                    tree->GetNode(line_index-1)->SetEW(ew);
                    tree->GetNode(line_index-1)->SetNW(nw);
                    tree->GetNode(line_index-1)->SetId(line_index-1);
                    
                    if(parent > 0){
                        tree->GetNode(parent-1)->AddChild(tree->GetNode(line_index-1));
                    }
                    
                    if(parent==0){
                        tree->SetRootId(line_index-1);
                    }
                    line_index++;
                }
            }
            
        }
    }while(OpenFile.good() );
    
    OpenFile.close();
}

void parse_tree(const char *filename,int * N ,int **prnts,double **nwghts,double **ewghts, double **mswghts){
    
    ifstream OpenFile(filename);
    char begin;
    char cur_char;
    int line_index = 0;
    bool nodes_cnt_read = false;
    unsigned int nb_of_nodes=0;
    string line;
    
    do
    {
        /*skip commentary lines*/
        begin = OpenFile.peek();
        if(OpenFile.good() ){
            
            if(begin=='%'){
                do{
                    cur_char = OpenFile.get();
                }while(cur_char != '\n' && OpenFile.good());
            }
            else{
                if(!nodes_cnt_read){
                    /* get the number of nodes and skip last trailing character*/
                    while(getline(OpenFile,line)){++nb_of_nodes;}
                    OpenFile.clear();
                    OpenFile.seekg(0,ios::beg);
                    *N = nb_of_nodes;
                    
                    //OpenFile.get(cur_char);
                    nodes_cnt_read = true;
                    /*allocate space for nodes*/
                    *prnts = new int[nb_of_nodes+1];
                    *nwghts = new double[nb_of_nodes+1];
                    *ewghts = new double[nb_of_nodes+1];
                    *mswghts = new double[nb_of_nodes+1];
                }
                else{
                    /*parse actual nodes*/
                    unsigned int id;
                    unsigned int parent;
                    double ew, nw, msw;
                    
                    OpenFile >> id >> parent >> nw >> msw >> ew;
                    //OpenFile >> id >> parent >> nw >> ew >> msw;
                    do{cur_char = OpenFile.get();}while(cur_char != '\n' && OpenFile.good());
                    parent = nb_of_nodes - parent + 1;//root has the largest id in the txt file
                    (*prnts)[nb_of_nodes-line_index] = parent;
                    (*nwghts)[nb_of_nodes-line_index] = nw;
                    (*ewghts)[nb_of_nodes-line_index] = ew;
                    (*mswghts)[nb_of_nodes-line_index] = msw;
                    
                    line_index++;
                }
            }
            
        }
    }while(OpenFile.good() );
    (*prnts)[1]=0;
    
    OpenFile.close();
}
// from he.
void parse_tree_Seq(const char *filename, int *N, int **prnts, double **nwghts, double **ewghts, double **mswghts)
{
    ifstream OpenFile(filename);
    char begin;
    char cur_char;
    int line_index = 0;
    bool nodes_cnt_read = false;
    unsigned int nb_of_nodes = 0;
    string line;

    do
    {
        /*skip commentary lines*/
        begin = OpenFile.peek();
        if (OpenFile.good())
        {

            if (begin == '%')
            {
                do
                {
                    cur_char = OpenFile.get();
                } while (cur_char != '\n' && OpenFile.good());
            }
            else
            {
                if (!nodes_cnt_read)
                {
                    /* get the number of nodes and skip last trailing character*/
                    while (getline(OpenFile, line))
                    {
                        ++nb_of_nodes;
                    }
                    OpenFile.clear();
                    OpenFile.seekg(0, ios::beg);
                    *N = nb_of_nodes;

                    //OpenFile.get(cur_char);
                    nodes_cnt_read = true;
                    /*allocate space for nodes*/
                    *prnts = new int[nb_of_nodes + 1];
                    *nwghts = new double[nb_of_nodes + 1];
                    *ewghts = new double[nb_of_nodes + 1];
                    *mswghts = new double[nb_of_nodes + 1];
                }
                else
                {
                    /*parse actual nodes*/
                    unsigned int id;
                    unsigned int parent;
                    double ew, nw, msw;
                    line_index++;

                    OpenFile >> id >> parent >> nw >> msw >> ew;
                    
                    //OpenFile >> id >> parent >> nw >> ew >> msw;
                    do
                    {
                        cur_char = OpenFile.get();
                    } while (cur_char != '\n' && OpenFile.good());

                    //正序输入，根节点是1
                    (*prnts)[line_index] = parent;
                    (*nwghts)[line_index] = nw;
                    (*mswghts)[line_index] = msw;
                    (*ewghts)[line_index] = ew;
                    
                    // cout<<id<<" "<< (*prnts)[line_index]<<" "<<(*nwghts)[line_index] <<" " <<(*mswghts)[line_index]<<" "<<(*ewghts)[line_index]<<endl;
                    
                }
            }
        }
    } while (OpenFile.good());
    (*prnts)[1] = 0;

    OpenFile.close();
}

//void ConvertToLiu(const Ctree * tree_us, Ctree * tree_liu) {
//
//    tree_liu->AllocateNodes(2*tree_us->GetNodes()->size());
//
//
//    tree_liu->SetRootId(tree_us->GetRootId());
//
//
//
//    unsigned int last_node = tree_us->GetNodes()->size()+1;
//
//    for(unsigned int node_index=1;node_index< tree_us->GetNodes()->size()+1;node_index++){
//        double ew, nw ;
//
//        Cnode * us_node = tree_us->GetNode(node_index);
//        Cnode * liu_edgenode = tree_liu->GetNode(node_index);
//        Cnode * liu_peaknode = tree_liu->GetNode(last_node);
//
//        liu_edgenode->SetId(node_index);
//        liu_peaknode->SetId(last_node);
//
//        nw = us_node->GetNW();
//        ew = us_node->GetEW();
//
//        liu_peaknode->SetParentId(liu_edgenode->GetId());
//        liu_peaknode->SetEW(0);
//        double lnw = ew+nw;
//        /*copy and compute weight of children*/
//        for(unsigned int j =0;j<us_node->GetChildren()->size();j++){
//            lnw+= us_node->GetChild(j)->GetEW();
//            liu_peaknode->AddChild(tree_liu->GetNode(us_node->GetChild(j)->GetId()));
//            tree_liu->GetNode(us_node->GetChild(j)->GetId())->SetParentId(liu_peaknode->GetId());
//        }
//        liu_peaknode->SetNW(lnw);
//
//        liu_edgenode->SetParentId(us_node->GetParentId());
//        liu_edgenode->AddChild(liu_peaknode);
//        liu_edgenode->SetEW(0);
//        liu_edgenode->SetNW(ew);
//
//        last_node++;
//    }
//}


//void ConvertToLiu(const int * oldprnts,const double * oldnwghts,const double * oldewghts, int N,const int* chstart,const int * children, int ** pprnts, double ** pnwghts, double ** pewghts) {
//
//
//    *pprnts = new int[2*N+1];
//    memcpy(*pprnts,oldprnts,(N+1)*sizeof(int));
//    *pnwghts = new double[2*N+1];
//    *pewghts = new double[2*N+1];
//
//    //double max_weight = 0;
//    unsigned int last_node = N+1;
//    for(int node_index=1;node_index< N + 1;node_index++){
//        double ew, nw ;
//
//        int liu_enode = node_index;
//        int liu_pnode = last_node;
//
//
//        nw = oldnwghts[node_index];
//        ew = oldewghts[node_index];
//
//        (*pprnts)[liu_pnode] = liu_enode;
//
//        double lnw = ew+nw;
//        /*copy and compute weight of children*/
//        for(int j = chstart[node_index]; j<chstart[node_index+1];j++){
//            int ch = children[j];
//            lnw+=oldewghts[ch];
//            (*pprnts)[ch] = liu_pnode;
//        }
//
//        (*pewghts)[liu_pnode] = 0;
//        (*pnwghts)[liu_pnode] = lnw;
//        //max_weight=max(max_weight,lnw);
//        (*pewghts)[liu_enode] = 0;
//        (*pnwghts)[liu_enode] = ew;
//
//        //if(liu_enode==308){
//        //cerr<<liu_pnode<<"("<< (*pnwghts)[liu_pnode]<<") >> "<<liu_enode<<"("<<(*pnwghts)[liu_enode]<<") children = {";
//        //        for(int j = chstart[node_index]; j<chstart[node_index+1];j++){
//        //            int ch = children[j];
//        //            cerr<<ch<<"("<<oldewghts[ch]<<") ";
//        //        }
//        //cerr<<"}"<<endl;
//        //}
//        last_node++;
//    }
//    //cerr<<"max "<<max_weight<<endl;
//}

void poaux(const int * chstart, const int * children, int N, int r, int * por, int *label){
    int * stack = new int[N+2];
    memset(stack,0,(N+2)*sizeof(*stack));
    
    int push = 1;
    stack[push]=r;
    push++;
    int * chtmp = (int *)malloc((N+2)*sizeof(int));
    memcpy(chtmp,chstart,(N+2)*sizeof(*chstart));
    
    
    while(push>1){
        int top = stack[push-1];
        if(chstart[top+1] - chtmp[top]>0){
            int ch = children[chtmp[top]];
            chtmp[top]++;
            stack[push] = ch;
            push++;
        }
        else
        {
            por[*label] = top;
            (*label)++;
            push--;
        }
    }
    
    delete[] chtmp;
    delete[] stack;
}
//(*chstart)[j]=i   (*children)[(*chstart)[j]]=i      j 是i的父结点
void po_construct(const int N, const int * prnts, int **chstart,int **chend,int **children, int * root){
    *chend = new int[N+2];
    *chstart = new int[N+2];
    *children = new int[N+1];
    memset ( (void *) *chstart, 0, (N+2)*sizeof(**chstart) );
    memset ( (void *) *children, 0, (N+1)*sizeof(**children) );
    
    *root = -1;
    // (*chstart)[j]=i  结点j的孩子个数为i
    for(int ii=1;ii<N+1;ii++){
        if(prnts[ii] > 0){
            (*chstart)[prnts[ii]]++;
        }
        else{
            *root = ii;
        }
    }
//    //from he
//    for(int ii=1;ii<N+1;ii++){
//	cout<<(*chstart)[ii]<<"  ";
//    }
//    cout<<endl<<endl<<endl<<endl<<endl<<endl;
    /*compute cumsum*/
    //(*chstart)[j]=i  结点j之前的结点数总和为i（不含结点j）
    int cum_val =1;
    for(int ii=1;ii<N+2;ii++){
        int val = cum_val;
        cum_val+= (*chstart)[ii];
        (*chstart)[ii] = val;
    }
    
//    for(int ii=1;ii<N+1;ii++){
//	cout<<(*chstart)[ii]<<"  ";
//    }    
//    cout<<endl<<endl<<endl<<endl<<endl;
    memcpy(*chend,*chstart,(N+2)*sizeof(**chstart));
    
    for(int ii=1;ii<N+1;ii++){
        if(prnts[ii] > 0){
	    //(*children)[j]=i  结点i之前的结点总数为j
            (*children)[(*chend)[prnts[ii]]] = ii;
	    //(*chend)[j]=i  结点j之前的结点数总和为i（含结点j）
            (*chend)[prnts[ii]]++;
        }
    }

   //from he
//    for(int ii=1;ii<10;ii++){
// 	cout<<  (*chstart)[ii] <<":"<<(*chend)[ii]<<":"<<(*children)[(*chend)[ii]]<<"  ";
//    }
//    cout<<endl;
//    //from he
//    for(int ii=1;ii<10;ii++){
//	cout<<(*chend)[ii]<<"  ";
//    }
//    cout<<endl;
}


double IOCounter(Ctree & tree, schedule_t & sub_schedule, double available_memory,bool divisible,int quiet){
    double memory_occupation = 0;
    double io_volume = 0;
    io_map unloaded_nodes;
    schedule_t loaded_nodes;
    
    /*iterates through the given permutation (schedule)*/
    for(schedule_t::iterator cur_task_id = sub_schedule.begin();cur_task_id!=sub_schedule.end();cur_task_id++){
        Cnode * cur_node = tree.GetNode(*cur_task_id);
        
        /*if the node was unloaded*/
        if(unloaded_nodes.find(*cur_task_id)!=unloaded_nodes.end()){
            if(!quiet){cerr<<"Loading "<< unloaded_nodes[*cur_task_id] <<"of "<<*cur_task_id<<" (IO)"<<endl;}
            
            memory_occupation += unloaded_nodes[*cur_task_id];
            unloaded_nodes.erase(*cur_task_id);
        }
        
        double data_to_unload = memory_occupation + cur_node->GetCost() - cur_node->GetEW() - available_memory;
        if(!quiet){cerr<<"min data to unload "<<data_to_unload<<endl;}
        if(data_to_unload > 0 )
        {
            /*if we dont have enough room, unload files and update both io and occupation*/
            double unloaded_data = 0;
            /*unload furthest non unloaded node which is NOT in current_node children first*/
            schedule_t::reverse_iterator far_node_id = loaded_nodes.rbegin();
            while((far_node_id!=loaded_nodes.rend()) && (unloaded_data < data_to_unload) ){
                /*try to unload this node*/
                bool is_already_unloaded = false;
                double remaining_loaded_data = tree.GetNode(*far_node_id)->GetEW();
                if(unloaded_nodes.find(*far_node_id)!=unloaded_nodes.end()){
                    is_already_unloaded = true;
                    remaining_loaded_data = max(remaining_loaded_data - unloaded_nodes[*far_node_id],0.0);
                }
                
                double local_data_to_unload;
                if(divisible){
                    local_data_to_unload = min(remaining_loaded_data, data_to_unload);
                }
                else{
                    local_data_to_unload = remaining_loaded_data;
                }
                
                if(!quiet){cerr<<"unloading (IO) "<<local_data_to_unload<<" of "<<*far_node_id<<endl;}
                unloaded_data += local_data_to_unload;
                if(is_already_unloaded){
                    unloaded_nodes[*far_node_id]+=local_data_to_unload;
                }
                else{
                    unloaded_nodes[*far_node_id] = local_data_to_unload;
                }
                
                if(remaining_loaded_data  == local_data_to_unload){
                    loaded_nodes.remove(*far_node_id);
                }
                
                far_node_id++;
            }
            
            io_volume+=unloaded_data;
            memory_occupation -= unloaded_data;
        }
        
        if(!quiet){cerr<<"occupation before processing "<<memory_occupation<<endl;}
        /*if we have enough memory to process the node, update occupation*/
        memory_occupation += cur_node->GetCost() - 2*cur_node->GetEW() - cur_node->GetNW();
        if(!quiet){cerr<<"processing "<<*cur_task_id<<endl;}
        if(!quiet){cerr<<"unloading "<<*cur_task_id<<endl;}
        loaded_nodes.remove(*cur_task_id);
        
        if(!quiet){cerr<<"loading ";}
        for(vector<Cnode*>::iterator child = cur_node->GetChildren()->begin();child!=cur_node->GetChildren()->end();child++){
            if(!quiet){cerr<<(*child)->GetId()<<" ";}
            loaded_nodes.push_back((*child)->GetId());
        }
        if(!quiet){cerr<<endl;}
        if(!quiet){cerr<<"New occupation after processing "<<memory_occupation<<endl;}
    }
    
    return io_volume;
    //    cerr<<"IO Volume "<<io_volume<<endl;
}

double unload_largest_first_fit(Ctree* tree, vector<unsigned int> & unloaded_nodes,list<node_ew> & loaded_nodes,const double data_to_unload,double *ewghts){
    double unloaded_data = 0.0;
    
    /*loaded_nodes already sorted, unload largest nodes till there is enough space*/
    list<node_ew>::iterator largest_node = loaded_nodes.begin();
    while ((largest_node!=loaded_nodes.end())&&(unloaded_data<data_to_unload)) {
        largest_node = loaded_nodes.begin();
        tree->GetNode(largest_node->first)->BreakEdge();//break this edge;
        //cout<<"******Largest first: break edge "<<largest_node->first<<endl;
        double edge_weight=largest_node->second;
        unloaded_data+=edge_weight;
        unloaded_nodes.push_back(largest_node->first);
        loaded_nodes.pop_front();
    }
    
    return unloaded_data;
}


double unload_furthest_nodes(Ctree* tree, vector<unsigned int> & unloaded_nodes, list<node_sche> & loaded_nodes,const double data_to_unload, double * ewghts, bool divisible){
    double unloaded_data = 0.0;
    
    /*unload furthest non unloaded node which is NOT in current_node children first*/
    list<node_sche>::iterator far_node = loaded_nodes.begin();
    while((far_node!=loaded_nodes.end()) && (unloaded_data < data_to_unload) ){
        /*try to unload this node*/
        far_node=loaded_nodes.begin();
        double remaining_loaded_data = ewghts[(*far_node).first];
        double local_data_to_unload;
        if(divisible){
            local_data_to_unload = min(remaining_loaded_data, data_to_unload);
        }
        else{
            local_data_to_unload = remaining_loaded_data;
        }
#if VERBOSE
        //cerr<<"unloading (IO) "<<local_data_to_unload<<" of "<<*far_node<<endl;
#endif
        unloaded_data += local_data_to_unload;
        
        unloaded_nodes.push_back((*far_node).first);
        
        ////cout<<"-------LSNF remove "<<local_data_to_unload<<endl;
        ////cout<<"break edge "<<far_node->first<<endl;
        tree->GetNode(far_node->first)->BreakEdge();//break this edge
        loaded_nodes.pop_front();
    }
    
    return unloaded_data;
}

double unload_furthest_first_fit(Ctree* tree, vector<unsigned int> & unloaded_nodes, list<node_sche> & loaded_nodes,const double data_to_unload, double * ewghts, bool divisible){
    double unloaded_data = 0.0;
    /*unload furthest non unloaded node which is NOT in current_node children first*/
    list<node_sche>::iterator far_node = loaded_nodes.begin();
    while((far_node!=loaded_nodes.end()) && (unloaded_data < data_to_unload) ){
        /*try to unload this node*/
        
        double remaining_loaded_data = ewghts[far_node->first];
        
        double local_data_to_unload;
        if(divisible){
            // divisible  内存是否可切割？？？可切割的话，就只移除所需内存的大小
            //local_data_to_unload 要卸载的数据大小 和 所需数据的  最小值 
            local_data_to_unload = min(remaining_loaded_data, data_to_unload);
        }
        else{
            local_data_to_unload = remaining_loaded_data;
        }
#if VERBOSE
        //cerr<<"unloading (IO) "<<local_data_to_unload<<" of "<<*far_node<<endl;
#endif
        //如果该还在内存中的数据 大于 还需要的数据大小
        /*if it "fits", that is if the amount of data is lower than what we need to unload*/
        if(local_data_to_unload >= data_to_unload){
            //cout<<"-------First fit success: ";
            //cout<<"break edge "<<far_node->first<<endl;
            unloaded_data += local_data_to_unload;
            unloaded_nodes.push_back(far_node->first);
            tree->GetNode(far_node->first)->BreakEdge();//break this edge
            loaded_nodes.erase(far_node);
            return unloaded_data;
        }
        far_node++;
    }
    
    /*if we did not unloaded enough data, call furthest node*/
    ////cout<<"-------First fit failed, go to LSNF: "<<endl;
    unloaded_data = unload_furthest_nodes(tree, unloaded_nodes, loaded_nodes,data_to_unload, ewghts, divisible);
    
    return unloaded_data;
}

double unload_furthest_best_fit(io_map & unloaded_nodes, schedule_t & loaded_nodes,const double data_to_unload, double * ewghts, bool divisible){
    double unloaded_data = 0.0;
    /*unload furthest non unloaded node which is NOT in current_node children first*/
    while(unloaded_data < data_to_unload){
        bool is_best_already_unloaded = false;
        double best_remaining_loaded_data = 0.0;
        unsigned int best_candidate = -1;
        double best_candi_score = -1.0;
        
        schedule_t::reverse_iterator far_node_id = loaded_nodes.rbegin();
        while((far_node_id!=loaded_nodes.rend())){
            /*try to unload this node*/
            bool is_already_unloaded = false;
            double remaining_loaded_data = ewghts[*far_node_id];
            if(unloaded_nodes.find(*far_node_id)!=unloaded_nodes.end()){
                is_already_unloaded = true;
                remaining_loaded_data = max(remaining_loaded_data - unloaded_nodes[*far_node_id],0.0);
            }
            double local_data_to_unload;
            if(divisible){
                local_data_to_unload = min(remaining_loaded_data, data_to_unload);
            }
            else{
                local_data_to_unload = remaining_loaded_data;
            }
#if VERBOSE
            cerr<<"unloading (IO) "<<local_data_to_unload<<" of "<<*far_node_id<<endl;
#endif
            
            /*if it "fits", that is if the amount of data is lower than what we need to unload*/
            if(local_data_to_unload <= data_to_unload - unloaded_data){
                if(local_data_to_unload > best_candi_score){
                    best_candi_score = local_data_to_unload;
                    best_candidate = *far_node_id;
                    is_best_already_unloaded = is_already_unloaded;
                    best_remaining_loaded_data = remaining_loaded_data;
                }
            }
            
            far_node_id++;
        }
        
        
        if(best_candi_score != -1){
            //cerr<<"best found :"<<best_candidate<<endl;
            unloaded_data += best_candi_score;
            if(is_best_already_unloaded){
                unloaded_nodes[best_candidate]+=best_candi_score;
            }
            else{
                unloaded_nodes[best_candidate] = best_candi_score;
            }
            if(best_remaining_loaded_data  == best_candi_score){
                loaded_nodes.remove(best_candidate);
            }
        }
        else{
            break;
        }
    }
    /*if we did not unloaded enough data, call furthest node*/
    if(unloaded_data<data_to_unload){
        //unloaded_data += unload_furthest_nodes(unloaded_nodes, loaded_nodes,data_to_unload - unloaded_data, ewghts, divisible);
    }
    
    
    
    return unloaded_data;
}






double unload_furthest_first_fit_abs(io_map & unloaded_nodes, schedule_t & loaded_nodes,const double data_to_unload, double * ewghts, bool divisible){
    double unloaded_data = 0.0;
    /*unload furthest non unloaded node which is NOT in current_node children first*/
    schedule_t::reverse_iterator far_node_id = loaded_nodes.rbegin();
    while((far_node_id!=loaded_nodes.rend()) && (unloaded_data < data_to_unload) ){
        /*try to unload this node*/
        bool is_already_unloaded = false;
        double remaining_loaded_data = ewghts[*far_node_id];
        if(unloaded_nodes.find(*far_node_id)!=unloaded_nodes.end()){
            is_already_unloaded = true;
            remaining_loaded_data = max(remaining_loaded_data - unloaded_nodes[*far_node_id],0.0);
        }
        double local_data_to_unload;
        if(divisible){
            local_data_to_unload = min(remaining_loaded_data, data_to_unload);
        }
        else{
            local_data_to_unload = remaining_loaded_data;
        }
#if VERBOSE
        cerr<<"unloading (IO) "<<local_data_to_unload<<" of "<<*far_node_id<<endl;
#endif
        
        /*if it "fits", that is if the amount of data is lower than what we need to unload*/
        if(local_data_to_unload >= data_to_unload - unloaded_data){
            unloaded_data += local_data_to_unload;
            if(is_already_unloaded){
                unloaded_nodes[*far_node_id]+=local_data_to_unload;
            }
            else{
                unloaded_nodes[*far_node_id] = local_data_to_unload;
            }
            if(remaining_loaded_data  == local_data_to_unload){
                loaded_nodes.remove(*far_node_id);
            }
        }
        
        far_node_id++;
    }
    
    return unloaded_data;
}



double unload_furthest_best_fit_abs(io_map & unloaded_nodes, schedule_t & loaded_nodes,const double data_to_unload, double * ewghts, bool divisible){
    double unloaded_data = 0.0;
    /*unload furthest non unloaded node which is NOT in current_node children first*/
    while(unloaded_data < data_to_unload){
        bool is_best_already_unloaded = false;
        double best_remaining_loaded_data = 0.0;
        unsigned int best_candidate = -1;
        double best_candi_score = -1.0;
        
        schedule_t::reverse_iterator far_node_id = loaded_nodes.rbegin();
        while((far_node_id!=loaded_nodes.rend())){
            /*try to unload this node*/
            bool is_already_unloaded = false;
            double remaining_loaded_data = ewghts[*far_node_id];
            if(unloaded_nodes.find(*far_node_id)!=unloaded_nodes.end()){
                is_already_unloaded = true;
                remaining_loaded_data = max(remaining_loaded_data - unloaded_nodes[*far_node_id],0.0);
            }
            double local_data_to_unload;
            if(divisible){
                local_data_to_unload = min(remaining_loaded_data, data_to_unload);
            }
            else{
                local_data_to_unload = remaining_loaded_data;
            }
#if VERBOSE
            cerr<<"unloading (IO) "<<local_data_to_unload<<" of "<<*far_node_id<<endl;
#endif
            
            /*if it "fits", that is if the amount of data is lower than what we need to unload*/
            if(abs(data_to_unload-unloaded_data-local_data_to_unload)<abs(data_to_unload-unloaded_data-best_candi_score)){
                best_candi_score = local_data_to_unload;
                best_candidate = *far_node_id;
                is_best_already_unloaded = is_already_unloaded;
                best_remaining_loaded_data = remaining_loaded_data;
            }
            
            far_node_id++;
        }
        
        
        if(best_candi_score != -1){
            //cerr<<"best found :"<<best_candidate<<endl;
            unloaded_data += best_candi_score;
            if(is_best_already_unloaded){
                unloaded_nodes[best_candidate]+=best_candi_score;
            }
            else{
                unloaded_nodes[best_candidate] = best_candi_score;
            }
            if(best_remaining_loaded_data  == best_candi_score){
                loaded_nodes.remove(best_candidate);
            }
        }
        else{
            break;
        }
    }
    
    return unloaded_data;
}











/*
 next_comb(int comb[], int k, int n)
 Generates the next combination of n elements as k after comb
 
 comb => the previous combination ( use (0, 1, 2, ..., k) for first)
 k => the size of the subsets to generate
 n => the size of the original set
 
 Returns: 1 if a valid combination was found
 0, otherwise
 */
int next_comb(int comb[], int k, int n) {
    int i = k - 1;
    ++comb[i];
    
    while ((i > 0) && (comb[i] >= n - k + 1 + i)) {
        --i;
        ++comb[i];
    }
    
    if (comb[0] > n - k) /* Combination (n-k, n-k+1, ..., n) reached */
        return 0; /* No more combinations can be generated */
    
    /* comb now looks like (..., x, n, n, n, ..., n).
     *  Turn it into (..., x, x + 1, x + 2, ...) */
    for (i = i + 1; i < k; ++i)
        comb[i] = comb[i - 1] + 1;
    
    return 1;
}



double unload_best_increasing_combi(io_map & unloaded_nodes, schedule_t & loaded_nodes,const double data_to_unload, double * ewghts, bool divisible, unsigned int init_combi_size, bool quiet){
    assert(!divisible);
    vector<unsigned int> candidates;
    double unloaded_data = 0.0;
    vector<unsigned int> best_combi;
    double best_combi_score = numeric_limits<double>::infinity();
    while(unloaded_data<data_to_unload){
        /*unload furthest non unloaded node which is NOT in current_node children first*/
        candidates.clear();
        schedule_t::reverse_iterator far_node_id = loaded_nodes.rbegin();
        while((far_node_id!=loaded_nodes.rend()) && (candidates.size() <= init_combi_size) ){
            /*try to unload this node*/
            if(unloaded_nodes.find(*far_node_id)==unloaded_nodes.end()){
                candidates.push_back(*far_node_id);
            }
            
            far_node_id++;
        }
        /*now we got our max_candidates candidates, compute every permutations*/
        
        vector<unsigned int> cur_combi;
        
        int n = candidates.size(); /* The size of the set; for {1, 2, 3, 4} it's 4 */
        for(int k = 1;k<=n;k++){
            /* k is the size of the subsets; for {1, 2}, {1, 3}, ... it's 2 */
            int comb[k]; /* comb[i] is the index of the i-th element in the
                          combination */
            
            /* Setup comb for the initial combination */
            for (int i = 0; i < k; ++i)
                comb[i] = i;
            /*compute score*/
            cur_combi.clear();
            double cur_unloaded_data = 0.0;
            if(!quiet){cerr<<"cur_combi is [";}
            for (int i = 0; i < k; ++i){
                unsigned int cur_node_id = candidates[comb[i]];
                cur_unloaded_data+= ewghts[cur_node_id];
                cur_combi.push_back(cur_node_id);
                if(!quiet){cerr<<" "<<cur_node_id;}
            }
            if(!quiet){cerr<<"], its score is "<<cur_unloaded_data<<endl;}
            
            
            if(abs(data_to_unload-unloaded_data-cur_unloaded_data)<abs(data_to_unload-unloaded_data-best_combi_score)){
                if(!quiet){cerr<<"cur_combi is better than prev best "<<cur_unloaded_data<<" vs "<<best_combi_score<<endl;}
                best_combi.assign(cur_combi.begin(),cur_combi.end());
                best_combi_score = cur_unloaded_data;
            }
            
            
            /* Generate and print all the other combinations */
            while (next_comb(comb, k, n)){
                /*compute score*/
                cur_combi.clear();
                double cur_unloaded_data = 0.0;
                if(!quiet){cerr<<"cur_combi is [";}
                for (int i = 0; i < k; ++i){
                    unsigned int cur_node_id = candidates[comb[i]];
                    cur_unloaded_data+= ewghts[cur_node_id];
                    cur_combi.push_back(cur_node_id);
                    if(!quiet){cerr<<" "<<cur_node_id;}
                }
                if(!quiet){cerr<<"], its score is "<<cur_unloaded_data<<endl;}
                
                if(abs(data_to_unload-unloaded_data-cur_unloaded_data)<abs(data_to_unload-unloaded_data-best_combi_score)){
                    if(!quiet){cerr<<"cur_combi is better than prev best "<<cur_unloaded_data<<" vs "<<best_combi_score<<endl;}
                    best_combi.assign(cur_combi.begin(),cur_combi.end());
                    best_combi_score = cur_unloaded_data;
                }
            }
        }
        
        unloaded_data = best_combi_score;
        
        if(unloaded_data < data_to_unload){
            init_combi_size*=2;
        }
        
        
    }
    if(!quiet){cerr<<"best_combi is [";}
    for(unsigned int i=0;i<best_combi.size();i++){
        unsigned int cur_id = best_combi[i];
        unloaded_nodes[cur_id] = ewghts[cur_id];
        loaded_nodes.remove(cur_id);
        if(!quiet){cerr<<" "<<cur_id;}
    }
    if(!quiet){cerr<<"], its score is "<<best_combi_score<<endl;}
    
    if(!quiet){cerr<<"We've unloaded "<<unloaded_data<<" out of "<<data_to_unload<<" so far"<<endl;}
    unloaded_data = best_combi_score;
    return unloaded_data;
}







double unload_best_furthest_nodes(io_map & unloaded_nodes, schedule_t & loaded_nodes,const double data_to_unload, double * ewghts, bool divisible, unsigned int max_candidates,bool quiet){
    assert(!divisible);
    vector<unsigned int> candidates;
    double unloaded_data = 0.0;
    while(unloaded_data<data_to_unload){
        /*unload furthest non unloaded node which is NOT in current_node children first*/
        candidates.clear();
        schedule_t::reverse_iterator far_node_id = loaded_nodes.rbegin();
        while((far_node_id!=loaded_nodes.rend()) && (candidates.size() <= max_candidates) ){
            /*try to unload this node*/
            if(unloaded_nodes.find(*far_node_id)==unloaded_nodes.end()){
                candidates.push_back(*far_node_id);
            }
            
            far_node_id++;
        }
        /*now we got our max_candidates candidates, compute every permutations*/
        
        
        
        vector<unsigned int> best_combi;
        vector<unsigned int> cur_combi;
        double best_combi_score = numeric_limits<double>::infinity();
        
        int n = candidates.size(); /* The size of the set; for {1, 2, 3, 4} it's 4 */
        for(int k = 1;k<=n;k++){
            /* k is the size of the subsets; for {1, 2}, {1, 3}, ... it's 2 */
            int comb[k]; /* comb[i] is the index of the i-th element in the
                          combination */
            
            /* Setup comb for the initial combination */
            for (int i = 0; i < k; ++i)
                comb[i] = i;
            /*compute score*/
            cur_combi.clear();
            double cur_unloaded_data = 0.0;
            if(!quiet){cerr<<"cur_combi is [";}
            for (int i = 0; i < k; ++i){
                unsigned int cur_node_id = candidates[comb[i]];
                cur_unloaded_data+= ewghts[cur_node_id];
                cur_combi.push_back(cur_node_id);
                if(!quiet){cerr<<" "<<cur_node_id;}
            }
            if(!quiet){cerr<<"], its score is "<<cur_unloaded_data<<endl;}
            
            
            if((double)abs(data_to_unload-unloaded_data-cur_unloaded_data)<(double)abs(data_to_unload-unloaded_data-best_combi_score)){
                if(!quiet){cerr<<"cur_combi is better than prev best "<<cur_unloaded_data<<" vs "<<best_combi_score<<endl;}
                best_combi.assign(cur_combi.begin(),cur_combi.end());
                best_combi_score = cur_unloaded_data;
            }
            
            
            /* Generate and print all the other combinations */
            while (next_comb(comb, k, n)){
                /*compute score*/
                cur_combi.clear();
                double cur_unloaded_data = 0.0;
                if(!quiet){cerr<<"cur_combi is [";}
                for ( int i = 0; i < k; ++i){
                    unsigned int cur_node_id = candidates[comb[i]];
                    cur_unloaded_data+= ewghts[cur_node_id];
                    cur_combi.push_back(cur_node_id);
                    if(!quiet){cerr<<" "<<cur_node_id;}
                }
                if(!quiet){cerr<<"], its score is "<<cur_unloaded_data<<endl;}
                
                if(abs(data_to_unload-unloaded_data-cur_unloaded_data)<abs(data_to_unload-unloaded_data-best_combi_score)){
                    if(!quiet){cerr<<"cur_combi is better than prev best "<<cur_unloaded_data<<" vs "<<best_combi_score<<endl;}
                    best_combi.assign(cur_combi.begin(),cur_combi.end());
                    best_combi_score = cur_unloaded_data;
                }
            }
        }
        
        
        if(!quiet){cerr<<"best_combi is [";}
        for(unsigned int i=0;i<best_combi.size();i++){
            unsigned int cur_id = best_combi[i];
            unloaded_nodes[cur_id] = ewghts[cur_id];
            loaded_nodes.remove(cur_id);
            if(!quiet){cerr<<" "<<cur_id;}
        }
        if(!quiet){cerr<<"], its score is "<<best_combi_score<<endl;}
        
        unloaded_data += best_combi_score;
        if(!quiet){cerr<<"We've unloaded "<<unloaded_data<<" out of "<<data_to_unload<<" so far"<<endl;}
        
    }
    
    return unloaded_data;
}

//  N 子树大小+1   spacewghts = nwght 结点权重      divisible =  内存是否可切割？？？   com_freq = 切割树的次数    brokenEdges
double IOCounter(Ctree* tree, int N, double * nwghts, double * ewghts, int * chstart,int * children, int * schedule, double available_memory,bool divisible,int quiet,unsigned int & com_freq, vector<unsigned int>* brokenEdges, io_method_t method){
    //memory_occupation第一个结点的输入文件大小
    double memory_occupation = ewghts[schedule[N-1]];
    double io_volume = 0;
    vector<unsigned int> unloaded_nodes;
    // typedef pair<unsigned int, unsigned int> node_sche
    list<node_sche> loaded_nodes;
    // typedef pair<unsigned int, double> node_ew;
    list<node_ew> loaded_nodes_ew;

    vector<int> schedule_vec(schedule,schedule+N);//使用schedule 来初始化schedule_vec   schedule=schedule_vec
    int cur_task_id;
    vector<unsigned int>::iterator unloaded;
    list<unsigned int> temp;
    list<unsigned int> queue;
    unsigned int child_start, child_end;
    unsigned int subtree_size;
    list<unsigned int>::iterator iter;
    double maxoutD, memory_required, IO_sub=0;
    uint64_t count;
    schedule_t * schedule_f = new schedule_t();
    list<int>::iterator ite_sche;
    int rootid;
    vector<unsigned int> subtreeBrokenEdges;
    
    /*iterates through the given permutation (schedule)*/
    for(int rank = N-1;rank>=1;rank--){
        cur_task_id = schedule[rank];
        
        // cout<<"current task id: "<<cur_task_id<<endl;
        if (cur_task_id!=0) {//0 means this node is on other subtrees
            /*if the node was unloaded*/
            unloaded = find(unloaded_nodes.begin() , unloaded_nodes.end(), cur_task_id);//  find()在指定范围内查找和目标元素值相等的第一个元素
            //结点已经切割，那么将以该结点为根结点的子树的所有结点从 原来的树上 移除，并递归地对该子树进行内存检测
            if (unloaded!=unloaded_nodes.end()) {//find node cur_task_id unloaded
                //cout<<", (break) "<<endl;
                brokenEdges->push_back(tree->GetNode(cur_task_id)->GetothersideID());
                ++com_freq;
                unloaded_nodes.erase(unloaded);
                temp.clear();
                queue.push_back(cur_task_id);
                //cout<<"children "<<endl;

             //将以cur_task_id为根的子树上的所有结点入队temp
                do {
                    child_start=*(chstart+queue.front());
                    child_end=*(chstart+queue.front()+1);
                    temp.push_back(queue.front());
                    queue.pop_front();
                    for (unsigned int i=child_start; i<child_end; ++i) {
                        //cout<<*(children+i)<<" ";
                        queue.push_back(*(children+i));
                    }
                } while (!queue.empty());
                //cout<<endl;
                
             //将属于子树的结点的遍历顺序置为0   作用 ：？？？该子树上的结点 已经从原来的树上切割出来
                subtree_size=temp.size();//just an approximation
                for (long i=rank-1; i>=0;i-- ) {
                    iter=find(temp.begin(), temp.end(), schedule[i]);
                    if (iter!=temp.end()) {
                        schedule[i]=0;//IO counter will pass 0;
                        temp.erase(iter);
                    }
                    if (temp.size()==1) {
                        break;
                    }
                }
                
                double *ewghtssub, *timewghtssub, *spacewghtssub;
                int *prntssub;
                Ctree* subtree = BuildSubtree(tree, tree->GetNode(cur_task_id), subtree_size, &prntssub, &ewghtssub, &timewghtssub, &spacewghtssub, chstart, children);
                
                subtree_size=subtree->GetNodes()->size();
                
                int * schedule_copy = new int[subtree_size+1];
                maxoutD = MaxOutDegree(subtree, true);
                schedule_f->clear();
                count=0;
                MinMem(subtree, maxoutD, memory_required, *schedule_f, true, count);
                ite_sche = schedule_f->begin();
                for (unsigned int i=subtree_size; i>=1; --i) {
                    schedule_copy[i]=*ite_sche;
                    advance(ite_sche, 1);
                }
                schedule_copy[0]=subtree_size+1;
                int *chstartsub, *chendsub, *childrensub;
                po_construct(subtree_size, prntssub, &chstartsub,&chendsub,&childrensub, &rootid);
                
                if (memory_required>available_memory) {
                    //cout<<"memory required "<<memory_required<<", is larger than what is available "<<available_memory<<endl;
                    //cout<<"----------------------Processing subtree!"<<endl;
                    IO_sub = IOCounter(subtree,subtree_size+1, spacewghtssub, ewghtssub, chstartsub,childrensub, schedule_copy, available_memory, divisible, quiet, com_freq, &subtreeBrokenEdges, method);
                    
                    for (vector<unsigned int>::iterator iter=subtreeBrokenEdges.begin(); iter!=subtreeBrokenEdges.end(); ++iter) {
                        brokenEdges->push_back(tree->GetNode(*iter)->GetothersideID());
                    }
                    //cout<<"----------------------Out of Processing subtree!"<<endl;
                }
                
                delete [] ewghtssub;
                delete [] timewghtssub;
                delete [] spacewghtssub;
                delete [] prntssub;
                delete [] chstartsub;
                delete [] chendsub;
                delete [] childrensub;
                delete [] schedule_copy;
                delete subtree;
               
                io_volume+=IO_sub;
                
            }else{// cur_task_id     was  loaded
                // 如果结点未切割，那就从判断执行该结点内存是否足够


            // 求结点cur_task_id  的权重加上孩子结点边权之和， 即node_cost为 结点的输入文件大小+输出文件大小+执行文件大小
            double node_cost = ewghts[cur_task_id] + nwghts[cur_task_id];
            //(*children)[(*chstart)[org]]=j     pnt[j]=org  
            for(int j = chstart[cur_task_id];j<chstart[cur_task_id+1];j++){  
                node_cost += ewghts[children[j]];
            }

            // 还需要的内存空间 [2011 on optimal Tree traversals ]   p563   [2020 partition ]  5.2.1
            // MemReq(i)=fi+ni+ 所有孩子结点j的输入文件fj
            // IOReq(i) =Mavail-(MemReq(i)-fi)
            // ?  是否可以考虑检测所有的不需要再用的文件并移除 ？  不可
            double data_to_unload = memory_occupation + node_cost - ewghts[cur_task_id] - available_memory;
                
            if(!quiet){cerr<<"min data to unload "<<data_to_unload<<endl;}
                    
                if(data_to_unload > 0 )
                {
                        //cerr<<"We must commit I/O in order to process node "<<cur_task_id<<" which requires "<< memory_occupation + node_cost - ewghts[cur_task_id]<< " but has "<<available_memory<<"available"<<endl;
                        /*if we dont have enough room, unload files and update both io and occupation*/

                    // 排序
                    switch (method) {
                        //从剩余文件中选择大小最少为IOReq /data_to_unload的第一个文件
                        //从剩于文件中选择最新执行的文件
                        case FIRST_FIT:
                        /**
                         * schedule_vec() 越后面越先执行
                         * end()-loc(ch)    ch遍历的顺序越前，loc(ch) 越大， end()-loc(ch) 就越小
                         * loaded_nodes 降序排列  则选择移除文件时，选择最新latest（最后？？？）执行的文件
                         * the data that will be used for processing the latest 
                         * select the first data from I until their total size exceeds or equals Need(j)
                         * */
                            loaded_nodes.remove(make_pair(cur_task_id, schedule_vec.end()-find(schedule_vec.begin(), schedule_vec.end(), cur_task_id)));
                            loaded_nodes.sort(sort_sche);//descending schedule order
                            break;
                        // 选择最大的数据文件
                        case LARGEST_FIT:
                            loaded_nodes_ew.remove(make_pair(cur_task_id, ewghts[cur_task_id]));
                            loaded_nodes_ew.sort(sort_ew);
                            break;
                            
                        default:
                            break;
                    }
                    
                    // 从内存中移除数据文件
                    double unloaded_data = 0.0;
                    switch (method){
                        case FIRST_FIT:
                            unloaded_data = unload_furthest_first_fit(tree,unloaded_nodes,loaded_nodes,data_to_unload, ewghts, divisible);
                            break;
                        case LARGEST_FIT:
                            unloaded_data = unload_largest_first_fit(tree,unloaded_nodes,loaded_nodes_ew,data_to_unload,ewghts);
                            break;
                        default:
                            unloaded_data = unload_furthest_first_fit(tree,unloaded_nodes, loaded_nodes,data_to_unload, ewghts, divisible);
                            break;
                    }
                    io_volume+=unloaded_data;
                    //double data_to_unload = memory_occupation + node_cost - ewghts[cur_task_id] - available_memory;
                    memory_occupation -= unloaded_data;
                }
                
                if(!quiet){cerr<<"occupation before processing "<<memory_occupation<<endl;}
                    /*if we have enough memory to process the node, update occupation*/
                    //为什么是两倍的ewghts[cur_task_id]？？？
                    /**
                     * memory_occupation已占用的空间
                     * memory_occupation+node_cost - 2*ewghts[cur_task_id]  - nwghts[cur_task_id]
                     * =memory_occupation+ewghts[cur_task_id]+ nwghts[cur_task_id]+fj- 2*ewghts[cur_task_id]  - nwghts[cur_task_id]
                     * =memory_occupation-ewghts[cur_task_id]+fj
                     * 
                     * memory_occupation = memory_occupation-ewghts[cur_task_id]+fj //占用的内存-当前结点的输入文件的大小+当前结点的孩子结点所有输入文件的大小=当前需要占用的内存
                     * 
                     * memory_occupation=memory_occupation- unloaded_data//当前需要占用的内存-移除的数据大小
                     * 
                     * */
                memory_occupation += node_cost - 2*ewghts[cur_task_id] - nwghts[cur_task_id];
                memory_occupation = max(0.0,memory_occupation);
                    
                if(!quiet){cerr<<"processing "<<cur_task_id<<endl;cerr<<"unloading "<<cur_task_id<<endl;cerr<<"loading ";}
                    
                switch (method) {
                    case FIRST_FIT:
                    
                        loaded_nodes.remove(make_pair(cur_task_id, schedule_vec.end()-find(schedule_vec.begin(), schedule_vec.end(), cur_task_id)));
                        break;
                    case LARGEST_FIT:
                        loaded_nodes_ew.remove(make_pair(cur_task_id, ewghts[cur_task_id]));
                        break;
                            
                    default:
                        break;
                }
                
                // 将cur_task_id的孩子结点加入loaded_nodes
                for(int j = chstart[cur_task_id];j<chstart[cur_task_id+1];j++){
                    int ch = children[j];
                    if(!quiet){cerr<<ch<<" ";}
                    switch (method) {
                        case FIRST_FIT:
                            loaded_nodes.push_back(make_pair(ch, schedule_vec.end()-find(schedule_vec.begin(), schedule_vec.end(), ch)));
                            break;
                        case LARGEST_FIT:
                            loaded_nodes_ew.push_back(make_pair(ch, ewghts[ch]));
                            break;
                                
                        default:
                            break;
                    }
                }
                    
                if(!quiet){cerr<<endl;}
                if(!quiet){cerr<<"New occupation after processing "<<memory_occupation<<endl;}
                }
                }
                }
    delete schedule_f;
                    
                return io_volume;
                    //    cerr<<"IO Volume "<<io_volume<<endl;
}

//double IOCounter(Ctree* tree, int N, int * prnts, double * nwghts, double * ewghts, int * schedule, double available_memory,int divisible, unsigned int & com_freq, io_method_t method){
//    bool div = (divisible==0);
//    
//    int * chstart,*chend,*children;
//    int root;
//    
//    po_construct(N, prnts, &chstart,&chend,&children, &root);
//    double io_volume = IOCounter(tree,N, nwghts, ewghts, chstart,children,schedule, available_memory,div,true,com_freq, method);
//    
//    
//    delete[] chstart;
//    delete[] chend;
//    delete[] children;
//    
//    return io_volume;
//    //    cerr<<"IO Volume "<<io_volume<<endl;
//}


bool check_schedule(int * prnts,int * sched,int N){
    bool valid = true;
    
    int * rev_sched = new int[N+1];
    
    for(int i = 1; i<N+1;i++){
        rev_sched[sched[i-1]] = i;
        
        //        cerr<<"Task "<<sched[i-1]<<" is scheduled at step "<<i<<endl;
    }
    
    for(int i = 1; i<N+1;i++){
        if(prnts[i]>0){
            //            cerr<<"T"<<i<<" at "<<rev_sched[i]<<" || T"<<prnts[i]<<" at "<<rev_sched[prnts[i]]<<endl;
            valid = valid & (rev_sched[prnts[i]]>rev_sched[i]);
            if(rev_sched[prnts[i]]<rev_sched[i]){
                cerr<<"Task "<<prnts[i]<<" is before Task "<<i<<endl;
            }
        }
    }
    
    delete[] rev_sched;
    
    return valid;
}

double MaxOutDegree(Ctree * tree,int quiet){
    double max_out =0;
    double max_j = 0;
    // cout<<tree->GetNodes()->size()<<endl;
    for(unsigned int j=1; j<=tree->GetNodes()->size(); j++) {
        double cur_out = tree->GetNode(j)->GetCost();
        if(cur_out>=max_out){
            max_out = cur_out;
            max_j = tree->GetNode(j)->GetId();
        }
    }
    if(!quiet){cerr<<"Max out degree "<<max_out<<" at "<<max_j<<endl;}
    return max_out;
}

double MaxOutDegree(int N, double * nwghts, double * ewghts, int * chstart,int * children){
    double max_out =0;
    //    double leaves =0;
    //    double max_j = 0;
    for(int j=1; j<N+1; j++) {
        double cur_out = nwghts[j]+ewghts[j];
        //         if(chstart[j]==chstart[j+1]){
        //            leaves++;
        //         }
        
        for(int ch = chstart[j];ch<chstart[j+1];ch++){
            cur_out+=ewghts[children[ch]];
        }
        
        if(cur_out>=max_out){
            max_out = cur_out;
            //            max_j = j;
        }
    }
    ////cout<<leaves<<" leaves"<<endl;
    ////cout<<max_j<<" is the max node"<<endl;
    return max_out;
}

double MaxOutDegree(int N, int * prnts, double * nwghts, double * ewghts){
    int * chstart,*chend,*children;
    int root;
    
    po_construct(N, prnts, &chstart,&chend,&children, &root);
    
    double max_out = MaxOutDegree(N, nwghts, ewghts, chstart, children);
    
    delete[] chstart;
    delete[] chend;
    delete[] children;
    
    return max_out;
}

Ctree* SubtreeRooted(Cnode* node){
    Ctree* subtree = new Ctree();
    
    subtree->SetRootId(1);
    subtree->SetTreeId(node->GetId());
    subtree->AddNode(node);
    
    vector<Cnode*> visit_next;
    vector<Cnode*>::iterator first_node;
    Cnode* end_node;
    if (node->IsLeaf()) {
        return subtree;
    } else {
        visit_next=*(node->GetChildren());
        while (!visit_next.empty()) {
            if (!visit_next.back()->IsBorken()) {// this child has not been cut
                subtree->AddNode(visit_next.back());
                end_node=visit_next.back();
                visit_next.pop_back();
                if (!end_node->IsLeaf()) {
                    visit_next.insert(visit_next.end(), end_node->GetChildren()->begin(), end_node->GetChildren()->end());
                }
            }else{visit_next.pop_back();}
        }
        return subtree;
    }
}

Ctree* BuildSubtree(Ctree* tree, Cnode* SubtreeRoot, unsigned int new_tree_size, int** prnts, double** ewghts, double** timewghts, double** spacewghts, int * chstart, int * children){
    *prnts = new int[new_tree_size+1];
    *ewghts = new double[new_tree_size+1];
    *timewghts = new double[new_tree_size+1];
    *spacewghts = new double[new_tree_size+1];
    unsigned int *originalIDs = new unsigned int [new_tree_size+1];
    
    unsigned int real_tree_size=0;
    unsigned int nodeID=1;
    //初始化根结点
    (*prnts)[1] = 0;
    originalIDs[1] = SubtreeRoot->GetId();
    (*ewghts)[1] = SubtreeRoot->GetEW();
    (*timewghts)[1] = SubtreeRoot->GetMSW();
    (*spacewghts)[1] = SubtreeRoot->GetNW();
    
    Cnode* currentNode;
    list<unsigned int> que;
    unsigned int originalID=SubtreeRoot->GetId();
    que.push_back(originalID);
    unsigned int parentID;
    unsigned int tempid=SubtreeRoot->GetothersideID();
    SubtreeRoot->SetothersideID(1);
    
    while (!que.empty()) {
        originalID=que.front();
        que.pop_front();
        parentID=tree->GetNode(originalID)->GetothersideID();
        for(int j = chstart[originalID];j<chstart[originalID+1];j++){
            currentNode=tree->GetNode(children[j]);
            if (!currentNode->IsBorken()) {// broken edge means node is on aother subtree
                nodeID++;
                originalIDs[nodeID]=children[j];
                (*prnts)[nodeID]=parentID;
                (*ewghts)[nodeID]=currentNode->GetEW();
                (*timewghts)[nodeID]=currentNode->GetMSW();
                (*spacewghts)[nodeID]=currentNode->GetNW();
                que.push_back(children[j]);
                tree->GetNode(children[j])->SetothersideID(nodeID);
            }
        }
    }
    real_tree_size=nodeID;
    
    SubtreeRoot->SetothersideID(tempid);
    
    //    for (int i=1; i<=real_tree_size; ++i) {
    //        //cout<<i<<" "<<prnts[i]<<" "<<timewghts[i]<<" "<<ewghts[i]<<endl;
    //    }
    
    Ctree *treeobj = new Ctree(real_tree_size,*prnts,*spacewghts,*ewghts,*timewghts);
    
    for (unsigned int i=1; i<=real_tree_size; i++) {
        treeobj->GetNode(i)->SetothersideID(originalIDs[i]);//corresponding to the original tree's id
    }
    
    delete [] originalIDs;
    
    return treeobj;
}
double *getRandom(double E, double S, int num)
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); //每次程序运行产生不同的随机数
    std::default_random_engine gen(seed);
    std::normal_distribution<double> n(E, S); //均值为0，标准差为1
    double *number = new double[num];
    double rnd;
    for (int i = 0; i < num; i++)
    {
        while ((rnd = n(gen)) < 0)
            ;
        number[i] = rnd;
        sleep(0.1);
    }
    return number;
}

void swap(double *prower, double *mems, int low, int high)
{
    // double temp1,temp2;
    swap(prower[low], prower[high]);
    swap(mems[low], mems[high]);
}
//快速排序（从大到小）
void quickSort(double *prower, double *mems, int left, int right)
{
    int low, high, m;
    double pivot, temp;
    low = left;
    high = right;

    m = (low + high) / 2;

    // if ((prower[low] / mems[low]) > prower[high] / mems[high])
    // {
    //     swap(prower, mems, low, high);
    // }
    // if ((prower[m] / mems[m]) > (prower[high] / mems[high]))
    // {
    //     swap(prower, mems, m, high);
    // }
    // if ((prower[m] / mems[m]) > (prower[low] / mems[low]))
    // {
    //     swap(prower, mems, m, low);
    // }
    //<>
    pivot = prower[low] / mems[low];
    while (low < high)
    {
        while (low < high && (prower[high] / mems[high]) <= pivot)
            high--;
        while (low < high && (prower[low] / mems[low]) >= pivot)
            low++;
        if (low < high)
            swap(prower, mems, low, high);
    }
    if (left < low)
    {
        swap(prower, mems, left, low);
        quickSort(prower, mems, left, low - 1);
        quickSort(prower, mems, low + 1, right);
    }
}