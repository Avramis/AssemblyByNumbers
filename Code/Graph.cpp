//
//  Graph.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "Graph.hpp"
void Graph::BFSUtil(int s, std::vector<bool> *visited, std::vector<int> *path, int q, std::vector<int> *network){
    std::vector<int> queue;
    //std::vector<int> nvi;
    queue.push_back(s);
    int n = q;
    // CHECK
    // https://eli.thegreenplace.net/2015/directed-graph-traversal-orderings-and-applications-to-data-flow-analysis/
    // http://www.geeksforgeeks.org/eulerian-path-and-circuit/
    // http://www.geeksforgeeks.org/detect-cycle-in-a-graph/
    // http://www.geeksforgeeks.org/detect-negative-cycle-graph-bellman-ford/
    // http://www.geeksforgeeks.org/detect-cycle-direct-graph-using-colors/
    // http://www.geeksforgeeks.org/depth-first-traversal-for-a-graph/
    //
    //
    // https://cs.stackexchange.com/questions/69226/complexity-of-finding-hamiltonian-path-using-dfs
    // http://www.geeksforgeeks.org/backtracking-set-7-hamiltonian-cycle/
    while(!queue.empty()){
        int  i = queue[0];
        queue.erase(queue.begin());
        if((*visited)[i] == false){
            (*visited)[i] = true;
            (*network)[i] = n;
            //nvi.push_back(i);
            (*path).push_back(i);
            for(int j = 0; j < (int)adjnodes[i].size(); j++){
                if((*visited)[adjnodes[i][j].first] == false){
                    queue.push_back(adjnodes[i][j].first);
                }
                else{
                    //for (int z = 0; z < (int)nvi.size(); z++){
                    for (int z = 0; z < (int)(*path).size(); z++){
                        if((*path)[z] == adjnodes[i][j].first){
                            //if(z != (int)nvi.size()-1){
                            if(z != (int)(*path).size()-1){
                                double minval = (double)std::numeric_limits<int>::max();
                                int midx = 0;
                                for (int y = 0; y < adjnodes[(*path)[z]].size(); y++){
                                    if(minval > adjnodes[(*path)[z]][y].second){
                                        minval = adjnodes[(*path)[z]][y].second;
                                        midx = adjnodes[(*path)[z]][y].first;
                                    }
                                }
                                //if(nvi[z] == adjnodes[i][j].first){
                                //if(adjnodes[i][j].second > adjnodes[(*path)[z]][0].second){
                                //if(adjnodes[i][j].second > adjnodes[(*path)[z]][0].second){
                                //if(adjnodes[i][j].second > adjnodes[(*path)[z]][midx].second){
                                if(adjnodes[i][j].second > minval){
                                    
                                    //nvi.push_back(nvi[z]);
                                    //nvi.erase(nvi.begin()+z);
                                    //int iiii;
                                    //iiii=0;
                                    (*path).push_back((*path)[z]);
                                    (*path).erase((*path).begin()+z);
                                }
                            }
                            break;
                        }
                    }
                    n = (*network)[adjnodes[i][j].first];
                }
            }
        }
    }
    
    if(n!=q){
        //for (int i = 0; i < (int)nvi.size(); i++){
        for (int i = 0; i < (int)(*path).size(); i++){
            //(*network)[nvi[i]] = n;
            (*network)[(*path)[i]] = n;
        }
    }
    std::vector<int>().swap(queue);
    //std::vector<int>().swap(nvi);
    
};

Graph::Graph(){};

Graph::Graph(int v){
    assingV(v);
};

void Graph::assingV(int v){
    V = v;
    std::vector<std::vector<std::pair<int, double>>>(V, std::vector<std::pair<int, double>>(0, std::make_pair(std::numeric_limits<int>::infinity(), std::numeric_limits<double>::infinity()))).swap(adjnodes);
};

void Graph::addEdge(int i, int j, double s){
    adjnodes[i].push_back(std::make_pair(j, s));
};

void Graph::BFS(std::vector<std::vector<int>>  *pathlist){
    std::vector<bool> visited(V, false);
    std::vector<int> network(V, -1);
    for (int i = 0; i < V; i++){
        if (visited[i] == false){
            std::vector<int> path;
            BFSUtil(i, &visited,&path, (int)(*pathlist).size(), &network);
            if(network[i] >=  (int)(*pathlist).size()){
                pathlist->push_back(path);
            }
            else{
                (*pathlist)[network[i]].insert( (*pathlist)[network[i]].begin(),path.begin(), path.end());
                
            }
        }
    }

   /*
    std::vector<bool>((int)visited.size(), false).swap(visited);
    
    int repfound = 0;
    for (int i = 0; i < (int)(*pathlist).size(); i++){
        for (int j = 0; j < (int)(*pathlist)[i].size(); j++){
            if(visited[(*pathlist)[i][j]] == true){
                repfound ++ ;
            }
            else{
                visited[(*pathlist)[i][j]] = false;
            }
        }
    }
    */
    //std::vector<bool>((int)visited.size(), false).swap(visited);
    
    //std::vector<int>cP((int)(*pathlist).size(), 0);
    //for (int i = 0; i < (int)(*pathlist).size(); i++){
    //    countPaths((*pathlist)[i][0], (*pathlist)[i][((int)(*pathlist)[i].size() - 1)], &visited, &cP[i]);
    //    if(cP[i] > 1){
    //        int testcp = 0;
    //        testcp = 1;
    //    }
    //}

    //int test2= 0;
    //test2 = 1;
    
    
};




void Graph::countPaths(int i, int d, std::vector<bool> *visited, int *pathCount){
    (*visited)[i] = true;
    if(i == d){
        (*pathCount)++;
    }
    else{
        for (int j = 0; j < adjnodes[i].size(); j++){
            int p = adjnodes[i][j].first;
            if((*visited)[adjnodes[i][j].first] == false){
                countPaths(adjnodes[i][j].first, d, visited, pathCount);
            }
        }
    }
    (*visited)[i] = false;
};


bool Graph::isCyclic(int v, std::vector<bool> *visited, std::vector<bool> *recStack){
    if((*visited)[v] == false){
        (*visited)[v] = true;
        (*recStack)[v] = true;
        for(int i = 0; i < (int)adjnodes[v].size(); i++){
            if(v != adjnodes[v][i].first){
                if((*visited)[adjnodes[v][i].first] == false && isCyclic(adjnodes[v][i].first, visited, recStack)){
                    return true;
                }
                else if((*recStack)[adjnodes[v][i].first] == true){
                    return true;
                }
            }
        }
    }
    (*recStack)[v] = false;
    return false;
};


bool Graph::isCyclicFull(){
    
    std::vector<bool>visited((int)adjnodes.size(), false);
    std::vector<bool>recStack((int)adjnodes.size(), false);
    
    for(int i = 0; i < (int)adjnodes.size(); i++){
        if (isCyclic(i, &visited, &recStack) == true){
            return true;
        }
    }
        return false;
    
};

