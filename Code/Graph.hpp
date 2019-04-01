//
//  Graph.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//
// CHECK
// https://eli.thegreenplace.net/2015/directed-graph-traversal-orderings-and-applications-to-data-flow-analysis/

#ifndef Graph_hpp
#define Graph_hpp

#include <stdio.h>
#include <vector>
#include <limits>

class Graph{
private:
    int V;
    std::vector<std::vector<std::pair<int, double>>> adjnodes;
    void BFSUtil(int s, std::vector<bool> *visited, std::vector<int> *path, int q, std::vector<int> *network);
    
public:
    Graph();
    Graph(int v);
    void assingV(int v);
    void addEdge(int i, int j, double s);
    void BFS(std::vector<std::vector<int>>  *pathlist);
    void countPaths(int i, int d, std::vector<bool> *visited, int *pathCount);
    bool isCyclic(int v, std::vector<bool> *visited, std::vector<bool> *recStack);
    
    
    bool isCyclicFull();
};
#endif /* Graph_hpp */

