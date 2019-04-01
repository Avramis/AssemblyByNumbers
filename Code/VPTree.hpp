//
//  VPTree.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//
#ifndef VPTree_hpp
#define VPTree_hpp

#include <stdio.h>
#include <vector>
#include <map>

#include "Node.hpp"
#include "NGS.hpp"
class VPTree{
private:

    Node vp;
public:
    
    VPTree();
    void TreeGeneration(std::vector<NGS> *NanoReads);
    void buildTree(Node *vp, int n, std::vector<NGS> *dataPoints);
    
    void KnnSearch(NGS *q, int k, std::map<double, std::vector<NGS>, std::greater<double>> *KNN, int *sc, int *absc);
    void FullKnnSearch(NGS *q, int k, std::map<double, std::vector<NGS>, std::greater<double>> *KNN, int *sc, int *absc);

};
#endif /* VPTree_hpp */
