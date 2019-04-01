//
//  VPTree.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "VPTree.hpp"

VPTree::VPTree(){};

void VPTree::TreeGeneration(std::vector<NGS> *dataPoints){
    buildTree(&vp, (int)dataPoints->size(), dataPoints);
};

void VPTree::buildTree(Node *vp, int n, std::vector<NGS> *dataPoints){
    vp->generateNode(dataPoints, 0, n, true);
};

void VPTree::KnnSearch(NGS *q, int k, std::map<double, std::vector<NGS>, std::greater<double>> *KNN, int *sc, int *absc){
    vp.KnnSearch(q, k, KNN, sc, absc);
};


void VPTree::FullKnnSearch(NGS *q, int k,std::map<double, std::vector<NGS>, std::greater<double>> *KNN, int *sc, int *absc){
    vp.initiateSearch(q, k, KNN, sc, absc);
};

