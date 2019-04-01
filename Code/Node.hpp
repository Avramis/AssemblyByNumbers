//
//  Node.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef Node_hpp
#define Node_hpp

#include <stdio.h>
#include <random>
#include <stdlib.h>
#include <memory>
#include <map>
#include "NGS.hpp"
#include "DistanceMethods.hpp"

class Node{
private:
    Node *lNode;
    Node *rNode;
    double tau, thresh, maxdis;
    bool Ntype = false;
    int pp=1;
    //srand (time(nullptr));
    NGS *vantagepoint;
    DistanceMethods DM;
    
    //double distance(T t1, T t2){
    //    return (double) sqrt(pow(t1-t2,2));
    //};
    
    //double distance(NGS &nr1, NGS &nr2, bool dir);
    
    void setNtype(bool b);
    
    void setINode(double t, double c, double m);
    
    void QuickSort(std::vector<NGS> *dataPoints, int leftmost, int rightmost);
    //void QuickSort(TreePoints *tp, int leftmost, int rightmost);
public:
    Node();
    int randominitiation(std::vector<NGS> *dataPoints, int m, int n);
    //int randominitiation(TreePoints *tP, int m, int n);
    
    //void generateNode(TreePoints *tP, int m, int n, bool r);
    void generateNode(std::vector<NGS> *dataPoints, int m, int n, bool r);
    

    void initiateSearch(NGS *q, int k, std::map<double, std::vector<NGS>, std::greater<double>> *KNN, int *sc, int *absc);
    
    void KnnSearch(NGS *q, int k, std::map<double, std::vector<NGS>, std::greater<double>> *KNN, int *sc, int *absc);
    
    
};
#endif /* Node_hpp */
