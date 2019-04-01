//
//  Node.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "Node.hpp"
void Node::setNtype(bool b){
    
    Ntype = b;
};

void Node::setINode(double t, double c, double m){
    tau = t;
    thresh = c;
    maxdis = m;
};

void Node::QuickSort(std::vector<NGS> *dataPoints, int leftmost, int rightmost){
    int i = leftmost, j = rightmost;
    double  pivot = (*dataPoints)[(leftmost + rightmost)/2].returnAlscore();
    
    while (i <= j) {
        
        while ((*dataPoints)[i].returnAlscore() < pivot){
            
            i++;
        };
        
        while ((*dataPoints)[j].returnAlscore() > pivot){
            
            j--;
        };
        
        if (i <= j) {
            std::swap((*dataPoints)[i], (*dataPoints)[j]);
            i++;
            j--;
        };
    };
    
    // recursion //
    if (leftmost < j){
        QuickSort(dataPoints, leftmost, j);
    };
    if (i < rightmost){
        QuickSort(dataPoints, i, rightmost);
    };
    
};

Node::Node(){
    tau = 0;
    thresh = 0;
    Ntype = false;
    pp=1;
    
};

int Node::randominitiation(std::vector<NGS> *dataPoints, int m, int n){
    srand ((unsigned)time(NULL));
    int v = rand() % n;
    vantagepoint  = &dataPoints->at(v);
    for(int i = 0;  i < n; i++){
        //(*dataPoints)[i].setAlscore(DM.returnEuclidianNorm(vantagepoint->returnTra(), (*dataPoints)[i].returnTra()));
        (*dataPoints)[i].setAlscore(DM.returnDistance(vantagepoint->returnTra(), (*dataPoints)[i].returnTra()));
    };
    
    
    QuickSort(dataPoints, m, n-1);
    return n-1;
};

void Node::generateNode(std::vector<NGS> *dataPoints, int m, int n, bool r){
    setNtype(true);
    if(r == true){
        n  = randominitiation(dataPoints, m, n);
    };
    vantagepoint  = &(*dataPoints)[n];
    
    if(m < n){
        
        for(int i = m;  i < n; i++){
            //(*dataPoints)[i].setAlscore(DM.returnEuclidianNorm(vantagepoint->returnTra(), (*dataPoints)[i].returnTra()));
            (*dataPoints)[i].setAlscore(DM.returnDistance(vantagepoint->returnTra(), (*dataPoints)[i].returnTra()));
        };
        QuickSort(dataPoints, m, n-1);
        
        int med = m +(((n-m))/2);
        for (int i = med+1; i < n; i++){
            if((*dataPoints)[i].returnAlscore() <= (*dataPoints)[med].returnAlscore()){
                med++;
            }
            else{
                break;
            }
        }
        n--;
        double tau = (*dataPoints)[med].returnAlscore();
        double thre =0;
        
        int q2 = med+((n-med)/2);//(3/4);
        q2 = ((n-m)*3/4);
        if (q2 > n){
            q2 = n;
        }
        if((*dataPoints)[q2].returnAlscore() >= tau){
            thre =  (*dataPoints)[q2].returnAlscore() - tau;
        }
        else{
            thre = tau - (*dataPoints)[q2].returnAlscore();
        }
        setINode(tau, thre, (*dataPoints)[n].returnAlscore());
        if(m < med){
            lNode = new Node;
            lNode->generateNode(dataPoints, m, med-1, false);
        }
        else{
            lNode = NULL;
        };
        
        if(med <= n){
            rNode = new Node;
            rNode->generateNode(dataPoints, med, n, false);
        }
        else{
            rNode = NULL;
        };
    }
    else{
        lNode = NULL;
        rNode = NULL;
    }
};

void Node::KnnSearch(NGS *q, int k, std::map<double, std::vector<NGS>, std::greater<double>> *KNN, int *sc, int *absc){
    //double dist = DM.returnEuclidianNorm(q->returnTra(), vantagepoint->returnTra());
    double dist = DM.returnDistance(q->returnTra(), vantagepoint->returnTra());
    double ftau = thresh;
    //typename std::map<double, std::vector<NGS>, std::greater<double>>::iterator it;
    std::map<double, std::vector<NGS>, std::greater<double>>::iterator it;
    
    if ((int)vantagepoint->returnGID().size()  ==  1 && q->returnSRid() == *vantagepoint->returnGID().at(0)){
        
    }
    else{
        if((int)KNN->size() < k){
            (*KNN)[dist].push_back(*vantagepoint);
        }
        else{
            it = KNN->begin();
            if(it->first >= dist){
                (*KNN)[dist].push_back(*vantagepoint);
                if((int)KNN->size() > k){
                    it = KNN->begin();
                    (*KNN).erase (it);
                }
            }
        }
    }
    
    
    it = KNN->begin();
    for (int i = 0; i < ((int)KNN->size())-1; i++){
        it++;
    }
    
    if(ftau > it->first){
        ftau = it->first; if(ftau > dist){ftau = dist;}
    }
    //ftau = ftau*(1.0/4.0);
    ftau = ftau*(2.0/4.0);
    
    if(dist<= maxdis+ftau){
        (*sc)++;
        if(dist-ftau <= tau){
            //if(dist < tau){
            if(lNode == NULL){}
            else{
                lNode->KnnSearch(q, k, KNN, sc, absc);
            }
        }
        if(dist + ftau >= tau){
            //if(dist >= tau){
            if(rNode == NULL){}
            else{
                rNode->KnnSearch(q, k, KNN, sc, absc);
            }
        }
    }
    else{
        (*absc)++;
    }
};


void Node::initiateSearch(NGS *q, int k, std::map<double, std::vector<NGS>, std::greater<double>> *KNN, int *sc, int *absc){
    KnnSearch(q, k, KNN, sc, absc);
};
