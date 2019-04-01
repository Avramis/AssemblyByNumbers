//
//  DistanceMethods.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "DistanceMethods.hpp"

DistanceMethods::DistanceMethods(){};
//double DistanceMethods::returnRepDistance(NGS *n, NGS *m){
//    return returnDistance( n->returnRep(), m->returnRep());
//};

//double DistanceMethods::returnTraDistance(NGS *n, NGS *m){
    //std::cout << "n->name: " << *n->returnGID().at(0) << " || " << *m->returnGID().at(0) << "\n";
//    return returnDistance( n->returnTra(), m->returnTra());
//};



double DistanceMethods::returnRepDistance(NGS *n, NGS *m, bool t){
    if( t == true){
        return returnEuclidianNorm(n->returnRep(), m->returnRep());
    }
    return returnDistance(n->returnRep(), m->returnRep());
};

double DistanceMethods::returnTraDistance(NGS *n, NGS *m, bool t){
    if( t == true){
        return returnEuclidianNorm( n->returnTra(), m->returnTra());
    }
     return returnDistance( n->returnTra(), m->returnTra());
};




double DistanceMethods::returnDistance(std::vector<std::vector<double>> *n, std::vector<std::vector<double>> *m){
    double d = 0.0;
    
    int y = (int)(*n).size();
    int x = (int)(*n)[0].size();
    
    if(x > (int)(*m)[0].size()){
        x = (int)(*m)[0].size();
    }
    
    for (int j = 0; j < y; j++){
        for (int i = 0; i < x; i++){
            d += ((((*n)[j][i]) - ((*m)[j][i])) * (((*n)[j][i]) - ((*m)[j][i])));
        }
    }
    
    d = sqrt(d);
    return d;
    
};



double DistanceMethods::returnEuclidianNorm(std::vector<std::vector<double>> *n, std::vector<std::vector<double>> *m){
    int y = (int)(*n).size();
    int x = (int)(*n)[0].size();
    
    if(x > (int)(*m)[0].size()){
        x = (int)(*m)[0].size();
    }
    
    std::vector<std::vector<double>> a;
    if(y >= x){
        std::vector<std::vector<double>>(y, std::vector<double>(x, 0.00)).swap(a);
        for (int i = 0; i < y; i++){
            for(int j = 0; j < x; j++){
                a[i][j] = ((*n)[i][j]) - ((*m)[i][j]);
            }
        }
    }
    else{
        std::vector<std::vector<double>>(x, std::vector<double>(y, 0.00)).swap(a);
        for (int i = 0; i < y; i++){
            for(int j = 0; j < x; j++){
                a[j][i] = ((*n)[i][j]) - ((*m)[i][j]);
            }
        }
    }
    
    
    return svdc.SVD(&a)[0];
    
};


double DistanceMethods::returnRepDTW(NGS *n, NGS *m, int w){
    return returnDTW((*n).returnRep(), (*m).returnRep(), w);
};

double DistanceMethods::returnDTW(std::vector<std::vector<double>> *n, std::vector<std::vector<double>> *m, int w){
    double d = 0.0;
    int wi = w/2;
    if( wi > (int)(*n)[0].size()){
        wi = (int)(*n)[0].size()/2;
    }
    else if( wi > (int)(*n)[0].size()/2){
        wi = (int)(*n)[0].size()/2;
    }
    
    int yDim = (int)(*n).size();
    int xDim = (int)(*n)[0].size();
    
    if(xDim > (int)(*m)[0].size()){
        xDim = (int)(*m)[0].size();
    }
    
    for (int i = 0; i < yDim; i++){
        for(int j = 0 ; j < xDim; j++){
            int sw = j - wi;
            int ew = j + wi;
            if(sw < 0){
                sw = 0;
            }
            if(ew >= xDim){
                ew = xDim-1;
            }
            double sv = (double)std::numeric_limits<int>::max();
            double lv = (double)std::numeric_limits<int>::min();
            
            for (int z = sw; z <= ew; z++){
                if(sv > (*n)[i][z]){
                    sv = (*n)[i][z];
                }
                
                if(lv < (*n)[i][z]){
                    lv = (*n)[i][z];
                }
            }
            
            if(sv <= (*m)[i][j] && lv >= (*m)[i][j]){
            }
            else{
                double cd = (sv < (*m)[i][j]) * ((*n)[i][j] - (*m)[i][j]) + (lv > (*m)[i][j]) * ((*n)[i][j] - (*m)[i][j]);
                d += (cd * cd);
            }
        }
    }
    return sqrt(d);
};
