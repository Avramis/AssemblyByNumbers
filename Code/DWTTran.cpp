//
//  DWTTran.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "DWTTran.hpp"

void DWTTran::generateTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int st, int len, int l){
    int X = st+len;
    if(X > (int)rep->at(0).size()){
        X = (int)rep->at(0).size();
    }
    int complev = (int)(ceil(log((double)len)/log(2.0))-ceil(log((double) l)/log(2.0)));
    if(complev <= 0){
        NoTran NT(rep, tran);
    }
    else{
        int inip = (int)(pow(2,ceil(log((double)(X-st))/log(2.0))));
        
        int inip2 = (int) ceil(log((double)(X-st))/log(2.0));
        int finp2 = (int) ceil(log((double)l)/log(2.0));
        std::vector<std::vector<double>>((int)rep->size(),std::vector<double>((int)pow(2,finp2), 0.0)).swap(*tran);
        //int y = (int)rep->size();
        
        for (int i = 0; i < (int)rep->size(); i++){
            std::vector<double> temparray(inip,0.0);
            int c = 0;
            for (int j = st; j < X; j++){
                temparray[c] = (*rep)[i][j];
                c++;
            }
            c = (int)NULL;
            for(int j = 0; j < (inip2-(finp2)); j++){
                for (int z = 0; z < (int) pow(2,inip2 - (j)); z+=2){
                    temparray[z/2] = (temparray[z] + temparray[z+1])/2;
                }
                int yy , xx;
                yy = 0;
                xx = 10;
                yy+=xx;
            }
            
            for (int j = 0; j < (int)pow(2,finp2); j++){
                (*tran)[i][j] = temparray[j];
            }
            std::vector<double>().swap(temparray);
        }
    }
};


DWTTran::DWTTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int l){
    generateTran(rep, tran, 0, (int)rep->at(0).size(), l);
};


DWTTran::DWTTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int st, int len, int l){
    generateTran(rep, tran, st, len, l);
};

DWTTran::~DWTTran(){};
