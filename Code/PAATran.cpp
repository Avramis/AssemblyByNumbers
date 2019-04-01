//
//  PAATran.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "PAATran.hpp"

void PAATran::generateTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int st, int len, int l){
    int X = len;
    if(X+st > (int)rep->at(0).size()){
        X = (int)rep->at(0).size() - st;
    }
    if(X <= l){
        NoTran NT(rep, tran, st, len);
    }
    else{
        for(int i = 0; i < (int)rep->size(); i++){
            tran->push_back(std::vector<double>(l, 0.0));
            double paawave = 0.0, remain = 0.0;
            int count1 = -1;
            for(int j = 0; j < l; j++){
                count1++;
                int count2 = 1;
                //Calculate the jth paa coefficient
                while((count2 * l)+ (int)remain < X)
                {
                    paawave += ((*rep)[i][count1+st] * (double)l);
                    count2++;
                    count1++;
                }
                
                remain = ((count2 * (double)l) + remain) - (double)X;
                paawave += ((*rep)[i][count1+st] * ((double)l - remain));
                (*tran)[i][j] = (paawave/ (double) X);
                paawave = 0.0;
                paawave += ((*rep)[i][count1+st] * (double)remain);
            }
        }
    }
};


PAATran::PAATran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int l){
    generateTran(rep, tran, 0, (int)rep->at(0).size(), l);
};


PAATran::PAATran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int st, int len, int l){
    generateTran(rep, tran, st, len, l);
};

PAATran::~PAATran(){};
