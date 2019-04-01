//
//  NoTran.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "NoTran.hpp"
void NoTran::generateTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int st, int len){
    std::vector<std::vector<double>>(rep->size(), std::vector<double>(len, 0.0));
    
    for (int i = 0; i < (int) rep->size(); i++){
        int c = 0;
        for (int j = st; j < st+len; j++){
            (*tran)[i][c] = (*rep)[i][j];
            c++;
        }
    }
};


NoTran::NoTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran){
    generateTran(rep, tran, 0, (int)rep->at(0).size());
};

NoTran::NoTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int st, int len){
    generateTran(rep, tran, 0, (int)rep->at(0).size());
};

NoTran::~NoTran(){};
