//
//  SAXTran.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//
#include "SAXTran.hpp"
void SAXTran::generateTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int st, int len, int l, std::vector<double> *cutpoints){
    PAATran PAA(rep, tran, st, len, l);
    for (int i = 0; i < (int)tran->size(); i++){
        for (int j = 0; j < (int)tran->at(i).size(); j++){
            int z;
            for (z = 0; z < (int)cutpoints->size(); z++){
                if((*tran)[i][j] <= (*cutpoints)[z]){
                    break;
                }
            }
            (*tran)[i][j] = (double)z;
        }
    }
}
SAXTran::SAXTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int l, std::vector<double> *cutpoints){
    generateTran(rep, tran, 0, (int)rep->at(0).size(), l, cutpoints);
};


SAXTran::SAXTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int st, int len, int l, std::vector<double> *cutpoints){
    generateTran(rep, tran, st, len, l, cutpoints);
};

SAXTran::~SAXTran(){};
