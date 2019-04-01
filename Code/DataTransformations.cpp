//
//  DataTransformations.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "DataTransformations.hpp"

DataTransformations::DataTransformations(){};

DataTransformations::DataTransformations(std::string method, std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int l){
    createTransformation(method, rep, tran, l);
};

DataTransformations::DataTransformations(std::string method, std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int l, int s){
    createTransformation(method, rep, tran, l, s);
};

void DataTransformations::createTransformation(std::string method, std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int l){
    createTransformation(method, rep, tran, 0, (int)rep->at(0).size(), l);
};

void DataTransformations::createTransformation(std::string method, std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int l, int s){
    setCutPoints(s);
    createTransformation(method, rep, tran, 0, (int)rep->at(0).size(), l);
};


void DataTransformations::createTransformation(std::string method, std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int sp, int x, int l){
    if(method == "DFT"){
        DFTTran DFT(rep, tran, l);
    }
    else if(method == "DWT"){
        DWTTran DWT(rep, tran, sp, x, l);
    }
    else if(method == "PAA"){
        PAATran PAA(rep, tran, sp, x, l);
    }
    else if(method == "SAX"){
        DataProcessing DT;
        DT.repZnormalisation(rep);
        SAXTran SAX(rep, tran, sp, x,  l, &cutpoints);
    }
    else{
        NoTran NT(rep, tran);
    }
};
void DataTransformations::createTransformation(std::string method, std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int sp, int x, int l, int s){
    setCutPoints(s);
    createTransformation(method, rep, tran, sp, x, l);
};


void DataTransformations::setCutPoints(int A){
    deleteCutPonts();
    NormalDistribution ND(0, 1);
    ND.invcdfsequence(&cutpoints, A);
};

//Realease meory of cutpoints array
void DataTransformations::deleteCutPonts(){
    std::vector<double>().swap(cutpoints);
};

DataTransformations::~DataTransformations(){
    deleteCutPonts();
};
