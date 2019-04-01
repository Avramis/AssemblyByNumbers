//
//  DataTransformations.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//
#ifndef DataTransformations_hpp
#define DataTransformations_hpp

#include <stdio.h>
#include <stdio.h>
#include <vector>
#include <string>
#include "DFTTran.hpp"
#include "DWTTran.hpp"
#include "PAATran.hpp"
#include "SAXTran.hpp"
#include "NoTran.hpp"
#include "NormalDistribution.hpp"
class DataTransformations{
private:
    std::vector <double> cutpoints;
public:
    DataTransformations();
    DataTransformations(std::string method, std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int l);
    DataTransformations(std::string method, std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int l, int s);
    
    void createTransformation(std::string method, std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int l);
    void createTransformation(std::string method, std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int l, int s);
    
    void createTransformation(std::string method, std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int sp, int x, int l);
    void createTransformation(std::string method, std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int sp, int x, int l, int s);
    
    void setCutPoints(int A);
    
    //Realease meory of cutpoints array
    void deleteCutPonts();
    
    ~DataTransformations();
};
#endif /* DataTransformations_hpp */
