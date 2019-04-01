//
//  SAXTran.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef SAXTran_hpp
#define SAXTran_hpp

#include <stdio.h>
#include <vector>

#include "PAATran.hpp"
class SAXTran{
private:
    void generateTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int st, int len, int l, std::vector<double> *cutpoints);
public:
    SAXTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int l, std::vector<double> *cutpoints);
    SAXTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int st, int len, int l, std::vector<double> *cutpoints);
    ~SAXTran();
};
#endif /* SAXTran_hpp */
