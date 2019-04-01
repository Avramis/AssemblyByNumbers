//
//  DWTTran.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef DWTTran_hpp
#define DWTTran_hpp

#include <stdio.h>
#include <vector>
#include <math.h>

#include "NoTran.hpp"
class DWTTran{
private:
    void generateTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int st, int len, int l);
public:
    DWTTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int l);
    DWTTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int st, int len, int l);
    ~DWTTran();
};
#endif /* DWTTran_hpp */
