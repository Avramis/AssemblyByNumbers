//
//  PAATran.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef PAATran_hpp
#define PAATran_hpp

#include <stdio.h>
#include <vector>

#include "NoTran.hpp"

class PAATran{
private:
    void generateTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int st, int len, int l);
public:
    PAATran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int l);
    PAATran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int st, int len, int l);
    ~PAATran();
};
#endif /* PAATran_hpp */
