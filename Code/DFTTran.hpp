//
//  DFTTran.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef DFTTran_hpp
#define DFTTran_hpp

#include <stdio.h>
#include <vector>
#include <math.h>

class DFTTran{
private:
    void fourier(std::vector<double> *data, int nn);
    void generateTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int st, int len, int l);
    
public:
    DFTTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int l);
    DFTTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int st, int len, int l);
    ~DFTTran();
};
#endif /* DFTTran_hpp */
