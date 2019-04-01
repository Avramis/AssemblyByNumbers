//
//  NoTran.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef NoTran_hpp
#define NoTran_hpp

#include <stdio.h>
#include <vector>

class NoTran{
private:
    void generateTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int st, int len);
public:
    NoTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran);
    NoTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int st, int len);
    ~NoTran();
};
#endif /* NoTran_hpp */
