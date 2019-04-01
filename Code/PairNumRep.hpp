//
//  PairNumRep.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef PairNumRep_hpp
#define PairNumRep_hpp

#include <stdio.h>
#include <vector>
#include <string>
class PairNumRep{
private:
public:
    PairNumRep(std::string read , std::vector<std::vector<double>> *rep);
    ~PairNumRep();
};
#endif /* PairNumRep_hpp */
