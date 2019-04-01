//
//  RealNumRep.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef RealNumRep_hpp
#define RealNumRep_hpp

#include <stdio.h>
#include <vector>
#include <string>
class RealNumRep{
private:
public:
    RealNumRep(std::string read , std::vector<std::vector<double>> *rep);
    ~RealNumRep();
};
#endif /* RealNumRep_hpp */
