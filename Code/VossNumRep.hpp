//
//  VossNumRep.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//
#ifndef VossNumRep_hpp
#define VossNumRep_hpp

#include <stdio.h>
#include <vector>
#include <string>

class VossNumRep{
private:
public:
    VossNumRep(std::string read , std::vector<std::vector<double>> *rep);
    ~VossNumRep();
};
#endif /* VossNumRep_hpp */
