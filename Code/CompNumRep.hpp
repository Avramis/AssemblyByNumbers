//
//  CompNumRep.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef CompNumRep_hpp
#define CompNumRep_hpp

#include <stdio.h>
#include <vector>
#include <string>
class CompNumRep{
private:
public:
    CompNumRep(std::string read , std::vector<std::vector<double>> *rep);
    ~CompNumRep();
};
#endif /* CompNumRep_hpp */
