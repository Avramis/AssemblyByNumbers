//
//  TetrNumRep.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef TetrNumRep_hpp
#define TetrNumRep_hpp

#include <stdio.h>
#include <vector>
#include <string>
#include <math.h>

class TetrNumRep{
private:
public:
    TetrNumRep(std::string read , std::vector<std::vector<double>> *rep);
    ~TetrNumRep();
};
#endif /* TetrNumRep_hpp */
