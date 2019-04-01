//
//  ZCurNumRep.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef ZCurNumRep_hpp
#define ZCurNumRep_hpp

#include <stdio.h>
#include <vector>
#include <string>

class ZCurNumRep{
private:
public:
    ZCurNumRep(std::string read , std::vector<std::vector<double>> *rep);
    ~ZCurNumRep();
};
#endif /* ZCurNumRep_hpp */
