//
//  EIIPNumRep.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef EIIPNumRep_hpp
#define EIIPNumRep_hpp

#include <stdio.h>
#include <vector>
#include <string>
class EIIPNumRep{
private:
public:
    EIIPNumRep(std::string read , std::vector<std::vector<double>> *rep);
    ~EIIPNumRep();
};
#endif /* EIIPNumRep_hpp */
