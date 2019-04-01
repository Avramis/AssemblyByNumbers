//
//  AtomNumRep.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef AtomNumRep_hpp
#define AtomNumRep_hpp

#include <stdio.h>
#include <vector>
#include <string>
class AtomNumRep{
private:
public:
    AtomNumRep(std::string read , std::vector<std::vector<double>> *rep);
    ~AtomNumRep();
};
#endif /* AtomNumRep_hpp */
