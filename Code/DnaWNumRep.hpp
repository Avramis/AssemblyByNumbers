//
//  DnaWNumRep.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef DnaWNumRep_hpp
#define DnaWNumRep_hpp

#include <stdio.h>
#include <vector>
#include <string>
class DnaWNumRep{
private:public:
    DnaWNumRep(std::string read , std::vector<std::vector<double>> *rep);
    ~DnaWNumRep();
};
#endif /* DnaWNumRep_hpp */
