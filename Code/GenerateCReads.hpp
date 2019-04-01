//
//  GenerateCReads.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef GenerateCReads_hpp
#define GenerateCReads_hpp

#include <stdio.h>
#include <vector>
#include <map>

#include "NucRepresentations.hpp"
#include "DataTransformations.hpp"
#include "NGS.hpp"
class GenerateCReads{
private:
public:
    GenerateCReads();
    GenerateCReads(std::vector<std::pair<std::string, std::string>> *Contigs, std::string repmeth, int k, std::string tranmeth, int t, std::vector<NGS> *contigReads);
    ~GenerateCReads();
};
#endif /* GenerateCReads_hpp */
