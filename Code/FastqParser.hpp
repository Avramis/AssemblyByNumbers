//
//  FastqParser.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef FastqParser_hpp
#define FastqParser_hpp

#include <stdio.h>
#include <iostream>
#include<vector>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include "NGS.hpp"
#include "NucRepresentations.hpp"
#include "DataTransformations.hpp"
#include "DataProcessing.hpp"

class FastqParser{
private:
    NucRepresentations rep;
    DataTransformations tra;
    DataProcessing DaPr;
    
    void processConactonatedReads(std::string filepath, std::string repmeth, int refkmerlen, std::string tranmeth, int tranlev, std::vector<NGS> *NuSequences1, std::vector<NGS> *NuSequences2, int sa, int ss);
    void processseparatedReads(std::string filepath, std::string repmeth, int refkmerlen, std::string tranmeth, int tranlev, std::vector<NGS> *NuSequences1, std::vector<NGS> *NuSequences2, int sa, int ss);
public:
    FastqParser(std::string filepath, std::string repmeth, int refkmerlen, std::string tranmeth, int tranlev, std::vector<NGS> *NuSequences, int sa, int ss);
    FastqParser(std::string filepath, std::string repmeth, int refkmerlen, std::string tranmeth, int tranlev, std::vector<NGS> *NuSequences1, std::vector<NGS> *NuSequences2, int sa, int ss, bool pe);
    //FastqParser(std::string filepath, std::string repmeth, int refkmerlen, std::string tranmeth, int tranlev, std::map<std::string, std::vector<NGS>> *NGSbins, int sa, int ss, bool bin);
    //FastqParser(std::string filepath, std::string repmeth, int refkmerlen, std::string tranmeth, int tranlev, std::vector<NGS> *NuSequences, std::map<std::string, std::vector<NGS*>> *NGSbins, int sa, int ss, bool bin);
    //FastqParser(std::string filepath, std::string repmeth, int refkmerlen, std::string tranmeth, int tranlev, std::vector<NGS> *NuSequences, std::map<std::string, std::vector<NGS*>> *NGSbins, int sa, int ss, bool bin);
    ~FastqParser();
};
#endif /* FastqParser_hpp */
