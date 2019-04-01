//
//  GenerateCReads.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "GenerateCReads.hpp"
GenerateCReads::GenerateCReads(){};

GenerateCReads::GenerateCReads(std::vector<std::pair<std::string, std::string>> *Contigs, std::string repmeth, int k, std::string tranmeth, int t, std::vector<NGS> *contigReads){
    NucRepresentations NR;
    DataTransformations DT;
    for(int i = 0; i < (int)(*Contigs).size(); i++){
        NGS n((*Contigs)[i].first, (*Contigs)[i].second, "");
        NR.createRepresentation(repmeth, (*Contigs)[i].second, n.returnRep());
        int r = k;
        if(r > (int)(*Contigs)[i].second.size()){
            r = (int)(*Contigs)[i].second.size();
        }
        DT.createTransformation(tranmeth, n.returnRep(), n.returnTra(), 0, r,  t);
        contigReads->push_back(n);
    }
};

GenerateCReads::~GenerateCReads(){};
