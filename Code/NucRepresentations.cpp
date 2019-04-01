//
//  NucRepresentations.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "NucRepresentations.hpp"
NucRepresentations::NucRepresentations(){};

NucRepresentations::NucRepresentations(std::string method, std::string read, std::vector<std::vector<double>> *rep){
     createRepresentation(method, read, rep);
};

NucRepresentations::NucRepresentations(std::string method, std::string read, std::vector<std::vector<double>> *rep, bool accu, bool norm){
    createRepresentation(method, read, rep, accu, norm);
};

void NucRepresentations::createRepresentation(std::string method, std::string read, std::vector<std::vector<double>> *rep){
    
    if(method == "Atomic_Numbers"){
        AtomNumRep Atomic_Numbers (read , rep);;
    }
    else if(method == "Complex_Numbers"){
        CompNumRep Complex_Numbers(read , rep);
    }
    else if(method == "Dna_Walk"){
        DnaWNumRep Dna_Walk(read , rep);
    }
    else if(method == "EIIP_numbers"){
        EIIPNumRep EIIP_numbers(read , rep);
    }
    else if(method == "Integer_numbers"){
        IntNumRep Integer_numbers(read , rep);
    }
    else if(method == "Pair_numbers"){
        PairNumRep Pair_numbers(read , rep);
    }
    else if(method == "Real_numbers"){
        RealNumRep Real_numbers(read , rep);
    }
    else if(method == "Tetrahedron"){
        TetrNumRep Tetrahedron(read , rep);
    }
    else if(method == "Voss_indicators"){
        VossNumRep Voss_indicators(read , rep);
    }
    else if(method == "Z_curve"){
        ZCurNumRep Z_curve(read , rep);
    }
};

//Generate the appropriate accumulated/normilised representation according to the representation method
void NucRepresentations::createRepresentation(std::string method, std::string read, std::vector<std::vector<double>> *rep, bool accu, bool norm){
    createRepresentation(method, read, rep);
    DataProcessing DT;
    if(accu == true){
        DT.repAccumulation(rep);
    }
    if(norm == true){
        DT.repZnormalisation(rep);
    }
};

NucRepresentations::~NucRepresentations(){};
