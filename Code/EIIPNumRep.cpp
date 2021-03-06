//
//  EIIPNumRep.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "EIIPNumRep.hpp"
EIIPNumRep::EIIPNumRep(std::string read , std::vector<std::vector<double>> *rep){
    
    std::vector<std::vector<double>>( 1 ,std::vector<double>(read.size(),0.0)).swap((*rep));
    //rep->push_back(std::vector<double>(read.size(),0));
    
    for (int i = 0; i < (int)read.size(); i++){
        if(read.substr(i,1) == "A" || read.substr(i,1) == "a"){
            (*rep)[0][i] = 0.1260;
        }
        else if(read.substr(i,1) == "C" || read.substr(i,1) == "c"){
            (*rep)[0][i] = 0.1340;
        }
        else if(read.substr(i,1) == "G" || read.substr(i,1) == "g"){
            (*rep)[0][i] = 0.0806;
        }
        else if(read.substr(i,1) == "T" || read.substr(i,1) == "t"){
            (*rep)[0][i] = 0.1335;
        }
        else if(read.substr(i,1) == "U" || read.substr(i,1) == "u"){
            (*rep)[0][i] = 0.1335;
        }
        else if(read.substr(i,1) == "N" || read.substr(i,1) == "n"){
            (*rep)[0][i] = 0.1185;
        }
        else if(read.substr(i,1) == "M" || read.substr(i,1) == "m"){
            (*rep)[0][i] = (0.1260 + 0.1340)/2;
        }
        else if(read.substr(i,1) == "R" || read.substr(i,1) == "r"){
            (*rep)[0][i] = (0.1260 + 0.0806)/2;
        }
        else if(read.substr(i,1) == "W" || read.substr(i,1) == "w"){
            (*rep)[0][i] = (0.1260 + 0.1335)/2;
        }
        else if(read.substr(i,1) == "S" || read.substr(i,1) == "s"){
            (*rep)[0][i] = (0.1340 + 0.0806)/2;
        }
        else if(read.substr(i,1) == "Y" || read.substr(i,1) == "y"){
            (*rep)[0][i] = (0.1340 + 0.1335)/2;
        }
        else if(read.substr(i,1) == "K" || read.substr(i,1) == "k"){
            (*rep)[0][i] = (0.0806 + 0.1335)/2;
        }
        else if(read.substr(i,1) == "V" || read.substr(i,1) == "v"){
            (*rep)[0][i] = (0.1260 + 0.1340 + 0.0806)/3;
        }
        else if(read.substr(i,1) == "H" || read.substr(i,1) == "h"){
            (*rep)[0][i] = (0.1260 + 0.1340 + 0.1335)/3;
        }
        else if(read.substr(i,1) == "D" || read.substr(i,1) == "d"){
            (*rep)[0][i] = (0.1260 + 0.0806 + 0.1335)/3;
        }
        else if(read.substr(i,1) == "B" || read.substr(i,1) == "b"){
            (*rep)[0][i] = (0.1340 + 0.0806 + 0.1335)/3.0;
        }
        else {
            (*rep)[0][i] = (0.1260 + 0.1340 + 0.0806 + 0.1335)/4.0;
        }
    }
};

EIIPNumRep::~EIIPNumRep(){};
