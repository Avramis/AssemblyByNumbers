//
//  RealNumRep.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "RealNumRep.hpp"
RealNumRep::RealNumRep(std::string read , std::vector<std::vector<double>> *rep){
    std::vector<std::vector<double>>(1 ,std::vector<double>(read.size(),0.0)).swap((*rep));
    //rep->push_back(std::vector<double>(read.size(),0));
    
    for (int i = 0; i < (int)read.size(); i++){
        if(read.substr(i,1) == "A" || read.substr(i,1) == "a"){
            (*rep)[0][1] = 1.5;
        }
        else if(read.substr(i,1) == "C" || read.substr(i,1) == "c"){
            (*rep)[0][1] = -0.5;
        }
        else if(read.substr(i,1) == "G" || read.substr(i,1) == "g"){
            (*rep)[0][1] = 0.5;
        }
        else if(read.substr(i,1) == "T" || read.substr(i,1) == "t"){
            (*rep)[0][1] = -1.5;
        }
        else if(read.substr(i,1) == "U" || read.substr(i,1) == "u"){
            (*rep)[0][1] = -1.5;
        }
        else if(read.substr(i,1) == "N" || read.substr(i,1) == "n"){
            (*rep)[0][1] = 0.0;
        }
        else if(read.substr(i,1) == "M" || read.substr(i,1) == "m"){
            (*rep)[0][1] = (1.5 + -0.5)/2.0;
        }
        else if(read.substr(i,1) == "R" || read.substr(i,1) == "r"){
            (*rep)[0][1] = (1.5 + 0.5)/2.0;
        }
        else if(read.substr(i,1) == "W" || read.substr(i,1) == "w"){
            (*rep)[0][1] = (1.5 + -1.5)/2.0;
        }
        else if(read.substr(i,1) == "S" || read.substr(i,1) == "s"){
            (*rep)[0][1] = (-0.5 + 0.5)/2.0;
        }
        else if(read.substr(i,1) == "Y" || read.substr(i,1) == "y"){
            (*rep)[0][1] = (-0.5 + -1.5)/2.0;
        }
        else if(read.substr(i,1) == "K" || read.substr(i,1) == "k"){
            (*rep)[0][1] = (0.5 + -1.5)/2.0;
        }
        else if(read.substr(i,1) == "V" || read.substr(i,1) == "v"){
            (*rep)[0][1] = (1.5 + -0.5 + -0.5)/3.0;
        }
        else if(read.substr(i,1) == "H" || read.substr(i,1) == "h"){
            (*rep)[0][1] =  (1.5 + -0.5 + -1.5)/3.0;
        }
        else if(read.substr(i,1) == "D" || read.substr(i,1) == "d"){
            (*rep)[0][1] = (1.5 + 0.5 + -1.5)/3.0;
        }
        else if(read.substr(i,1) == "B" || read.substr(i,1) == "b"){
            (*rep)[0][1] = (-0.5 + 0.5 + -1.5)/3.0;
        }
        else {
            (*rep)[0][1] = (1.5 + -0.5 + 0.5 + -1.5)/4.0;
        }
    }
};

RealNumRep::~RealNumRep(){};
