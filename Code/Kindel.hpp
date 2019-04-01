//
//  Kindel.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef Kindel_hpp
#define Kindel_hpp

#include <stdio.h>
#include <string>
#include <map>
#include <vector>
#include <ctype.h>

#include "SideDetails.hpp"
class Kindel{
private:
    std::string cName;
    std::string contiq;
    int cSize;
    //int countsize;
    std::vector<SideDetails> contMap;
public:
    Kindel();
    
    Kindel(std::string id, int s);
    
    void generateKindle(std::string id, int s);
    
    void addRead(int x , std::string c, std::string r);
    //void addRead(int x , std::string r);
    
    // void processReads();
    
    void generateConsensus();
    
    void generateConsensus(bool b);
    
    void generateConsensus(std::vector<int> *dp);
    
    void generateConsensus(bool b, std::vector<int> *dp);
    
    std::string getConsensus();
    
    std::string getContigName();
    
    void decomposeCigar(std::string c, std::map<size_t, std::string, std::less<size_t>> *cdec);
    
    void deleteRefMap();
    
    ~Kindel();
};
#endif /* Kindel_hpp */
