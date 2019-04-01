//
//  SideDetails.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef SideDetails_hpp
#define SideDetails_hpp

#include <stdio.h>
#include <stdio.h>
#include <map>
#include <string>
class SideDetails{
private:
    // c = coverage || ic = icoverage
    double c = 0, ic = 0;
    std::map<std::string, double> NucDetails;
    std::map<std::string, double> Insertions;
    std::map<std::string, double> SClip;
    std::map<std::string, double> HClip;
    std::string cons = "";
public:
    SideDetails();
    void addValue(std::string s);
    void addInsertion(std::string s);
    void addSoftClip(std::string s);
    void addHardClip(std::string s);
    void addDeletion();
    void generateConcensus();
    void generateConcensus(bool b);
    
    void generateConcensus(int *dp);
    void generateConcensus(bool b, int *dp);
    
    
    
    void deleteVariable();
    std::string getConsensus();
    ~SideDetails();
};
#endif /* SideDetails_hpp */
