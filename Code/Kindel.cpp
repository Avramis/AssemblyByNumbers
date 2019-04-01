//
//  Kindel.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "Kindel.hpp"
Kindel::Kindel(){};

Kindel::Kindel(std::string id, int s){
    generateKindle(id, s);
};

void Kindel::generateKindle(std::string id, int s){
    cName = id;
    std::vector<SideDetails>().swap(contMap);
    std::vector<SideDetails>(s+200, SideDetails()).swap(contMap);
};

void Kindel::addRead(int x, std::string c, std::string r){
    std::map<size_t, std::string, std::less<size_t>> cmap;
    decomposeCigar(c, &cmap);
    int csp = -1;
    int rsp = 0;
    int mps = 0;
    
    for(std::map<size_t, std::string, std::less<size_t>>::iterator cit = cmap.begin(); cit != cmap.end(); cit++){
        
        int len  = stoi(c.substr(csp+1, (cit->first - csp)-1));
        //std::cout << c.substr(csp+1, (cit->first - csp)-1) << "\n";
        if(cit->second == "M"){
            for (int i = 0; i < len; i++){
                contMap[x+mps].addValue(r.substr(rsp , 1));
                rsp++;
                mps++;
            }
        }
        else if(cit->second == "I"){
            contMap[x+mps].addInsertion(r.substr(rsp+1, len));
            rsp += len;
        }
        else if(cit->second == "D"){
            
            for (int i = 0; i < len; i++){
                
                contMap[x+mps].addDeletion();
                mps++;
            }
        }
        else if(cit->second == "S"){
            contMap[x+mps].addSoftClip(r.substr(rsp+1, len));
            rsp += len;
        }
        else if(cit->second == "H"){
            contMap[x+mps].addHardClip(r.substr(rsp+1, len));
            rsp += len;
        }
        csp = (int)cit->first;
    }
    std::map<size_t, std::string, std::less<size_t>>().swap( cmap);
};

void Kindel::generateConsensus(){
    for (int i = 0; i < (int)contMap.size(); i++){
        contMap[i].generateConcensus();
        contiq+=contMap[i].getConsensus();
    }
    cSize = (int)contiq.size();
};


void Kindel::generateConsensus(bool b){
    
    for (int i = 0; i < (int)contMap.size(); i++){
        contMap[i].generateConcensus(b);
        contiq+=contMap[i].getConsensus();
    }
    cSize = (int)contiq.size();
};

void Kindel::generateConsensus(std::vector<int> *dp){
    std::vector<int>((int)contMap.size(), 0).swap((*dp));
    contMap[0].generateConcensus(&(*dp)[0]);
    contiq+=contMap[0].getConsensus();
    
    for (int i = 1; i < (int)contMap.size(); i++){
        contMap[i].generateConcensus(false, &(*dp)[i]);
        contiq+=contMap[i].getConsensus();
        (*dp)[i] += (*dp)[i-1];
    }
    cSize = (int)contiq.size();
};

void Kindel::generateConsensus(bool b, std::vector<int> *dp){
    std::vector<int>((int)contMap.size(), 0).swap((*dp));
    contMap[0].generateConcensus(b, &(*dp)[0]);
    contiq+=contMap[0].getConsensus();
    
    for (int i = 1; i < (int)contMap.size(); i++){
        (*dp)[i] += (*dp)[i-1];
        contMap[i].generateConcensus(b, &(*dp)[i]);
        contiq+=contMap[i].getConsensus();
        
    }
    cSize = (int)contiq.size();
};

std::string Kindel::getConsensus(){
    return contiq;
};


std::string Kindel::getContigName(){
    return cName;
};

void Kindel::decomposeCigar(std::string c, std::map<size_t, std::string, std::less<size_t>> *cdec){
    size_t pos =c.find("M",0);
    while(pos != std::string::npos){
        (*cdec)[pos] = "M";
        pos = c.find("M",pos+1);
    };
    
    pos =c.find("I",0);
    while(pos != std::string::npos){
        (*cdec)[pos] = "I";
        pos = c.find("I",pos+1);
    };
    
    pos =c.find("D",0);
    while(pos != std::string::npos){
        (*cdec)[pos] = "D";
        pos = c.find("D",pos+1);
    };
    
    pos =c.find("S",0);
    while(pos != std::string::npos){
        (*cdec)[pos] = "S";
        pos = c.find("S",pos+1);
    };
    
    pos =c.find("H",0);
    while(pos != std::string::npos){
        (*cdec)[pos] = "H";
        pos = c.find("H",pos+1);
    };
};

void Kindel::deleteRefMap(){
    for (int i = 0; i < (int)contMap.size(); i++){
        contMap[i].deleteVariable();
    };
    std::vector<SideDetails>().swap(contMap);
    cName = "";
    contiq = "";
};

Kindel::~Kindel(){};
