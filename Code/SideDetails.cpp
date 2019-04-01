//
//  SideDetails.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "SideDetails.hpp"

SideDetails::SideDetails(){
    NucDetails["A"] = 0.0;
    NucDetails["C"] = 0.0;
    NucDetails["G"] = 0.0;
    NucDetails["T"] = 0.0;
    // NucDetails["R"] = 0.0;
    // NucDetails["Y"] = 0.0;
    // NucDetails["K"] = 0.0;
    // NucDetails["M"] = 0.0;
    // NucDetails["S"] = 0.0;
    // NucDetails["W"] = 0.0;
    // NucDetails["B"] = 0.0;
    // NucDetails["D"] = 0.0;
    // NucDetails["H"] = 0.0;
    // NucDetails["V"] = 0.0;
    NucDetails["N"] = 0.0;
    NucDetails["Del"] = 0.0;
    NucDetails["Indel"] = 0.0;
    NucDetails["SC"] = 0.0;
    NucDetails["HC"] = 0.0;
};

void SideDetails::addValue(std::string s){
    std::map<std::string, double>::iterator it = NucDetails.find(s);
    if(it != NucDetails.end()){
        //Exist
        NucDetails[s]++;
        c++;
    };
};

void SideDetails::addInsertion(std::string s){
    NucDetails["Indel"]++;
    if(Insertions.find(s) != Insertions.end()){
        Insertions[s]++;
    }
    else{
        Insertions[s] = 1.0;
    }
};

void SideDetails::addSoftClip(std::string s){
    NucDetails["SC"]++;
    if(SClip.find(s) != SClip.end()){
        SClip[s]++;
    }
    else{
        SClip[s] = 1.0;
    }
};

void SideDetails::addHardClip(std::string s){
    NucDetails["HC"]++;
    if(HClip.find(s) != HClip.end()){
        HClip[s]++;
    }
    else{
        HClip[s] = 1.0;
    }
};

void SideDetails::addDeletion(){
    NucDetails["Del"]++;
};

void SideDetails::generateConcensus(){
    int i = 0;
    generateConcensus(&i);
};

void SideDetails::generateConcensus(bool b){
    int i = 0;
    generateConcensus(b, &i);
};

void SideDetails::deleteVariable(){
    c = 0;
    ic = 0;
    std::map<std::string, double>().swap(NucDetails);
    std::map<std::string, double>().swap(Insertions);
    std::map<std::string, double>().swap(SClip);
    std::map<std::string, double>().swap(HClip);
    cons = "";
};

std::string SideDetails::getConsensus(){
    return cons;
};

SideDetails::~SideDetails(){};

void SideDetails::generateConcensus(int *dp){
    cons = "";
    double Acov = 0;
    double Ccov = 0;
    double Gcov = 0;
    double Tcov = 0;
    if(c > 0){
        Acov = NucDetails["A"]/c;
        Ccov = NucDetails["C"]/c;
        Gcov = NucDetails["G"]/c;
        Tcov = NucDetails["T"]/c;
    }
    if((c/2.0) < NucDetails["Indel"]){
        std::string indel="";
        double maxins = 0;
        for (std::map<std::string, double>::iterator it = Insertions.begin(); it != Insertions.end(); it++){
            if(maxins <= it->second){
                if(maxins == it->second){
                    indel += it->first;
                }
                else{
                    maxins = it->second;
                    indel = "";
                    indel = it->first;
                }
            }
        }
        cons = indel;
    }
    
    if(NucDetails["Del"] < c/2.0 ){
        if(Acov > Ccov  && Acov > Gcov && Acov > Tcov){
            cons+="A";
        }
        else if(Ccov > Acov  && Ccov > Gcov && Ccov > Tcov){
            cons+="C";
        }
        else if(Gcov > Acov  && Gcov > Ccov && Gcov > Tcov){
            cons+="G";
        }
        else if(Tcov > Acov  && Tcov > Ccov && Tcov > Gcov){
            cons+="T";
        }
        else{
            cons+="N";
        }
    }
    else{
        (*dp)++;
    }
};

void SideDetails::generateConcensus(bool b, int *dp){
    generateConcensus(dp);
    if(b == true){
        std::map<std::string, double>().swap(NucDetails);
        std::map<std::string, double>().swap(Insertions);
        std::map<std::string, double>().swap(SClip);
        std::map<std::string, double>().swap(HClip);
    }
};
