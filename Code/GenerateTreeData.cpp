//
//  GenerateTreeData.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "GenerateTreeData.hpp"
GenerateTreeData::GenerateTreeData(){
    
};

GenerateTreeData::GenerateTreeData(std::vector<NGS> *ShortReads, std::vector<NGS> *ReferenveListCom, int k, std::string tranmeth, int l, int sa){
    DataTransformations DT;
    DT.setCutPoints(sa);
    std::map<std::string, int> dubReads;
    int c = 0;
    
    for (int i = 0; i < (int)ShortReads->size(); i++){
        std::string read = ShortReads->at(i).returnRead();
        int x = (int) read.size();
        if(x > k){
            for (int j = 0; j < (x - k) + 1; j++){
                std::map<std::string, int>::iterator it = dubReads.find(read.substr(j,k));
                if(it != dubReads.end()){
                    ReferenveListCom->at(it->second).addCoordinates(i, j);
                    ReferenveListCom->at(it->second).addGnp(ShortReads->at(i).returnUIP());
                }
                else{
                    NGS n(ShortReads->at(i).returnUIP(), ShortReads->at(i).returnSRidP(), i, j, k);
                    DT.createTransformation(tranmeth, ShortReads->at(i).returnRep(), n.returnTra(), j, k,  l);
                    ReferenveListCom->push_back(n);
                    dubReads[read.substr(j, k)] = c;
                    c++;
                }
            }
        }
        else{
            std::map<std::string, int>::iterator it = dubReads.find(read);
            if(it != dubReads.end()){
                ReferenveListCom->at(it->second).addCoordinates(i, 0);
                ReferenveListCom->at(it->second).addGnp(ShortReads->at(i).returnUIP());
            }
            else{
                NGS n(ShortReads->at(i).returnUIP(), ShortReads->at(i).returnSRidP(), i, 0, (int)read.size());
                DT.createTransformation(tranmeth, ShortReads->at(i).returnRep(), n.returnTra(), 0, k,  l);
                ReferenveListCom->push_back(n);
                //dubReads[read.substr(i, 0)] = c;
                dubReads[read] = c;
                c++;
            }
        }
        ShortReads->at(i).clearRep();
    }
    
    std::map<std::string, int>().swap(dubReads);
    DT.deleteCutPonts();
    DT.~DataTransformations();
};


GenerateTreeData::GenerateTreeData(std::vector<NGS*> *ShortReads, std::vector<NGS> *ReferenveListCom, int k, std::string tranmeth, int l, int sa){
    DataTransformations DT;
    DT.setCutPoints(sa);
    std::map<std::string, int> dubReads;
    int c = 0;
    
    for (int i = 0; i < (int)ShortReads->size(); i++){
        std::string read = (*(ShortReads->at(i))).returnRead();
        int x = (int) read.size();
        if(x > k){
            for (int j = 0; j < (x - k) + 1; j++){
                std::map<std::string, int>::iterator it = dubReads.find(read.substr(j,k));
                if(it != dubReads.end()){
                    ReferenveListCom->at(it->second).addCoordinates(i, j);
                    ReferenveListCom->at(it->second).addGnp((*(ShortReads->at(i))).returnUIP());
                }
                else{
                    NGS n((*(ShortReads->at(i))).returnUIP(), (*(ShortReads->at(i))).returnSRidP(), i, j, k);
                    DT.createTransformation(tranmeth, (*(ShortReads->at(i))).returnRep(), n.returnTra(), j, k,  l);
                    ReferenveListCom->push_back(n);
                    dubReads[read.substr(j, k)] = c;
                    c++;
                }
            }
        }
        else{
            std::map<std::string, int>::iterator it = dubReads.find(read);
            if(it != dubReads.end()){
                ReferenveListCom->at(it->second).addCoordinates(i, 0);
                ReferenveListCom->at(it->second).addGnp((*(ShortReads->at(i))).returnUIP());
            }
            else{
                NGS n((*(ShortReads->at(i))).returnUIP(), (*(ShortReads->at(i))).returnSRidP(), i, 0, (int)read.size());
                DT.createTransformation(tranmeth, (*(ShortReads->at(i))).returnRep(), n.returnTra(), 0, k,  l);
                ReferenveListCom->push_back(n);
                dubReads[read] = c;
                c++;
            }
        }
       (*(ShortReads->at(i))).clearRep();
    }
    
    std::map<std::string, int>().swap(dubReads);
    DT.deleteCutPonts();
    DT.~DataTransformations();
};

GenerateTreeData::~GenerateTreeData(){};


