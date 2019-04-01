//
//  IniAssemby.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "IniAssemby.hpp"

void IniAssemby::printLabels(int i, int j){
    if(i == 0){
        printLabels(10, 0);
        std::cout << "> Assembly initialisation.\n";
        printLabels(10, 0);
        
    }
    else if(i == 1){
        std::cout << "> Generate data bins.\n";
    }
    else if(i == 2){
        std::cout << "\n";
        printLabels(10, 0);
    }
    else if(i == 3){
        std::cout << "> Bin: " << j << "\n";
    }
    else if(i == 4){
        std::cout << "> Bin size: " << j << "\n";
    }
    else if(i == 5){
        if(j == 0){
            std::cout << "> Prepare data for VP-tree.\n";
        }
        else{
            std::cout << "> Prepare data for VP-trees.\n";
        }
    }
    else if(i == 6){
        std::cout << "======================================\n";
    }
    else if(i == 7){
        printLabels(10, 0);
        std::cout << "\n";
    }
    else if (i == 8){
        printLabels(10, 0);
        std::cout << " > Further assembly improvement attempt_" << j << "\n";
        printLabels(10, 0);
    }
    else if(i == 9){
        std::cout << "> Create contig reads.\n";
    }
    else if(i == 10){
        std::cout << "===========================================================\n";
    }
    else if(i == 11){
        if( j == 1){
             std::cout << " > " << j << " contig was identified.\n";
        }
        else{
            std::cout << " > " << j << " contigs were identified.\n";
        }
        std::cout << "> Start building fasta and sam files\n";
    }
    else if(i == 12){
        printLabels(10, 0);
        std::cout << "\n";
    }
    else if(i == 13){
        printLabels(10, 0);
        std::cout << "> Assembly completed.\n";
        printLabels(10, 0);
        std::cout << "\n";
    }
    else if(i == 14){
        std::cout << "===========================================================\n";
        std::cout << ">                No contigs were identified                \n";
        std::cout << ">                   Terminating assembly                   \n";
        std::cout << "===========================================================\n\n";
    }
    else if(i == 15){
        if(j == 0){
            std::cout << "> Importing reads from file.\n";
        }
        else{
            std::cout << "> Importing reads from file " << j <<".\n";
        }
    }
    else if(i == 16){
        if(j == 0){
            std::cout << "> Preparing mate pair 1 reads.\n";
        }
        else{
            std::cout << "> Preparing mate pair 2 reads.\n";
        }
    }
    
};

//void IniAssemby::importReads(std::string inputfile, std::string repmeth, int k, std::string tranmeth, int t, std::map<std::string, std::vector<NGS>> *NGSbins, int sa, int ss, bool bin)
void IniAssemby::importReads(std::string inputfile, std::string repmeth, int k, std::string tranmeth, int t, std::vector<NGS> *NuSequences, int sa, int ss){
    FastqParser (inputfile, repmeth, k, tranmeth, t, NuSequences, sa, ss);
    //FastqParser(std::string filepath, std::string repmeth, int refkmerlen, std::string tranmeth, int tranlev, std::vector<NGS> *NuSequences, int sa, int ss)
};

void IniAssemby::importReads(std::string p1inputfile, std::string p2inputfile, std::string repmeth, int k, std::string tranmeth, int t, std::vector<NGS> *NuSequences1, std::vector<NGS> *NuSequences2, int sa, int ss){
    FastqParser (p1inputfile, repmeth, k, tranmeth, t, NuSequences1, sa, ss);
    FastqParser (p2inputfile, repmeth, k, tranmeth, t, NuSequences2, sa, ss);
};

void IniAssemby::importReads(std::string inputfile, std::string repmeth, int k, std::string tranmeth, int t, std::vector<NGS> *NuSequences1, std::vector<NGS> *NuSequences2, int sa, int ss, bool pe){
    FastqParser(inputfile, repmeth, k, tranmeth, t, NuSequences1, NuSequences2, sa, ss, pe);
    
    //FastqParser(std::string filepath, std::string repmeth, int refkmerlen, std::string tranmeth, int tranlev, std::vector<NGS> *NuSequences, int sa, int ss)
};



void IniAssemby::ReferenceList(std::vector<NGS> *ShortReads, std::vector<NGS> *ReferenveListCom, int k, std::string tranmeth, int l, int sa){
    //std::cout << "> Prepare data for VP-tree.\n";
    GenerateTreeData GTD(ShortReads, ReferenveListCom, k, tranmeth, l, sa);
    GTD.~GenerateTreeData();
};

void IniAssemby::ReferenceList(std::vector<NGS*> *ShortReads, std::vector<NGS> *ReferenveListCom, int k, std::string tranmeth, int l, int sa){
    GenerateTreeData GTD(ShortReads, ReferenveListCom, k, tranmeth, l, sa);
    GTD.~GenerateTreeData();
};



void IniAssemby::ReferenceList(std::vector<NGS> *ShortReads1, std::vector<NGS> *ReferenveListCom1, std::vector<NGS> *ShortReads2, std::vector<NGS> *ReferenveListCom2, int k, std::string tranmeth, int l, int sa){
    //std::cout << "> Prepare data for VP-tree.\n";
    GenerateTreeData GTD1(ShortReads1, ReferenveListCom1, k, tranmeth, l, sa);
    GTD1.~GenerateTreeData();
    GenerateTreeData GTD2(ShortReads2, ReferenveListCom2, k, tranmeth, l, sa);
    GTD2.~GenerateTreeData();
};




void IniAssemby::processReads(std::vector<NGS> *ShortReads, std::vector<NGS> *ReferenveList, std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::vector<int>> *gVaridx, int knn, std::string rmeth, int k, double als){
    ProcessReads Pr(ShortReads, ReferenveList, Contigs, gVaridx, knn, rmeth, k, als);
    Pr.~ProcessReads();
};

void IniAssemby::processReads(std::vector<NGS*> *ShortReads, std::vector<NGS> *ReferenveList, std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::vector<int>> *gVaridx, int knn, std::string rmeth, int k, double als){
    ProcessReads Pr(ShortReads, ReferenveList, Contigs, gVaridx, knn, rmeth, k, als);
    Pr.~ProcessReads();
};


void IniAssemby::processReads(std::vector<NGS> *ShortReads1, std::vector<NGS> *ReferenveList1, std::vector<std::pair<std::string, std::string>> *Contigs1, std::vector<std::vector<int>> *gVaridx1, std::vector<NGS> *ShortReads2, std::vector<NGS> *ReferenveList2, std::vector<std::pair<std::string, std::string>> *Contigs2, std::vector<std::vector<int>> *gVaridx2, int knn, std::string rmeth, int k, double als, bool missmatch){
    ProcessReads Pr(ShortReads1, ReferenveList1, Contigs1, gVaridx1, ShortReads2, ReferenveList2, Contigs2, gVaridx2, knn, rmeth, k, als, missmatch);
    Pr.~ProcessReads();
    
};

void IniAssemby::generateCReads(std::vector<std::pair<std::string, std::string>> *Contigs, std::string repmeth, int k, std::string tranmeth, int t, std::vector<NGS> *NuSequences){
    GenerateCReads GCR(Contigs, repmeth, k, tranmeth, t,NuSequences);
};


void IniAssemby::generateFiles(std::vector<NGS> *ShortReads, std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::vector<int>> *varIdx, std::string outputsam, std::string outputcontigs){
    FileGeneration FG(ShortReads, Contigs, varIdx, outputsam, outputcontigs);
};

void IniAssemby::generateFiles(std::vector<NGS> *ShortReads1, std::vector<std::pair<std::string, std::string>> *Contigs1, std::vector<std::vector<int>> *varIdx1, std::vector<NGS> *ShortReads2, std::vector<std::pair<std::string, std::string>> *Contigs2, std::vector<std::vector<int>> *varIdx2, std::string outputsam, std::string outputcontigs){
    FileGeneration FG(ShortReads1, Contigs1, varIdx1, ShortReads2, Contigs2, varIdx2, outputsam, outputcontigs);
};

void IniAssemby::removedunmatched(std::vector<std::pair<std::string, std::string>> *tContigs, std::vector<std::vector<int>> *gVaridx, std::vector<std::vector<int>> *fgVaridx, std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<NGS> *NuSequences){
    int cidx = ((int)(*tContigs).size()) - 1;
    while(!(*gVaridx).empty()){
        if((int)(*gVaridx)[cidx].size() > 1){
            (*fgVaridx).push_back((*gVaridx)[cidx]);
            (*Contigs).push_back((*tContigs)[cidx]);
        }
        else{
            (*NuSequences)[(*gVaridx)[cidx][0]].clearAldetails();
        }
        (*gVaridx).erase((*gVaridx).begin()+cidx);
        (*tContigs).erase((*tContigs).begin()+cidx);
        cidx--;
    }
    std::vector<std::vector<int>>().swap((*gVaridx));
    std::vector<std::pair<std::string, std::string>>().swap((*tContigs));
};

void IniAssemby::removedunmatched(std::vector<NGS> *NuSequences, std::vector<std::vector<int>> *tgVaridx, std::vector<std::vector<int>> *gVaridx, std::vector<std::vector<int>> *fgVaridx, std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::pair<std::string, std::string>> *tContigs, std::vector<NGS> *ContigReads){
    std::vector<std::vector<int>>().swap((*gVaridx));
    (*gVaridx).swap((*fgVaridx));
    std::vector<std::vector<int>>().swap((*fgVaridx));
    int ts = (int)(*tgVaridx).size();
    ts--;
    int fvadd = 0;
    while (!(*tgVaridx).empty()){
        if((int)(*tgVaridx)[ts].size() == 1 &&  (int)(*gVaridx)[(*tgVaridx)[ts][0]].size() == 1){
            (*NuSequences)[(*gVaridx)[(*tgVaridx)[ts][0]][0]].clearAldetails();
        }
        else{
            (*fgVaridx).push_back((*gVaridx)[(*tgVaridx)[ts][0]]);
            (*Contigs).push_back((*tContigs)[ts]);
            for (int ti = 1; ti < (int)(*tgVaridx)[ts].size(); ti++){
                int ccoor = (*ContigReads)[(*tgVaridx)[ts][ti]].returnCoordinates()[0].second;
                (*fgVaridx)[fvadd].insert((*fgVaridx)[fvadd].end(), (*gVaridx)[(*tgVaridx)[ts][ti]].begin(),  (*gVaridx)[(*tgVaridx)[ts][ti]].end());
                for (int k = 0; k < (int)(*gVaridx)[(*tgVaridx)[ts][ti]].size(); k++){
                    int ncoor = (*NuSequences)[(*gVaridx)[(*tgVaridx)[ts][ti]][k]].returnCoordinates()[0].second;
                    (*NuSequences)[(*gVaridx)[(*tgVaridx)[ts][ti]][k]].editCoordinates(0,ccoor+ncoor);
                }
            }
            fvadd++;
        }
        (*tContigs).erase((*tContigs).begin()+ts);
        (*tgVaridx).erase((*tgVaridx).begin()+ts);
        ts--;
    }
    std::vector<std::pair<std::string, std::string>>().swap((*tContigs));
    std::vector<std::vector<int>>().swap((*tgVaridx));
};

void IniAssemby::readsBinning(std::vector<NGS> *NuSequences, std::map<std::string, std::vector<NGS*>> *NGSbins, int sa, int ss){
    NucRepresentations rep;
    DataTransformations tra;
    tra.setCutPoints(sa);
    DataProcessing DaPr;
    printLabels(1, 0);
    std::vector<std::vector<double>> saxrep ,saxtra;
    std::vector<std::vector<double>>().swap(saxrep);
    std::vector<std::vector<double>>().swap(saxtra);
    
    for(int ni = 0; ni < (int)(*NuSequences).size(); ni++){
        std::string binkey ="";
        saxrep = *((*NuSequences)[ni].returnRep());
        DaPr.repAccumulation(&saxrep);
        tra.createTransformation("SAX", &saxrep, &saxtra ,ss);
        for (int i = 0; i < (int)saxtra.size(); i++){
            for (int j = 0; j < (int)saxtra.at(0).size(); j++){
                binkey += std::to_string((int)saxtra[i][j]);
            };
        };
        (*NGSbins)[binkey].push_back(&(*NuSequences)[ni]);
        std::vector<std::vector<double>>().swap(saxrep);
        std::vector<std::vector<double>>().swap(saxtra);
    }
};

void IniAssemby::readsBinning(std::vector<NGS> *NuSequences1, std::vector<NGS> *NuSequences2, std::map<std::string, std::vector<NGS*>> *NGSbins1, std::map<std::string, std::vector<NGS*>> *NGSbins2, int sa, int ss){
    NucRepresentations rep;
    DataTransformations tra;
    tra.setCutPoints(sa);
    DataProcessing DaPr;
    printLabels(1, 0);
    std::vector<std::vector<double>> saxrep1 ,saxtra1;
    std::vector<std::vector<double>>().swap(saxrep1);
    std::vector<std::vector<double>>().swap(saxtra1);
    
    std::vector<std::vector<double>> saxrep2 ,saxtra2;
    std::vector<std::vector<double>>().swap(saxrep2);
    std::vector<std::vector<double>>().swap(saxtra2);
    
    for(int ni = 0; ni < (int)(*NuSequences1).size(); ni++){
        std::string binkey1 ="";
        std::string binkey2 ="";
        
        saxrep1 = *((*NuSequences1)[ni].returnRep());
        saxrep2 = *((*NuSequences2)[ni].returnRep());
        
        DaPr.repAccumulation(&saxrep1);
        DaPr.repAccumulation(&saxrep2);
        
        tra.createTransformation("SAX", &saxrep1, &saxtra1 ,ss);
        tra.createTransformation("SAX", &saxrep2, &saxtra2 ,ss);
        
        for (int i = 0; i < (int)saxtra1.size(); i++){
            for (int j = 0; j < (int)saxtra1.at(0).size(); j++){
                binkey1 += std::to_string((int)saxtra1[i][j]);
                binkey2 += std::to_string((int)saxtra2[i][j]);
            };
        };
        
        (*NGSbins1)[binkey1].push_back(&(*NuSequences1)[ni]);
        (*NGSbins2)[binkey2].push_back(&(*NuSequences2)[ni]);
        
        std::vector<std::vector<double>>().swap(saxrep1);
        std::vector<std::vector<double>>().swap(saxtra1);
        
        std::vector<std::vector<double>>().swap(saxrep2);
        std::vector<std::vector<double>>().swap(saxtra2);
    }
    
}

void IniAssemby::processBins(std::map<std::string, std::vector<NGS*>> *NGSbins, std::vector<NGS> *ReferenveListCom, std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::vector<int>> *fgVaridx, int k, std::string repmeth, std::string tranmeth, int t, int sa, int knn, double als){
    int bi = 0;
    while(!(*NGSbins).empty()){
        printLabels(2, 0);
        printLabels(3, bi);
        bi++;
        
        std::map<std::string, std::vector<NGS*>>::iterator it  = (*NGSbins).begin();
        std::vector<NGS> tempReferenveListCom;
        std::vector<std::pair<std::string, std::string>> tempContigs;
        printLabels(4, (int)it->second.size());
        printLabels(5, 0);
        ReferenceList(&(*NGSbins)[it->first], ReferenveListCom, k, tranmeth, t, sa);
        std::vector<std::vector<int>> tempgVaridx, tempfgVaridx;
        
        processReads(&(*NGSbins)[it->first], ReferenveListCom, &tempContigs, &tempfgVaridx, knn, repmeth, k, als);
        //int ii = (int)tempfgVaridx.size();
        (*NGSbins).erase(it);
        (*Contigs).insert((*Contigs).end(), tempContigs.begin(), tempContigs.end());
        (*fgVaridx).insert((*fgVaridx).end(), tempfgVaridx.begin(), tempfgVaridx.end());
        std::vector<std::pair<std::string, std::string>>().swap(tempContigs);
        std::vector<std::vector<int>>().swap(tempfgVaridx);
        std::vector<NGS>().swap((*ReferenveListCom));
        printLabels(6, 0);
    }
};



void IniAssemby:: processBins(std::map<std::string, std::vector<NGS*>> *NGSbins1, std::vector<NGS> *ReferenveListCom1, std::vector<std::pair<std::string, std::string>> *Contigs1, std::vector<std::vector<int>> *fgVaridx1, std::map<std::string, std::vector<NGS*>> *NGSbins2, std::vector<NGS> *ReferenveListCom2, std::vector<std::pair<std::string, std::string>> *Contigs2, std::vector<std::vector<int>> *fgVaridx2, int k, std::string repmeth, std::string tranmeth, int t, int sa, int knn, double als){
        int bi = 0;
    if((*NGSbins1).size() > (*NGSbins1).size()){
        while(!(*NGSbins1).empty()){
            printLabels(2, 0);
            printLabels(3, bi);
            bi++;
            std::map<std::string, std::vector<NGS*>>::iterator it1  = (*NGSbins1).begin();
            std::vector<NGS> tempReferenveListCom1;
            std::vector<std::pair<std::string, std::string>> tempContigs1;
            printLabels(4, (int)it1->second.size());
            printLabels(5, 0);
            ReferenceList(&(*NGSbins1)[it1->first], ReferenveListCom1, k, tranmeth, t, sa);
            std::vector<std::vector<int>> tempgVaridx1, tempfgVaridx1;
            processReads(&(*NGSbins1)[it1->first], ReferenveListCom1, &tempContigs1, &tempfgVaridx1, knn, repmeth, k, als);
            //int ii = (int)tempfgVaridx.size();
            (*NGSbins1).erase(it1);
            (*Contigs1).insert((*Contigs1).end(), tempContigs1.begin(), tempContigs1.end());
            (*fgVaridx1).insert((*fgVaridx1).end(), tempfgVaridx1.begin(), tempfgVaridx1.end());
            std::vector<std::pair<std::string, std::string>>().swap(tempContigs1);
            std::vector<std::vector<int>>().swap(tempfgVaridx1);
            std::vector<NGS>().swap((*ReferenveListCom1));
            printLabels(6, 0);
            
            
            
            if(!(*NGSbins2).empty()){
                std::map<std::string, std::vector<NGS*>>::iterator it2  = (*NGSbins2).begin();
                std::vector<NGS> tempReferenveListCom2;
                std::vector<std::pair<std::string, std::string>> tempContigs2;
                ReferenceList(&(*NGSbins2)[it2->first], ReferenveListCom2, k, tranmeth, t, sa);
                std::vector<std::vector<int>> tempgVaridx2, tempfgVaridx2;
                processReads(&(*NGSbins2)[it2->first], ReferenveListCom2, &tempContigs2, &tempfgVaridx2, knn, repmeth, k, als);
                (*NGSbins2).erase(it2);
                (*Contigs2).insert((*Contigs2).end(), tempContigs2.begin(), tempContigs2.end());
                (*fgVaridx2).insert((*fgVaridx2).end(), tempfgVaridx2.begin(), tempfgVaridx2.end());
                std::vector<std::pair<std::string, std::string>>().swap(tempContigs2);
                std::vector<std::vector<int>>().swap(tempfgVaridx2);
                std::vector<NGS>().swap((*ReferenveListCom2));
            }
        }
    }
    else{
        while(!(*NGSbins2).empty()){
            printLabels(2, 0);
            printLabels(3, bi);
            bi++;
            std::map<std::string, std::vector<NGS*>>::iterator it2  = (*NGSbins2).begin();
            std::vector<NGS> tempReferenveListCom2;
            std::vector<std::pair<std::string, std::string>> tempContigs2;
            printLabels(4, (int)it2->second.size());
            printLabels(5, 0);
            ReferenceList(&(*NGSbins2)[it2->first], ReferenveListCom2, k, tranmeth, t, sa);
            std::vector<std::vector<int>> tempgVaridx2, tempfgVaridx2;
            processReads(&(*NGSbins2)[it2->first], ReferenveListCom2, &tempContigs2, &tempfgVaridx2, knn, repmeth, k, als);
            (*NGSbins2).erase(it2);
            (*Contigs2).insert((*Contigs2).end(), tempContigs2.begin(), tempContigs2.end());
            (*fgVaridx2).insert((*fgVaridx2).end(), tempfgVaridx2.begin(), tempfgVaridx2.end());
            std::vector<std::pair<std::string, std::string>>().swap(tempContigs2);
            std::vector<std::vector<int>>().swap(tempfgVaridx2);
            std::vector<NGS>().swap((*ReferenveListCom2));
            printLabels(6, 0);
            
            if(!(*NGSbins1).empty()){
                std::map<std::string, std::vector<NGS*>>::iterator it1  = (*NGSbins1).begin();
                std::vector<NGS> tempReferenveListCom1;
                std::vector<std::pair<std::string, std::string>> tempContigs1;
                ReferenceList(&(*NGSbins1)[it1->first], ReferenveListCom1, k, tranmeth, t, sa);
                std::vector<std::vector<int>> tempgVaridx1, tempfgVaridx1;
                processReads(&(*NGSbins1)[it1->first], ReferenveListCom1, &tempContigs1, &tempfgVaridx1, knn, repmeth, k, als);
                (*NGSbins1).erase(it1);
                (*Contigs1).insert((*Contigs1).end(), tempContigs1.begin(), tempContigs1.end());
                (*fgVaridx1).insert((*fgVaridx1).end(), tempfgVaridx1.begin(), tempfgVaridx1.end());
                std::vector<std::pair<std::string, std::string>>().swap(tempContigs1);
                std::vector<std::vector<int>>().swap(tempfgVaridx1);
                std::vector<NGS>().swap((*ReferenveListCom1));
            }
        }
    }
};

void IniAssemby::joinContigs(std::vector<NGS> *NuSequences, std::vector<NGS> *ContigReads, std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::pair<std::string, std::string>> *tContigs, std::vector<std::vector<int>> *gVaridx, std::vector<std::vector<int>> *fgVaridx, std::vector<std::vector<int>> *tgVaridx){
    std::vector<std::pair<std::string, std::string>>().swap((*Contigs));// tContigs;
    (*tContigs).swap((*Contigs));
    std::vector<std::pair<std::string, std::string>>().swap((*tContigs));// tContigs;
    std::vector<std::vector<int>>().swap((*gVaridx));
    (*gVaridx).swap((*fgVaridx));
    std::vector<std::vector<int>>().swap((*fgVaridx));
    int ts = (int)(*tgVaridx).size();
    std::vector<std::vector<int>>(ts, std::vector<int>()).swap((*fgVaridx));
    ts--;
    int cv = 0;
    while (!(*tgVaridx).empty()){
        int cc = 0;
        for (int ti = 0; ti < (int)(*tgVaridx)[ts].size(); ti++){
            cc += (int)(*gVaridx)[(*tgVaridx)[ts][ti]].size();
        }
        std::vector<int>(cc, 0).swap((*fgVaridx)[ts]);
        int kk = 0;
        for (int k = 0; k < (int)(*gVaridx)[(*tgVaridx)[ts][0]].size(); k++){
            (*fgVaridx)[ts][kk] = (*gVaridx)[(*tgVaridx)[ts][0]][k];
            kk++;
        }
        for (int ti = 1; ti < (int)(*tgVaridx)[ts].size(); ti++){
            int ccoor = (*ContigReads)[(*tgVaridx)[ts][ti]].returnCoordinates()[0].second;
            for (int k = 0; k < (int)(*gVaridx)[(*tgVaridx)[ts][ti]].size(); k++){
                int ncoor = (*NuSequences)[(*gVaridx)[(*tgVaridx)[ts][ti]][k]].returnCoordinates()[0].second;
                (*NuSequences)[(*gVaridx)[(*tgVaridx)[ts][ti]][k]].editCoordinates(0,ccoor+ncoor);
                (*fgVaridx)[ts][kk] = (*gVaridx)[(*tgVaridx)[ts][ti]][k];
                kk++;
            }
        }
        (*tgVaridx).erase((*tgVaridx).end()-1);
        cv++;
        ts--;
    }
    
};

IniAssemby::IniAssemby(std::string inputfile, std::string outputsam, std::string outputcontigs, std::string repmeth, int k, std::string tranmeth, int t, int ss, int sa, int binmet, bool bin, int knn, double als, bool isen, double  siv){
    std::vector<NGS> NuSequences;
    std::map<std::string, std::vector<NGS*>> NGSbins;
    bool discardunmatch = true;
    printLabels(0, 0);
    printLabels(15, 0);
    importReads(inputfile, repmeth, k, tranmeth, t, &NuSequences, sa, ss);
    std::vector<NGS> ReferenveListCom;
    std::vector<std::pair<std::string, std::string>> Contigs;
    std::vector<std::vector<int>> gVaridx, fgVaridx;
    //for (std::map<std::string, std::vector<NGS>>::iterator it = NGSbins.begin(); it != NGSbins.end(); it++){
    if(bin == false){
        std::vector<std::pair<std::string, std::string>> tContigs;
        printLabels(5, 0);
        ReferenceList(&NuSequences, &ReferenveListCom, k, tranmeth, t, sa);
        processReads(&NuSequences, &ReferenveListCom, &tContigs, &gVaridx, knn, repmeth, k, als);
        removedunmatched(&tContigs, &gVaridx, &fgVaridx, &Contigs, &NuSequences);
        discardunmatch = false;
    }
    else{
        readsBinning(&NuSequences, &NGSbins, sa, ss);
        processBins(&NGSbins, &ReferenveListCom, &Contigs, &fgVaridx, k, repmeth, tranmeth, t,  sa, knn, als);
    }
    std::map<std::string, std::vector<NGS*>>().swap(NGSbins);
    std::vector<NGS>().swap(ReferenveListCom);
    printLabels(7, 0);
    if(Contigs.size() < NuSequences.size()){
        std::vector<std::pair<std::string, std::string>> tempContigs;
        int a  = 1;
        while(Contigs.size() != tempContigs.size()){
            printLabels(8, a);
            a++;
            std::vector<NGS> ContigReads;
            printLabels(9, 0);
            std::vector<std::pair<std::string, std::string>>().swap(tempContigs);
            tempContigs.swap(Contigs);
            std::vector<std::pair<std::string, std::string>>().swap(Contigs);
            generateCReads(&tempContigs, repmeth, k, tranmeth, t, &ContigReads);
            std::vector<NGS>().swap(ReferenveListCom);
            ReferenceList(&ContigReads, &ReferenveListCom, k, tranmeth, t, sa);
            std::vector<std::vector<int>> tgVaridx;
            if(isen == true && als < 2.00){
                als +=siv;
                if(als > 2.00){
                    als = 2.00;
                }
            }
            std::vector<std::pair<std::string, std::string>> tContigs;
            processReads(&ContigReads, &ReferenveListCom, &tContigs, &tgVaridx, knn, repmeth, k, als);
            std::vector<NGS>().swap(ReferenveListCom);
            if(fgVaridx.size() != tgVaridx.size()){
                if( discardunmatch == true){
                    removedunmatched(&NuSequences, &tgVaridx, &gVaridx, &fgVaridx, &Contigs, &tContigs, &ContigReads);
                    discardunmatch = false;
                }
                else{
                    joinContigs(&NuSequences, &ContigReads, &Contigs, &tContigs, &gVaridx, &fgVaridx, &tgVaridx);
                }
                std::vector<std::vector<int>>().swap(tgVaridx);
                std::vector<std::vector<int>>().swap(gVaridx);
                
            }
            else{
                tContigs.swap(Contigs);
                std::vector<std::pair<std::string, std::string>>().swap( tContigs);
            }
            printLabels(7, 0);
        }
        printLabels(10, 0);
        printLabels(11, (int)Contigs.size());
        printLabels(12, 0);
        generateFiles(&NuSequences, &Contigs, &fgVaridx, outputsam, outputcontigs);
        printLabels(13, 0);
    }
    else{
        printLabels(14, 0);
    }
};

IniAssemby::IniAssemby(std::string inputfile1, std::string inputfile2, std::string outputsam, std::string outputcontigs, std::string repmeth, int k, std::string tranmeth, int t, int ss, int sa, int binmet, bool bin, int knn, double als, bool isen, double siv, bool discardunmatch){
    std::vector<NGS> NuSequences1, NuSequences2;
    //std::map<std::string, std::vector<NGS*>> NGSbins1, NGSbins2;
    printLabels(0, 0);
    printLabels(15, 1);
    importReads(inputfile1, repmeth, k, tranmeth, t, &NuSequences1, sa, ss);
    printLabels(15, 2);
    importReads(inputfile2, repmeth, k, tranmeth, t, &NuSequences2, sa, ss);
    processPEReads (&NuSequences1, &NuSequences2, outputsam, outputcontigs, repmeth, k, tranmeth,t, ss, sa, binmet, bin, knn, als, isen, siv, discardunmatch);
};

IniAssemby::IniAssemby(std::string inputfile, std::string outputsam, std::string outputcontigs, std::string repmeth, int k, std::string tranmeth, int t, int ss, int sa, int binmet, bool bin, int knn, double als, bool isen, double siv, bool pe, bool dealmissmatch){
    std::vector<NGS> NuSequences1, NuSequences2;
    printLabels(0, 0);
    printLabels(15, 0);
    importReads(inputfile, repmeth, k, tranmeth, t, &NuSequences1, &NuSequences2, sa, ss, pe);
    processPEReads (&NuSequences1, &NuSequences2, outputsam, outputcontigs, repmeth, k, tranmeth,t, ss, sa, binmet, bin, knn, als, isen, siv, dealmissmatch);
};



void IniAssemby::processPEReads (std::vector<NGS> *NuSequences1, std::vector<NGS> *NuSequences2, std::string outputsam, std::string outputcontigs, std::string repmeth, int k, std::string tranmeth, int t, int ss, int sa, int binmet, bool bin, int knn, double als, bool isen, double siv, bool dealmissmatch){
    std::map<std::string, std::vector<NGS*>> NGSbins1, NGSbins2;
    bool discardunmatch = true;
    std::vector<NGS> ReferenveListCom1;
    std::vector<NGS> ReferenveListCom2;
    
    std::vector<std::pair<std::string, std::string>> Contigs1;
    std::vector<std::vector<int>> gVaridx1, fgVaridx1;
    std::vector<std::pair<std::string, std::string>> Contigs2;
    std::vector<std::vector<int>> gVaridx2, fgVaridx2;
    
    if(bin == false){
        std::vector<std::pair<std::string, std::string>> tContigs1;
        std::vector<std::pair<std::string, std::string>> tContigs2;
        
        printLabels(5, 1);
        printLabels(16, 0);
        ReferenceList(NuSequences1, &ReferenveListCom1, k, tranmeth, t, sa);
        printLabels(16, 1);
        ReferenceList(NuSequences2, &ReferenveListCom2, k, tranmeth, t, sa);
        processReads(NuSequences1, &ReferenveListCom1, &tContigs1, &gVaridx1, NuSequences2, &ReferenveListCom2, &tContigs2, &gVaridx2, knn, repmeth, k, als, dealmissmatch);
        removedunmatched(&tContigs1, &gVaridx1, &fgVaridx1, &Contigs1, NuSequences1);
        removedunmatched(&tContigs2, &gVaridx2, &fgVaridx2, &Contigs2, NuSequences2);
        discardunmatch = false;
    }
    else{
        readsBinning(NuSequences1, NuSequences2,&NGSbins1, &NGSbins2, sa, ss);
        processBins(&NGSbins1, &ReferenveListCom1, &Contigs1, &fgVaridx1, &NGSbins2, &ReferenveListCom2, &Contigs2, &fgVaridx2, k, repmeth, tranmeth, t, sa, knn, als);
    }
    std::map<std::string, std::vector<NGS*>>().swap(NGSbins1);
    std::map<std::string, std::vector<NGS*>>().swap(NGSbins2);
    std::vector<NGS>().swap(ReferenveListCom1);
    std::vector<NGS>().swap(ReferenveListCom2);
    printLabels(7, 0);
    
    
    bool c1 = Contigs1.size() < (*NuSequences1).size();
    bool c2 = Contigs2.size() < (*NuSequences2).size();
    if(c1 == true || c2 == true){
        int a  = 1;
        std::vector<std::pair<std::string, std::string>> tempContigs1, tempContigs2;
        while(Contigs1.size() != tempContigs1.size() || Contigs2.size() != tempContigs2.size()){
            printLabels(8, a);
            a++;
            if(Contigs1.size() != tempContigs1.size() && Contigs2.size() != tempContigs2.size()){
                std::vector<NGS> ContigReads1, ContigReads2;
                printLabels(9, 0);
                std::vector<std::pair<std::string, std::string>>().swap(tempContigs1);
                std::vector<std::pair<std::string, std::string>>().swap(tempContigs2);
                tempContigs1.swap(Contigs1);
                std::vector<std::pair<std::string, std::string>>().swap(Contigs1);
                generateCReads(&tempContigs1, repmeth, k, tranmeth, t, &ContigReads1);
                std::vector<NGS>().swap(ReferenveListCom1);
                ReferenceList(&ContigReads1, &ReferenveListCom1, k, tranmeth, t, sa);
                
                tempContigs2.swap(Contigs2);
                std::vector<std::pair<std::string, std::string>>().swap(Contigs2);
                generateCReads(&tempContigs2, repmeth, k, tranmeth, t, &ContigReads2);
                std::vector<NGS>().swap(ReferenveListCom2);
                ReferenceList(&ContigReads2, &ReferenveListCom2, k, tranmeth, t, sa);
                std::vector<std::vector<int>> tgVaridx1;
                std::vector<std::vector<int>> tgVaridx2;
                if(isen == true && als < 2.00){
                    als +=siv;
                    if(als > 2.00){
                        als = 2.00;
                    }
                }
                std::vector<std::pair<std::string, std::string>> tContigs1;
                std::vector<std::pair<std::string, std::string>> tContigs2;
                processReads(&ContigReads1, &ReferenveListCom1, &tContigs1, &tgVaridx1, &ContigReads2, &ReferenveListCom2, &tContigs2, &tgVaridx2, knn, repmeth, k, als, false);
                std::vector<NGS>().swap(ReferenveListCom1);
                std::vector<NGS>().swap(ReferenveListCom2);
                if( discardunmatch == true){
                    removedunmatched(NuSequences1, &tgVaridx1, &gVaridx1, &fgVaridx1, &Contigs1, &tContigs1, &ContigReads1);
                    removedunmatched(NuSequences2, &tgVaridx2, &gVaridx2, &fgVaridx2, &Contigs2, &tContigs2, &ContigReads2);
                    discardunmatch = false;
                }
                else{
                    joinContigs(NuSequences1, &ContigReads1, &Contigs1, &tContigs1, &gVaridx1, &fgVaridx1, &tgVaridx1);
                    joinContigs(NuSequences2, &ContigReads2, &Contigs2, &tContigs2, &gVaridx2, &fgVaridx2, &tgVaridx2);
                }
                std::vector<std::vector<int>>().swap(tgVaridx1);
                std::vector<std::vector<int>>().swap(gVaridx1);
                
                std::vector<std::vector<int>>().swap(tgVaridx2);
                std::vector<std::vector<int>>().swap(gVaridx2);
            }
            else if(Contigs1.size() != tempContigs1.size()){
                std::vector<NGS> ContigReads1;
                printLabels(9, 0);
                std::vector<std::pair<std::string, std::string>>().swap(tempContigs1);
                
                tempContigs1.swap(Contigs1);
                std::vector<std::pair<std::string, std::string>>().swap(Contigs1);
                generateCReads(&tempContigs1, repmeth, k, tranmeth, t, &ContigReads1);
                std::vector<NGS>().swap(ReferenveListCom1);
                ReferenceList(&ContigReads1, &ReferenveListCom1, k, tranmeth, t, sa);
                
                std::vector<std::vector<int>> tgVaridx1;
                if(isen == true && als < 2.00){
                    als +=siv;
                    if(als > 2.00){
                        als = 2.00;
                    }
                }
                std::vector<std::pair<std::string, std::string>> tContigs1;
                processReads(&ContigReads1, &ReferenveListCom1, &tContigs1, &tgVaridx1, knn, repmeth, k, als);
                std::vector<NGS>().swap(ReferenveListCom1);
                if( discardunmatch == true){
                    removedunmatched(NuSequences1, &tgVaridx1, &gVaridx1, &fgVaridx1, &Contigs1, &tContigs1, &ContigReads1);
                    discardunmatch = false;
                }
                else{
                    joinContigs(NuSequences1, &ContigReads1, &Contigs1, &tContigs1, &gVaridx1, &fgVaridx1, &tgVaridx1);
                }
                std::vector<std::vector<int>>().swap(tgVaridx1);
                std::vector<std::vector<int>>().swap(gVaridx1);

            }
            else if(Contigs2.size() != tempContigs2.size()){
                std::vector<NGS> ContigReads1, ContigReads2;
                printLabels(9, 0);
                std::vector<std::pair<std::string, std::string>>().swap(tempContigs2);
                
                tempContigs2.swap(Contigs2);
                std::vector<std::pair<std::string, std::string>>().swap(Contigs2);
                generateCReads(&tempContigs2, repmeth, k, tranmeth, t, &ContigReads2);
                std::vector<NGS>().swap(ReferenveListCom2);
                ReferenceList(&ContigReads2, &ReferenveListCom2, k, tranmeth, t, sa);
                std::vector<std::vector<int>> tgVaridx2;
                if(isen == true && als < 2.00){
                    als +=siv;
                    if(als > 2.00){
                        als = 2.00;
                    }
                }
                std::vector<std::pair<std::string, std::string>> tContigs2;
                processReads(&ContigReads2, &ReferenveListCom2, &tContigs2, &tgVaridx2, knn, repmeth, k, als);
                std::vector<NGS>().swap(ReferenveListCom2);
                if( discardunmatch == true){
                    removedunmatched(NuSequences2, &tgVaridx2, &gVaridx2, &fgVaridx2, &Contigs2, &tContigs2, &ContigReads2);
                    discardunmatch = false;
                }
                else{
                    joinContigs(NuSequences2, &ContigReads2, &Contigs2, &tContigs2, &gVaridx2, &fgVaridx2, &tgVaridx2);
                }
                
                std::vector<std::vector<int>>().swap(tgVaridx2);
                std::vector<std::vector<int>>().swap(gVaridx2);
                
            }
        }
        
        printLabels(10, 0);
        
        if(c1 == true && c2 == true){
            printLabels(11, (int)Contigs1.size());
            printLabels(11, (int)Contigs2.size());
            printLabels(12, 0);
            generateFiles(NuSequences1, &Contigs1, &fgVaridx1, NuSequences2, &Contigs2, &fgVaridx2, outputsam, outputcontigs);
        }
        else if(c1 == true){
            printLabels(11, (int)Contigs1.size());
            printLabels(12, 0);
            generateFiles(NuSequences1, &Contigs1, &fgVaridx1, outputsam, outputcontigs);
        }
        else if(c2 == true){
            printLabels(11, (int)Contigs2.size());
            printLabels(12, 0);
            generateFiles(NuSequences2, &Contigs2, &fgVaridx2, outputsam, outputcontigs);
        }
        printLabels(13, 0);
    }
    else{
        printLabels(14, 0);
    }
    
};


IniAssemby::~IniAssemby(){};
