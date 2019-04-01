//
//  FileGeneration.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "FileGeneration.hpp"
FileGeneration::FileGeneration(){};

FileGeneration::FileGeneration(std::vector<NGS> *ShortReads, std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::vector<int>> *gVaridx, std::string outputsam, std::string outputcontigs){
    generatefile(ShortReads, Contigs, gVaridx, outputsam, outputcontigs);
    
};

void FileGeneration::generatefile(std::vector<NGS> *ShortReads, std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::vector<int>> *gVaridx, std::string outputsam, std::string outputcontigs){
    //AssemblyByNumbers
    std::ofstream fout(outputcontigs);
    std::ofstream sout(outputsam);
    sout << "@HD" << "\t" << "VN:1.0" << "\t" << "SO:unknown" << "\n";
    int cus = 0;
    int ccs = 0;
    for(int i = 0; i < (int)(*Contigs).size(); i++){
        std::string curcontigname = "";
        if((*gVaridx)[i].size() == 1){
            
            curcontigname = "Unassembled_Read_" + std::to_string(cus);
            cus++;
        }
        else{
            curcontigname = "Contig_" + std::to_string(ccs);
            ccs++;
        }
        (*Contigs)[i].first = curcontigname;
        
        sout << "@SQ" << "\t" << "SN:" << (*Contigs)[i].first  << "\t" << "LN:"<< std::to_string((int)(*Contigs)[i].second.size())<< "\n";
        fout << ">" << (*Contigs)[i].first <<"\n" << (*Contigs)[i].second << "\n";
        
    }
    fout.close();
    sout << "@PG" <<"\t" <<"ID:AssemblyByNumbers" << "\t" << "PN:AssemblyByNumbers" << "\t" << "VN:1" << "\t" << "CL:N/A" << "\n";
    for(int i = 0; i < (int)(*gVaridx).size(); i++){
        for (int j =0; j < (int)(*gVaridx)[i].size(); j++){
            sout << (*ShortReads)[(*gVaridx)[i][j]].returnSRid() << "\t";
            if((int)(*ShortReads)[(*gVaridx)[i][j]].returnCoordinates().size() == 0){
                sout << "4\t";
                sout << "*\t";
                sout << "0\t";
                sout << "0\t";
                sout << "*\t";
                sout << "*\t";
                sout << "0\t";
                sout << "0\t";
                sout << (*ShortReads)[(*gVaridx)[i][j]].returnRead() << "\t";
                sout << (*ShortReads)[(*gVaridx)[i][j]].returnQual() << "\n";
            }
            else{
                sout << "0\t";
                sout << (*Contigs)[i].first << "\t";
                sout << (*ShortReads)[(*gVaridx)[i][j]].returnCoordinates()[0].second + 1 << "\t";
                sout << "40\t";
                sout << (*ShortReads)[(*gVaridx)[i][j]].returnCigar()[0] << "\t";
                sout << "*" << "\t";
                sout << "0" << "\t";
                sout << "0" << "\t";
                sout << (*ShortReads)[(*gVaridx)[i][j]].returnRead() << "\t";
                sout << (*ShortReads)[(*gVaridx)[i][j]].returnQual() << "\n";
        }
        }
    }
    sout.close();
};





FileGeneration::FileGeneration(std::vector<NGS> *ShortReads1, std::vector<std::pair<std::string, std::string>> *Contigs1, std::vector<std::vector<int>> *gVaridx1, std::vector<NGS> *ShortReads2, std::vector<std::pair<std::string, std::string>> *Contigs2, std::vector<std::vector<int>> *gVaridx2, std::string outputsam, std::string outputcontigs){
    std::ofstream fout(outputcontigs);
    std::ofstream sout(outputsam);
    sout << "@HD" << "\t" << "VN:1.0" << "\t" << "SO:unknown" << "\n";
    int cus = 0;
    int ccs = 0;
    for(int i = 0; i < (int)(*Contigs1).size(); i++){
        std::string curcontigname = "";
        if((*gVaridx1)[i].size() == 1){
            
            curcontigname = "Pair_1_Unassembled_Read" + std::to_string(cus);
            cus++;
        }
        else{
            curcontigname = "Pair_1_Contig_" + std::to_string(ccs);
            ccs++;
        }
        (*Contigs1)[i].first = curcontigname;
        
        sout << "@SQ" << "\t" << "SN:" << (*Contigs1)[i].first  << "\t" << "LN:"<< std::to_string((int)(*Contigs1)[i].second.size())<< "\n";
        fout << ">" << (*Contigs1)[i].first <<"\n" << (*Contigs1)[i].second << "\n";
    }
    
    cus = 0;
    ccs = 0;
    for(int i = 0; i < (int)(*Contigs2).size(); i++){
        std::string curcontigname = "";
        if((*gVaridx2)[i].size() == 1){
            
            curcontigname = "Pair_2_Unassembled_Read_" + std::to_string(cus);
            cus++;
        }
        else{
            curcontigname = "Pair_2_Contig_" + std::to_string(ccs);
            ccs++;
        }
        (*Contigs2)[i].first = curcontigname;
        
        sout << "@SQ" << "\t" << "SN:" << (*Contigs2)[i].first  << "\t" << "LN:"<< std::to_string((int)(*Contigs2)[i].second.size())<< "\n";
        fout << ">" << (*Contigs2)[i].first <<"\n" << (*Contigs2)[i].second << "\n";
        
    }
    fout.close();
    sout << "@PG" <<"\t" <<"ID:AssemblyByNumbers" << "\t" << "PN:AssemblyByNumbers" << "\t" << "VN:1" << "\t" << "CL:N/A" << "\n";
    
    for(int i = 0; i < (int)(*gVaridx1).size(); i++){
        for (int j =0; j < (int)(*gVaridx1)[i].size(); j++){
            sout << (*ShortReads1)[(*gVaridx1)[i][j]].returnSRid() << "\t";
            if((int)(*ShortReads1)[(*gVaridx1)[i][j]].returnCoordinates().size() == 0){
                sout << "4\t";
                sout << "*\t";
                sout << "0\t";
                sout << "0\t";
                sout << "*\t";
                sout << "*\t";
                sout << "0\t";
                sout << "0\t";
                sout << (*ShortReads1)[(*gVaridx1)[i][j]].returnRead() << "\t";
                sout << (*ShortReads1)[(*gVaridx1)[i][j]].returnQual() << "\n";
            }
            else{
                sout << "0\t";
                sout << (*Contigs1)[i].first << "\t";
                sout << (*ShortReads1)[(*gVaridx1)[i][j]].returnCoordinates()[0].second + 1 << "\t";
                sout << "40\t";
                sout << (*ShortReads1)[(*gVaridx1)[i][j]].returnCigar()[0] << "\t";
                sout << "*" << "\t";
                sout << "0" << "\t";
                sout << "0" << "\t";
                sout << (*ShortReads1)[(*gVaridx1)[i][j]].returnRead() << "\t";
                sout << (*ShortReads1)[(*gVaridx1)[i][j]].returnQual() << "\n";
            }
        }
    }
    
    for(int i = 0; i < (int)(*gVaridx2).size(); i++){
        for (int j =0; j < (int)(*gVaridx2)[i].size(); j++){
            sout << (*ShortReads2)[(*gVaridx2)[i][j]].returnSRid() << "\t";
            if((int)(*ShortReads2)[(*gVaridx2)[i][j]].returnCoordinates().size() == 0){
                sout << "4\t";
                sout << "*\t";
                sout << "0\t";
                sout << "0\t";
                sout << "*\t";
                sout << "*\t";
                sout << "0\t";
                sout << "0\t";
                sout << (*ShortReads2)[(*gVaridx2)[i][j]].returnRead() << "\t";
                sout << (*ShortReads2)[(*gVaridx2)[i][j]].returnQual() << "\n";
            }
            else{
                sout << "0\t";
                sout << (*Contigs2)[i].first << "\t";
                sout << (*ShortReads2)[(*gVaridx2)[i][j]].returnCoordinates()[0].second + 1 << "\t";
                sout << "40\t";
                sout << (*ShortReads2)[(*gVaridx2)[i][j]].returnCigar()[0] << "\t";
                sout << "*" << "\t";
                sout << "0" << "\t";
                sout << "0" << "\t";
                sout << (*ShortReads2)[(*gVaridx2)[i][j]].returnRead() << "\t";
                sout << (*ShortReads2)[(*gVaridx2)[i][j]].returnQual() << "\n";
            }
        }
    }
    sout.close();
};


FileGeneration::~FileGeneration(){};
