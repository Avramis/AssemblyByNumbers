//
//  FileGeneration.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef FileGeneration_hpp
#define FileGeneration_hpp

#include <stdio.h>
#include <vector>
#include <string>
#include <fstream>

#include "NGS.hpp"
#include "ReverseCompliment.hpp"
class FileGeneration{
private:

    void generatefile(std::vector<NGS> *ShortReads, std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::vector<int>> *gVaridx, std::string outputsam, std::string outputcontigs);
public:
    FileGeneration();
    
    FileGeneration(std::vector<NGS> *ShortReads, std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::vector<int>> *gVaridx, std::string outputsam, std::string outputcontigs);
    
    
    FileGeneration(std::vector<NGS> *ShortReads1, std::vector<std::pair<std::string, std::string>> *Contigs1, std::vector<std::vector<int>> *gVaridx1, std::vector<NGS> *ShortReads2, std::vector<std::pair<std::string, std::string>> *Contigs2, std::vector<std::vector<int>> *gVaridx2, std::string outputsam, std::string outputcontigs);
    
    
    ~FileGeneration();
};
#endif /* FileGeneration_hpp */
