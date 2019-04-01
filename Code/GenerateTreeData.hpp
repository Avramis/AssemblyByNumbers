//
//  GenerateTreeData.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef GenerateTreeData_hpp
#define GenerateTreeData_hpp

#include <stdio.h>
#include <vector>
#include <map>

#include "NGS.hpp"
#include "DataTransformations.hpp"

class GenerateTreeData{
private:
public:
    GenerateTreeData();
    GenerateTreeData(std::vector<NGS> *ShortReads, std::vector<NGS> *ReferenveListCom, int k, std::string tranmeth, int l, int sa);
    GenerateTreeData(std::vector<NGS*> *ShortReads, std::vector<NGS> *ReferenveListCom, int k, std::string tranmeth, int l, int sa);
    ~GenerateTreeData();
};
#endif /* GenerateTreeData_hpp */
