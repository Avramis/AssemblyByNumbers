//
//  DecompCigar.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//


#ifndef DecompCigar_hpp
#define DecompCigar_hpp

#include <stdio.h>
#include <map>
#include <string>
#include <iostream>
//#include "DecompCigar.hpp"
class  DecompCigar{
private:
    std::string cigar;
    int ap;
    
    
public:
    DecompCigar();
    DecompCigar(std::string tc, std::string c, int a);
    void decomposeCigars(std::string c, std::string *d);
    void CorrectCigarstep(std::string tc, std::string c, int a);
    void decomposePipe(std::string *pipe);
    std::string returnEditedstring();
    int returnEditedalPos();
    
};

#endif /* DecompCigar_hpp */
