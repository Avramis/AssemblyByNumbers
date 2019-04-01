//
//  ReverseCompliment.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef ReverseCompliment_hpp
#define ReverseCompliment_hpp

#include <stdio.h>
#include <stdio.h>
#include <string>
#include <algorithm>

class ReverseCompliment{
private:
    //Variable
    std::string rvread;
    
    //Compares nucleotide letters and return the reverse compliment
    char NucleotideComparison(char x);
    
public:
    //Destructor
    ~ReverseCompliment();
    
    //Constructors:
    //Accepts a read and generates the reverse compliment
    ReverseCompliment(std::string read);
    
    //Empty constractor
    ReverseCompliment();
    
    
    //Generate the reverse complimen of the nucleotide read
    void setReverese(std::string read);
    
    //Return the reverse compliment nucleotide
    std::string getReverse();
    
};
#endif /* ReverseCompliment_hpp */
