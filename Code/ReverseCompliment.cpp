//
//  ReverseCompliment.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "ReverseCompliment.hpp"
//Compares nucleotide letters and return the reverse compliment
char ReverseCompliment::NucleotideComparison(char x){
    if(toupper(x) == 'A' )
    {
        x='T';
    }
    else if(toupper(x) == 'T' || toupper(x) == 'U' )
    {
        x='A';
    }
    else if(toupper(x) == 'C' )
    {
        x='G';
    }
    else if(toupper(x) == 'G' )
    {
        x='C';
    }
    else if(toupper(x) == 'N' )
    {
        x='N';
    }
    else if(toupper(x) == 'M' )
    {
        x='K';
    }
    else if(toupper(x) == 'R' )
    {
        x='Y';
    }
    else if(toupper(x) == 'W' )
    {
        x='W';
    }
    else if(toupper(x) == 'S' )
    {
        x='S';
    }
    else if(toupper(x) == 'Y' )
    {
        x='R';
    }
    else if(toupper(x) == 'K' )
    {
        x='M';
    }
    else if(toupper(x) == 'V' )
    {
        x='B';
    }
    else if(toupper(x) == 'H' )
    {
        x='D';
    }
    else if(toupper(x) == 'D' )
    {
        x='H';
    }
    else if(toupper(x) == 'B' )
    {
        x='V';
    }
    return x;
}

//Destructor
ReverseCompliment::~ReverseCompliment(){};

//Constructors
//Accepts a read and generates the reverse compliment
ReverseCompliment::ReverseCompliment(std::string read){
    setReverese(read);
};

//Empty constractor
ReverseCompliment::ReverseCompliment(){};


//Generate the reverse complimen of the nucleotide read
void ReverseCompliment::setReverese(std::string read){
    rvread = read;
    std::reverse(rvread.begin(),rvread.end());
    for(int i = 0; i < (int) rvread.size(); i++){
        rvread[i]=NucleotideComparison(rvread[i]);
    }
    
};



//Return the reverse compliment nucleotide
std::string ReverseCompliment::getReverse(){
    return rvread;
}
