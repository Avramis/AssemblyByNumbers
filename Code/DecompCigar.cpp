//
//  DecompCigar.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "DecompCigar.hpp"
void DecompCigar::decomposeCigars(std::string c, std::string *d){
    std::map<size_t, std::string, std::less<size_t>> cidec;
    size_t pos;
    
    pos = c.find("M",0);
    while(pos != std::string::npos){
        cidec[pos] = "M";
        pos = c.find("M",pos+1);
    }
    
    pos = c.find("I",0);
    while(pos != std::string::npos){
        cidec[pos] = "I";
        pos = c.find("I",pos+1);
    }
    
    pos = c.find("D",0);
    while(pos != std::string::npos){
        cidec[pos] = "D";
        pos = c.find("D",pos+1);
    }
    
    pos = c.find("S",0);
    while(pos != std::string::npos){
        cidec[pos] = "S";
        pos = c.find("S",pos+1);
    }
    
    pos = c.find("H",0);
    while(pos != std::string::npos){
        cidec[pos] = "H";
        pos = c.find("H",pos+1);
    }
    
    std::map<size_t, std::string, std::less<size_t>>::iterator fit;
    pos = 0;
    for (fit = cidec.begin(); fit !=cidec.end(); fit++){
        for (int i = 0; i < stoi(c.substr(pos, fit->first-pos)); i++){
            (*d)+=c.substr(fit->first, 1);
        }
        pos = fit->first+1;
    }
    
    std::map<size_t, std::string, std::less<size_t>>().swap(cidec);
};


DecompCigar::DecompCigar(){
};

DecompCigar::DecompCigar(std::string tc, std::string c, int a){
    CorrectCigarstep(tc, c, a);
};




std::string DecompCigar::returnEditedstring(){
    return cigar;
};

int DecompCigar::returnEditedalPos(){
    return ap;
};



void DecompCigar::CorrectCigarstep(std::string tc, std::string c, int a){
    std::string pipe = "";
    cigar = "";
    ap = 0;
    std::string fullt, fullc;
    
    fullt = "";
    decomposeCigars(tc, &fullt);
    
    fullc = "";
    //decomposeCigars(tc, &fullt);
    decomposeCigars(c, &fullc);
    
    
    int tap = a;
    while(tap > 0 && fullt.size() > 0){
        //if(fullt[0] != 'D'){
        if(fullt[0] == 'M'){
            tap--;
            
        }
        else if(fullt[0] == 'I'){
            tap--;
            ap--;
        }
        else{
            ap++;
        }
        fullt = fullt.substr(1);
    }
    
    //bool change=true, incc = false, inct = false;
    std::string value = "";
    //int len = 0;
    int temcoun = 0, cigcoun = 0;

    while(temcoun < (int)fullt.size() && cigcoun < (int)fullc.size()){
    
        if(fullt[temcoun] == 'M'){
            if(fullc[cigcoun] == 'M'){
                pipe+='M';
                cigcoun++;
                temcoun++;
            }
            //case cigar INSERTIONS
            else if(fullc[cigcoun] == 'I'){
                pipe+='I';
                cigcoun++;
                //temcoun++;
            }
            //case cigar DELETIONS
            else if(fullc[cigcoun] == 'D'){
                pipe+='D';
                cigcoun++;
                temcoun++;
            }
            //case cigar Substitution
            else if(fullc[cigcoun] == 'S'){
                pipe+='S';
                cigcoun++;
                temcoun++;
            }
        }
        //case template INSERTIONS
        else if(fullt[temcoun] == 'I'){
            //case cigar MATCH
            if(fullc[cigcoun] == 'M'){
                pipe+='I';
                cigcoun++;
                temcoun++;
            }
            //case cigar INSERTIONS
            else if(fullc[cigcoun] == 'I'){
                pipe+='I';
                cigcoun++;
            }
            //case cigar DELETIONS
            else if(fullc[cigcoun] == 'D'){
                cigcoun++;
                temcoun++;
                
            }
            //case cigar Substitution
            else if(fullc[cigcoun] == 'S'){
                pipe+='S';
                cigcoun++;
                temcoun++;
            }
        }
        //case template Seletions
        else if(fullt[temcoun] == 'D'){
            //case cigar MATCH
            if(fullc[cigcoun] == 'M'){
                pipe+='D';
                //cigcoun++;
                temcoun++;
            }
            //case cigar INSERTIONS
            else if(fullc[cigcoun] == 'I'){
                pipe+='M';
                cigcoun++;
                temcoun++;
            }
            //case cigar DELETIONS
            else if(fullc[cigcoun] == 'D'){
                pipe+='D';
                cigcoun++;
                //temcoun++;
            }
            //case cigar Substitution
            else if(fullc[cigcoun] == 'S'){
                pipe+='S';
                cigcoun++;
                temcoun++;
            }
        }
        
        
        //case template Seletions
        else if(fullt[temcoun] == 'S'){
            temcoun++;
        }
        
    }
    
    while(cigcoun < (int)fullc.size()){
        pipe +=fullc[cigcoun];
        cigcoun++;
    }
    decomposePipe( &pipe);
};




void DecompCigar::decomposePipe(std::string *pipe){
    
    std::string m = "";
    int n = 0;
    if ((*pipe).substr(0,1) == "D"){
        m = (*pipe).substr(0,1);
        while(m == (*pipe).substr(0,1)){
            ap++;
            (*pipe) = (*pipe).substr(1, (*pipe).size()-1);
        }
        m= "";
    }
    
    while ((*pipe).size() > 0){
        if(m == ""){
            m = (*pipe).substr(0,1);
            n = 1;
            (*pipe) = (*pipe).substr(1, (*pipe).size()-1);
        }
        else{
            if(m == (*pipe).substr(0,1)){
                n++;
                (*pipe) = (*pipe).substr(1,(*pipe).size()-1);
            }
            else{
                cigar+=(std::to_string(n)+ m);
                m = "";
            }
        }
        
    }
    
    if(m != ""){
        cigar+=(std::to_string(n)+ m);
    }
};

