//
//  FastqParser.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "FastqParser.hpp"

//FastqParser::FastqParser(std::string filepath, std::string repmeth, std::string tranmeth, int tranlev, std::vector<NGS> *NuSequences, int refkmerlen, int alphasize)
FastqParser::FastqParser(std::string filepath, std::string repmeth, int refkmerlen, std::string tranmeth, int tranlev, std::vector<NGS> *NuSequences, int sa, int ss){
    std::ifstream fastqFile(filepath);
    std::string lineContents;
    std::string read, uid, qualstr;
    int count = 0, count1 = 0, readlength = 0;
    std::vector<std::vector<double>> saxrep ,saxtra;
    int  c = 0;
    while(!fastqFile.eof()){
        getline(fastqFile, lineContents);
        std::stringstream readstream(lineContents);
        if(count == 0 && lineContents[0] == '@'){
            //If line starts with @ char and count = 0 then asinge line content to uid
            uid = *new std::string;
            //uid = "";
            //read = "";
            //qualstr = "";
            //readlength = 0;
            std::getline (readstream, uid, '\n');
            //readstream >> uid;
            uid.erase(0,1);
            count++;
        }
        else{
            if(count == 1)
            {
                readlength = *new int;
                read = *new std::string;
                readstream >> read;
                readlength = (int)read.size();
                count++;
            }
            else if(count == 2){
                //If line doesn't starts with @ char and count = 2  do nothing
                count++;
            }
            else if(count == 3){
                qualstr = *new std::string;
                readstream >> qualstr;
                int reqlen = readlength;
                if(reqlen > 0)
                {
                    if(reqlen >= refkmerlen){
                        reqlen = refkmerlen;
                    }
                }
                
                // NGS n(uid, read, qualstr);
                // Declare read variable
                NGS n(uid, read, qualstr);
                n.setIdx(c);
                rep.createRepresentation(repmeth, read, n.returnRep());
                tra.createTransformation(tranmeth, n.returnRep(), n.returnTra(), 0, reqlen,  tranlev);
                (*NuSequences).push_back(n);
                c++;
                count = 0;
                count1++;
                if((count1%100000) == 0){
                    std::cout << "> Fastq read: " << count1 << "\n";
                }
            }
        }
    };
};

FastqParser::FastqParser(std::string filepath, std::string repmeth, int refkmerlen, std::string tranmeth, int tranlev, std::vector<NGS> *NuSequences1, std::vector<NGS> *NuSequences2, int sa, int ss, bool pe){
    //bool concatoantedreads = false;
    std::ifstream testpecase(filepath);
    std::string lineContents;
    int rf = 0;
    std::string name1 = "", name2 = "";
    while(!testpecase.eof() && rf < 2 ){
        getline(testpecase, lineContents);
        std::stringstream readstream(lineContents);
        if(lineContents[0] == '@'){
            if(rf == 0){
                name1 = lineContents.substr(0, lineContents.size()-2);
            }
            else{
                name2 = lineContents.substr(0, lineContents.size()-2);
            }
            rf++;
        }
    }
    bool conper = (name1 != name2);
    
    std::ifstream fastqFile(filepath);
    lineContents = "";
    std::string read, uid, qualstr;
    int count = 0, count1 = 0, readlength = 0;
    std::vector<std::vector<double>> saxrep ,saxtra;
    int  c = 0;
    bool p1 = true;
    while(!fastqFile.eof()){
        getline(fastqFile, lineContents);
        std::stringstream readstream(lineContents);
        if(count == 0 && lineContents[0] == '@'){
            //If line starts with @ char and count = 0 then asinge line content to uid
            uid = *new std::string;
            //uid = "";
            //read = "";
            //qualstr = "";
            //readlength = 0;
            std::getline (readstream, uid, '\n');
            //readstream >> uid;
            uid.erase(0,1);
            count++;
        }
        else{
            if(count == 1)
            {
                readlength = *new int;
                read = *new std::string;
                readstream >> read;
                readlength = (int)read.size();
                count++;
            }
            else if(count == 2){
                //If line doesn't starts with @ char and count = 2  do nothing
                count++;
            }
            else if(count == 3){
                qualstr = *new std::string;
                readstream >> qualstr;
                
                // NGS n(uid, read, qualstr);
                // Declare read variable
                
                
                if(conper == true){
                    int rsize = (int)read.size()/2;
                    
                    int reqlen = rsize;
                    if(reqlen > 0)
                    {
                        if(reqlen >= refkmerlen){
                            reqlen = refkmerlen;
                        }
                    }
                    NGS n1(uid+"/1", read.substr(0, rsize), qualstr.substr(0, rsize));
                    NGS n2(uid+"/2", read.substr(rsize, rsize), qualstr.substr(rsize, rsize));
                    n1.setIdx(c);
                    n2.setIdx(c);
                    rep.createRepresentation(repmeth, read.substr(0, rsize), n1.returnRep());
                    rep.createRepresentation(repmeth, read.substr(rsize, rsize), n2.returnRep());
                    tra.createTransformation(tranmeth, n1.returnRep(), n1.returnTra(), 0, reqlen,  tranlev);
                    tra.createTransformation(tranmeth, n2.returnRep(), n2.returnTra(), 0, reqlen,  tranlev);
                    (*NuSequences1).push_back(n1);
                    (*NuSequences2).push_back(n2);
                }
                else{
                    int reqlen = readlength;
                    if(reqlen > 0)
                    {
                        if(reqlen >= refkmerlen){
                            reqlen = refkmerlen;
                        }
                    }
                    
                    NGS n(uid, read, qualstr);
                    n.setIdx(c);
                    rep.createRepresentation(repmeth, read, n.returnRep());
                    tra.createTransformation(tranmeth, n.returnRep(), n.returnTra(), 0, reqlen,  tranlev);

                    if(p1 == true){
                        p1 = false;
                        (*NuSequences1).push_back(n);
                    }
                    else{
                        p1 = true;
                        (*NuSequences2).push_back(n);
                    }
                }
                c++;
                count = 0;
                count1++;
                if((count1%100000) == 0){
                    std::cout << "> Fastq read: " << count1 << "\n";
                }
            }
        }
    }
    
    
    
    
    
};


void FastqParser::processseparatedReads(std::string filepath, std::string repmeth, int refkmerlen, std::string tranmeth, int tranlev, std::vector<NGS> *NuSequences1, std::vector<NGS> *NuSequences2, int sa, int ss){
    
};


void FastqParser::processConactonatedReads(std::string filepath, std::string repmeth, int refkmerlen, std::string tranmeth, int tranlev, std::vector<NGS> *NuSequences1, std::vector<NGS> *NuSequences2, int sa, int ss){};

FastqParser::~FastqParser(){
    rep.~NucRepresentations();
    tra.~DataTransformations();
};
