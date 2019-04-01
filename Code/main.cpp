//
//  main.cpp
//  AsembleByNumberSamCode
//
//  Created by Avraam Tapinos on 28/06/2017.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include <iostream>
#include <map>
#include <sys/stat.h>
#include <sys/types.h>
#include "IniAssemby.hpp"
#include "InstructionClass.hpp"

int main(int argc, const char * argv[]) {
    
    std::string filepath, filepath1, filepath2, repmeth, tranmeth, soutdir, foutdir;
    //std::string repmeth, tranmeth, outdir;
    int kmer, tranlevl, blvl, alpha_size, knn;
    double als  = 1.0, siv = 0.5;;
    repmeth = "Tetrahedron";
    tranmeth="DWT";
    kmer = 100;
    tranlevl = 4;
    blvl = 1;
    alpha_size = 1;
    knn = 5;
    bool fqc = false, fqc1 = false, fqc2 = false, soutp = false, foutp = false, isen = false, pe = false, twofiles = false, dealmismatches = false;
    int civ = (argc - 1) % 2;
    int binmet = 1;
    bool bin = false;
    
    if (civ || argc - 1 == 0) {
        for (int i = 1; i < argc; i++) {
            if ((std::string) (argv[i]) == "-h"
                || (std::string) (argv[i]) == "-help"
                || (std::string) (argv[i]) == "-H"
                || (std::string) (argv[i]) == "-Help"
                || (std::string) (argv[i]) == "-HELP") {
                InstructionClass(true);
                return 1;
            }
        }
        InstructionClass instreuctions(false);
        return 1;
    }
    else {
        for (int i = 1; i < argc; i += 2) {
            if ((std::string) (argv[i]) == "-h"
                || (std::string) (argv[i]) == "-help"
                || (std::string) (argv[i]) == "-H"
                || (std::string) (argv[i]) == "-Help"
                || (std::string) (argv[i]) == "-HELP") {
                InstructionClass instreuctions(true);
                return 1;
            }
            else if ((std::string) (argv[i]) == "-fq") {
                filepath = (std::string) (argv[i + 1]);
                if (FILE *file = fopen(filepath.c_str(), "r")) {
                    fclose(file);
                    fqc = true;
                    pe = false;
                    twofiles = false;
                }
            }
            else if ((std::string) (argv[i]) == "-fq1") {
                filepath1 = (std::string) (argv[i + 1]);
                if (FILE *file = fopen(filepath1.c_str(), "r")) {
                    fclose(file);
                    fqc1 = true;
                    pe = true;
                    twofiles = true;
                }
            }
            else if ((std::string) (argv[i]) == "-fq2") {
                filepath2 = (std::string) (argv[i + 1]);
                if (FILE *file = fopen(filepath2.c_str(), "r")) {
                    fclose(file);
                    fqc2 = true;
                    pe = true;
                    twofiles = true;
                }
            }
            else if((std::string) (argv[i]) == "-pe") {
                std::cout <<"\n\nValue in pe is " << ((std::string) argv[i + 1]) <<"\n\n";
                if((std::string) argv[i + 1] == "true"){
                    pe = true;
                }
            }
            else if((std::string) (argv[i]) == "-so"){
                soutdir = (std::string) (argv[i + 1]);
                size_t sidx = soutdir.find_last_of("/");
                struct stat sb;
                if (!(stat(soutdir.substr(0,sidx+1).c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))){
                    std::cout << "Couldn't find Folder " << soutdir.substr(0,sidx+1) <<"\n";
                }
                else{
                    soutp = true;
                }
            }
            else if((std::string) (argv[i]) == "-fo"){
                foutdir = (std::string) (argv[i + 1]);
                size_t sidx = foutdir.find_last_of("/");
                struct stat sb;
                if (!(stat(foutdir.substr(0,sidx+1).c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))){
                    std::cout << "Couldn't find Folder " << foutdir.substr(0,sidx+1) <<"\n";
                }
                else{
                    foutp = true;
                }
            }
            else if((std::string) (argv[i]) == "-ovlp"){
                kmer = atoi(argv[i + 1]);
            }
            else if((std::string) (argv[i]) == "-rep"){
                repmeth = (std::string) (argv[i + 1]);
            }
            else if((std::string) (argv[i]) == "-tra"){
                tranmeth = (std::string) (argv[i + 1]);
            }
            else if((std::string) (argv[i]) == "-clvl"){
                tranlevl = atoi(argv[i + 1]);
            }
            else if((std::string) (argv[i]) == "-knn"){
                knn = atoi(argv[i + 1]);
            }
            else if ((std::string) (argv[i]) == "-bin") {
                if( (std::string) argv[i + 1] == "true"){
                    alpha_size = 2;
                    bin = true;
                }
            }
            else if ((std::string) (argv[i]) == "-binmet") {
                if( atoi(argv[i + 1])  >= 1 && atoi(argv[i + 1]) <= 2 ){
                    binmet = atoi(argv[i + 1]);
                }
            }
 
            else if ((std::string) (argv[i]) == "-siv") {
                siv = atof(argv[i + 1]);
            }
            else if((std::string) (argv[i]) == "-ssi") {
                if((std::string) argv[i + 1] == "true"){
                    isen = true;
                }
                else if((std::string) argv[i + 1] == "false"){
                    isen = false;
                }
                else{
                    std::cout << "Wrong -sei parameter provided. ["<< (std::string) argv[i + 1] << "instead of true or false value]\n";
                    InstructionClass instreuctions(false);
                    return -1;
                }
                
            }
            else if((std::string) (argv[i]) == "-sen"){
                if(-2.0 <= atof(argv[i + 1]) && atof(argv[i + 1]) <= 2.0 ){
                    als = atof(argv[i + 1]);
                }
                else{
                    std::cout << "Non valid sensitivity value provided. Value of 1.0 will be used as sensitivity value.\n";
                }
            }
        }
    }
    
    if(fqc == false){
        if(fqc1 == false || fqc2 == false){
            std::cout << "Please provide (a) valid fq file(s).\n";
            InstructionClass instreuctions(false);
            return 1;
        }
    }

    if(pe == true){
        if(fqc1 == false || fqc2 == false){
            std::cout << "Please provide (a) valid fq file(s).\n";
            InstructionClass instreuctions(false);
            return 1;
        }   
    }
    
    
    if(soutp == false){
        size_t sidx;
        if(foutp == true){
            sidx = foutdir.find_last_of("/");
            soutdir = foutdir.substr(0, sidx+1) + "_Approximate_OLC_Alignment.sam";
            std::cout <<"Sam a file directory has not provided.\n";
            std::cout << "File will be saved in " << soutdir << "\n";
        }
        else{
            if(twofiles == true){
                filepath = filepath1;
            }
            sidx = filepath.find_last_of(".");
            soutdir = filepath.substr(0, sidx) + "_Approximate_OLC_Alignment.sam";
            std::cout <<"Sam a file directory has not provided.\n";
            std::cout << "File will be saved in " << soutdir << "\n";
            
        }
        std::cout <<  "Sam file will be save in the following directory: " << soutdir<<"\n";
    }
    
    if(foutp == false){
        if(soutp == true){
            size_t sidx = soutdir.find_last_of("/");
            foutdir = soutdir.substr(0, sidx+1) + "Contigs.fasta";
            std::cout <<"Fasta a file directory has not provided.\n";
            std::cout << "File will be saved in " << foutdir << "\n";
        }
        else{
            size_t sidx = filepath.find_last_of(".");
            foutdir = filepath.substr(0, sidx) + "Contigs.fasta";
            std::cout <<"Fasta a file directory has not provided.\n";
            std::cout << "File will be saved in " << foutdir << "\n";
        }
        std::cout <<  "Fatsa file will be save in the following directory: " << foutdir <<"\n";
    }
    
    //ReadsBining(filepath, repmeth, kmer, tranmeth, tranlevl, blvl, alpha_size, knn, soutdir, foutdir, als, isen, siv, binmet, bin);

    if( pe == false){
        IniAssemby InA(filepath, soutdir, foutdir, repmeth, kmer, tranmeth, tranlevl, blvl, alpha_size, binmet, bin, knn, als, isen, siv);
    }
    else{
        if(twofiles == false){
            IniAssemby InA(filepath, soutdir, foutdir, repmeth, kmer, tranmeth, tranlevl, blvl, alpha_size, binmet, bin, knn, als, isen, siv, pe, dealmismatches);
        }
        else{
            IniAssemby InA(filepath1, filepath2, soutdir, foutdir, repmeth, kmer, tranmeth, tranlevl, blvl, alpha_size, binmet, bin, knn, als, isen, siv, dealmismatches);
        }
    }
    
    return 0;
};

