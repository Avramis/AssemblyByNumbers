//
//  NGS.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "NGS.hpp"
NGS::NGS(){};
//For short reads
NGS::NGS(std::string u, std::string r, std::string q){
    setSRid(u);
    setSread(r);
    setQual(q);
};

NGS::NGS(std::string u, std::string r, std::string q, int k){
    setSRid(u);
    setSread(r);
    setQual(q);
    setKmer(k);
};

//For reference dubliacted reads

NGS::NGS(int g, int s){
    
};

NGS::NGS( int g, int s, int k){
    
};

NGS::NGS(std::string *u, int g, int s){
    addGnp(u);
    addCoordinates(g, s);
    
};
NGS::NGS(std::string *u, int g, int s, int k){
    addGnp(u);
    addCoordinates(g, s);
    setKmer(k);
};


NGS::NGS(int g, int s, std::string *r){
    setRefreafp(r);
    addCoordinates(g, s);
};

NGS::NGS(int g, int s, int k, std::string *r){
    setRefreafp(r);
    addCoordinates(g, s);
    setKmer(k);
};

NGS::NGS(std::string *u, std::string *r, int g, int s){
    addGnp(u);
    setRefreafp(r);
    addCoordinates(g, s);
};

NGS::NGS(std::string *u, std::string *r, int g, int s, int k){
    addGnp(u);
    setRefreafp(r);
    addCoordinates(g, s);
    setKmer(k);
};

void NGS::setSRid(std::string s){
    rid = s;
};

void NGS::setSread(std::string s){
    read = s;
    //fulllen = (int)read.size();
};

void NGS::setQual(std::string s){
    qual = s;
};

void NGS::setKmer(int i){
    kmer = i;
};

void NGS::setIdx(int i){
    idx = i;
};

void NGS::setAldetails(int i, int j, std::string s, bool b){
    addCoordinates(i, j);
    cigars.push_back(s);
    alDir.push_back(b);
};

void NGS::setAldetails(int i, int j, std::string s, std::string r){
    addCoordinates(i, j);
    cigars.push_back(s);
    alRef.push_back(r);
};

void NGS::addGnp(std::string *s){
    gid.push_back(s);
};

void NGS::addCoordinates(int i, int p){
    coordinates.push_back(std::make_pair(i, p));
};

void NGS::clearRep(){
    std::vector<std::vector<double>>().swap(nFwrep);
};

void NGS::clearTra(){
    std::vector<std::vector<double>>().swap(nFwtra);
};

void NGS::setRefreafp(std::string *r){
    refread = r;
    //gid.push_back(r);
};

void NGS::deleteRead(){
    read = "";
    //fulllen = 0;
    
};

void NGS::deleteUid(){
    rid = "";
};


void NGS::setAlscore(double i){
    alSc = i;
};



void NGS::setAlVecscore(double i){
    alScVec.push_back(i);
};


void NGS::clearAldetails(){
    alSc = 0;
    std::vector <double>().swap(alScVec);
    std::vector <std::string>().swap(alRef);
    //std::vector <std::string*>().swap(gid);
    std::vector<std::pair<int, int>>().swap(coordinates);
    std::vector<std::string>().swap(cigars);
};

void NGS::setCigar(std::string c){
    cigars.push_back(c);
};

void NGS::setAlRef(std::string s){
    alRef.push_back(s);
};

void NGS::editCigar(int i, std::string c){
    cigars.at(i) = c;
};

void NGS::editCoordinates(int i, int c){
    coordinates.at(i).second = c;
};

//void setRep(std::string m);
//void setRep(std::string m, int s);
//void setTra(std::string m, int l);

std::string NGS::returnRead(){
    return read;
};

std::string NGS::returnSRid(){
    return rid ;
};

std::string NGS::returnQual(){
    return qual;
};

int NGS::returnKmer(){
    return kmer;
};

std::vector<std::vector<double>> * NGS::returnRep(){
    return &nFwrep;
};

std::vector<std::vector<double>> * NGS::returnTra(){
    return &nFwtra;
};

std::string* NGS::returnrefReads(){
    return refread;
};

std::vector<std::string*> NGS::returnrefUID(){
    return gid;
};

std::string *NGS::returnSRidP(){
    return &read;
};

std::string *NGS::returnUIP(){
    return &rid;
};

double NGS::returnAlscore(){
    return alSc;
};

std::vector<double> NGS::returnAlVec(){
    return alScVec;
};

std::vector<std::pair<int, int>> NGS::returnCoordinates(){
    return  coordinates;
};


std::vector <std::string*> NGS::returnGID(){
    return gid;
};


std::vector <std::string> NGS::returnAlref(){
    return alRef;
};

std::vector <std::string> NGS::returnCigar(){
    return cigars;
};

int NGS::returnIdx(){
    return idx;
};

