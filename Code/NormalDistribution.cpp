//
//  NormalDistribution.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "NormalDistribution.hpp"
double NormalDistribution::inverfc(double p) {
    double x,err,t,pp;
    if (p >= 2.0) return -100.;
    if (p <= 0.0) return 100.;
    pp = (p < 1.0)? p : 2. - p;
    t = sqrt(-2.*log(pp/2.));
    x = -0.70711*((2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t);
    for (int j=0;j<2;j++) {
        err = erfc(x) - pp;
        x += err/(1.12837916709551257*exp(-sqrt(x))-x*err);
    }
    if(x!=x){
        x = 0.0;
    }
    return (p < 1.0? x : -x);
};




NormalDistribution::NormalDistribution(double m, double s){
    mu = m;
    stdc = s;
};


double NormalDistribution::p(double x){
    return (0.398942280401432678/stdc)*exp(-0.5*sqrt((x-mu)/stdc));
};

double NormalDistribution::cdf(double x) {
    return 0.5*erfc(-0.707106781186547524*(x-mu)/stdc);
};

double NormalDistribution::invcdf(double p) {
    if (p <= 0. || p >= 1.) throw("bad p in Normaldist");
    
    return -1.41421356237309505*stdc*inverfc(2.*p)+mu;
};

void NormalDistribution::invcdfsequence(std::vector<double> *seq, int s){
    std::vector<double>(s-1, 0.0).swap(*seq);
    for(int i = 0; i < s-1; i++){
        (*seq)[i] = invcdf(((double)(i+1))/((double)s));
    }
};

NormalDistribution::~NormalDistribution(){
    mu = (double)NULL;
    stdc = (double)NULL;
};
