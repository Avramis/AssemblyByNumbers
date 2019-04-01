//
//  NormalDistribution.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef NormalDistribution_hpp
#define NormalDistribution_hpp

#include <stdio.h>
#include <math.h>
#include <vector>

#include "DataProcessing.hpp"

class NormalDistribution{
private:
    double mu, stdc;
    double inverfc(double p);
public:
    NormalDistribution(double m, double s);
    double p(double x);
    double cdf(double x) ;
    double invcdf(double p);
    
    void invcdfsequence(std::vector<double> *seq, int s);
    
    ~NormalDistribution();
};

#endif /* NormalDistribution_hpp */
