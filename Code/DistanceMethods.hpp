//
//  DistanceMethods.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef DistanceMethods_hpp
#define DistanceMethods_hpp

#include <stdio.h>
#include <vector>
#include <math.h>

#include "NGS.hpp"
#include "SVDClass.hpp"

class DistanceMethods{
private:
    SVDClass svdc;
public:
    DistanceMethods();
    //double returnRepDistance(NGS *n, NGS *m);
    //double returnTraDistance(NGS *n, NGS *m);
    
    double returnRepDistance(NGS *n, NGS *m, bool t);
    double returnTraDistance(NGS *n, NGS *m, bool t);
    
    double returnDistance(std::vector<std::vector<double>> *n, std::vector<std::vector<double>> *m);
    double returnEuclidianNorm(std::vector<std::vector<double>> *n, std::vector<std::vector<double>> *m);
    
    double returnRepDTW(NGS *n, NGS *m, int w);
    double returnDTW(std::vector<std::vector<double>> *n, std::vector<std::vector<double>> *m, int w);
    
};
#endif /* DistanceMethods_hpp */
