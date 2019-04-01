//
//  ProcessReads.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef ProcessReads_hpp
#define ProcessReads_hpp

#include <stdio.h>
#include <vector>
#include <map>
#include <math.h>
#include <limits.h>

#include "NucRepresentations.hpp"
#include "DataTransformations.hpp"
#include "DistanceMethods.hpp"
#include "NGS.hpp"
#include "VPTree.hpp"
#include "DecompCigar.hpp"


#include "Graph.hpp"
#include "Kindel.hpp"
#include "ssw_cpp.h"

class ProcessReads{
private:
    static bool sortNGSAl(std::pair<int, NGS*> n, std::pair<int, NGS*> m);
    
    void fixCigar(std::string *cigar, int *sp);
    
    void adjustpath(std::map<int, int> *route, int x, NGS *n, bool search);
    
    void ProcessContigs(std::vector<NGS> *ShortReads, std::vector< std::vector<int>> *pathlist,  std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::vector<int>> *gVaridx);
    
    void ProcessContigs(std::vector<NGS*> *ShortReads, std::vector< std::vector<int>> *pathlist,  std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::vector<int>> *gVaridx);

    //void CigarsCorrection(std::vector<NGS> *ShortReads, std::vector<int> *pathlist, std::map<int, int> *keyalig, std::vector<std::pair<int, NGS*>> *grouping, std::string pathname, int *conlen);
    void CigarsCorrection(std::vector<NGS> *ShortReads, std::vector<int> *pathlist, std::vector<std::pair<int, NGS*>> *grouping, std::string pathname, int *conlen);
    
    void CigarsCorrection(std::vector<NGS*> *ShortReads, std::vector<int> *pathlist, std::vector<std::pair<int, NGS*>> *grouping, std::string pathname, int *conlen);
    
    void TreeSearch(VPTree *VPT, NGS *query, std::vector<NGS> *ShortReads, std::vector <NGS> *neighbours, int knn, int *sc, int *absc);
    
    void TreeSearch(VPTree *VPT, NGS *query, std::vector<NGS*> *ShortReads, std::vector <NGS> *neighbours, int knn, int *sc, int *absc);
    
    void EvalauteNeighbours(std::vector<NGS> *ShortReads, int i, std::vector <NGS> *neighbours,  std::vector<std::vector<double>> *qrep, std::vector<std::vector<double>> *knrep, std::string rmeth, NucRepresentations *NR, DistanceMethods *DM);
    
    void EvalauteNeighbours(std::vector<NGS*> *ShortReads, int i, std::vector <NGS> *neighbours,  std::vector<std::vector<double>> *qrep, std::vector<std::vector<double>> *knrep, std::string rmeth, NucRepresentations *NR, DistanceMethods *DM);
    
    void SmithWatermanEvalaution(std::vector<NGS> *ShortReads, int i, std::vector <NGS> *neighbours, StripedSmithWaterman::Aligner *aligner, StripedSmithWaterman::Filter *filter, StripedSmithWaterman::Alignment *alignment, double als);
    
    void SmithWatermanEvalaution(std::vector<NGS*> *ShortReads, int i, std::vector <NGS> *neighbours, StripedSmithWaterman::Aligner *aligner, StripedSmithWaterman::Filter *filter, StripedSmithWaterman::Alignment *alignment, double als);

public:
    ProcessReads();
    
    ProcessReads(std::vector<NGS> *ShortReads, std::vector<NGS> *ReferenveList, std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::vector<int>> *gVaridx, int knn, std::string rmeth, int k, double als);
    
    ProcessReads(std::vector<NGS*> *ShortReads, std::vector<NGS> *ReferenveList, std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::vector<int>> *gVaridx, int knn, std::string rmeth, int k, double als);
    
    ProcessReads(std::vector<NGS> *ShortReads1, std::vector<NGS> *ReferenveList1, std::vector<std::pair<std::string, std::string>> *Contigs1, std::vector<std::vector<int>> *gVaridx1, std::vector<NGS> *ShortReads2, std::vector<NGS> *ReferenveList2, std::vector<std::pair<std::string, std::string>> *Contigs2, std::vector<std::vector<int>> *gVaridx2, int knn, std::string rmeth, int k, double als, bool rvmissmatch);
    
    
    ProcessReads(std::vector<NGS*> *ShortReads1, std::vector<NGS> *ReferenveList1, std::vector<std::pair<std::string, std::string>> *Contigs1, std::vector<std::vector<int>> *gVaridx1, std::vector<NGS*> *ShortReads2, std::vector<NGS> *ReferenveList2, std::vector<std::pair<std::string, std::string>> *Contigs2, std::vector<std::vector<int>> *gVaridx2, int knn, std::string rmeth, int k, double als, bool rvmissmatch);
    
    
    ~ProcessReads();
};
#endif /* ProcessReads_hpp */
