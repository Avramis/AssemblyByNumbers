//
//  IniAssemby.hpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#ifndef IniAssemby_hpp
#define IniAssemby_hpp

#include <stdio.h>
#include <string>
#include <map>

#include "NGS.hpp"
#include "FastqParser.hpp"
#include "GenerateTreeData.hpp"
#include "ProcessReads.hpp"
#include "GenerateCReads.hpp"
#include "FileGeneration.hpp"

class IniAssemby{
private:
    void printLabels(int i, int j);
    
    //void importReads(std::string inputfile, std::string repmeth, int k, std::string tranmeth, int t, std::map<std::string, std::vector<NGS>> *NGSbins, int sa, int ss, bool bin);
    void importReads(std::string inputfile, std::string repmeth, int k, std::string tranmeth, int t, std::vector<NGS> *NuSequences, int sa, int ss);
    
    void importReads(std::string p1inputfile, std::string p2inputfile, std::string repmeth, int k, std::string tranmeth, int t, std::vector<NGS> *NuSequences1, std::vector<NGS> *NuSequences2, int sa, int ss);
    
    void importReads(std::string inputfile, std::string repmeth, int k, std::string tranmeth, int t, std::vector<NGS> *NuSequences1, std::vector<NGS> *NuSequences2, int sa, int ss, bool pe);
    
    void ReferenceList(std::vector<NGS> *ShortReads, std::vector<NGS> *ReferenveList, int k, std::string tranmeth, int l, int sa);
    
    void ReferenceList(std::vector<NGS*> *ShortReads, std::vector<NGS> *ReferenveListCom, int k, std::string tranmeth, int l, int sa);
    
    //For pairedendreads
    void ReferenceList(std::vector<NGS> *ShortReads1, std::vector<NGS> *ReferenveListCom1, std::vector<NGS> *ShortReads2, std::vector<NGS> *ReferenveListCom2, int k, std::string tranmeth, int l, int sa);
    
    void processReads(std::vector<NGS> *ShortReads, std::vector<NGS> *ReferenveList, std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::vector<int>> *gVaridx, int knn, std::string rmeth, int k, double als);
    
    void processReads(std::vector<NGS*> *ShortReads, std::vector<NGS> *ReferenveList, std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::vector<int>> *gVaridx, int knn, std::string rmeth, int k, double als);
    
    //For pairedendreads
    void processReads(std::vector<NGS> *ShortReads1, std::vector<NGS> *ReferenveList1, std::vector<std::pair<std::string, std::string>> *Contigs1, std::vector<std::vector<int>> *gVaridx1, std::vector<NGS> *ShortReads2, std::vector<NGS> *ReferenveList2, std::vector<std::pair<std::string, std::string>> *Contigs2, std::vector<std::vector<int>> *gVaridx2, int knn, std::string rmeth, int k, double als, bool missmatch);
    
    void processReads(std::vector<NGS*> *ShortReads1, std::vector<NGS> *ReferenveList1, std::vector<std::pair<std::string, std::string>> *Contigs1, std::vector<std::vector<int>> *gVaridx1, std::vector<NGS*> *ShortReads2, std::vector<NGS> *ReferenveList2, std::vector<std::pair<std::string, std::string>> *Contigs2, std::vector<std::vector<int>> *gVaridx2, int knn, std::string rmeth, int k, double als, bool missmatch);
    
    void generateCReads(std::vector<std::pair<std::string, std::string>> *Contigs, std::string repmeth, int k, std::string tranmeth, int t, std::vector<NGS> *NuSequences);
    
    void generateFiles(std::vector<NGS> *ShortReads, std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::vector<int>> *varIdx, std::string outputsam, std::string outputcontigs);
    
    
    void generateFiles(std::vector<NGS> *ShortReads1, std::vector<std::pair<std::string, std::string>> *Contigs1, std::vector<std::vector<int>> *varIdx1, std::vector<NGS> *ShortReads2, std::vector<std::pair<std::string, std::string>> *Contigs2, std::vector<std::vector<int>> *varIdx2, std::string outputsam, std::string outputcontigs);
    
    void processPEReads (std::vector<NGS> *NuSequences1, std::vector<NGS> *NuSequences2, std::string outputsam, std::string outputcontigs, std::string repmeth, int k, std::string tranmeth, int t, int ss, int sa, int binmet, bool bin, int knn, double als, bool isen, double siv, bool missmatch);
    
    void removedunmatched(std::vector<std::pair<std::string, std::string>> *tContigs, std::vector<std::vector<int>> *gVaridx, std::vector<std::vector<int>> *fgVaridx, std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<NGS> *NuSequences);
    
    void removedunmatched(std::vector<NGS> *NuSequences, std::vector<std::vector<int>> *tgVaridx, std::vector<std::vector<int>> *gVaridx, std::vector<std::vector<int>> *fgVaridx, std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::pair<std::string, std::string>> *tContigs, std::vector<NGS> *ContigReads);
    
    void readsBinning(std::vector<NGS> *NuSequences, std::map<std::string, std::vector<NGS*>> *NGSbins, int sa, int ss);
    
    void readsBinning(std::vector<NGS> *NuSequences1, std::vector<NGS> *NuSequences2, std::map<std::string, std::vector<NGS*>> *NGSbins1, std::map<std::string, std::vector<NGS*>> *NGSbins2, int sa, int ss);
    
    void processBins(std::map<std::string, std::vector<NGS*>> *NGSbins, std::vector<NGS> *ReferenveListCom, std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::vector<int>> *fgVaridx, int k, std::string repmeth, std::string tranmeth, int t, int sa, int knn, double als);
    
    
    void processBins(std::map<std::string, std::vector<NGS*>> *NGSbins1, std::vector<NGS> *ReferenveListCom1, std::vector<std::pair<std::string, std::string>> *Contigs1, std::vector<std::vector<int>> *fgVaridx1, std::map<std::string, std::vector<NGS*>> *NGSbins2, std::vector<NGS> *ReferenveListCom2, std::vector<std::pair<std::string, std::string>> *Contigs2, std::vector<std::vector<int>> *fgVaridx2, int k, std::string repmeth, std::string tranmeth, int t, int sa, int knn, double als);
    
    void joinContigs(std::vector<NGS> *NuSequences, std::vector<NGS> *ContigReads, std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::pair<std::string, std::string>> *tContigs, std::vector<std::vector<int>> *gVaridx, std::vector<std::vector<int>> *fgVaridx, std::vector<std::vector<int>> *tgVaridx);
    
public:
    
    IniAssemby(std::string inputfile, std::string outputsam, std::string outputcontigs, std::string repmeth, int k, std::string tranmeth, int t, int ss, int sa, int binmet, bool bin, int knn, double als, bool isen, double siv);
    
    IniAssemby(std::string inputfile1, std::string inputfile2, std::string outputsam, std::string outputcontigs, std::string repmeth, int k, std::string tranmeth, int t, int ss, int sa, int binmet, bool bin, int knn, double als, bool isen, double siv, bool dealmissmatch);
    
    IniAssemby(std::string inputfile, std::string outputsam, std::string outputcontigs, std::string repmeth, int k, std::string tranmeth, int t, int ss, int sa, int binmet, bool bin, int knn, double als, bool isen, double siv, bool pe, bool dealmissmatch);
    
    ~IniAssemby();
};

#endif /* IniAssemby_hpp */

