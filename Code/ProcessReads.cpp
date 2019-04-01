//
//  ProcessReads.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "ProcessReads.hpp"

//ProcessReads::ProcessReads(){};

bool ProcessReads::sortNGSAl(std::pair<int, NGS*> n, std::pair<int, NGS*> m){
    return n.second->returnCoordinates()[0].second < m.second->returnCoordinates()[0].second;
};

void ProcessReads::adjustpath(std::map<int, int> *route, int x, NGS *n, bool search){
    int nidx = 0, pos = 0;
    double alSc = 0.0;
    std::string cigar = "", ref = "";
    if(search == true){
        bool exist = false;
        for (int i = 0; i < (int)(*n).returnCoordinates().size(); i++){
            if((*route).find((*n).returnCoordinates()[i].first)!= ((*route).end())){
                nidx = (*n).returnCoordinates()[i].first;
                pos = (*n).returnCoordinates()[i].second;
                cigar = (*n).returnCigar()[i];
                alSc = (*n).returnAlVec()[i];
                if((int)(*n).returnAlref().size() > i){
                    ref = (*n).returnAlref()[i];
                }
                exist = true;
                break;
            }
        }
        
        if(exist == true){
            (*n).clearAldetails();
            (*n).setAldetails(nidx, pos, cigar,ref);
            (*n).setAlVecscore(alSc);
        }
        else{
            nidx = x;
            pos = 0;
            alSc = 0;
            cigar = std::to_string((*n).returnRead().size()) + "M";
            (*n).clearAldetails();
            (*n).setAldetails(nidx, pos, cigar,ref);
            (*n).setAlVecscore(alSc);
        }
    }
    else{
        nidx = x;
        pos = 0;
        alSc = 0;
        cigar = std::to_string((*n).returnRead().size()) + "M";
        (*n).clearAldetails();
        (*n).setAldetails(nidx, pos, cigar,ref);
        (*n).setAlVecscore(alSc);
    }
    (*route)[x] = nidx;
};

void ProcessReads::fixCigar(std::string *cigar, int *sp){
    
    std::size_t cigM, cigS, cigI, cigD;
    cigM = (*cigar).find("M");
    cigS = (*cigar).find("S");
    cigI = (*cigar).find("I");
    cigD = (*cigar).find("D");
    if(cigS < cigM && cigS < cigI && cigS < cigI &&  cigS < cigD){
        (*sp) -= stoi((*cigar).substr(0,cigS));
        (*cigar) = (*cigar).substr(0,cigS)+ "M" + (*cigar).substr(cigS+1);
    }
    if((*cigar)[(*cigar).size()-1] == 'S'){
        (*cigar)[(*cigar).size()-1] = 'M';
    }
};

void ProcessReads::TreeSearch(VPTree *VPT, NGS *query, std::vector<NGS> *ShortReads, std::vector <NGS> *neighbours, int knn, int *sc, int *absc){
    std::map<double, std::vector<NGS>, std::greater<double>> dneighbours;
    (*VPT).FullKnnSearch(query, knn, &dneighbours, sc, absc);
    
    if((int)dneighbours.size() > 0){
        for (std::map<double, std::vector<NGS>, std::greater<double>>::iterator it = dneighbours.begin(); it != dneighbours.end(); it++){
            for(int j = 0; j < (int)dneighbours[it->first].size(); j++){
                for (int z = 0; z < (int)dneighbours[it->first].at(j).returnCoordinates().size(); z++){
                    if((*ShortReads)[dneighbours[it->first].at(j).returnCoordinates().at(z).first].returnSRid() != (*query).returnSRid()){
                        NGS n((*ShortReads)[dneighbours[it->first].at(j).returnCoordinates().at(z).first].returnUIP(), dneighbours[it->first].at(j).returnCoordinates().at(z).first, dneighbours[it->first].at(j).returnCoordinates().at(z).second, (*ShortReads)[dneighbours[it->first].at(j).returnCoordinates().at(z).first].returnKmer());
                        (*neighbours).push_back(n);
                    }
                }
            }
        }
    }
    std::map<double, std::vector<NGS>, std::greater<double>>().swap(dneighbours);
};



void ProcessReads::TreeSearch(VPTree *VPT, NGS *query, std::vector<NGS*> *ShortReads, std::vector <NGS> *neighbours, int knn, int *sc, int *absc){
    std::map<double, std::vector<NGS>, std::greater<double>> dneighbours;
    (*VPT).FullKnnSearch(query, knn, &dneighbours, sc, absc);
    
    if((int)dneighbours.size() > 0){
        for (std::map<double, std::vector<NGS>, std::greater<double>>::iterator it = dneighbours.begin(); it != dneighbours.end(); it++){
            for(int j = 0; j < (int)dneighbours[it->first].size(); j++){
                for (int z = 0; z < (int)dneighbours[it->first].at(j).returnCoordinates().size(); z++){
                    if((*(*ShortReads)[dneighbours[it->first].at(j).returnCoordinates().at(z).first]).returnSRid() != (*query).returnSRid()){
                        NGS n((*(*ShortReads)[dneighbours[it->first].at(j).returnCoordinates().at(z).first]).returnUIP(), dneighbours[it->first].at(j).returnCoordinates().at(z).first, dneighbours[it->first].at(j).returnCoordinates().at(z).second, (*(*ShortReads)[dneighbours[it->first].at(j).returnCoordinates().at(z).first]).returnKmer());
                        (*neighbours).push_back(n);
                    }
                }
            }
        }
    }
    std::map<double, std::vector<NGS>, std::greater<double>>().swap(dneighbours);
};


void ProcessReads::EvalauteNeighbours(std::vector<NGS> *ShortReads, int i, std::vector <NGS> *neighbours,  std::vector<std::vector<double>> *qrep, std::vector<std::vector<double>> *knrep, std::string rmeth, NucRepresentations *NR, DistanceMethods *DM){
    std::vector<std::vector<double>>().swap(*qrep);
    std::vector<std::vector<double>>().swap(*knrep);
    double bestdist = std::numeric_limits<double>::infinity();
    int sOVLP = std::numeric_limits<int>::max();
    (*NR).createRepresentation(rmeth, (*ShortReads)[i].returnRead().substr(0), qrep);
    std::vector<NGS> tempneighbours;
    tempneighbours.swap((*neighbours));
    std::vector<NGS>().swap((*neighbours));
    
    for (int j = 0; j < (int)tempneighbours.size(); j++){
        for (int z = 0; z < (int)tempneighbours[j].returnGID().size(); z++){
            if((*ShortReads)[tempneighbours[j].returnCoordinates()[z].first].returnSRid() != (*ShortReads)[i].returnSRid()){
                std::vector<std::vector<double>>().swap(*knrep);
                (*NR).createRepresentation(rmeth, (*ShortReads)[tempneighbours[j].returnCoordinates()[z].first].returnRead().substr(tempneighbours[j].returnCoordinates()[z].second), knrep);
                int tx = (int)(*qrep)[0].size();
                if(tx > (int)(*knrep)[0].size()){
                    tx = (int)(*knrep)[0].size();
                    
                }
                int w = 6;
                double tempdist = (*DM).returnDTW(qrep, knrep, w);
                if(bestdist >= tempdist){
                    int tOVLP = (int)(*ShortReads)[tempneighbours[j].returnCoordinates()[z].first].returnRead().substr(tempneighbours[j].returnCoordinates()[z].second).size();
                    if(bestdist > tempdist){
                        bestdist = tempdist;
                        std::vector<NGS >().swap((*neighbours));
                        // Testing length addition
                        sOVLP = tOVLP;
                        
                    }
                    if(sOVLP >= tOVLP){
                        if(sOVLP > tOVLP){
                            sOVLP = tOVLP;
                            std::vector<NGS >().swap((*neighbours));
                        }
                        NGS n((*ShortReads)[tempneighbours[j].returnCoordinates()[z].first].returnUIP(), tempneighbours[j].returnCoordinates()[z].first, tempneighbours[j].returnCoordinates()[z].second, (*ShortReads)[tempneighbours[j].returnCoordinates()[z].first].returnKmer());
                        (*neighbours).push_back(n);
                    }
                }
            }

        }
    }
    
    if((int)(*neighbours).size() <= 0){
        (*ShortReads)[i].clearAldetails();
    }
    std::vector<std::vector<double>>().swap(*qrep);
    std::vector<std::vector<double>>().swap(*knrep);
};


void ProcessReads::EvalauteNeighbours(std::vector<NGS*> *ShortReads, int i, std::vector <NGS> *neighbours,  std::vector<std::vector<double>> *qrep, std::vector<std::vector<double>> *knrep, std::string rmeth, NucRepresentations *NR, DistanceMethods *DM){
    std::vector<std::vector<double>>().swap(*qrep);
    std::vector<std::vector<double>>().swap(*knrep);
    double bestdist = std::numeric_limits<double>::infinity();
    int sOVLP = std::numeric_limits<int>::max();
    (*NR).createRepresentation(rmeth, (*(*ShortReads)[i]).returnRead().substr(0), qrep);
    std::vector<NGS> tempneighbours;
    tempneighbours.swap((*neighbours));
    std::vector<NGS>().swap((*neighbours));
    
    for (int j = 0; j < (int)tempneighbours.size(); j++){
        for (int z = 0; z < (int)tempneighbours[j].returnGID().size(); z++){
            if((*(*ShortReads)[tempneighbours[j].returnCoordinates()[z].first]).returnSRid() != (*(*ShortReads)[i]).returnSRid()){
                std::vector<std::vector<double>>().swap(*knrep);
                (*NR).createRepresentation(rmeth, (*(*ShortReads)[tempneighbours[j].returnCoordinates()[z].first]).returnRead().substr(tempneighbours[j].returnCoordinates()[z].second), knrep);
                int tx = (int)(*qrep)[0].size();
                if(tx > (int)(*knrep)[0].size()){
                    tx = (int)(*knrep)[0].size();
                    
                }
                int w = 6;
                double tempdist = (*DM).returnDTW(qrep, knrep, w);
                if(bestdist >= tempdist){
                    int tOVLP = (int)(*(*ShortReads)[tempneighbours[j].returnCoordinates()[z].first]).returnRead().substr(tempneighbours[j].returnCoordinates()[z].second).size();
                    if(bestdist > tempdist){
                        bestdist = tempdist;
                        std::vector<NGS >().swap((*neighbours));
                        // Testing length addition
                        sOVLP = tOVLP;
                    }
                    if(sOVLP >= tOVLP){
                        if(sOVLP > tOVLP){
                            sOVLP = tOVLP;
                            std::vector<NGS >().swap((*neighbours));
                        }
                        NGS n((*(*ShortReads)[tempneighbours[j].returnCoordinates()[z].first]).returnUIP(), tempneighbours[j].returnCoordinates()[z].first, tempneighbours[j].returnCoordinates()[z].second, (*(*ShortReads)[tempneighbours[j].returnCoordinates()[z].first]).returnKmer());
                        (*neighbours).push_back(n);
                    }
                }
            }
        }
    }
    
    if((int)(*neighbours).size() <= 0){
        (*(*ShortReads)[i]).clearAldetails();
    }
    std::vector<std::vector<double>>().swap(*qrep);
    std::vector<std::vector<double>>().swap(*knrep);
};

void ProcessReads::SmithWatermanEvalaution(std::vector<NGS> *ShortReads, int i, std::vector <NGS> *neighbours, StripedSmithWaterman::Aligner *aligner, StripedSmithWaterman::Filter *filter, StripedSmithWaterman::Alignment *alignment, double als){
    bool alas = false;
    std::map<int, int> dubalingment;
    std::string cigarstring;
    //uint16_t bestscore = 0;
    double bestscore = 0;
    int bwsize = 0;
    for (int j = 0; j < (int)(*neighbours).size(); j++){
        //int ap = (*neighbours)[j].returnCoordinates()[0].second-10;
        int ap = (*neighbours)[j].returnCoordinates()[0].second;
        if(ap < 0){
            ap = 0;
        }; //int ap = neighbours[j].returnCoordinates()[0].second;
        int pidx = (*neighbours)[j].returnCoordinates()[0].first;
        //int wsize = (int)(*ShortReads)[i].returnRead().size()+20; //int wsize = (int)(*ShortReads)[i].returnRead().size();
        int wsize = (int)(*ShortReads)[i].returnRead().size(); //int wsize = (int)(*ShortReads)[i].returnRead().size();
        if(wsize > (int)(*ShortReads)[pidx].returnRead().substr(ap).size()){
            wsize = (int)(*ShortReads)[pidx].returnRead().substr(ap).size();
        }
        (*aligner).Align((*ShortReads)[i].returnRead().c_str(),(*ShortReads)[(*neighbours)[j].returnCoordinates()[0].first].returnRead().substr(ap).c_str(), wsize, (*filter), alignment);
        if(((double)(*alignment).sw_score) /((double) wsize) >= als){
            alas = true;
            if(bestscore <= ((double)(*alignment).sw_score)/((double) wsize)){
                if(bestscore < ((double)(*alignment).sw_score)/((double) wsize)){
                    bestscore = ((double)(*alignment).sw_score)/((double) wsize);
            //if(bestscore <= (*alignment).sw_score){
                //if(bestscore < (*alignment).sw_score){
                    //bestscore = (*alignment).sw_score;
                    (*ShortReads)[i].clearAldetails();
                    //(*ShortReads)[i].setAlscore(bestscore);
                    //(*ShortReads)[i].setAlscore((double)(*alignment).sw_score);
                    bwsize = 0;
                    bwsize = wsize;
                    std::map<int, int>().swap(dubalingment);
                }
                std::map<int, int>::iterator dit = dubalingment.find((*neighbours)[j].returnCoordinates()[0].first);
                if(dit != dubalingment.end()){
                    if(dit->second != ap + (*alignment).ref_begin){
                        //(*ShortReads)[i].setAldetails(neighbours[j].returnCoordinates()[0].first, neighbours[j].returnCoordinates()[0].second + alignment.ref_begin, alignment.cigar_string, *neighbours[j].returnGID()[0]);
                        dubalingment[(*neighbours)[j].returnCoordinates()[0].first] = ap + (*alignment).ref_begin;
                        (*ShortReads)[i].setAldetails((*neighbours)[j].returnCoordinates()[0].first, ap + (*alignment).ref_begin, (*alignment).cigar_string, (*(*neighbours)[j].returnGID()[0]));
                        //(*ShortReads)[i].setAlVecscore((double)(*alignment).sw_score);
                        (*ShortReads)[i].setAlVecscore(-1);
                    }
                }
                else{
                    dubalingment[(*neighbours)[j].returnCoordinates()[0].first] = ap + (*alignment).ref_begin;
                    (*ShortReads)[i].setAldetails((*neighbours)[j].returnCoordinates()[0].first, ap + (*alignment).ref_begin, (*alignment).cigar_string, *((*neighbours)[j].returnGID()[0]));
                    //(*ShortReads)[i].setAlscore((double)(*alignment).sw_score);
                    (*ShortReads)[i].setAlVecscore(bestscore);
                }
            }
        }
    }
    std::map<int, int>().swap(dubalingment);
    std::vector<NGS>().swap((*neighbours));
    if(alas == false){
        (*ShortReads)[i].clearAldetails();
        (*ShortReads)[i].setAldetails(i, 0, (std::to_string((*ShortReads)[i].returnRead().size())+ "M"), (*ShortReads)[i].returnUIP());
        //(*ShortReads)[i].setAlscore(0.0);
        //(*ShortReads)[i].setAlVecscore((double)(*alignment).sw_score);
        (*ShortReads)[i].setAlVecscore(-1);
    }
    else{
        for (int j = 0; j < (int)(*ShortReads)[i].returnCoordinates().size(); j++){
            
            std::string cigar = (*ShortReads)[i].returnCigar().at(j);
            int sp = (*ShortReads)[i].returnCoordinates().at(j).second;
            fixCigar(&cigar, &sp);
            (*ShortReads)[i].editCoordinates(j,sp);
            (*ShortReads)[i].editCigar(j, cigar);
            //(*G).addEdge(i, (*ShortReads)[i].returnCoordinates()[j].first, (*ShortReads)[i].returnAlscore());
        }
    }
};


void ProcessReads::SmithWatermanEvalaution(std::vector<NGS*> *ShortReads, int i, std::vector <NGS> *neighbours, StripedSmithWaterman::Aligner *aligner, StripedSmithWaterman::Filter *filter, StripedSmithWaterman::Alignment *alignment, double als){
    bool alas = false;
    std::map<int, int> dubalingment;
    std::string cigarstring;
    //uint16_t bestscore = 0;
    double bestscore = 0;
    int bwsize = 0;
    for (int j = 0; j < (int)(*neighbours).size(); j++){
        //int ap = (*neighbours)[j].returnCoordinates()[0].second-10;
        int ap = (*neighbours)[j].returnCoordinates()[0].second;
        if(ap < 0){
            ap = 0;
        }; //int ap = neighbours[j].returnCoordinates()[0].second;
        int pidx = (*neighbours)[j].returnCoordinates()[0].first;
        //int wsize = (int)(*(*ShortReads)[i]).returnRead().size()+20; //int wsize = (int)(*ShortReads)[i].returnRead().size();
        int wsize = (int)(*(*ShortReads)[i]).returnRead().size(); //int wsize = (int)(*ShortReads)[i].returnRead().size();
        if(wsize > (int)(*(*ShortReads)[pidx]).returnRead().substr(ap).size()){
            wsize = (int)(*(*ShortReads)[pidx]).returnRead().substr(ap).size();
        }
        (*aligner).Align((*(*ShortReads)[i]).returnRead().c_str(),(*(*ShortReads)[(*neighbours)[j].returnCoordinates()[0].first]).returnRead().substr(ap).c_str(), wsize, (*filter), alignment);
        if(((double)(*alignment).sw_score) /((double) wsize) >= als){
            alas = true;
            if(bestscore <= ((double)(*alignment).sw_score)/((double) wsize)){
                if(bestscore < ((double)(*alignment).sw_score)/((double) wsize)){
                    bestscore = ((double)(*alignment).sw_score)/((double) wsize);
            //if(bestscore <= (*alignment).sw_score){
                //if(bestscore < (*alignment).sw_score){
                    //bestscore = (*alignment).sw_score;
                    (*(*ShortReads)[i]).clearAldetails();
                    //(*(*ShortReads)[i]).setAlscore(bestscore);
                    //(*(*ShortReads)[i]).setAlscore((double)(*alignment).sw_score);
                    bwsize = 0;
                    bwsize = wsize;
                    std::map<int, int>().swap(dubalingment);
                }
                std::map<int, int>::iterator dit = dubalingment.find((*neighbours)[j].returnCoordinates()[0].first);
                if(dit != dubalingment.end()){
                    if(dit->second != ap + (*alignment).ref_begin){
                        //(*ShortReads)[i].setAldetails(neighbours[j].returnCoordinates()[0].first, neighbours[j].returnCoordinates()[0].second + alignment.ref_begin, alignment.cigar_string, *neighbours[j].returnGID()[0]);
                        dubalingment[(*neighbours)[j].returnCoordinates()[0].first] = ap + (*alignment).ref_begin;
                        (*(*ShortReads)[i]).setAldetails((*neighbours)[j].returnCoordinates()[0].first, ap + (*alignment).ref_begin, (*alignment).cigar_string, (*(*neighbours)[j].returnGID()[0]));
                        (*(*ShortReads)[i]).setAlscore((*alignment).sw_score);
                    }
                }
                else{
                    dubalingment[(*neighbours)[j].returnCoordinates()[0].first] = ap + (*alignment).ref_begin;
                    (*(*ShortReads)[i]).setAldetails((*neighbours)[j].returnCoordinates()[0].first, ap + (*alignment).ref_begin, (*alignment).cigar_string, *((*neighbours)[j].returnGID()[0]));
                    (*(*ShortReads)[i]).setAlscore((*alignment).sw_score);
                }
            }
        }
    }
    std::map<int, int>().swap(dubalingment);
    std::vector<NGS>().swap((*neighbours));
    if(alas == false){
        (*(*ShortReads)[i]).clearAldetails();
        (*(*ShortReads)[i]).setAldetails(i, 0, (std::to_string((*(*ShortReads)[i]).returnRead().size())+ "M"), (*(*ShortReads)[i]).returnUIP());
        (*(*ShortReads)[i]).setAlscore(0.0);
    }
    else{
        for (int j = 0; j < (int)(*(*ShortReads)[i]).returnCoordinates().size(); j++){
            
            std::string cigar = (*(*ShortReads)[i]).returnCigar().at(j);
            int sp = (*(*ShortReads)[i]).returnCoordinates().at(j).second;
            fixCigar(&cigar, &sp);
            (*(*ShortReads)[i]).editCoordinates(j,sp);
            (*(*ShortReads)[i]).editCigar(j, cigar);
            //(*G).addEdge(i, (*ShortReads)[i].returnCoordinates()[j].first, (*ShortReads)[i].returnAlscore());
        }
    }
};


ProcessReads::ProcessReads(std::vector<NGS> *ShortReads1, std::vector<NGS> *ReferenveList1, std::vector<std::pair<std::string, std::string>> *Contigs1, std::vector<std::vector<int>> *gVaridx1, std::vector<NGS> *ShortReads2, std::vector<NGS> *ReferenveList2, std::vector<std::pair<std::string, std::string>> *Contigs2, std::vector<std::vector<int>> *gVaridx2, int knn, std::string rmeth, int k, double als, bool rvmissmatch){
    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;
    DistanceMethods DM;
    NucRepresentations NR;
    std::vector<std::vector<double>> qrep, knrep;
    VPTree VPT1;
    VPTree VPT2;
    std::cout << "> Build VP-trees.\n";
    std::cout << "> Generate VPTree for mate pair 1 reads.\n";
    VPT1.TreeGeneration(ReferenveList1);
    std::cout << "> Generate VPTree for mate pair 2 reads.\n";
    VPT2.TreeGeneration(ReferenveList2);
    
    Graph G1((int) (*ShortReads1).size());
    Graph G2((int) (*ShortReads2).size());
    
    int pv = 4;
    
    int datasize = (int) (*ShortReads1).size();
    if(datasize < (int) (*ShortReads2).size()){
        datasize = (int) (*ShortReads2).size();
        rvmissmatch = false;
    }
    if(datasize  < 100){
        pv = 2;
    }
    if(datasize  < 50){
        pv = 1;
    }
    int per = 100/pv;
    int cv = (int)ceil(((double)datasize/(double)pv));
    std::cout<< "> [] 0%\n";
    for (int i = 0; i < datasize; i++){
        if((i%cv == 0) && i != 0){
            std::cout<< "> [" << std::string( (per * i /(cv))/2, '*')<<"] " << per * i /(cv) <<"%\n";
        }
        
        int sc = 0;
        int absc = 0;
        
        if(i < (int) (*ShortReads1).size()){
            std::vector <NGS> neighbours1;
            TreeSearch(&VPT1, &(*ShortReads1).at(i), ShortReads1, &neighbours1, knn, &sc, &absc);
            if((int)neighbours1.size() > 0){
                EvalauteNeighbours(ShortReads1, i, &neighbours1,  &qrep, &knrep, rmeth, &NR, &DM);
                if((int)neighbours1.size() > 0){
                    SmithWatermanEvalaution(ShortReads1, i, &neighbours1, &aligner, &filter, &alignment, als);
                }
            }
        }
        if(i < (int) (*ShortReads2).size()){
            std::vector <NGS> neighbours2;
            sc = 0;
            absc = 0;
            TreeSearch(&VPT2, &(*ShortReads2).at(i), ShortReads2, &neighbours2, knn, &sc, &absc);
            if((int)neighbours2.size() > 0){
                EvalauteNeighbours(ShortReads2, i, &neighbours2,  &qrep, &knrep, rmeth, &NR, &DM);
                if((int)neighbours2.size() > 0){
                    SmithWatermanEvalaution(ShortReads2, i, &neighbours2, &aligner, &filter, &alignment, als);
                }
            }
        }
        if(rvmissmatch == true){
            if((int)(*ShortReads1)[i].returnCoordinates().size() ==  0 || (int)(*ShortReads2)[i].returnCoordinates().size() ==  0){
                (*ShortReads1)[i].clearAldetails();
                (*ShortReads2)[i].clearAldetails();
            }
        }
        
        if(i < (int) (*ShortReads1).size()){
            for (int j = 0; j < (int)(*ShortReads1)[i].returnCoordinates().size(); j++){
                //G1.addEdge(i,(*ShortReads1)[i].returnCoordinates()[j].first, (*ShortReads1)[i].returnAlscore()[j]);
                G1.addEdge(i,(*ShortReads1)[i].returnCoordinates()[j].first, (*ShortReads1)[i].returnAlVec()[j]);
            }
            (*ShortReads1)[i].clearTra();
            (*ShortReads1)[i].clearRep();
        }
        
        if(i < (int) (*ShortReads2).size()){
            for (int j = 0; j < (int)(*ShortReads2)[i].returnCoordinates().size(); j++){
                G2.addEdge(i,(*ShortReads2)[i].returnCoordinates()[j].first, (*ShortReads2)[i].returnAlVec()[j]);
            }
            (*ShortReads2)[i].clearTra();
            (*ShortReads2)[i].clearRep();
        }
        
    }
    std::cout<< "> [" << std::string(50, '*')<<"] 100%\n";
    std::vector<NGS>().swap(*ReferenveList1);
    std::vector<NGS>().swap(*ReferenveList2);
    std::vector<std::vector<int>> pathlist1;
    std::vector<std::vector<int>> pathlist2;
    std::cout << "> Identify paths in graphs.\n";
    
    G1.BFS(&pathlist1);
    for (int i = 0; i < (int)pathlist1.size(); i++){
        std::map<int, int> route;
        adjustpath(&route, pathlist1[i][pathlist1[i].size()-1], &(*ShortReads1)[pathlist1[i][pathlist1[i].size()-1]], false);
        for (int j = (int)pathlist1[i].size()-2; j >= 0; j--){
            adjustpath(&route, pathlist1[i][j], &(*ShortReads1)[pathlist1[i][j]], true);
        }
        std::map<int, int>().swap(route);
    }
    G2.BFS(&pathlist2);
    for (int i = 0; i < (int)pathlist2.size(); i++){
        std::map<int, int> route;
        adjustpath(&route, pathlist2[i][pathlist2[i].size()-1], &(*ShortReads2)[pathlist2[i][pathlist2[i].size()-1]], false);
        for (int j = (int)pathlist2[i].size()-2; j >= 0; j--){
            adjustpath(&route, pathlist2[i][j], &(*ShortReads2)[pathlist2[i][j]], true);
        }
        std::map<int, int>().swap(route);
    }
    
    std::cout <<"> Generate " << (int)pathlist1.size() << " contig(s) from " << (*ShortReads1).size() << " mate pair 1 reads datapoint(s).\n";
    std::cout <<"> Generate " << (int)pathlist2.size() << " contig(s) from " << (*ShortReads2).size() << " mate pair 2 reads datapoint(s).\n";
    ProcessContigs(ShortReads1, &pathlist1, Contigs1, gVaridx1);
    ProcessContigs(ShortReads2, &pathlist2, Contigs2, gVaridx2);
};

ProcessReads::ProcessReads(std::vector<NGS> *ShortReads, std::vector<NGS> *ReferenveList, std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::vector<int>> *gVaridx, int knn, std::string rmeth, int k, double als){
    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;
    DistanceMethods DM;
    NucRepresentations NR;
    VPTree VPT;
    std::vector<std::vector<double>> qrep, knrep;
    std::cout << "> Build VP-tree.\n";
    VPT.TreeGeneration(ReferenveList);
    
    Graph G((int) ShortReads->size());

    std::cout << "> Link reads and build graph.\n";
    int pv = 4;
    
    if(ShortReads->size() < 100){
        pv = 2;
    }
    if(ShortReads->size() < 50){
        pv = 1;
    }
    
    int per = 100/pv;
    
    int cv = (int)ceil(((double)ShortReads->size()/(double)pv));
    std::cout<< "> [] 0%\n";
    for (int i = 0; i < (int) ShortReads->size(); i++){
        if((i%cv == 0) && i != 0){
            std::cout<< "> [" << std::string( (per * i /(cv))/2, '*')<<"] " << per * i /(cv) <<"%\n";
        }
        int sc = 0;
        int absc = 0;
        std::vector <NGS> neighbours;
        TreeSearch(&VPT, &ShortReads->at(i), ShortReads, &neighbours, knn, &sc, &absc);
        if((int)neighbours.size() > 0){
            EvalauteNeighbours(ShortReads, i, &neighbours,  &qrep, &knrep, rmeth, &NR, &DM);
            if((int)neighbours.size() > 0){
                SmithWatermanEvalaution(ShortReads, i, &neighbours, &aligner, &filter, &alignment, als);
            }
        }
        for (int j = 0; j < (*ShortReads)[i].returnCoordinates().size(); j++){
            G.addEdge(i, (*ShortReads)[i].returnCoordinates()[j].first, (*ShortReads)[i].returnAlVec()[j]);
        }
        
        ShortReads->at(i).clearTra();
        ShortReads->at(i).clearRep();
    }
    std::cout<< "> [" << std::string(50, '*')<<"] 100%\n";
    std::vector<NGS>().swap( *ReferenveList);
    std::vector<std::vector<int>> pathlist;
    std::cout << "> Identify paths in graph.\n";
    bool check = G.isCyclicFull();
    G.BFS(&pathlist);

    for (int i = 0; i < (int)pathlist.size(); i++){
        std::map<int, int> route;
        adjustpath(&route, pathlist[i][pathlist[i].size()-1], &(*ShortReads)[pathlist[i][pathlist[i].size()-1]], false);
        for (int j = (int)pathlist[i].size()-2; j >= 0; j--){
            adjustpath(&route, pathlist[i][j], &(*ShortReads)[pathlist[i][j]], true);
        }
        std::map<int, int>().swap(route);
    }
    
    std::cout <<"> Generate " << (int)pathlist.size() << " contig(s) from " << (*ShortReads).size() <<" datapoint(s).\n";
    ProcessContigs(ShortReads, &pathlist, Contigs, gVaridx);
};

void ProcessReads::ProcessContigs(std::vector<NGS> *ShortReads, std::vector< std::vector<int>> *pathlist,  std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::vector<int>> *gVaridx){
    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;
    std::vector<std::pair<std::string, int>> contlengths;
    for (int i = 0; i < (*pathlist).size(); i++){
        int conlen = std::numeric_limits<int>::min();
        std::vector<std::pair<int, NGS*>> grouping;
        std::string UID = "Contig_"+std::to_string((*Contigs).size());
        
        //if((int)(*pathlist)[i].size() > 1){
            CigarsCorrection(ShortReads , &(*pathlist)[i], &grouping, UID, &conlen);
            Kindel Ki(UID, conlen);
            
            for (int j = 0; j < grouping.size(); j++){
                Ki.addRead((*grouping[j].second).returnCoordinates()[0].second, (*grouping[j].second).returnCigar()[0], (*grouping[j].second).returnRead());
            }
            
            std::vector<int> dp;
            Ki.generateConsensus(true, &dp);
            
            std::string iniContig = Ki.getConsensus();
            std::string Cname = Ki.getContigName();
            Ki.deleteRefMap();
            for (int j = 0; j < grouping.size(); j++){
                //
                int  sp = (*grouping[j].second).returnCoordinates()[0].second - dp[(*grouping[j].second).returnCoordinates()[0].second];
                
                if(sp - 10 > 0 ){
                    sp -= 10;
                }
                else{
                    sp = 0;
                }
                int ep = sp + (int)(*grouping[j].second).returnRead().size();
                if(ep + 10 < (int)iniContig.size()){
                    ep += 10;
                }
                else{
                    ep = (int)iniContig.size();
                }
                
                int wsize = ep - sp;
                aligner.Align((*grouping[j].second).returnRead().c_str(),iniContig.substr(sp, wsize).c_str(), wsize, filter, &alignment);
                
                std::string tcigar = alignment.cigar_string;
                fixCigar(&tcigar, &sp);
                (*grouping[j].second).clearAldetails();
                (*grouping[j].second).setAldetails(0, sp + alignment.ref_begin, tcigar, Cname);
            }
            sort(grouping.begin(), grouping.end(), sortNGSAl);
            if((*grouping[0].second).returnCoordinates()[0].second < 0){
                conlen = std::numeric_limits<int>::min();
                int incval = abs((*grouping[0].second).returnCoordinates()[0].second)+1;
                for (int s = 0; s < grouping.size(); s++){
                    (*grouping[s].second).editCoordinates(0, (*grouping[s].second).returnCoordinates()[0].second + incval);
                    if(conlen < (*grouping[s].second).returnCoordinates()[0].second +  (int)(*grouping[s].second).returnRead().size()){
                        conlen = (*grouping[s].second).returnCoordinates()[0].second +  (int)(*grouping[s].second).returnRead().size();
                    }
                }
            }
            Ki.deleteRefMap();
            Ki.generateKindle(UID, conlen);
            
            int gidx = (int)grouping.size() - 1;
            std::vector<int>idx(gidx+1, 0);
            while(!grouping.empty()){
                Ki.addRead((*grouping[gidx].second).returnCoordinates()[0].second, (*grouping[gidx].second).returnCigar()[0], (*grouping[gidx].second).returnRead());
                idx[gidx] = grouping[gidx].first;
                gidx--;
                grouping.pop_back();
            }
            Ki.generateConsensus(true);
            Contigs->push_back( std::make_pair(Ki.getContigName(), Ki.getConsensus()));
            (*gVaridx).push_back(idx);
            std::vector<int>().swap(dp);
            std::vector<int>().swap(idx);
            std::vector<std::pair<int, NGS*>>().swap(grouping);
        //}
        //else{
        //    Contigs->push_back( std::make_pair(UID, (*ShortReads)[(*pathlist)[i][0]].returnRead()));
        //}
    }
};

void ProcessReads::CigarsCorrection(std::vector<NGS> *ShortReads, std::vector<int> *pathlist, std::vector<std::pair<int, NGS*>> *grouping, std::string pathname, int *conlen){
    DecompCigar DC;
    int ps = (int)(*pathlist).size();
    int sav = 0;
    int oav = 0;
    int nidx = (*pathlist)[ps - 1];
    (*ShortReads)[nidx].clearAldetails();
    (*ShortReads)[nidx].setAldetails(nidx, 0, (std::to_string((int)(*ShortReads)[nidx].returnRead().size()) + "M"), pathname);
    //(*keyalig)[nidx] =  0;
    (*grouping).push_back(std::make_pair(nidx, &(*ShortReads)[nidx]));
    
    if((*conlen) < (*ShortReads)[nidx].returnCoordinates()[0].second + (int)(*ShortReads)[nidx].returnRead().size()){
        (*conlen) = (*ShortReads)[nidx].returnCoordinates()[0].second + (int)(*ShortReads)[nidx].returnRead().size();
    }
    
    //std::cout << "ps: " << ps << "\n";
    for (int i = ps-2; i >= 0; i--){
        nidx = (*pathlist)[i];
        int ap1 = (*ShortReads)[nidx].returnCoordinates()[0].second;
        DC.CorrectCigarstep((*ShortReads)[(*ShortReads)[nidx].returnCoordinates()[0].first].returnCigar()[0],(*ShortReads)[nidx].returnCigar()[0], (*ShortReads)[nidx].returnCoordinates()[0].second);
        int cc = (*ShortReads)[nidx].returnCoordinates()[0].first;
        NGS nnidx = (*ShortReads)[nidx];
        (*ShortReads)[nidx].clearAldetails();
        //(*ShortReads)[nidx].returnCoordinates()[0].first
        NGS ncc = (*ShortReads)[cc];
        if((int)(*ShortReads)[cc].returnCoordinates().size() > 0){
        ap1 += (*ShortReads)[cc].returnCoordinates()[0].second + DC.returnEditedalPos();
        }
        else{
            ap1 += DC.returnEditedalPos();
        }
        (*ShortReads)[nidx].setAldetails(cc, ap1, DC.returnEditedstring(), pathname);
        (*grouping).push_back(std::make_pair(nidx, &(*ShortReads)[nidx]));
        
        if((*conlen) < (*ShortReads)[nidx].returnCoordinates()[0].second + (int)(*ShortReads)[nidx].returnRead().size()){
            (*conlen) = (*ShortReads)[nidx].returnCoordinates()[0].second + (int)(*ShortReads)[nidx].returnRead().size();
        }
        
        if(ap1 < 0){
            if(sav > ap1){
                sav = ap1;
            }
        }
    }

    sort((*grouping).begin(), (*grouping).end(), sortNGSAl);
    
    if(sav < 0){
        sav = abs(sav);
        (*conlen) = 0;
        for (int i = 0; i < (*grouping).size(); i++){
            int tc = (*(*grouping)[i].second).returnCoordinates()[0].second;
            (*(*grouping)[i].second).editCoordinates(0, tc + sav);
            if((*conlen) < (*ShortReads)[(*grouping)[i].first].returnCoordinates()[0].second + (int)(*ShortReads)[(*grouping)[i].first].returnRead().size()){
                (*conlen) = (*ShortReads)[(*grouping)[i].first].returnCoordinates()[0].second + (int)(*ShortReads)[(*grouping)[i].first].returnRead().size();
            }
            
        }
    }
    
};

ProcessReads::~ProcessReads(){};

ProcessReads::ProcessReads(std::vector<NGS*> *ShortReads1, std::vector<NGS> *ReferenveList1, std::vector<std::pair<std::string, std::string>> *Contigs1, std::vector<std::vector<int>> *gVaridx1, std::vector<NGS*> *ShortReads2, std::vector<NGS> *ReferenveList2, std::vector<std::pair<std::string, std::string>> *Contigs2, std::vector<std::vector<int>> *gVaridx2, int knn, std::string rmeth, int k, double als, bool rvmissmatch){
    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;
    DistanceMethods DM;
    NucRepresentations NR;
    std::vector<std::vector<double>> qrep, knrep;
    VPTree VPT1;
    VPTree VPT2;
    std::cout << "> Build VP-trees.\n";
    VPT1.TreeGeneration(ReferenveList1);
    VPT2.TreeGeneration(ReferenveList2);
    
    Graph G1((int) (*ShortReads1).size());
    Graph G2((int) (*ShortReads2).size());
    
    int pv = 4;
    
    int datasize = (int) (*ShortReads1).size();
    if(datasize < (int) (*ShortReads2).size()){
        datasize = (int) (*ShortReads2).size();
        rvmissmatch = false;
    }
    if(datasize  < 100){
        pv = 2;
    }
    if(datasize  < 50){
        pv = 1;
    }
    int per = 100/pv;
    int cv = (int)ceil(((double)datasize/(double)pv));
    std::cout<< "> [] 0%\n";
    
    for (int i = 0; i < datasize; i++){
        if((i%cv == 0) && i != 0){
            std::cout<< "> [" << std::string( (per * i /(cv))/2, '*')<<"] " << per * i /(cv) <<"%\n";
        }
        
        int sc = 0;
        int absc = 0;
        
        if(i < (int) (*ShortReads1).size()){
            std::vector <NGS> neighbours1;
            TreeSearch(&VPT1, (*ShortReads1).at(i), ShortReads1, &neighbours1, knn, &sc, &absc);
            if((int)neighbours1.size() > 0){
                EvalauteNeighbours(ShortReads1, i, &neighbours1,  &qrep, &knrep, rmeth, &NR, &DM);
                if((int)neighbours1.size() > 0){
                    SmithWatermanEvalaution(ShortReads1, i, &neighbours1, &aligner, &filter, &alignment, als);
                }
            }
        }
        if(i < (int) (*ShortReads2).size()){
            std::vector <NGS> neighbours2;
            sc = 0;
            absc = 0;
            TreeSearch(&VPT2, (*ShortReads2).at(i), ShortReads2, &neighbours2, knn, &sc, &absc);
            if((int)neighbours2.size() > 0){
                EvalauteNeighbours(ShortReads2, i, &neighbours2,  &qrep, &knrep, rmeth, &NR, &DM);
                if((int)neighbours2.size() > 0){
                    SmithWatermanEvalaution(ShortReads2, i, &neighbours2, &aligner, &filter, &alignment, als);
                }
            }
        }
        if(rvmissmatch == true){
            if((int)(*(*ShortReads1)[i]).returnCoordinates().size() ==  0 || (int)(*(*ShortReads2)[i]).returnCoordinates().size() ==  0){
                (*(*ShortReads1)[i]).clearAldetails();
                (*(*ShortReads2)[i]).clearAldetails();
            }
        }
        
        if(i < (int) (*ShortReads1).size()){
            for (int j = 0; j < (int)(*(*ShortReads1)[i]).returnCoordinates().size(); j++){
                G1.addEdge(i,(*(*ShortReads1)[i]).returnCoordinates()[j].first, (*(*ShortReads1)[i]).returnAlVec()[j]);
            }
            (*(*ShortReads1)[i]).clearTra();
            (*(*ShortReads1)[i]).clearRep();
        }
        
        if(i < (int) (*ShortReads2).size()){
            for (int j = 0; j < (int)(*(*ShortReads2)[i]).returnCoordinates().size(); j++){
                G2.addEdge(i,(*(*ShortReads2)[i]).returnCoordinates()[j].first, (*(*ShortReads2)[i]).returnAlVec()[j]);
            }
            (*(*ShortReads2)[i]).clearTra();
            (*(*ShortReads2)[i]).clearRep();
        }
    }
    
};


ProcessReads::ProcessReads(std::vector<NGS*> *ShortReads, std::vector<NGS> *ReferenveList, std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::vector<int>> *gVaridx, int knn, std::string rmeth, int k, double als){
    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;
    DistanceMethods DM;
    NucRepresentations NR;
    VPTree VPT;
    std::vector<std::vector<double>> qrep, knrep;
    std::cout << "> Build VP-tree.\n";
    VPT.TreeGeneration(ReferenveList);
    
    Graph G((int) ShortReads->size());
    
    std::cout << "> Link reads and build graph.\n";
    int pv = 4;
    
    if(ShortReads->size() < 100){
        pv = 2;
    }
    if(ShortReads->size() < 50){
        pv = 1;
    }
    int per = 100/pv;
    int cv = (int)ceil(((double)ShortReads->size()/(double)pv));
    std::cout<< "> [] 0%\n";
    
    for (int i = 0; i < (int) ShortReads->size(); i++){
        if((i%cv == 0) && i != 0){
            std::cout<< "> [" << std::string( (per * i /(cv))/2, '*')<<"] " << per * i /(cv) <<"%\n";
        }
        int sc = 0;
        int absc = 0;
        std::vector <NGS> neighbours;
        TreeSearch(&VPT, ShortReads->at(i), ShortReads, &neighbours, knn, &sc, &absc);
        if((int)neighbours.size() > 0){
            EvalauteNeighbours(ShortReads, i, &neighbours,  &qrep, &knrep, rmeth, &NR, &DM);
            if((int)neighbours.size() > 0){
                SmithWatermanEvalaution(ShortReads, i, &neighbours, &aligner, &filter, &alignment, als);
            }
        }
        for (int j = 0; j < (*(*ShortReads)[i]).returnCoordinates().size(); j++){
            G.addEdge(i, (*(*ShortReads)[i]).returnCoordinates()[j].first, (*(*ShortReads)[i]).returnAlVec()[j]);
        }
        (*(*ShortReads)[i]).clearTra();
        (*(*ShortReads)[i]).clearRep();
    }
    std::cout<< "> [" << std::string(50, '*')<<"] 100%\n";
    std::vector<NGS>().swap( *ReferenveList);
    std::vector<std::vector<int>> pathlist;
    std::cout << "> Identify paths in graph.\n";
    G.BFS(&pathlist);
    for (int i = 0; i < (int)pathlist.size(); i++){
        std::map<int, int> route;
        adjustpath(&route, pathlist[i][pathlist[i].size()-1], (*ShortReads)[pathlist[i][pathlist[i].size()-1]], false);
        for (int j = (int)pathlist[i].size()-2; j >= 0; j--){
            adjustpath(&route, pathlist[i][j], (*ShortReads)[pathlist[i][j]], true);
        }
        std::map<int, int>().swap(route);
    }
    
    std::vector<std::vector<double>>().swap(qrep);
    std::vector<std::vector<double>>().swap(knrep);
    std::cout <<"> Generate " << (int)pathlist.size() << " contig(s) from " << (*ShortReads).size() <<" datapoint(s).\n";
    ProcessContigs(ShortReads, &pathlist, Contigs, gVaridx);
};

void ProcessReads::ProcessContigs(std::vector<NGS*> *ShortReads, std::vector< std::vector<int>> *pathlist,  std::vector<std::pair<std::string, std::string>> *Contigs, std::vector<std::vector<int>> *gVaridx){
    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;
    std::vector<std::pair<std::string, int>> contlengths;
    for (int i = 0; i < (*pathlist).size(); i++){
        
        int conlen = std::numeric_limits<int>::min();
        std::vector<std::pair<int, NGS*>> grouping;
        //std::map<int, int> keyalig;
        std::string UID = "Contig_"+std::to_string((*Contigs).size());
        //CigarsCorrection(ShortReads , &(*pathlist)[i], &keyalig, &grouping, UID, &conlen);
        CigarsCorrection(ShortReads , &(*pathlist)[i], &grouping, UID, &conlen);
        Kindel Ki(UID, conlen);
        
        for (int j = 0; j < grouping.size(); j++){
            Ki.addRead((*grouping[j].second).returnCoordinates()[0].second, (*grouping[j].second).returnCigar()[0], (*grouping[j].second).returnRead());
        }
        
        std::vector<int> dp;
        Ki.generateConsensus(true, &dp);
        
        std::string iniContig = Ki.getConsensus();
        std::string Cname = Ki.getContigName();
        Ki.deleteRefMap();

        for (int j = 0; j < grouping.size(); j++){
            //
            int  sp = (*grouping[j].second).returnCoordinates()[0].second - dp[(*grouping[j].second).returnCoordinates()[0].second];
            
            if(sp - 10 > 0 ){
                sp -= 10;
            }
            else{
                sp = 0;
            }
            int ep = sp + (int)(*grouping[j].second).returnRead().size();
            if(ep + 10 < (int)iniContig.size()){
                ep += 10;
            }
            else{
                ep = (int)iniContig.size();
            }
            
            int wsize = ep - sp;
            aligner.Align((*grouping[j].second).returnRead().c_str(),iniContig.substr(sp, wsize).c_str(), wsize, filter, &alignment);
            /*std::string tcigar = alignment.cigar_string;
            std::size_t cigM, cigS, cigI, cigD;
            cigM = tcigar.find("M");
            cigS = tcigar.find("S");
            cigI = tcigar.find("I");
            cigD = tcigar.find("D");
            if(cigS < cigM && cigS < cigI && cigS < cigI &&  cigS < cigD){
                sp -= stoi(tcigar.substr(0,cigS));
                tcigar = tcigar.substr(0,cigS)+ "M" + tcigar.substr(cigS+1);
            }
            
            if(tcigar[tcigar.size()-1] == 'S'){
                tcigar[tcigar.size()-1] = 'M';
            }*/
            (*grouping[j].second).clearAldetails();
            (*grouping[j].second).setAldetails(0, sp + alignment.ref_begin,alignment.cigar_string, Cname);
        }
        
        sort(grouping.begin(), grouping.end(), sortNGSAl);
        
        if((*grouping[0].second).returnCoordinates()[0].second < 0){
            conlen = std::numeric_limits<int>::min();
            int incval = abs((*grouping[0].second).returnCoordinates()[0].second)+1;
            for (int s = 0; s < grouping.size(); s++){
                (*grouping[s].second).editCoordinates(0, (*grouping[s].second).returnCoordinates()[0].second + incval);
                if(conlen < (*grouping[s].second).returnCoordinates()[0].second +  (int)(*grouping[s].second).returnRead().size()){
                    conlen = (*grouping[s].second).returnCoordinates()[0].second +  (int)(*grouping[s].second).returnRead().size();
                }
            }
        }
        
        Ki.deleteRefMap();
        Ki.generateKindle(UID, conlen);
        
        int gidx = (int)grouping.size() - 1;
        std::vector<int>idx(gidx+1, 0);
        while(!grouping.empty()){
            Ki.addRead((*grouping[gidx].second).returnCoordinates()[0].second, (*grouping[gidx].second).returnCigar()[0], (*grouping[gidx].second).returnRead());
            idx[gidx] = (*(*ShortReads)[grouping[gidx].first]).returnIdx();
            gidx--;
            grouping.pop_back();
        }
        Ki.generateConsensus(true);
        Contigs->push_back( std::make_pair(Ki.getContigName(), Ki.getConsensus()));
        (*gVaridx).push_back(idx);
        std::vector<int>().swap(dp);
        std::vector<int>().swap(idx);
        std::vector<std::pair<int, NGS*>>().swap(grouping);
    }
};

void ProcessReads::CigarsCorrection(std::vector<NGS*> *ShortReads, std::vector<int> *pathlist, std::vector<std::pair<int, NGS*>> *grouping, std::string pathname, int *conlen){
    DecompCigar DC;
    int ps = (int)(*pathlist).size();
    int sav = 0;
    int nidx = (*pathlist)[ps - 1];
    (*(*ShortReads)[nidx]).clearAldetails();
    (*(*ShortReads)[nidx]).setAldetails(nidx, 0, (std::to_string((int)(*(*ShortReads)[nidx]).returnRead().size()) + "M"), pathname);
    (*grouping).push_back(std::make_pair(nidx, (*ShortReads)[nidx]));
    
    if((*conlen) < (*(*ShortReads)[nidx]).returnCoordinates()[0].second + (int)(*(*ShortReads)[nidx]).returnRead().size()){
        (*conlen) = (*(*ShortReads)[nidx]).returnCoordinates()[0].second + (int)(*(*ShortReads)[nidx]).returnRead().size();
    }

    for (int i = ps-2; i >= 0; i--){
        nidx = (*pathlist)[i];
        int ap1 = (*(*ShortReads)[nidx]).returnCoordinates()[0].second;
        DC.CorrectCigarstep((*(*ShortReads)[(*(*ShortReads)[nidx]).returnCoordinates()[0].first]).returnCigar()[0],(*(*ShortReads)[nidx]).returnCigar()[0], (*(*ShortReads)[nidx]).returnCoordinates()[0].second);
        int cc = (*(*ShortReads)[nidx]).returnCoordinates()[0].first;
        (*(*ShortReads)[nidx]).clearAldetails();
        if((int)(*(*ShortReads)[cc]).returnCoordinates().size() > 0){
            ap1 +=(*(*ShortReads)[cc]).returnCoordinates()[0].second + DC.returnEditedalPos();
        }
        else{
            ap1 += DC.returnEditedalPos();
        }
        (*(*ShortReads)[nidx]).setAldetails(cc, ap1, DC.returnEditedstring(), pathname);
        (*grouping).push_back(std::make_pair(nidx, (*ShortReads)[nidx]));
        
        if((*conlen) < (*(*ShortReads)[nidx]).returnCoordinates()[0].second + (int)(*(*ShortReads)[nidx]).returnRead().size()){
            (*conlen) = (*(*ShortReads)[nidx]).returnCoordinates()[0].second + (int)(*(*ShortReads)[nidx]).returnRead().size();
        }
        if(ap1 < 0){
            if(sav > ap1){
                sav = ap1;
            }
        }
    }
    
    sort((*grouping).begin(), (*grouping).end(), sortNGSAl);
    
    if(sav < 0){
        sav = abs(sav);
        (*conlen) = 0;
        for (int i = 0; i < (*grouping).size(); i++){
            int tc = (*(*grouping)[i].second).returnCoordinates()[0].second;
            (*(*grouping)[i].second).editCoordinates(0, tc + sav);
           
            if((*conlen) < (*(*ShortReads)[(*grouping)[i].first]).returnCoordinates()[0].second + (int)(*(*ShortReads)[(*grouping)[i].first]).returnRead().size()){
                (*conlen) = (*(*ShortReads)[(*grouping)[i].first]).returnCoordinates()[0].second + (int)(*(*ShortReads)[(*grouping)[i].first]).returnRead().size();
            }
            
        }
    }
    
};
