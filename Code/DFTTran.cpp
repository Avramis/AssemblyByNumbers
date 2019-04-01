//
//  DFTTran.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "DFTTran.hpp"
void DFTTran::fourier(std::vector<double> *data, int nn){
    int i, n, m, j, mmax, istep, isign;
    isign = -1;
    double wtemp, wr, wpr, wpi, wi, theta, tempr, tempi;
    n = nn<<1;
    j = 1;
    for (i=1; i < n; i+=2) {
        if (j>i) {
            double tempswapj, tempswapi;
            tempswapj = (*data)[j-1];
            tempswapi = (*data)[i-1];
            (*data)[i-1] = tempswapj;
            (*data)[j-1] = tempswapi;
            tempswapj = (*data)[j];
            tempswapi = (*data)[i];
            (*data)[i] = tempswapj;
            (*data)[j] = tempswapi;
        }
        m = nn;
        while (m>=2 && j>m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    
    mmax=2;
    while (n>mmax) {
        istep = mmax<<1;
        theta = (isign)*(2.0*M_PI/((double)mmax));
        wtemp = sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m=1; m < mmax; m += 2) {
            for (i=m; i <= n; i += istep) {
                j=i+mmax;
                tempr = wr * (*data)[j-1] - wi * (*data)[j];
                tempi = wr * (*data)[j] + wi * (*data)[j-1];
                
                
                (*data)[j-1] = (*data)[i-1] - tempr;
                (*data)[j] = (*data)[i] - tempi;
                (*data)[i-1] += tempr;
                (*data)[i] += tempi;
            }
            wtemp=wr;
            wr += wr*wpr - wi*wpi;
            wi += wi*wpr + wtemp*wpi;
        }
        mmax=istep;
    }
};



void DFTTran::generateTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int st, int len, int l){
    int X = st+len;
    if(X > (int)rep->at(0).size()){
        X = (int)rep->at(0).size();
    }
    if(len <= l){
        l = X - st;
        if(l < 1 ){
            l  = 1;
        }
    }
    std::vector<std::vector<double>>((int)rep->size(),std::vector<double>(l * 2, 0.0)).swap(*tran);
    for (int i = 0; i < (int)rep->size(); i++){
        //tran->push_back(std::vector<double>(l * 2, 0.00));
        int np2,
        nn;
        nn = X - st;
        double snn = sqrt((double)nn);
        np2 = (int) (ceil(log2((double)nn)));
        if (nn !=  pow(2,np2 )){
            nn = (int)(pow(2,np2 ));
        }
        std::vector<double> ffttran(nn*2, 0.0);
        int c = 0;
        for (int j = st; j < X; j++){
            if(j < X){
                ffttran[c*2] = ((*rep)[i][j]) / snn;
            }
            else{
                ffttran[c*2] = (0)/snn;
            }
            ffttran[(c*2)+1] = (0)/snn;
            c++;
        }
        DFTTran::fourier(&ffttran,nn);
        //tran->push_back(std::vector<double>(l*2,0.0));
        for (int j = 0; j < l*2; j++){
            (*tran)[i][j] = ffttran[j];
        }
        std::vector<double>().swap(ffttran);
    }
};



DFTTran::DFTTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int l){
    generateTran(rep, tran, 0, (int)rep->at(0).size(), l);
};

DFTTran::~DFTTran(){};





DFTTran::DFTTran(std::vector<std::vector<double>> *rep, std::vector<std::vector<double>> *tran, int st, int len, int l){
    generateTran(rep, tran, st, len, l);
};

