#ifndef CDSET_H
#define	CDSET_H

#include "main.h"
#include "CParam.h"
#include "CInterferogram.h"
#include "CImage.h"

class CSet
{
public:
    float TYPE; //0-range offsets or insar, 1- azimuth offsets
    float AZIMUTH;
    float INCIDENCE;
    float S[3]; // in order - north, east, up
    string TIME;
    string NAME;
    float *TIME_MATRIX;

    float *LV_THETA;
    float *LV_PHI;
    
    vector <CInterferogram> InSAR;
    vector <CImage> SLC;
    int SETN;
      
    CSet(CParam *par, const int setn);
    void ReadInSARParFile();
    void MakeSLC();
    int RemoveDisconnectedInSAR();
    void ApplyBoundaryCorrectionInSAR(const vector <CImage> &slc);

    void MakeTimeMatrix(const vector <CImage> &slc);
    
private:
    CParam *PAR;
};

#endif	/* CDSET_H */

     