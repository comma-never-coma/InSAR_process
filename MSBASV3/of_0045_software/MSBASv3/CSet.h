#ifndef CDSET_H
#define	CDSET_H

#include "main.h"
#include "CParam.h"
#include "CInterferogram.h"
#include "CImage.h"

class CSet
{
public:
    float AZIMUTH;
    float INCIDENCE;
    float S[3]; // in order - north, east, up
    string TIME;
    string NAME;
    float *TIME_MATRIX;
    
    vector <CInterferogram> InSAR;
    vector <CImage> SLC;
    int SNUM;
       
    CSet(CParam *par, const int snum);
    void MakeInSARfromFile();
    void MakeSLC();
    void FilterInSAR(const vector <CImage> &slc);
    void MakeTimeMatrix1D(const vector <CImage> &slc);
    void MakeTimeMatrix2D(const vector <CImage> &slc);
     
private:
    CParam *PAR;
};

#endif	/* CDSET_H */
