#ifndef CINTERFEROGRAM_H
#define	CINTERFEROGRAM_H

#include "main.h"
#include "CParam.h"
#include "CImage.h"

class CInterferogram
{
public:
    CImage *pM;
    CImage *pS;

    float BASELINE;
    float SPAN;
    float BNDCOR;
    float OFFSET;
    float *DATA;
    float MEAN;
    float STD;    

    string NAME;

    
    CInterferogram(CParam *par, const int setn, const string time, const string filein, const string baseline, const string master, const string slave);
    void ReadInterferogram();
    void CalibrateInterferogram();
    void ComputeStatistics();
    int SETN;
    
private:
    CParam *PAR;
};

#endif	/* CINTERFEROGRAM_H */

