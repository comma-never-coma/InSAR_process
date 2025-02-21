#ifndef CINTERFEROGRAM_H
#define	CINTERFEROGRAM_H

#include "main.h"
#include "CParam.h"
#include "CImage.h"

class CInterferogram
{
public:
    string NAME;
    float BASELINE;
    CImage *pM;
    CImage *pS;
    float SPAN;
    float BNDCOR;
    float OFFSET;
    float *DATA;
    float MEAN;
    float STD;    

    CInterferogram(CParam *par, const int snum, const string time, const string filein, const string baseline, const string master, const string slave);
    void ReadInterferogram();
    void CalibrateInterferogram();
    void ComputeStatistics();
    int SNUM;
    
private:
    CParam *PAR;
};

#endif	/* CINTERFEROGRAM_H */

