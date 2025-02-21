#ifndef CAREA_H
#define	CAREA_H

#include "main.h"
#include "CParam.h"
#include "CSet.h"
#include "CImage.h"
#include "CInterferogram.h"

class CArea
{
public:
    vector <CSet> DSET;
    vector <CImage> SLC;
    vector <CInterferogram *> pInSAR;
    
    int DIM;
    int ROWS;
    int COLUMNS;
    float *TIME_MATRIX;
    
    float *RATE_LOS;
    float *RATE_EW;
    float *RATE_UD;

    float *RATE_ERROR_LOS;
    float *RATE_ERROR_EW;
    float *RATE_ERROR_UD;

    float *TOPO_CORRECTION;
    
    float *NORM_X;
    float *NORM_AXY;
    
    float *ZSCORE_MASK;
    
    CArea(CParam *par);
    void ComputeDIM();
    void MakeSLC();
    void MakeTimeMatrix1D();
    void MakeTimeMatrix2D();
    void ReadInterferograms();
    void CalibrateInterferograms();
    void ComputeZScoreMask();
    void Inversion();
    void Inversion_old();
    void PostProcessing();
    void TemporalGaussianFilter();
    void WriteResultsToDisk();
    void ComputeLinearRate();
    void InteractiveMode();
    int sgelss(float *a, int m, int n, float *s, int lwork, float *work, float *b, int nrhs);
    void linreg(float *x, float *y, int length, float &aa, float &bb, float &sigmaa, float &sigmab, float &chi2);
    
private:
    CParam *PAR;

};

#endif	/* CAREA_H */

