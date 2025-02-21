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
    
    int LWORK;
     
    bool *MASK;

    int DIM;
    int ROWS;
    int COLUMNS;
    float *TIME_MATRIX;

    float *RATE_LOS;
    float *RATE_NS;
    float *RATE_EW;
    float *RATE_UD;
    float *RATE_UD1;
    
    float *RATE_STD_LOS;
    float *RATE_STD_NS;
    float *RATE_STD_EW;
    float *RATE_STD_UD;
    float *RATE_STD_UD1;
    
    float *RATE_R2_LOS;
    float *RATE_R2_NS;
    float *RATE_R2_EW;
    float *RATE_R2_UD;
    float *RATE_R2_UD1;

    float *NORM_X;
    float *NORM_AXY;
    float *RANK;
    float *COND_NUM;
    
    float *ZSCORE_MASK;
    
    CArea(CParam *par);
    void ComputeDIM();
    void MakeTimeMatrix();
    
    void ReadInterferograms();
    void CalibrateInterferograms();
    void ComputeZScoreMask();
    void Inversion();
    void InversionByLine();
    void InversionByPixel();
    void PostProcessing();
    void WriteResultsToDisk();
    void ComputeLinearRate();
    void InteractiveMode();
    void sgelss_query_lwork(int m, int n, int nrhs);
    inline int sgelss(float *a, int m, int n, float *b, int nrhs, int &rank_num, float &cond_num);
    void linreg(float *x, float *y, int length, float &aa, float &bb, float &sigmaa, float &sigmab, float &r2);
    void write_float_data(string name, float *data);
    void ReadBinaryFile(string filename, float *dataout);
        
private:
    CParam *PAR;

};

#endif	/* CAREA_H */