#ifndef CPARAM_H
#define	CPARAM_H

#include "main.h"

class CParam
{
public:
    string NAME;
    int FORMAT;
    int FORMAT_INT_FLAG;
    int FWIDTH, FLENGTH;
    int CSTART, CSTOP, LSTART, LSTOP, WIDTH, LENGTH;
    int C_FLAG, CPOS[10], LPOS[10], CSIZE, LSIZE;
    int R_FLAG;
    float R_LAMBDA;
    bool T_FLAG;
    float TAV_FLAG;
    int I_FLAG;
    string INAME;
    GDALDataset *poSrcDS;
    
    struct SSet
    {
        string TIME;
        string AZIMUTH;
        string INCIDENCE;
        string NAME;
        SSet(string time, string azimuth, string incidence, string name)
        {
            TIME = time;
            AZIMUTH = azimuth;
            INCIDENCE = incidence;
            NAME = name;
        }
    };
    
    vector <SSet> DSET;
    
    CParam(string name);
    ~CParam(); 
    void ReadParFile();
    void ProcessPar(string parameter, string value);
    void ValidatePar();
    void read_float_data(string name, float *data);
    void write_float_data(string name, float *data);
    
private:

};

#endif	/* CPARAM_H */

