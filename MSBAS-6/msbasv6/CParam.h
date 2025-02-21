#ifndef CPARAM_H
#define	CPARAM_H

#include "main.h"

class CParam
{
public:
    string NAME;
    int FORMAT;
    int FORMAT_INT_RADIUS;
    int FWIDTH, FLENGTH;
    int CSTART, CSTOP, LSTART, LSTOP, WIDTH, LENGTH;
    int C_FLAG, CPOS[10], LPOS[10], CSIZE, LSIZE;
    int R_FLAG;
    float R_LAMBDA;
    int I_FLAG;
    string I_FLAG_FILE;
    int V_FLAG;
    int D_FLAG;
    
    float *DEM, *DDNS, *DDEW;
    string DEM_FILE, DDNS_FILE, DDEW_FILE;
    
    GDALDataset *poSrcDS;
    
    struct SSet
    {
        string TYPE; //0-range offsets or insar, 1- azimuth offsets
        string TIME;
        string AZIMUTH;
        string INCIDENCE;
        string NAME;
        string LV_THETA_FILE;
        string LV_PHI_FILE;

        SSet(string type, string time, string azimuth, string incidence, string name, string lv_theta_file, string lv_phi_file)
        {
            TYPE=type;
            TIME = time;
            AZIMUTH = azimuth;
            INCIDENCE = incidence;
            NAME = name;
            LV_THETA_FILE=lv_theta_file;
            LV_PHI_FILE=lv_phi_file;
        }
    };
    
    vector <SSet> DSET;
    
    CParam();
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

