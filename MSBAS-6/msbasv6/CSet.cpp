#include "CSet.h"

CSet::CSet(CParam *par, const int snum)
{
    PAR = par;
    SETN = snum;
    
    TIME = PAR->DSET[SETN].TIME;
    TYPE = atof(PAR->DSET[SETN].TYPE.c_str());
    AZIMUTH = atof(PAR->DSET[SETN].AZIMUTH.c_str());
    INCIDENCE = atof(PAR->DSET[SETN].INCIDENCE.c_str());
       
    
    NAME = PAR->DSET[SETN].NAME;
    TIME_MATRIX = NULL;
    LV_THETA = NULL;
    LV_PHI = NULL;
    
    if (PAR->DSET[SETN].LV_THETA_FILE.size() != 0) 
    {
        float_mem_alloc(LV_THETA, (PAR->WIDTH)*(PAR->LENGTH));
        par->read_float_data(PAR->DSET[SETN].LV_THETA_FILE, LV_THETA);
    }

    if (PAR->DSET[SETN].LV_PHI_FILE.size() != 0) 
    {
        float_mem_alloc(LV_PHI, (PAR->WIDTH)*(PAR->LENGTH));
        par->read_float_data(PAR->DSET[SETN].LV_PHI_FILE, LV_PHI);
    }
    
    // range
    if (TYPE == 0)
    {
        S[0] = sin(AZIMUTH * PI / 180.) * sin(INCIDENCE * PI / 180.);
        S[1] = -1 * cos(AZIMUTH * PI / 180.) * sin(INCIDENCE * PI / 180.);
        S[2] = cos(INCIDENCE * PI / 180.);
    }
    // azimuth
    else if (TYPE == 1)
    {
        S[0] = cos(AZIMUTH * PI / 180.);
        S[1] = sin(AZIMUTH * PI / 180.);
        S[2] = 0;
    }
    else 
        throw std::invalid_argument(ERROR("detected incorrect set type"));
   
    WriteLog("\nreading set " + i2s(SETN) + " file " + NAME + " directional cosines [" + f2s(S[0]) + ", " + f2s(S[1]) + ", " + f2s(S[2]) + "]:\n");
    
    ReadInSARParFile();
};

void CSet::ReadInSARParFile()
{
    InSAR.clear();
    
    FILE * pFile;
    pFile = fopen(NAME.c_str(), "r");
    if (pFile == NULL)
        throw std::invalid_argument(ERROR("can not open set file " + NAME));
    
    char v1[MAX_LENGTH_FILENAME], v2[MAX_LENGTH_FILENAME], v3[MAX_LENGTH_FILENAME], v4[MAX_LENGTH_FILENAME];
    int k=0;
    while (fscanf(pFile, "%s %s %s %s", v1, v2, v3, v4) == 4)
    {
        string val1(v1), val2(v2), val3(v3), val4(v4);
        WriteLog(i2sf(k) + ": " + val1 + " " + val2 + " " + val3 + " " + val4);
        InSAR.push_back(CInterferogram(PAR, SETN, TIME, val1, val2, val3, val4));
        k++;
    }
    fclose(pFile);

    if(InSAR.empty()) throw std::invalid_argument(ERROR("no interferograms are available for processing"));
    
    // recursively remove disconnected interferograms
    do MakeSLC(); while (RemoveDisconnectedInSAR()>0);

};

void CSet::MakeSLC()
{
    if(InSAR.empty()) throw std::invalid_argument(ERROR("no interferograms are available for processing"));
    
    SLC.clear();
    
    for (int i = 0; i < InSAR.size(); i++)
    {
        bool used = false;
        for (int j = 0; j < SLC.size(); j++)
        {
            if (InSAR[i].pM->DATE == SLC[j].DATE)
            {
                SLC[j].MCOUNT++;
                used = true;
            }
        }
        if (!used)
        {
            SLC.push_back(CImage(*InSAR[i].pM));
            SLC[SLC.size()-1].MCOUNT++;
        }    
    }
    
    for (int i = 0; i < InSAR.size(); i++)
    {
        bool used = false;
        for (int j = 0; j < SLC.size(); j++)
        {
            if (InSAR[i].pS->DATE == SLC[j].DATE)
            {
                SLC[j].SCOUNT++;
                used = true;
            }
        }
        if (!used)
        {
            SLC.push_back(CImage(*InSAR[i].pS));
            SLC[SLC.size()-1].SCOUNT++;
        }    
    }
 
    sort(SLC.begin(), SLC.end());
    SLC.erase(unique(SLC.begin(), SLC.end()), SLC.end());
    
    if (SLC.empty()) throw std::invalid_argument(ERROR("no SLCs are available for processing"));

    WriteLog("\ndetected " + i2s(SLC.size()) + " SLCs [NAME DATE SETN MCOUNT SCOUNT]:\n");
    for (int i=0; i< SLC.size(); i++) WriteLog(i2sf(i) + ": " + SLC[i].NAME + " " + f2s(SLC[i].DATE) + " " + i2sf(SLC[i].SETN) + " " + i2sf(SLC[i].MCOUNT) + " " + i2sf(SLC[i].SCOUNT));    
    
};

int CSet::RemoveDisconnectedInSAR()
{
    int removed = 0;
    
    // uncomment for removing disconnected interferograms, excluding first and last SLC
    /*
    for(int i = InSAR.size() - 1; i >= 0; i--)
        for (int j = 1; j < (SLC.size() - 1); j++)
            if (SLC[j].MCOUNT == 0 && InSAR[i].pS->DATE == SLC[j].DATE)
            {
                if (removed == 0) WriteLog("\nremoving disconnected interferograms:");
                WriteLog(i2sf(i) + ": " + InSAR[i].NAME);
                InSAR.erase(InSAR.begin() + i);
                removed++;
            }
            else if (SLC[j].SCOUNT == 0 && InSAR[i].pM->DATE == SLC[j].DATE)
            {
                if (removed == 0) WriteLog("\nremoving disconnected interferograms:");
                WriteLog(i2sf(i) + ": " + InSAR[i].NAME);
                InSAR.erase(InSAR.begin() + i);
                removed++;
            }
    
    if(SLC.empty()) throw std::invalid_argument(ERROR("no SLCs are available for processing"));
    if(InSAR.empty()) throw std::invalid_argument(ERROR("no interferograms are available for processing"));
    */
    return removed;
}

void CSet::ApplyBoundaryCorrectionInSAR(const vector <CImage> &slc)
{
    // vector is sorted
    for (int i=1; i<slc.size(); i++) if(slc[i-1].DATE >= slc[i].DATE) throw std::invalid_argument(ERROR("slc[i-1].DATE >= slc[i].DATE"));
    
    // exclude interferograms outside of the common window
    WriteLog("\nremoving interferograms outside of common time span (if required):");
    for(int i = InSAR.size() - 1; i >= 0; i--)
        if (InSAR[i].pS->DATE <= slc[0].DATE || InSAR[i].pM->DATE >= slc[slc.size()-1].DATE)
        {
            WriteLog(i2sf(i) + ": " + InSAR[i].NAME);
            InSAR.erase(InSAR.begin() + i); 
        }
    
    if (InSAR.empty()) throw std::invalid_argument(ERROR("no interferograms are available for processing"));
 
    // compute boundary correction
    WriteLog("\ncomputing and applying boundary correction (if required):");
    for (int i=0; i<InSAR.size(); i++)
        {
            InSAR[i].BNDCOR = 1; // reinitialize in case of repeated calls
            if (InSAR[i].pM->DATE < slc[0].DATE)
            {
                InSAR[i].BNDCOR = InSAR[i].BNDCOR*(InSAR[i].pS->DATE - slc[0].DATE) / (InSAR[i].pS->DATE - InSAR[i].pM->DATE);
                InSAR[i].pM->DATE = slc[0].DATE;
                InSAR[i].pM->NAME = slc[0].NAME;
                WriteLog(i2sf(i) + ": " + InSAR[i].NAME + " " + f2s(InSAR[i].BNDCOR) + " " + InSAR[i].pM->NAME);
            }
            if (InSAR[i].pS->DATE > slc[slc.size()-1].DATE)
            {
                InSAR[i].BNDCOR = InSAR[i].BNDCOR*(slc[slc.size()-1].DATE - InSAR[i].pM->DATE) / (InSAR[i].pS->DATE - InSAR[i].pM->DATE);
                InSAR[i].pS->DATE = slc[slc.size()-1].DATE;
                InSAR[i].pS->NAME = slc[slc.size()-1].NAME;
                WriteLog(i2sf(i) + ": " + InSAR[i].NAME + " " + f2s(InSAR[i].BNDCOR) + " " + InSAR[i].pS->NAME);
            }
        }
    
    MakeSLC();
};

void CSet::MakeTimeMatrix(const vector <CImage> &slc)
{
    if(slc.empty()) throw std::invalid_argument(ERROR("no SLCs are available for processing"));
    if(InSAR.empty()) throw std::invalid_argument(ERROR("no interferograms are available for processing"));
    
    float_mem_alloc(TIME_MATRIX, (slc.size()-1)*InSAR.size());

    for (int j = 0; j < InSAR.size(); j++)
        for (int i = 0; i < (slc.size() - 1); i++)
            if ((slc[i].DATE >= InSAR[j].pM->DATE)&&(slc[i].DATE < InSAR[j].pS->DATE))
                TIME_MATRIX[ i + (slc.size() - 1) * j] = (slc[i+1].DATE - slc[i].DATE);
};