#include "CSet.h"

CSet::CSet(CParam *par, const int snum)
{
    PAR = par;
    SNUM = snum;
    
    TIME = PAR->DSET[SNUM].TIME;
    AZIMUTH = atof(PAR->DSET[SNUM].AZIMUTH.c_str());
    INCIDENCE = atof(PAR->DSET[SNUM].INCIDENCE.c_str());
    NAME = PAR->DSET[SNUM].NAME;
    TIME_MATRIX = NULL;
    
    S[0] = sin(AZIMUTH * PI / 180.) * sin(INCIDENCE * PI / 180.);
    S[1] = -1 * cos(AZIMUTH * PI / 180.) * sin(INCIDENCE * PI / 180.);
    S[2] = cos(INCIDENCE * PI / 180.);
    
    WriteLog("\nreading set " + i2s(SNUM) + " file " + NAME + " directional cosines [" + f2s(S[0]) + ", " + f2s(S[1]) + ", " + f2s(S[2]) + "]:\n");
    
    MakeInSARfromFile();
};

void CSet::MakeInSARfromFile()
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
        WriteLog(i2s(k) + ": " + val1 + " " + val2 + " " + val3 + " " + val4);
        InSAR.push_back(CInterferogram(PAR, SNUM, TIME, val1, val2, val3, val4));
        k++;
    }
    fclose(pFile);

    if(InSAR.empty()) throw std::invalid_argument(ERROR("no interferograms are available for processing"));
    
    MakeSLC();
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

    WriteLog("\ndetected " + i2s(SLC.size()) + " SLCs [NAME DATE SNUM MCOUNT SCOUNT]:");
    for (int i=0; i< SLC.size(); i++) WriteLog(i2sf(i) + ": " + SLC[i].NAME + " " + f2s(SLC[i].DATE) + " " + i2sf(SLC[i].SNUM) + " " + i2sf(SLC[i].MCOUNT) + " " + i2sf(SLC[i].SCOUNT));    

};

void CSet::FilterInSAR(const vector <CImage> &slc)
{
    if(slc.empty()) throw std::invalid_argument(ERROR("no SLCs are available for processing"));
    
    // is sorted
    for (int i=1; i<slc.size(); i++) if(slc[i-1].DATE >= slc[i].DATE) throw std::invalid_argument(ERROR("slc[i-1].DATE >= slc[i].DATE"));
    
    if(InSAR.empty()) throw std::invalid_argument(ERROR("no interferograms are available for processing"));

    // exclude interferograms outside of common span
    WriteLog("excluded interferograms:");
    for(int i = InSAR.size() - 1; i >= 0; i--)
        if (InSAR[i].pS->DATE <= slc[0].DATE || InSAR[i].pM->DATE >= slc[slc.size()-1].DATE)
        {
            WriteLog(i2sf(i) + ": " + InSAR[i].NAME);
            InSAR.erase(InSAR.begin() + i); 
        }

    if (InSAR.empty()) throw std::invalid_argument(ERROR("no interferograms are available for processing"));
 
    // compute boundary correction
    WriteLog("computed boundary correction:");
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

void CSet::MakeTimeMatrix1D(const vector <CImage> &slc)
{
    if(slc.empty()) throw std::invalid_argument(ERROR("no SLCs are available for processing"));
    if(InSAR.empty()) throw std::invalid_argument(ERROR("no interferograms are available for processing"));
    
    float_mem_alloc(TIME_MATRIX, (slc.size()-1)*InSAR.size());

    for (int j = 0; j < InSAR.size(); j++)
        for (int i = 0; i < (slc.size() - 1); i++)
            if ((slc[i].DATE >= InSAR[j].pM->DATE)&&(slc[i].DATE < InSAR[j].pS->DATE))
                TIME_MATRIX[ i + (slc.size() - 1) * j] = (slc[i+1].DATE - slc[i].DATE);
};

void CSet::MakeTimeMatrix2D(const vector <CImage> &slc)
{
    if(slc.empty()) throw std::invalid_argument(ERROR("no interferograms are available for processing"));
    if(InSAR.empty()) throw std::invalid_argument(ERROR("no SLCs are available for processing"));
    
    float_mem_alloc(TIME_MATRIX, 2*(slc.size()-1)*InSAR.size());

    for (int j = 0; j < InSAR.size(); j++)
        for (int i = 0; i < (slc.size()-1); i++)
            if ((slc[i].DATE >= InSAR[j].pM->DATE)&&(slc[i].DATE < InSAR[j].pS->DATE))
            {
                TIME_MATRIX[2 * i + 2 * (slc.size() - 1) * j] = (slc[i+1].DATE - slc[i].DATE) * S[1];
                TIME_MATRIX[2 * i + 1 + 2 * (slc.size() - 1) * j] = (slc[i+1].DATE - slc[i].DATE) * S[2];
            }
};