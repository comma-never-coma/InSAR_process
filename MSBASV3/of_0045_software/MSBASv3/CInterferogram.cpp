#include "CInterferogram.h"

CInterferogram::CInterferogram(CParam *par, const int snum, const string time, const string name, const string baseline, const string m, const string s)
{
    PAR = par;
    SNUM = snum;
    NAME = name;
    BASELINE = atof(baseline.c_str());

    pM = new CImage(PAR, SNUM, m, time);
    pS = new CImage(PAR, SNUM, s, time);

    SPAN = pS->DATE - pM->DATE;
    BNDCOR = 1;
    OFFSET = 0;
    DATA = NULL; // do not allocate memory now because some interferograms will be deleted later
    MEAN = 0;
    STD = 0;
};

void CInterferogram::ReadInterferogram()
{
    // allocate memory for interferogram
    float_mem_alloc(DATA, PAR->FWIDTH * PAR->LENGTH);
    
    PAR->read_float_data(NAME, DATA);
    
    for (int i = 0; i < PAR->WIDTH; i++)
        for (int j = 0; j < PAR->LENGTH; j++)
            DATA[i + PAR->WIDTH * j] = DATA[i + PAR->WIDTH * j] * BNDCOR;

    ComputeStatistics();
};

void CInterferogram::CalibrateInterferogram()
{
    float ave = 0;
    int num = 0, i = 0, j = 0, k = 0;

    if (PAR->C_FLAG == 0)
    {
        //  do nothing because calibrated data is provided
    }
    else if ((PAR->C_FLAG > 0)&&(PAR->C_FLAG < 10)) // set average of reference regions to zero, use up to 9 reference regions
    {
        for (k = 0; k < PAR->C_FLAG; k++)
        {
            for (i = (PAR->CPOS[k] - PAR->CSIZE); i <= (PAR->CPOS[k] + PAR->CSIZE); i++)
                for (j = (PAR->LPOS[k] - PAR->LSIZE); j <= (PAR->LPOS[k] + PAR->LSIZE); j++)
                {
                    if ((i >= 0)&&(j >= 0)&&(i < PAR->WIDTH)&&(j < PAR->LENGTH) && DATA[i + j * PAR->WIDTH] != 0)
                    {
                        ave = ave + DATA[i + j * PAR->WIDTH];
                        num++;
                    }
                }
        }
        if (num > 0)
            ave = ave / num;
        else
            throw std::invalid_argument(ERROR("no valid points for calibration, check location of reference region(s) in C_FLAG"));

        for (i = 0; i < PAR->WIDTH; i++)
            for (j = 0; j < PAR->LENGTH; j++)
                if (DATA[i + j * PAR->WIDTH] != 0) DATA[i + j * PAR->WIDTH] = DATA[i + j * PAR->WIDTH] - ave;
    }
    else if (PAR->C_FLAG == 10) // set average of interferogram to zero
    {
        for (i = 0; i < PAR->WIDTH; i++)
            for (j = 0; j < PAR->LENGTH; j++)
                if (DATA[i + j * PAR->WIDTH] != 0)
                {
                    ave = ave + DATA[i + j * PAR->WIDTH];
                    num++;
                }

        if (num > 0)
            ave = ave / num;
        else
            throw std::invalid_argument(ERROR("no valid points for calibration, check location of reference region(s) in C_FLAG"));

        for (i = 0; i < PAR->WIDTH; i++)
            for (j = 0; j < PAR->LENGTH; j++)
                if (DATA[i + j * PAR->WIDTH] != 0) DATA[i + j * PAR->WIDTH] = DATA[i + j * PAR->WIDTH] - ave;
    }
    else
        throw std::invalid_argument(ERROR("incorrect value in C_FLAG"));
    
    OFFSET = ave;
};

void CInterferogram::ComputeStatistics()
{
    int num = 0, i = 0, j = 0, k = 0;

    
    for (i = 0; i < PAR->WIDTH; i++)
        for (j = 0; j < PAR->LENGTH; j++)
            if (DATA[i + j * PAR->WIDTH] != 0)
            {
                MEAN = MEAN + DATA[i + j * PAR->WIDTH];
                num++;
            }

    if (num > 0)
            MEAN = MEAN/num;
        else
            
            throw std::invalid_argument(ERROR("failed computing mean value - no valid pixels found"));

    num = 0;
    
    for (i = 0; i < PAR->WIDTH; i++)
        for (j = 0; j < PAR->LENGTH; j++)
            if (DATA[i + j * PAR->WIDTH] != 0)
            {
                STD = STD + (DATA[i + j * PAR->WIDTH]-MEAN)*(DATA[i + j * PAR->WIDTH]-MEAN);
                num++;
            }

    if (num > 0)
        STD = sqrt(STD/(num-1));
    else
        throw std::invalid_argument(ERROR("failed computing standard deviation value - no valid pixels found"));
};


