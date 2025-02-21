#include "CArea.h"

CArea::CArea(CParam *par)
{
    PAR = par; 
    DIM = 0;
    ROWS = COLUMNS = 0;
    TIME_MATRIX = NULL;
    RATE_LOS = NULL;
    RATE_EW = NULL;
    RATE_UD = NULL;
    RATE_ERROR_LOS = NULL;
    RATE_ERROR_EW = NULL;
    RATE_ERROR_UD = NULL;
    TOPO_CORRECTION = NULL;
    NORM_X = NULL;
    NORM_AXY = NULL;
    ZSCORE_MASK = NULL;

    for (int i = 0; i < PAR->DSET.size(); i++)
        DSET.push_back(CSet(PAR, i));
};

void CArea::ComputeDIM()
{
    // compute DIM
    if (DSET.size() <= 0)
        throw std::invalid_argument(ERROR("incorrect number of sets"));
    else
    {
        int asc = 0, dsc = 0;
        for (int i = 0; i < DSET.size(); i++)
        {
            if (DSET[i].S[1] < 0)
                asc++;
            else if (DSET[i].S[1] > 0)
                dsc++;
            else
                throw std::invalid_argument(ERROR("azimuth cannot be equal to zero"));
        }
        if (asc > 0 && dsc > 0)
        {
            WriteLog("\ndetected ascending and descending data: applying MSBAS processing...");
            DIM = 2;
        }
        else if (asc > 0)
        {
            WriteLog("\ndetected ascending data: applying SBAS processing...");
            DIM = 1;
        }
        else if (dsc > 0)
        {
            WriteLog("\ndetected descending data: applying SBAS processing...");
            DIM = 1;
        }
    }
};

void CArea::MakeSLC()
{
    SLC.clear();

    if (DIM == 1)
    {
        for (int i = 0; i < DSET.size(); i++)
            for (int j = 0; j < DSET[i].SLC.size(); j++)
                SLC.push_back(CImage(DSET[i].SLC[j]));

    }
    else if (DIM == 2)
    {
        double start_asc = LARGE_NUM_FLOAT;
        double start_dsc = LARGE_NUM_FLOAT;

        double stop_asc = 0;
        double stop_dsc = 0;

        for (int i = 0; i < DSET.size(); i++)
        {
            if (DSET[i].S[1] < 0)
            {
                if (DSET[i].SLC[0].DATE < start_asc) start_asc = DSET[i].SLC[0].DATE;
                if (DSET[i].SLC[DSET[i].SLC.size() - 1].DATE > stop_asc) stop_asc = DSET[i].SLC[DSET[i].SLC.size() - 1].DATE;
            }
            else if (DSET[i].S[1] > 0)
            {
                if (DSET[i].SLC[0].DATE < start_dsc) start_dsc = DSET[i].SLC[0].DATE;
                if (DSET[i].SLC[DSET[i].SLC.size() - 1].DATE > stop_dsc) stop_dsc = DSET[i].SLC[DSET[i].SLC.size() - 1].DATE;
            }
        }

        double datestart = MAXIMUM(start_asc, start_dsc);
        double datestop = MINIMUM(stop_asc, stop_dsc);

        for (int i = 0; i < DSET.size(); i++)
            for (int j = 0; j < DSET[i].SLC.size(); j++)
                if (DSET[i].SLC[j].DATE >= datestart && DSET[i].SLC[j].DATE <= datestop) SLC.push_back(CImage(DSET[i].SLC[j]));

    }

    sort(SLC.begin(), SLC.end());
    SLC.erase(unique(SLC.begin(), SLC.end()), SLC.end());

    if (SLC.empty()) throw std::invalid_argument(ERROR("no valid SLCs are available"));

    WriteLog("\ncomputed overall temporal coverage: " + SLC[0].NAME + "-" + SLC[SLC.size() - 1].NAME + " or " + f2s(SLC[0].DATE) + "-" + f2s(SLC[SLC.size() - 1].DATE));

    WriteLog("\nselected " + i2s(SLC.size()) + " SLCs for SVD processing - NAME DATE SNUM MCOUNT SCOUNT: ");
    for (int i = 0; i < SLC.size(); i++) WriteLog(i2sf(i) + ": " + SLC[i].NAME + " " + f2s(SLC[i].DATE) + " " + i2sf(SLC[i].SNUM) + " " + i2sf(SLC[i].MCOUNT) + " " + i2sf(SLC[i].SCOUNT));

    // exclude interferograms outside of common span and compute boundary correction
    WriteLog("\nexcluding interferograms outside of common time span and computing boundary correction...");
    for (int i = 0; i < DSET.size(); i++)
    {
        WriteLog("\nset " + i2s(i) + ":");
        DSET[i].FilterInSAR(SLC);
    }
};

void CArea::MakeTimeMatrix1D()
{
    // topographic correction is only defined for R_FLAG 0 or 1
    if (PAR->T_FLAG && PAR->R_FLAG > 1) throw std::invalid_argument(ERROR("topographic correction is only defined for R_FLAG 0 or 1"));

    ROWS = 0;

    for (int i = 0; i < DSET.size(); i++)
    {
        DSET[i].MakeTimeMatrix1D(SLC);
        ROWS = ROWS + DSET[i].InSAR.size();
    }

    COLUMNS = (SLC.size() - 1);
    if (PAR->T_FLAG) COLUMNS += 1;

    //regularization
    if (PAR->R_FLAG == 1)
        ROWS = ROWS + COLUMNS;
    else if (PAR->R_FLAG == 2)
        ROWS = ROWS + COLUMNS - 1;
    else if (PAR->R_FLAG == 3)
        ROWS = ROWS + COLUMNS - 2;

    float_mem_alloc(TIME_MATRIX, COLUMNS * ROWS);

    int rows = 0;

    for (int k = 0; k < DSET.size(); k++)
    {
        for (int j = 0; j < DSET[k].InSAR.size(); j++)
        {
            for (int i = 0; i < (SLC.size() - 1); i++) TIME_MATRIX[i + rows * COLUMNS] = DSET[k].TIME_MATRIX[i + j * (SLC.size() - 1)];
            //add baseline
            if (PAR->T_FLAG) TIME_MATRIX[(COLUMNS - 1) + rows * COLUMNS] = DSET[k].InSAR[j].BASELINE;
            rows++;
        }
    }

    // add regularization lambda
    if (PAR->R_FLAG == 1)
        for (int i = 0; i < COLUMNS; i++)
        {
            TIME_MATRIX[(rows + i) * COLUMNS + i] = PAR->R_LAMBDA;
        }
    else if (PAR->R_FLAG == 2)
        for (int i = 0; i < (COLUMNS - 1); i++)
        {
            TIME_MATRIX[(rows + i) * COLUMNS + i] = -PAR->R_LAMBDA;
            TIME_MATRIX[(rows + i) * COLUMNS + i + 1] = PAR->R_LAMBDA;
        }
    else if (PAR->R_FLAG == 3)
        for (int i = 0; i < (COLUMNS - 2); i++)
        {
            TIME_MATRIX[(rows + i) * COLUMNS + i] = PAR->R_LAMBDA;
            TIME_MATRIX[(rows + i) * COLUMNS + i + 1] = -2 * PAR->R_LAMBDA;
            TIME_MATRIX[(rows + i) * COLUMNS + i + 2] = PAR->R_LAMBDA;
        }

    FILE * pFile;
    pFile = fopen("MSBAS_TIME_MATRIX.txt", "w");
    if (pFile == 0) throw std::invalid_argument(ERROR("cannot open file MSBAS_TIME_MATRIX.txt"));
    for (int j = 0; j < ROWS; j++)
    {
        for (int i = 0; i < COLUMNS; i++)
        {
            fprintf(pFile, "%11.6f ", TIME_MATRIX[i + COLUMNS * j]);
        }
        fprintf(pFile, "\n");
    }
    fclose(pFile);

    // create vector of pointers for InSAR, do not change afterwards the order of elements in the vector
    for (int l = 0; l < DSET.size(); l++) 
        for (int k = 0; k < DSET[l].InSAR.size(); k++)
            pInSAR.push_back(&DSET[l].InSAR[k]);
};

void CArea::MakeTimeMatrix2D()
{
    // topographic correction is only defined for R_FLAG 0 or 1
    if (PAR->T_FLAG && PAR->R_FLAG > 1) throw std::invalid_argument(ERROR("topographic correction is only defined for R_FLAG 0 or 1"));

    ROWS = 0;

    for (int i = 0; i < DSET.size(); i++)
    {
        DSET[i].MakeTimeMatrix2D(SLC);
        ROWS = ROWS + DSET[i].InSAR.size();
    }

    COLUMNS = 2 * (SLC.size() - 1);
    if (PAR->T_FLAG) COLUMNS += 1;

    //regularization
    if (PAR->R_FLAG == 1)
        ROWS = ROWS + COLUMNS;
    else if (PAR->R_FLAG == 2)
        ROWS = ROWS + COLUMNS - 2;
    else if (PAR->R_FLAG == 3)
        ROWS = ROWS + COLUMNS - 4;

    float_mem_alloc(TIME_MATRIX, COLUMNS * ROWS);

    int rows = 0;

    for (int k = 0; k < DSET.size(); k++)
    {
        for (int j = 0; j < DSET[k].InSAR.size(); j++)
        {
            for (int i = 0; i < 2 * (SLC.size() - 1); i++) TIME_MATRIX[i + rows * COLUMNS] = DSET[k].TIME_MATRIX[i + j * 2 * (SLC.size() - 1)];
            //add baseline
            if (PAR->T_FLAG) TIME_MATRIX[(COLUMNS - 1) + rows * COLUMNS] = DSET[k].InSAR[j].BASELINE;
            rows++;
        }
    }

    // add regularization lambda
    if (PAR->R_FLAG == 1)
        for (int i = 0; i < COLUMNS; i++)
        {
            TIME_MATRIX[(rows + i) * COLUMNS + i] = PAR->R_LAMBDA;
        }
    else if (PAR->R_FLAG == 2)
        for (int i = 0; i < (COLUMNS - 2); i++)
        {
            TIME_MATRIX[(rows + i) * COLUMNS + i] = -PAR->R_LAMBDA;
            TIME_MATRIX[(rows + i) * COLUMNS + i + 2] = PAR->R_LAMBDA;
        }
    else if (PAR->R_FLAG == 3)
        for (int i = 0; i < (COLUMNS - 4); i++)
        {
            TIME_MATRIX[(rows + i) * COLUMNS + i] = PAR->R_LAMBDA;
            TIME_MATRIX[(rows + i) * COLUMNS + i + 2] = -2 * PAR->R_LAMBDA;
            TIME_MATRIX[(rows + i) * COLUMNS + i + 4] = PAR->R_LAMBDA;
        }

    FILE * pFile;
    pFile = fopen("MSBAS_TIME_MATRIX.txt", "w");
    if (pFile == 0) throw std::invalid_argument(ERROR("cannot open file MSBAS_TIME_MATRIX.txt"));
    for (int j = 0; j < ROWS; j++)
    {
        for (int i = 0; i < COLUMNS; i++)
        {
            fprintf(pFile, "%11.6f ", TIME_MATRIX[i + COLUMNS * j]);
        }
        fprintf(pFile, "\n");
    }
    fclose(pFile);
    
    // create vector of pointers for InSAR, do not change afterwards the order of elements in the vector
    for (int l = 0; l < DSET.size(); l++) 
        for (int k = 0; k < DSET[l].InSAR.size(); k++)
            pInSAR.push_back(&DSET[l].InSAR[k]);
};

void CArea::ReadInterferograms()
{
    WriteLog("\nreading interferograms (order is managed by OpenMP and may appear irregular)...");
    
    for (int i = 0; i < DSET.size(); i++)
    {
        WriteLog("\nset " + i2s(i) + ":");
        
        #pragma omp parallel for
        for (int j = 0; j < DSET[i].InSAR.size(); j++)
        {
            DSET[i].InSAR[j].ReadInterferogram();
            WriteLog(i2sf(j) + ": " + DSET[i].InSAR[j].NAME);
        }
    }
};

void CArea::CalibrateInterferograms()
{
    WriteLog("\ncalibrating interferograms by computing and removing offset (order is managed by OpenMP and may appear irregular)...");
    
    // for C_FLAG=100 locate region with smallest ZScore
    if (PAR->C_FLAG == 100)
    {
        //reset C_FLAG=1  
        PAR->C_FLAG=1;

        WriteLog("\ncomputing location of minimal ZSCORE...");
        float zscore_min=1e9, fcov_min=0;
    
        // find coordinates of point with minimal ZSCORE
        #pragma omp parallel for collapse(2)
        for (int ii=0; ii<PAR->WIDTH; ii++)
            for (int jj=0; jj<PAR->LENGTH; jj++)
            {
                float zscore=0, fcov=0;
                int num=0, num_all=0;
                
                for (int i = (ii - PAR->CSIZE); i <= (ii + PAR->CSIZE); i++)
                    for (int j = (jj - PAR->LSIZE); j <= (jj + PAR->LSIZE); j++)
                    {
                        if ((i >= 0)&&(j >= 0)&&(i < PAR->WIDTH)&&(j < PAR->LENGTH) && ZSCORE_MASK[i+j*PAR->WIDTH] != 0)
                        {
                            zscore += ZSCORE_MASK[i+j*PAR->WIDTH];
                            num++;
                        }
                        num_all++;
                    }
                if (num > 0) zscore /= num;

                // fraction of the number of valid pixels in the region, consider regions with coverage greater than 0.75
                fcov=num/num_all;    
                
                #pragma omp critical
                {
                if(zscore != 0 && fcov > 0.5 && zscore < zscore_min)
                {
                    zscore_min=zscore;
                    PAR->CPOS[0] = ii;
                    PAR->LPOS[0] = jj;
                    fcov_min=fcov;
                }
                }
            }

        // print absolute coordinates even if WINDOW_SIZE is specified, but in processing use relative coordinates
        WriteLog("\nminimal ZSCORE " + f2s(zscore_min) +  " with coverage of " + i2s(int(100*fcov_min)) + "% is observed at column and row: " + i2sf(PAR->CPOS[0]+PAR->CSTART) + " " + i2sf(PAR->LPOS[0]+PAR->LSTART));
    }
        
    for (int i = 0; i < DSET.size(); i++)
    {
        WriteLog("\nSet " + i2s(i) + ":");

        #pragma omp parallel for
        for (int j = 0; j < DSET[i].InSAR.size(); j++)
        {
            DSET[i].InSAR[j].CalibrateInterferogram();
            WriteLog(i2sf(j) + ": " + DSET[i].InSAR[j].NAME + " " + f2s(DSET[i].InSAR[j].OFFSET));
        }
    }
};

void CArea::ComputeZScoreMask()
{
    WriteLog("\ncomputing ZSCORE_MASK...");
    
    float_mem_alloc(ZSCORE_MASK, PAR->WIDTH * PAR->LENGTH);

    #pragma omp parallel for collapse(2)
    for (int ii=0; ii<PAR->WIDTH; ii++)
        for (int jj=0; jj<PAR->LENGTH; jj++)
        {
            long num = 0;

            for (int i = 0; i < DSET.size(); i++)
                for (int j = 0; j < DSET[i].InSAR.size(); j++)
                {
                    if(num>=0 && DSET[i].InSAR[j].DATA[ii+jj*PAR->WIDTH] != 0)
                    {
                        ZSCORE_MASK[ii+jj*PAR->WIDTH] += abs(DSET[i].InSAR[j].DATA[ii+jj*PAR->WIDTH]-DSET[i].InSAR[j].MEAN)/DSET[i].InSAR[j].STD;
                        num++;
                    }
                    else
                    {
                        num=-1;
                        ZSCORE_MASK[ii+jj*PAR->WIDTH]=0;
                    }
                }
            
            if (num>0) ZSCORE_MASK[ii+jj*PAR->WIDTH] = ZSCORE_MASK[ii+jj*PAR->WIDTH]/num;
        }

    long count = 0;
    
    // compute percentage of coverage
    #pragma omp parallel for collapse(2) reduction (+:count)
    for (int ii=0; ii<PAR->WIDTH; ii++)
        for (int jj=0; jj<PAR->LENGTH; jj++)
            if(ZSCORE_MASK[ii+jj*PAR->WIDTH] != 0)
                count++;

    if (count == 0)
        throw std::invalid_argument(ERROR("no coherent pixels detected, processing cannot continue"));
    
    WriteLog("\nselected " + i2s(floor(100.0 * count / (PAR->WIDTH * PAR->LENGTH))) + "% coherent pixels for further processing");
    
    PAR->write_float_data("MSBAS_ZSCORE_MASK", ZSCORE_MASK);
};

void CArea::Inversion()
{
    WriteLog("\nstarting SVD inversion... ");
    WriteLog("inversion parameters for problem Ax=Y:");
    WriteLog("A: " + i2s(COLUMNS) + " x " + i2s(ROWS));
    WriteLog("x: " + i2s(PAR->WIDTH) + " x " + i2s(COLUMNS));
    WriteLog("Y: " + i2s(PAR->WIDTH) + " x " + i2s(ROWS));
    WriteLog("in summary " + i2s(COLUMNS) + " unknowns & " + i2s(ROWS) + " equations");

    if (COLUMNS > ROWS) throw std::invalid_argument(ERROR("COLUMNS > ROWS, under-determined problem"));

    if (PAR->T_FLAG) float_mem_alloc(TOPO_CORRECTION, PAR->WIDTH * PAR->LENGTH);

    float_mem_alloc(NORM_X, PAR->WIDTH * PAR->LENGTH);
    float_mem_alloc(NORM_AXY, PAR->WIDTH * PAR->LENGTH);

    for (int k = 0; k < SLC.size(); k++)
    {
        if (DIM == 1)
        {
            float_mem_alloc(SLC[k].DISP_LOS, PAR->LENGTH * PAR->WIDTH);
        }
        else if (DIM == 2)
        {
            float_mem_alloc(SLC[k].DISP_EW, PAR->LENGTH * PAR->WIDTH);
            float_mem_alloc(SLC[k].DISP_UD, PAR->LENGTH * PAR->WIDTH);
        }
    }
    
    long num = 0;
    if (DIM == 1)
        num = SLC.size() * PAR->WIDTH * PAR->LENGTH;
    else if (DIM == 2)
        num = 2 * SLC.size() * PAR->WIDTH * PAR->LENGTH;
    WriteLog("successfully allocated " + i2s(num * sizeof (float) / 1024000) + " MB of RAM for output data");

    int steps_completed = 0;
    #pragma omp parallel for
    for (int j = 0; j < PAR->LENGTH; j++)
    {
        // allocating memory required for inversion
        float *V;
        float_mem_alloc(V, ROWS * PAR->WIDTH);

        float *AT;
        float_mem_alloc(AT, COLUMNS * ROWS);

        float *WORK;
        long lwork = 10 * PAR->WIDTH * PAR->LENGTH;
        float_mem_alloc(WORK, lwork);

        float *S;
        float_mem_alloc(S, MINIMUM(COLUMNS, ROWS));

        // re-initialize V, column-major mode
        for (int i = 0; i < PAR->WIDTH; i++)
            for (int k = 0; k < pInSAR.size(); k++)
                V[k + i * ROWS] = pInSAR[k]->DATA[i + j * PAR->WIDTH];
        
        // re-initialize AT, column-major mode
        for (int ii = 0; ii < COLUMNS; ii++) for (int jj = 0; jj < ROWS; jj++) AT[jj + ii * ROWS] = TIME_MATRIX[ii + jj * COLUMNS];
        
        // perform inversion
        if (sgelss(AT, ROWS, COLUMNS, S, lwork, WORK, V, PAR->WIDTH) != 0) throw std::invalid_argument(ERROR("LAPACK sgelss function signaled error"));

        for (int i = 0; i < PAR->WIDTH; i++)
            if (ZSCORE_MASK[i + j * PAR->WIDTH] != 0)
            {
                // reconstruct time series
                for (int k = 0; k < (SLC.size() - 1); k++)
                {
                    double dt = (SLC[k + 1].DATE - SLC[k].DATE);
                    if (DIM == 1)
                    {
                        SLC[k + 1].DISP_LOS[i + j * PAR->WIDTH] = SLC[k].DISP_LOS[i + j * PAR->WIDTH] + V[k + i * ROWS] * dt;

                        // set to a small number to avoid conflict with NaN==0    
                        if (k == 0)
                            SLC[k].DISP_LOS[i + j * PAR->WIDTH]=SMALL_NUM_FLOAT;
                    }
                    else if (DIM == 2)
                    {
                        SLC[k + 1].DISP_EW[i + j * PAR->WIDTH] = SLC[k].DISP_EW[i + j * PAR->WIDTH] + V[2 * k + i * ROWS] * dt;
                        SLC[k + 1].DISP_UD[i + j * PAR->WIDTH] = SLC[k].DISP_UD[i + j * PAR->WIDTH] + V[2 * k + 1 + i * ROWS] * dt;
                        
                        // set to a small number to avoid conflict with NaN==0
                        if (k == 0)
                        {
                            SLC[k].DISP_EW[i + j * PAR->WIDTH]=SMALL_NUM_FLOAT;
                            SLC[k].DISP_UD[i + j * PAR->WIDTH]=SMALL_NUM_FLOAT;
                        }
                    }
                }
                
                if (PAR->T_FLAG && DIM == 1)
                    TOPO_CORRECTION[i + j * PAR->WIDTH] = V[(SLC.size() - 1) + i * ROWS];
                else if (PAR->T_FLAG && DIM == 2)
                    TOPO_CORRECTION[i + j * PAR->WIDTH] = V[2 * (SLC.size() - 1) + i * ROWS];

                // compute ||x||
                NORM_X[i + j * PAR->WIDTH] = 0;
                for (int k = 0; k < (SLC.size() - 1); k++)
                {
                    if (DIM == 1)
                    {
                        NORM_X[i + j * PAR->WIDTH] += V[k + i * ROWS] * V[k + i * ROWS];
                    }
                    else if (DIM == 2)
                    {
                        NORM_X[i + j * PAR->WIDTH] += V[2 * k + i * ROWS] * V[2 * k + i * ROWS] + V[2 * k + 1 + i * ROWS] * V[2 * k + 1 + i * ROWS];
                    }
                }
                if (PAR->T_FLAG && DIM == 1)
                    NORM_X[i + j * PAR->WIDTH] += V[(SLC.size() - 1) + i * ROWS] * V[(SLC.size() - 1) + i * ROWS];
                else if (PAR->T_FLAG && DIM == 2)
                    NORM_X[i + j * PAR->WIDTH] += V[2 * (SLC.size() - 1) + i * ROWS] * V[2 * (SLC.size() - 1) + i * ROWS];
                
                NORM_X[i + j * PAR->WIDTH] = sqrt(NORM_X[i + j * PAR->WIDTH]);

                // compute ||Ax-Y||
                NORM_AXY[i + j * PAR->WIDTH] = 0;
                float norm_tmp = 0; 
                for (int l = 0; l < pInSAR.size(); l++)
                {
                    norm_tmp = 0;
                    for (int k = 0; k < (SLC.size() - 1); k++)
                        if (DIM == 1)
                        {
                            norm_tmp += TIME_MATRIX[k + l * COLUMNS] * V[k + i * ROWS];
                        }
                        else if (DIM == 2)
                        {
                            norm_tmp += TIME_MATRIX[2 * k + l * COLUMNS] * V[2 * k + i * ROWS] + TIME_MATRIX[2 * k + 1 + l * COLUMNS] * V[2 * k + 1 + i * ROWS];
                        }

                    if (PAR->T_FLAG && DIM == 1)
                        norm_tmp += TIME_MATRIX[(SLC.size() - 1) + l * COLUMNS] * V[(SLC.size() - 1) + i * ROWS];  
                    else if (PAR->T_FLAG && DIM == 2)
                        norm_tmp += TIME_MATRIX[2*(SLC.size() - 1) + l * COLUMNS] * V[2 * (SLC.size() - 1) + i * ROWS];  
                                                
                    norm_tmp -= pInSAR[l]->DATA[i + j * PAR->WIDTH];
                    
                    NORM_AXY[i + j * PAR->WIDTH] += norm_tmp*norm_tmp;
                }
                NORM_AXY[i + j * PAR->WIDTH] = sqrt(NORM_AXY[i + j * PAR->WIDTH]);
            }
        free(AT);
        free(WORK);
        free(S);
        free(V);
        
        #pragma omp atomic
        ++steps_completed;
        
        if (steps_completed % int(round(PAR->LENGTH/10.0)) == 1)
        {
            #pragma omp critical
            printf("completed %i%% \n", (int) round(100.0*steps_completed/PAR->LENGTH));
        }
    }

    WriteLog("SVD inversion completed successfully");
    
    // compute average norm for the entire image
    double norm_x = 0, norm_axy = 0;
    long count = 0;
   
    for (int j = 0; j < PAR->LENGTH; j++)
        for (int i = 0; i < PAR->WIDTH; i++)
            if (ZSCORE_MASK[i + j * PAR->WIDTH] != 0)
            {
                norm_x += NORM_X[i + j * PAR->WIDTH];
                norm_axy += NORM_AXY[i + j * PAR->WIDTH];
                count ++;
            }
    
    if (count > 0)
    {
        norm_x = norm_x/count;
        norm_axy = norm_axy/count;
        
        WriteLog("computed ||x|| and ||Ax-Y|| norms: " + f2s(norm_x) + " " + f2s(norm_axy));
    }
    else
        throw std::invalid_argument(ERROR("no coherent in all interferograms pixels detected"));
        
    
};

void CArea::PostProcessing()
{
    WriteLog("\nstarting post processing... ");

    //TemporalGaussianFilter();

    ComputeLinearRate();

    WriteResultsToDisk();

    InteractiveMode();

};

void CArea::TemporalGaussianFilter()
{
    if (PAR->TAV_FLAG > 0)
    {
        WriteLog("applying Gaussian filter in time domain");

        float *w;
        float_mem_alloc(w, SLC.size());

        double dt2 = 0;
        float fsum = 0;

        float ts = 0;

        for (int k = 0; k < SLC.size(); k++)
        {
            fsum = 0;
            for (int kk = 0; kk < SLC.size(); kk++)
            {
                dt2 = (SLC[kk].DATE - SLC[k].DATE)*(SLC[kk].DATE - SLC[k].DATE);
                w[kk] = exp(-dt2 / (2 * PAR->TAV_FLAG * PAR->TAV_FLAG));
                fsum = fsum + w[kk];
            }
            for (int kk = 0; kk < SLC.size(); kk++) w[kk] = w[kk] / fsum;

            for (int j = 0; j < PAR->LENGTH; j++)
                for (int i = 0; i < PAR->WIDTH; i++)
                {
                    if (DIM == 1 && ZSCORE_MASK[i + j * PAR->WIDTH] != 0)
                    {
                        ts = 0;
                        for (int kk = 0; kk < SLC.size(); kk++) ts = ts + SLC[kk].DISP_LOS[i + j * PAR->WIDTH] * w[kk];
                        SLC[k].DISP_LOS[i + j * PAR->WIDTH] = ts;
                    }
                    else if (DIM == 2 && ZSCORE_MASK[i + j * PAR->WIDTH] != 0)
                    {
                        ts = 0;
                        for (int kk = 0; kk < SLC.size(); kk++) ts = ts + SLC[kk].DISP_EW[i + j * PAR->WIDTH] * w[kk];
                        SLC[k].DISP_EW[i + j * PAR->WIDTH] = ts;

                        ts = 0;
                        for (int kk = 0; kk < SLC.size(); kk++) ts = ts + SLC[kk].DISP_UD[i + j * PAR->WIDTH] * w[kk];
                        SLC[k].DISP_UD[i + j * PAR->WIDTH] = ts;
                    }

                }
        }

        free(w);
    }
};

void CArea::WriteResultsToDisk()
{
    WriteLog("writing results to a disk");

    // create file that will be used in post-processing for extracting time series from binary MSBAS output
    FILE * pFile0 = fopen("MSBAS_TSOUT.txt", "w");
    if (!pFile0) throw std::invalid_argument(ERROR("can not write to file MSBAS_TSOUT.txt"));
    if (PAR->FORMAT == 2)
        fprintf(pFile0, "%s=%i,%i\n", "FORMAT", PAR->FORMAT, PAR->FORMAT_INT_FLAG);
    else
        fprintf(pFile0, "%s=%i\n", "FORMAT", PAR->FORMAT);
    
    fprintf(pFile0, "%s=%i,%i\n", "FILE_SIZE", PAR->FWIDTH, PAR->LENGTH);
    fprintf(pFile0, "%s=%i,", "C_FLAG", PAR->C_FLAG);
    if ((PAR->C_FLAG > 0)&&(PAR->C_FLAG < 10))
        for (int i = 0; i < PAR->C_FLAG; i++)
            fprintf(pFile0, "%i,%i,", PAR->CPOS[i]+PAR->CSTART, PAR->LPOS[i]+PAR->LSTART);
    fprintf(pFile0, "%i,%i\n", PAR->CSIZE, PAR->LSIZE);
    
    for (int k = 0; k < SLC.size(); k++)
    {
        if (DIM == 1)
        {
            string f1="MSBAS_" + SLC[k].NAME + "_LOS";
            if (PAR->FORMAT == 0 || PAR->FORMAT == 1)
            {
                f1 += ".bin";
            }
            else if (PAR->FORMAT == 2)
            {
                f1 += ".tif";
            }
            fprintf(pFile0, "%s %f %s\n", SLC[k].NAME.c_str(), SLC[k].DATE, f1.c_str());
        }
        else if (DIM == 2)
        {
            string f1="MSBAS_" + SLC[k].NAME + "_EW";
            string f2="MSBAS_" + SLC[k].NAME + "_UD";
            
            if (PAR->FORMAT == 0 || PAR->FORMAT == 1)
            {
                f1 += ".bin";
                f2 += ".bin";
            }
            else if (PAR->FORMAT == 2)
            {
                f1 += ".tif";
                f2 += ".tif";
            }
            fprintf(pFile0, "%s %f %s %s\n", SLC[k].NAME.c_str(), SLC[k].DATE, f1.c_str(), f2.c_str());
        }
    }
    fclose(pFile0);
    
    #pragma omp parallel for
    for (int k = 0; k < SLC.size(); k++) SLC[k].WriteData();

    if (DIM == 1)
    {
        if (RATE_LOS) PAR->write_float_data("MSBAS_LINEAR_RATE_LOS", RATE_LOS);
        if (RATE_ERROR_LOS) PAR->write_float_data("MSBAS_LINEAR_RATE_ERROR_LOS", RATE_ERROR_LOS);
    }
    else if (DIM == 2)
    {
        if (RATE_EW) PAR->write_float_data("MSBAS_LINEAR_RATE_EW", RATE_EW);
        if (RATE_ERROR_EW) PAR->write_float_data("MSBAS_LINEAR_RATE_ERROR_EW", RATE_ERROR_EW);
        if (RATE_UD) PAR->write_float_data("MSBAS_LINEAR_RATE_UD", RATE_UD);
        if (RATE_ERROR_UD) PAR->write_float_data("MSBAS_LINEAR_RATE_ERROR_UD", RATE_ERROR_UD);
    }

    if (PAR->T_FLAG) PAR->write_float_data("MSBAS_TOPO_CORRECTION", TOPO_CORRECTION);

    PAR->write_float_data("MSBAS_NORM_X", NORM_X);
    PAR->write_float_data("MSBAS_NORM_AXY", NORM_AXY);
};

void CArea::ComputeLinearRate()
{
    WriteLog("computing linear rate(s)");

    if (DIM == 1)
    {
        float_mem_alloc(RATE_LOS, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_ERROR_LOS, PAR->WIDTH * PAR->LENGTH);
    }
    else if (DIM == 2)
    {
        float_mem_alloc(RATE_EW, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_ERROR_EW, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_UD, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_ERROR_UD, PAR->WIDTH * PAR->LENGTH);
    }

    #pragma omp parallel for collapse(2)
    for (int j = 0; j < PAR->LENGTH; j++)
        for (int i = 0; i < PAR->WIDTH; i++)
        {
            float *x, *y;
            float_mem_alloc(x, SLC.size());
            float_mem_alloc(y, SLC.size());

            float Xa, Xb, Xsigmaa, Xsigmab, chi2;

            if (DIM == 1 && ZSCORE_MASK[i + j * PAR->WIDTH] != 0)
            {
                for (int k = 0; k < SLC.size(); k++)
                {
                    x[k] = SLC[k].DATE;
                    y[k] = SLC[k].DISP_LOS[i + j * PAR->WIDTH];
                }

                linreg(x, y, SLC.size(), Xa, Xb, Xsigmaa, Xsigmab, chi2);
                RATE_LOS[i + j * PAR->WIDTH] = Xb;
                RATE_ERROR_LOS[i + j * PAR->WIDTH] = Xsigmab;
            }
            else if (DIM == 2 && ZSCORE_MASK[i + j * PAR->WIDTH] != 0)
            {
                for (int k = 0; k < SLC.size(); k++)
                {
                    x[k] = SLC[k].DATE;
                    y[k] = SLC[k].DISP_EW[i + j * PAR->WIDTH];
                }

                linreg(x, y, SLC.size(), Xa, Xb, Xsigmaa, Xsigmab, chi2);
                RATE_EW[i + j * PAR->WIDTH] = Xb;
                RATE_ERROR_EW[i + j * PAR->WIDTH] = Xsigmab;

                for (int k = 0; k < SLC.size(); k++)
                {
                    x[k] = SLC[k].DATE;
                    y[k] = SLC[k].DISP_UD[i + j * PAR->WIDTH];
                }

                linreg(x, y, SLC.size(), Xa, Xb, Xsigmaa, Xsigmab, chi2);
                RATE_UD[i + j * PAR->WIDTH] = Xb;
                RATE_ERROR_UD[i + j * PAR->WIDTH] = Xsigmab;
            }
            free(x);
            free(y);
        }
};

void CArea::InteractiveMode()
{
    int cref = 0;
    int lref = 0;
    int cradius = 0;
    int lradius = 0;
    int num = 0;
    float a1 = 0;
    float a2 = 0;
    float sd1 = 0;
    float sd2 = 0;

    if (PAR->I_FLAG == 1)
    {
        WriteLog("starting interactive mode");
        int flag = 1;
        while (flag == 1)
        {
            printf("\nenter five space delimited parameters - region location [column row], half-dimensions [column_radius row_radius] and repeat flag [0-no, 1-yes], e.g. for a region centered at [10 20] with dimensions of [2*5+1 2*8+1] enter '10 20 5 8 1': ");
            scanf("%i %i %i %i %i", &cref, &lref, &cradius, &lradius, &flag);
            string name = "MSBAS_" + i2s(cref) + "_" + i2s(lref) + "_" + i2s(cradius) + "_" + i2s(lradius) + ".txt";

            cref = cref - PAR->CSTART;
            lref = lref - PAR->LSTART;
            if ((cref >= 0)&&(cref < PAR->WIDTH)&&(lref >= 0)&&(lref < PAR->LENGTH))
            {
                FILE * pFile = fopen(name.c_str(), "w");
                if (!pFile) throw std::invalid_argument(ERROR("can not open file " + name));
                WriteLog("writing results (YYMMDDTHHMMSS YYYY.YYYYYY [for each dimension] DISP DISP_ERROR) to a file " + name);
                
                for (int k = 0; k < SLC.size(); k++)
                {
                    a1 = 0;
                    a2 = 0;
                    sd1 = 0;
                    sd2 = 0;

                    num = 0;

                    for (int ii = (cref - cradius); ii <= (cref + cradius); ii++)
                        for (int jj = (lref - lradius); jj <= (lref + lradius); jj++)
                            if ((ii >= 0)&&(jj >= 0)&&(ii < PAR->WIDTH)&&(jj < PAR->LENGTH))
                            {
                                if (DIM == 1)
                                {
                                    a1 = a1 + SLC[k].DISP_LOS[ii + jj * PAR->WIDTH];
                                    num++;
                                }
                                else if (DIM == 2)
                                {
                                    a1 = a1 + SLC[k].DISP_EW[ii + jj * PAR->WIDTH];
                                    a2 = a2 + SLC[k].DISP_UD[ii + jj * PAR->WIDTH];
                                    num++;
                                }
                            }

                    if (DIM == 1)
                    {
                        if (num > 0)
                        {
                            a1 = a1 / num;
                        }
                        else
                            a1 = 0;
                    }
                    else if (DIM == 2)
                    {
                        if (num > 0)
                        {
                            a1 = a1 / num;
                            a2 = a2 / num;
                        }
                        else
                            a1 = a2 = 0;
                    }

                    num = 0;

                    for (int ii = (cref - cradius); ii <= (cref + cradius); ii++)
                        for (int jj = (lref - lradius); jj <= (lref + lradius); jj++)
                            if ((ii >= 0)&&(jj >= 0)&&(ii < PAR->WIDTH)&&(jj < PAR->LENGTH))
                            {
                                if (DIM == 1)
                                {
                                    sd1 = sd1 + (SLC[k].DISP_LOS[ii + jj * PAR->WIDTH] - a1)*(SLC[k].DISP_LOS[ii + jj * PAR->WIDTH] - a1);
                                    num++;
                                }
                                else if (DIM == 2)
                                {
                                    sd1 = sd1 + (SLC[k].DISP_EW[ii + jj * PAR->WIDTH] - a1)*(SLC[k].DISP_EW[ii + jj * PAR->WIDTH] - a1);
                                    sd2 = sd2 + (SLC[k].DISP_UD[ii + jj * PAR->WIDTH] - a2)*(SLC[k].DISP_UD[ii + jj * PAR->WIDTH] - a2);
                                    num++;
                                }
                            }

                    if (DIM == 1)
                    {
                        if (num > 1)
                            sd1 = sqrt(sd1 / (num - 1));
                        else
                            sd1 = 0;

                        fprintf(pFile, "%s %f %f %f\n", SLC[k].NAME.c_str(), SLC[k].DATE, a1, sd1);
                    }
                    else if (DIM == 2)
                    {
                        if (num > 1)
                        {
                            sd1 = sqrt(sd1 / (num - 1));
                            sd2 = sqrt(sd2 / (num - 1));
                        }
                        else sd1 = sd2 = 0;

                        fprintf(pFile, "%s %f %f %f %f %f\n", SLC[k].NAME.c_str(), SLC[k].DATE, a1, sd1, a2, sd2);
                    }
                }
                fclose(pFile);
            }
            else
            {
                printf("warning: this point exceeds boundaries of the region, try again\n");
            }
        }
    }
        // read from INAME file
    else if (PAR->I_FLAG == 2)
    {
        FILE * pIFile = fopen(PAR->INAME.c_str(), "r");
        if (!pIFile) throw std::invalid_argument(ERROR("can not open file" + PAR->INAME));

        while (fscanf(pIFile, "%i %i %i %i", &cref, &lref, &cradius, &lradius) == 4)
        {
            string name = "MSBAS_" + i2s(cref) + "_" + i2s(lref) + "_" + i2s(cradius) + "_" + i2s(lradius) + ".txt";
            cref = cref - PAR->CSTART;
            lref = lref - PAR->LSTART;
            if ((cref >= 0)&&(cref < PAR->WIDTH)&&(lref >= 0)&&(lref < PAR->LENGTH))
            {
                FILE * pFile = fopen(name.c_str(), "w");
                if (!pFile) throw std::invalid_argument(ERROR("can not open file " + name));
                WriteLog("writing results (YYMMDDTHHMMSS YYYY.YYYYYY [for each dimension] DISP DISP_ERROR) to a file " + name);
                
                for (int k = 0; k < SLC.size(); k++)
                {
                    a1 = 0;
                    a2 = 0;
                    sd1 = 0;
                    sd2 = 0;
                            
                    num = 0;

                    for (int ii = (cref - cradius); ii <= (cref + cradius); ii++)
                        for (int jj = (lref - lradius); jj <= (lref + lradius); jj++)
                            if ((ii >= 0)&&(jj >= 0)&&(ii < PAR->WIDTH)&&(jj < PAR->LENGTH))
                            {
                                if (DIM == 1)
                                {
                                    a1 = a1 + SLC[k].DISP_LOS[ii + jj * PAR->WIDTH];
                                    num++;
                                }
                                else if (DIM == 2)
                                {
                                    a1 = a1 + SLC[k].DISP_EW[ii + jj * PAR->WIDTH];
                                    a2 = a2 + SLC[k].DISP_UD[ii + jj * PAR->WIDTH];
                                    num++;
                                }
                            }

                    if (DIM == 1)
                    {
                        if (num > 0)
                        {
                            a1 = a1 / num;
                        }
                        else
                            a1 = 0;
                    }
                    else if (DIM == 2)
                    {
                        if (num > 0)
                        {
                            a1 = a1 / num;
                            a2 = a2 / num;
                        }
                        else
                            a1 = a2 = 0;
                    }

                    num = 0;

                    for (int ii = (cref - cradius); ii <= (cref + cradius); ii++)
                        for (int jj = (lref - lradius); jj <= (lref + lradius); jj++)
                            if ((ii >= 0)&&(jj >= 0)&&(ii < PAR->WIDTH)&&(jj < PAR->LENGTH))
                            {
                                if (DIM == 1)
                                {
                                    sd1 = sd1 + (SLC[k].DISP_LOS[ii + jj * PAR->WIDTH] - a1)*(SLC[k].DISP_LOS[ii + jj * PAR->WIDTH] - a1);
                                    num++;
                                }
                                else if (DIM == 2)
                                {
                                    sd1 = sd1 + (SLC[k].DISP_EW[ii + jj * PAR->WIDTH] - a1)*(SLC[k].DISP_EW[ii + jj * PAR->WIDTH] - a1);
                                    sd2 = sd2 + (SLC[k].DISP_UD[ii + jj * PAR->WIDTH] - a2)*(SLC[k].DISP_UD[ii + jj * PAR->WIDTH] - a2);
                                    num++;
                                }
                            }

                    if (DIM == 1)
                    {
                        if (num > 1)
                            sd1 = sqrt(sd1 / (num - 1));
                        else
                            sd1 = 0;

                        fprintf(pFile, "%s %f %f %f\n", SLC[k].NAME.c_str(), SLC[k].DATE, a1, sd1);
                    }
                    else if (DIM == 2)
                    {
                        if (num > 1)
                        {
                            sd1 = sqrt(sd1 / (num - 1));
                            sd2 = sqrt(sd2 / (num - 1));
                        }
                        else sd1 = sd2 = 0;

                        fprintf(pFile, "%s %f %f %f %f %f\n", SLC[k].NAME.c_str(), SLC[k].DATE, a1, sd1, a2, sd2);
                    }
                }
                fclose(pFile);
            }
            else
            {
                printf("warning.: this point exceeds boundaries of the region, skipping\n");
            }
        }
        fclose(pIFile);
    }
        // print all to a single file
    else if (PAR->I_FLAG == 3)
    {
        FILE * pFile = fopen("MSBAS_TIME_SERIES.txt", "w");
        if (!pFile) throw std::invalid_argument(ERROR("can not open MSBAS_TIME_SERIES.txt file."));

        if (DIM == 1)
        {
            fprintf(pFile, "%11s %11s %11s %11s ", "X", "Y", "RATE_LOS", "RATE_ERR_LOS");
            for (int k = 0; k < SLC.size(); k++) fprintf(pFile, "%11.6f ", SLC[k].DATE);
            fprintf(pFile, "\n");
        }
        else if (DIM == 2)
        {
            fprintf(pFile, "%11s %11s %11s %11s %11s %11s ", "X", "Y", "RATE_EW", "RATE_ERR_EW", "RATE_UD", "RATE_ERR_UD");
            for (int k = 0; k < SLC.size(); k++) fprintf(pFile, "%11.6f %11.6f ", SLC[k].DATE, SLC[k].DATE);
            fprintf(pFile, "\n");
        }

        for (int i = 0; i < PAR->WIDTH; i++) for (int j = 0; j < PAR->LENGTH; j++) if (ZSCORE_MASK[i + j * PAR->WIDTH] != 0)
                {
                    cref = i - PAR->CSTART;
                    lref = j - PAR->LSTART;
                    if ((cref >= 0)&&(cref < PAR->WIDTH)&&(lref >= 0)&&(lref < PAR->LENGTH))
                    {
                        if (DIM == 1)
                        {
                            fprintf(pFile, "%11i %11i %11.6f %11.6f ", i, j, RATE_LOS[i + j * PAR->WIDTH], RATE_ERROR_LOS[i + j * PAR->WIDTH]);
                            for (int k = 0; k < SLC.size(); k++) fprintf(pFile, "%11.6f ", SLC[k].DISP_LOS[i + j * PAR->WIDTH]);
                            fprintf(pFile, "\n");
                        }
                        else if (DIM == 2)
                        {
                            fprintf(pFile, "%11i %11i %11.6f %11.6f %11.6f %11.6f ", i, j, RATE_EW[i + j * PAR->WIDTH], RATE_ERROR_EW[i + j * PAR->WIDTH], RATE_UD[i + j * PAR->WIDTH], RATE_ERROR_UD[i + j * PAR->WIDTH]);

                            for (int k = 0; k < SLC.size(); k++) fprintf(pFile, "%11.6f %11.6f ", SLC[k].DISP_EW[i + j * PAR->WIDTH], SLC[k].DISP_UD[i + j * PAR->WIDTH]);
                            fprintf(pFile, "\n");
                        }
                    }
                }
        fclose(pFile);
    }
};

int CArea::sgelss(float *a, int m, int n, float *s, int lwork, float *work, float *b, int nrhs)
{
    int info, rank, ldb;
    float rcond = 0.00001; // condition number ... choose a small value
    int lda = m; // leading dimension of A
    if (m > n) // leading dimension of B
        ldb = m;
    else
        ldb = n;

    sgelss_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, &info);
    //int sdim = MINIMUM(m,n);
    //printf("Rank of matrix A: %d\n", rank);
    //printf("Singular values of A: ");
    //for(int i = 0; i < sdim; i++) printf("%f\t", s[i]);
    //printf("\n");

    if (info < 0)
    {
        WriteLog("\nLapack routine sgelss returned error");
        WriteLog("invalid argument: " + f2s(info));
    }
    else if (info > 0)
    {
        WriteLog("\nSVD failed to converge and returned: " + f2s(info));
        WriteLog("see man page for details");
    }
    else
    {
        //  int sdim = MINIMUM(m,n);
        //  printf("Rank of matrix A: %d\n", rank);
        //  printf("Singular values of A: ");
        //  for(int i = 0; i < sdim; i++) printf("%f\t", s[i]);
        //  printf("\n");
        //  printf("Conditioning number of A: %f\n ", s[0]/s[sdim-1]);
        //  printf("SVD solution is: ");
        //  for(i = 0; i < sdim; i++) printf("%f\t", b[i]);
        //  printf("\n");
    }

    return info;
};

void CArea::linreg(float *x, float *y, int length, float &aa, float &bb, float &sigmaa, float &sigmab, float &chi2)
{
    if (length < 2) throw std::invalid_argument(ERROR("at least two points required for computing linear trend"));

    aa = bb = sigmaa = sigmab = chi2 = 0;

    float ss = 0, sx = 0, sy = 0, st2 = 0, t = 0, sxoss = 0, sigdat = 0;

    for (int i = 0; i < length; i++)
    {
        sx += x[i];
        sy += y[i];
    }

    ss = length;
    sxoss = sx / ss;

    for (int i = 0; i < length; i++)
    {
        t = x[i] - sxoss;
        st2 += t*t;
        bb += t * y[i];
    }

    bb /= st2;
    aa = (sy - sx * bb) / ss;
    sigmaa = sqrt((1.0 + sx * sx / (ss * st2)) / ss);
    sigmab = sqrt(1.0 / st2);

    for (int i = 0; i < length; i++) chi2 += (y[i] - aa - bb * x[i])*(y[i] - aa - bb * x[i]);

    if (length > 2)
    {
        sigdat = sqrt(chi2 / (length - 2));
        sigmaa *= sigdat;
        sigmab *= sigdat;
    }
};
