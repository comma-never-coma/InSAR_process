#include "CArea.h"

CArea::CArea(CParam *par)
{
    PAR = par;
  
    LWORK = 0;
    RANK = 0;
    COND_NUM = 0;
    
    DIM = 0;
    ROWS = COLUMNS = 0;
    TIME_MATRIX = NULL;
    RATE_LOS = NULL;
    RATE_NS = NULL;
    RATE_EW = NULL;
    RATE_UD = NULL;
    RATE_UD1 = NULL;
    RATE_STD_LOS = NULL;
    RATE_STD_NS = NULL;
    RATE_STD_EW = NULL;
    RATE_STD_UD = NULL;
    RATE_STD_UD1 = NULL;
    RATE_R2_LOS = NULL;
    RATE_R2_NS = NULL;
    RATE_R2_EW = NULL;
    RATE_R2_UD = NULL;
    RATE_R2_UD1 = NULL;
    NORM_X = NULL;
    NORM_AXY = NULL;
    RANK=NULL;
    COND_NUM=NULL;
    ZSCORE_MASK = NULL;

    for (int i = 0; i < PAR->DSET.size(); i++)
        DSET.push_back(CSet(PAR, i));
};

void CArea::ComputeDIM()
{
     SLC.clear();
    
    // compute DIM
    if (DSET.size() <= 0)
        throw std::invalid_argument(ERROR("incorrect number of sets"));
    else
    {
        int ra = 0, aa = 0, rd = 0, ad = 0, dem = 0; // range ascending, azimuth ascending, range descending, azimuth descending, DEM
        for (int i = 0; i < DSET.size(); i++)
        {
            if (cos(DSET[i].AZIMUTH * PI / 180.) > 0 && DSET[i].TYPE == 0)
                ra++;
            else if (cos(DSET[i].AZIMUTH * PI / 180.) > 0 && DSET[i].TYPE == 1)
                aa++;
            else if (cos(DSET[i].AZIMUTH * PI / 180.) < 0 && DSET[i].TYPE == 0)
                rd++;
            else if (cos(DSET[i].AZIMUTH * PI / 180.) < 0 && DSET[i].TYPE == 1)
                ad++;
            else
                throw std::invalid_argument(ERROR("incorrect set type and/or azimuth"));
        }

        // DEM 
        if (PAR->DDNS && PAR->DDEW) dem = 1;

        // This is not an exhaustive list!

        // just one kind of a set
        if (((aa > 0 ) &&  (ad + ra + rd) == 0) || ((ad > 0 ) &&  (aa + ra + rd) == 0) || ((ra > 0 ) &&  (aa + ad + rd) == 0) || ((rd > 0 ) &&  (aa + ad + ra) == 0))
        {
            DIM = 1;
            for (int i = 0; i < DSET.size(); i++) for (int j = 0; j < DSET[i].SLC.size(); j++) SLC.push_back(CImage(DSET[i].SLC[j]));
        }
        // no azimuth offsets, no DEM
        else if ((aa + ad) == 0 && dem == 0 && ra > 0 && rd > 0)
        {
            DIM = 2;
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
            for (int i = 0; i < DSET.size(); i++) for (int j = 0; j < DSET[i].SLC.size(); j++) if (DSET[i].SLC[j].DATE >= datestart && DSET[i].SLC[j].DATE <= datestop) SLC.push_back(CImage(DSET[i].SLC[j]));
        }
        // (DEM (SPF) or azimuth offsets) and range offsets
        else if ((dem > 0 || (aa + ad) > 0) && ra > 0 && rd > 0)
        {
            if (dem > 0 && (aa + ad) > 0 && PAR->D_FLAG == 1)
                DIM = 4;  
            else
                DIM = 3;
            double start_asc = LARGE_NUM_FLOAT;
            double start_dsc = LARGE_NUM_FLOAT;
            double stop_asc = 0;
            double stop_dsc = 0;
            for (int i = 0; i < DSET.size(); i++)
                if (DSET[i].TYPE == 0)
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
            for (int i = 0; i < DSET.size(); i++) for (int j = 0; j < DSET[i].SLC.size(); j++) if (DSET[i].SLC[j].DATE >= datestart && DSET[i].SLC[j].DATE <= datestop) SLC.push_back(CImage(DSET[i].SLC[j]));
        }
        // either ascending or descending azimuth and range and DEM (SPF)
        else if (((aa > 0 && ra > 0 && ad == 0 && rd == 0) || (ad > 0 && rd > 0 && aa == 0 && ra == 0)) && dem > 0)
        {
            DIM = 3;
            for (int i = 0; i < DSET.size(); i++) for (int j = 0; j < DSET[i].SLC.size(); j++) SLC.push_back(CImage(DSET[i].SLC[j]));
        }
        else
            throw std::invalid_argument(ERROR("cannot determine problem dimension, detected unsupported combination of input data"));
    }
    
    WriteLog("\napplying MSBAS-"+i2s(DIM)+"D processing...");
    sort(SLC.begin(), SLC.end());
    SLC.erase(unique(SLC.begin(), SLC.end()), SLC.end());

    if (SLC.empty()) throw std::invalid_argument(ERROR("no valid SLCs are available"));

    WriteLog("\ncomputed overall temporal coverage: " + SLC[0].NAME + "-" + SLC[SLC.size() - 1].NAME + " or " + f2s(SLC[0].DATE) + "-" + f2s(SLC[SLC.size() - 1].DATE));

    WriteLog("\nselected " + i2s(SLC.size()) + " SLCs for SVD processing - NAME DATE SETN MCOUNT SCOUNT: ");
    for (int i = 0; i < SLC.size(); i++) WriteLog(i2sf(i) + ": " + SLC[i].NAME + " " + f2s(SLC[i].DATE) + " " + i2sf(SLC[i].SETN) + " " + i2sf(SLC[i].MCOUNT) + " " + i2sf(SLC[i].SCOUNT));

    // exclude interferograms outside of common span and compute boundary correction
    WriteLog("\nexcluding interferograms outside of common time span and computing boundary correction...");
    for (int i = 0; i < DSET.size(); i++)
    {
        WriteLog("\nset " + i2s(i) + ":");
        DSET[i].ApplyBoundaryCorrectionInSAR(SLC);
    }     
};

void CArea::MakeTimeMatrix()
{
    // estimate matrix size

    ROWS = 0;
    COLUMNS = DIM * (SLC.size() - 1);
    
    for (int i = 0; i < DSET.size(); i++)
    {
        DSET[i].MakeTimeMatrix(SLC);
        ROWS = ROWS + DSET[i].InSAR.size();
    }

    // constraint for motion along slope
    if (PAR->DDNS && PAR->DDEW)
        ROWS = ROWS + (SLC.size() - 1); 

    // regularization
    if (PAR->R_FLAG == 1)
        ROWS = ROWS + COLUMNS;
    else if (PAR->R_FLAG == 2 && (COLUMNS - 1*DIM)>0)
        ROWS = ROWS + COLUMNS - 1*DIM;
    else if (PAR->R_FLAG == 3 && (COLUMNS - 2*DIM)>0)
        ROWS = ROWS + COLUMNS - 2*DIM;

    // allocate memory for matrix 
    float_mem_alloc(TIME_MATRIX, COLUMNS * ROWS);
    
    // populate matrix
    int rows = 0;

    for (int k = 0; k < DSET.size(); k++)
    {
        for (int j = 0; j < DSET[k].InSAR.size(); j++)
        {
            for (int i = 0; i < (SLC.size() - 1); i++)
                for (int l=0; l<DIM; l++)
                    TIME_MATRIX[l + i*DIM + rows * COLUMNS] = DSET[k].TIME_MATRIX[i + j *  (SLC.size() - 1)];
            rows++;
        }
    }

    // constraint for motion along slope
    if (PAR->DDNS && PAR->DDEW)
        rows = rows + (SLC.size() - 1);

    // regularization
    if (PAR->R_FLAG == 1)
        for (int i = 0; i < COLUMNS; i++)
        {
            TIME_MATRIX[(rows + i) * COLUMNS + i] = PAR->R_LAMBDA;
        }
    else if (PAR->R_FLAG == 2)
        for (int i = 0; i < (COLUMNS - 1*DIM); i++)
        {
            TIME_MATRIX[(rows + i) * COLUMNS + i] = -PAR->R_LAMBDA;
            TIME_MATRIX[(rows + i) * COLUMNS + i + 1*DIM] = PAR->R_LAMBDA;
        }
    else if (PAR->R_FLAG == 3)
        for (int i = 0; i < (COLUMNS - 2*DIM); i++)
        {
            TIME_MATRIX[(rows + i) * COLUMNS + i] = PAR->R_LAMBDA;
            TIME_MATRIX[(rows + i) * COLUMNS + i + 1*DIM] = -2 * PAR->R_LAMBDA;
            TIME_MATRIX[(rows + i) * COLUMNS + i + 2*DIM] = PAR->R_LAMBDA;
        }

    FILE * pFile;
    pFile = fopen("MSBAS_TIME_MATRIX.txt", "w");
    if (pFile == 0) throw std::invalid_argument("Error: CArea::MakeTimeMatrix3DSPF, pFile == 0. Cannot open file MSBAS_TIME_MATRIX.txt.");
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

    // for C_FLAG=100 locate region with smallest ZSCORE
    if (PAR->C_FLAG == 100)
    {

        WriteLog("\ncomputing location of minimal ZSCORE...");

        float zscore_min = 1e9, cov_min = 0;

        PAR->CSIZE++;
        PAR->LSIZE++;
        
        while (zscore_min == 1e9)
        {
            if (PAR->CSIZE > 0) PAR->CSIZE--;
            if (PAR->LSIZE > 0) PAR->LSIZE--;
            
            WriteLog("\ntrying region width " + i2s(PAR->CSIZE) + " and length " + i2s(PAR->LSIZE) + "...");

            #pragma omp parallel for collapse(2)
            for (int ii = 0; ii < PAR->WIDTH; ii++)
                for (int jj = 0; jj < PAR->LENGTH; jj++)
                {
                    float zscore = 0, cov = 0;
                    long num = 0, num_all = 0;

                    for (int i = (ii - PAR->CSIZE); i <= (ii + PAR->CSIZE); i++)
                        for (int j = (jj - PAR->LSIZE); j <= (jj + PAR->LSIZE); j++)
                        {
                            if ((i >= 0)&&(j >= 0)&&(i < PAR->WIDTH)&&(j < PAR->LENGTH) && ZSCORE_MASK[i + j * PAR->WIDTH] != 0)
                            {
                                zscore += ZSCORE_MASK[i + j * PAR->WIDTH];
                                num++;
                            }
                            num_all++;
                        }
                    if (num > 0 && num_all > 0)
                    {
                        zscore /= num;
                        cov = num / num_all;
                        #pragma omp critical
                        {
                            if (cov > 0.68 && zscore < zscore_min)
                            {
                                zscore_min = zscore;
                                cov_min = cov;
                                PAR->C_FLAG = 1;
                                PAR->CPOS[0] = ii;
                                PAR->LPOS[0] = jj;
                            }
                        }
                    }
                }
        }
        // print absolute coordinates even if WINDOW_SIZE is specified, but in processing use relative coordinates
        WriteLog("\nminimal ZSCORE " + f2s(zscore_min) + " with coverage of " + i2s(int(100 * cov_min)) + "% is found at column and row: " + i2sf(PAR->CPOS[0] + PAR->CSTART) + " " + i2sf(PAR->LPOS[0] + PAR->LSTART));
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
    for (int ii = 0; ii < PAR->WIDTH; ii++)
        for (int jj = 0; jj < PAR->LENGTH; jj++)
        {
            long num = 0;

            for (int i = 0; i < DSET.size(); i++)
                for (int j = 0; j < DSET[i].InSAR.size(); j++)
                {
                    if (num >= 0 && DSET[i].InSAR[j].DATA[ii + jj * PAR->WIDTH] != 0)
                    {
                        ZSCORE_MASK[ii + jj * PAR->WIDTH] += abs(DSET[i].InSAR[j].DATA[ii + jj * PAR->WIDTH] - DSET[i].InSAR[j].MEAN) / DSET[i].InSAR[j].STD;
                        num++;
                    }
                    else
                    {
                        num = -1;
                        ZSCORE_MASK[ii + jj * PAR->WIDTH] = 0;
                    }

                    // remove pixels without dem coverage    
                    if (PAR->DDNS && PAR->DDEW && (PAR->DDNS[ii + jj * PAR->WIDTH] == 0 || PAR->DDEW[ii + jj * PAR->WIDTH] == 0))
                    {
                        num = -1;
                        ZSCORE_MASK[ii + jj * PAR->WIDTH] = 0;
                    }
                }

            if (num > 0) ZSCORE_MASK[ii + jj * PAR->WIDTH] = ZSCORE_MASK[ii + jj * PAR->WIDTH] / num;
        }

    long count = 0;

    // compute percentage of coverage
    #pragma omp parallel for collapse(2) reduction (+:count)
    for (int ii = 0; ii < PAR->WIDTH; ii++)
        for (int jj = 0; jj < PAR->LENGTH; jj++)
            if (ZSCORE_MASK[ii + jj * PAR->WIDTH] != 0)
                count++;

    if (count == 0)
        throw std::invalid_argument(ERROR("no coherent pixels detected, processing cannot continue"));

    WriteLog("\nselected " + i2s(floor(100.0 * count / (PAR->WIDTH * PAR->LENGTH))) + "% coherent pixels for further processing");

    PAR->write_float_data("MSBAS_ZSCORE_MASK", ZSCORE_MASK);
};

void CArea::ReadBinaryFile(string filename, float *dataout)
{
    ifstream in(filename.c_str(), ifstream::in | ifstream::binary);
    if (in.fail()) throw std::invalid_argument("Error: CArea::ReadBinaryFile, in.fail(). Can not open parameter file " + filename + ".");

    char tmpstr[4];
    int i = 0, j = 0, num = 0;

    float *fdata;
    void *v;

    int floatwidth = 4 * PAR->FWIDTH;
    char* data = (char*) calloc(floatwidth, sizeof (char));

    if (!data) throw std::invalid_argument("Error: CArea::ReadBinaryFile, !data. Can not allocate space for data.");

    j = 0;
    while (!in.eof())
    {
        in.read(data, floatwidth);
        num = in.gcount();

        if (num == floatwidth)
        {
            for (i = 0; i < PAR->FWIDTH; i++)
            {
                if (PAR->FORMAT == 1)
                {
                    tmpstr[3] = data[4 * i + 0];
                    tmpstr[2] = data[4 * i + 1];
                    tmpstr[1] = data[4 * i + 2];
                    tmpstr[0] = data[4 * i + 3];
                }
                else if (PAR->FORMAT == 0)
                {
                    tmpstr[0] = data[4 * i + 0];
                    tmpstr[1] = data[4 * i + 1];
                    tmpstr[2] = data[4 * i + 2];
                    tmpstr[3] = data[4 * i + 3];
                }
                else
                    throw std::invalid_argument("Error: CArea::ReadBinaryFile, PAR->FORMAT != 0 or 1. Incorrect value in FORMAT.");

                v = (void*) tmpstr;
                fdata = reinterpret_cast<float*> (v);

                //Set NaN values to 0
                if (*fdata != *fdata) *fdata = 0;

                if ((i >= PAR->CSTART)&&(i <= PAR->CSTOP)&&(j >= PAR->LSTART)&&(j <= PAR->LSTOP))
                    dataout[(i - PAR->CSTART)+(j - PAR->LSTART)*(PAR->CSTOP - PAR->CSTART + 1)] = *fdata;
            }
            j++;
        }
        else if (num != 0)
            throw std::invalid_argument("Error: CArea::ReadBinaryFile, num != floatwidth, Incorrect number of columns in interferogram file detected, check FILE_SIZE.");
    }

    in.close();
    free(data);

};


void CArea::Inversion()
{
    // are there sets with variable geometry?
    bool vargeo = false;
    for (int k = 0; k < DSET.size(); k++) if (DSET[k].LV_THETA && DSET[k].LV_PHI) vargeo = true;
    
    // are there DEM directional derivatives?
    bool demdd = false; 
    if (PAR->DDNS && PAR->DDEW) demdd = true;
    
    if ((DIM == 1) || (DIM == 2 && vargeo == false) || (DIM == 3 && vargeo == false && demdd == false)) // demdd is always true for DIM=4 
        InversionByLine();
    else
        InversionByPixel();
}

void CArea::InversionByLine()
{
    WriteLog("\nstarting SVD inversion (by line) ... ");
    WriteLog("inversion parameters for problem Ax=Y:");
    WriteLog("A: " + i2s(COLUMNS) + " x " + i2s(ROWS));
    WriteLog("x: " + i2s(PAR->WIDTH) + " x " + i2s(COLUMNS));
    WriteLog("Y: " + i2s(PAR->WIDTH) + " x " + i2s(ROWS));
    WriteLog("in summary " + i2s(COLUMNS) + " unknowns & " + i2s(ROWS) + " equations");

    if (COLUMNS > ROWS) throw std::invalid_argument(ERROR("COLUMNS > ROWS, under-determined problem"));

    float_mem_alloc(NORM_X, PAR->WIDTH * PAR->LENGTH);
    float_mem_alloc(NORM_AXY, PAR->WIDTH * PAR->LENGTH);
    float_mem_alloc(RANK, PAR->WIDTH * PAR->LENGTH);
    float_mem_alloc(COND_NUM, PAR->WIDTH * PAR->LENGTH);

    long mem = 0;
    if (DIM == 1)
    {
        for (int k = 0; k < SLC.size(); k++)
        {
            float_mem_alloc(SLC[k].DISP_LOS, PAR->LENGTH * PAR->WIDTH);
        }
        mem = SLC.size() * PAR->WIDTH * PAR->LENGTH;
    }
    else if (DIM == 2)
    {
        for (int k = 0; k < SLC.size(); k++)
        {
            float_mem_alloc(SLC[k].DISP_EW, PAR->LENGTH * PAR->WIDTH);
            float_mem_alloc(SLC[k].DISP_UD, PAR->LENGTH * PAR->WIDTH);
        }
        mem = 2 * SLC.size() * PAR->WIDTH * PAR->LENGTH;
    }
    else if (DIM == 3)
    {
        for (int k = 0; k < SLC.size(); k++)
        {
            float_mem_alloc(SLC[k].DISP_NS, PAR->LENGTH * PAR->WIDTH);
            float_mem_alloc(SLC[k].DISP_EW, PAR->LENGTH * PAR->WIDTH);
            float_mem_alloc(SLC[k].DISP_UD, PAR->LENGTH * PAR->WIDTH);
        }
        mem = 3 * SLC.size() * PAR->WIDTH * PAR->LENGTH;
    }

    WriteLog("successfully allocated " + i2s(mem * sizeof (float) / 1024000) + " MB of RAM for output data");

    // populate time matrix
    
    float *TIME_MATRIX1;
    float_mem_alloc(TIME_MATRIX1, COLUMNS * ROWS);
    memcpy(TIME_MATRIX1, TIME_MATRIX, COLUMNS * ROWS * sizeof (float));

    int rows = 0;

    for (int k = 0; k < DSET.size(); k++)
    {
        // constant geometry
        float s_north = DSET[k].S[0], s_east = DSET[k].S[1], s_up = DSET[k].S[2];

        // populate time matrix
        for (int jj = 0; jj < DSET[k].InSAR.size(); jj++)
        {
            for (int ii = 0; ii < DIM * (SLC.size() - 1); ii = ii + DIM)
            {
                if (DIM == 1)
                {
                    TIME_MATRIX1[ii + 0 + rows * COLUMNS] = TIME_MATRIX[ii + 0 + rows * COLUMNS];

                }
                else if (DIM == 2)
                {
                    TIME_MATRIX1[ii + 0 + rows * COLUMNS] = TIME_MATRIX[ii + 0 + rows * COLUMNS] * s_east;
                    TIME_MATRIX1[ii + 1 + rows * COLUMNS] = TIME_MATRIX[ii + 1 + rows * COLUMNS] * s_up;
                }
                else if (DIM == 3)
                {

                    TIME_MATRIX1[ii + 0 + rows * COLUMNS] = TIME_MATRIX[ii + 0 + rows * COLUMNS] * s_north;
                    TIME_MATRIX1[ii + 1 + rows * COLUMNS] = TIME_MATRIX[ii + 1 + rows * COLUMNS] * s_east;
                    TIME_MATRIX1[ii + 2 + rows * COLUMNS] = TIME_MATRIX[ii + 2 + rows * COLUMNS] * s_up;
                }
            }
            rows++;
        }
    }  

    // determine optimal size of LWORK (from LAPACK)
    sgelss_query_lwork(ROWS, COLUMNS, PAR->WIDTH);
        
    int steps_completed = 0;

    omp_sched_t kind;
    int chunk_size;
    omp_get_schedule(&kind, &chunk_size);

    string schedule = "";
    switch(kind)
    {
        case omp_sched_static:
            schedule = "static";
            break;
        case omp_sched_dynamic:
            schedule = "dynamic";
            break;
        case omp_sched_guided:
            schedule = "guided";
            break;
        case omp_sched_auto:
            schedule = "auto";
            break;
        default:
            schedule = "implementatio nspecific";
            break;
    }

    WriteLog("default schedule type is " + schedule + " " + i2s(chunk_size) + ", changing to dynamic");

    #pragma omp parallel for schedule (dynamic)
    for (int j = 0; j < PAR->LENGTH; j++)
    {
        // allocating memory required for inversion
        float *V;
        float_mem_alloc(V, ROWS * PAR->WIDTH);

        float *AT;
        float_mem_alloc(AT, COLUMNS * ROWS);

        float rowmask = 0;
        
        // re-initialize V, column-major mode
        for (int i = 0; i < PAR->WIDTH; i++)
        {
            for (int k = 0; k < pInSAR.size(); k++) V[k + i * ROWS] = pInSAR[k]->DATA[i + j * PAR->WIDTH];
            rowmask = rowmask + ZSCORE_MASK[i + j * PAR->WIDTH];
        }
        
        // if there are valid pixels in this row than proceed
        if (rowmask != 0)
        {
            // re-initialize AT, column-major mode
            for (int ii = 0; ii < COLUMNS; ii++) for (int jj = 0; jj < ROWS; jj++) AT[jj + ii * ROWS] = TIME_MATRIX1[ii + jj * COLUMNS];

            // perform inversion
            int rank_num=0; float cond_num=0;
            if (sgelss(AT, ROWS, COLUMNS, V, PAR->WIDTH, rank_num, cond_num) != 0) throw std::invalid_argument(ERROR("LAPACK sgelss function signaled error"));

            for (int i = 0; i < PAR->WIDTH; i++)
                if (ZSCORE_MASK[i + j * PAR->WIDTH] != 0)
                {
                    NORM_X[i + j * PAR->WIDTH] = 0;
                    RANK[i + j * PAR->WIDTH] = float (rank_num);
                    COND_NUM[i + j * PAR->WIDTH] = cond_num;

                    // reconstruct time series and compute ||x||
                    for (int k = 0; k < (SLC.size() - 1); k++)
                    {
                        double dt = (SLC[k + 1].DATE - SLC[k].DATE);
                        if (DIM == 1)
                        {
                            if(PAR->V_FLAG == 0)
                            {
                                SLC[k + 1].DISP_LOS[i + j * PAR->WIDTH] = SLC[k].DISP_LOS[i + j * PAR->WIDTH] + V[k + i * ROWS] * dt;
                            }
                            else
                            {
                                SLC[k + 1].DISP_LOS[i + j * PAR->WIDTH] = V[k + i * ROWS];
                            }

                            NORM_X[i + j * PAR->WIDTH] += V[k + i * ROWS] * V[k + i * ROWS];

                            // set to a small number to avoid conflict with NaN==0    
                            if (k == 0)
                                SLC[k].DISP_LOS[i + j * PAR->WIDTH]=SMALL_NUM_FLOAT;
                        }
                        else if (DIM == 2)
                        {
                            if(PAR->V_FLAG == 0)
                            {
                                SLC[k + 1].DISP_EW[i + j * PAR->WIDTH] = SLC[k].DISP_EW[i + j * PAR->WIDTH] + V[2 * k + 0 + i * ROWS] * dt;
                                SLC[k + 1].DISP_UD[i + j * PAR->WIDTH] = SLC[k].DISP_UD[i + j * PAR->WIDTH] + V[2 * k + 1 + i * ROWS] * dt;
                            }
                            else
                            {
                                SLC[k + 1].DISP_EW[i + j * PAR->WIDTH] = V[2 * k + 0 + i * ROWS];
                                SLC[k + 1].DISP_UD[i + j * PAR->WIDTH] = V[2 * k + 1 + i * ROWS];
                            }
                            
                            NORM_X[i + j * PAR->WIDTH] += V[2 * k + i * ROWS] * V[2 * k + i * ROWS] + V[2 * k + 1 + i * ROWS] * V[2 * k + 1 + i * ROWS];

                            // set to a small number to avoid conflict with NaN==0
                            if (k == 0)
                            {
                                SLC[k].DISP_EW[i + j * PAR->WIDTH]=SMALL_NUM_FLOAT;
                                SLC[k].DISP_UD[i + j * PAR->WIDTH]=SMALL_NUM_FLOAT;
                            }
                        }
                        else if (DIM == 3)
                        {
                            if(PAR->V_FLAG == 0)
                            {
                                SLC[k + 1].DISP_NS[i + j * PAR->WIDTH] = SLC[k].DISP_NS[i + j * PAR->WIDTH] + V[3 * k + 0 + i * ROWS] * dt;
                                SLC[k + 1].DISP_EW[i + j * PAR->WIDTH] = SLC[k].DISP_EW[i + j * PAR->WIDTH] + V[3 * k + 1 + i * ROWS] * dt;
                                SLC[k + 1].DISP_UD[i + j * PAR->WIDTH] = SLC[k].DISP_UD[i + j * PAR->WIDTH] + V[3 * k + 2 + i * ROWS] * dt;
                            }
                            else
                            {
                                SLC[k + 1].DISP_NS[i + j * PAR->WIDTH] = V[3 * k + 0 + i * ROWS];
                                SLC[k + 1].DISP_EW[i + j * PAR->WIDTH] = V[3 * k + 1 + i * ROWS];
                                SLC[k + 1].DISP_UD[i + j * PAR->WIDTH] = V[3 * k + 2 + i * ROWS];
                            }
                            
                            NORM_X[i + j * PAR->WIDTH] += V[3 * k + i * ROWS] * V[3 * k + i * ROWS] + V[3 * k + 1 + i * ROWS] * V[3 * k + 1 + i * ROWS] + V[3 * k + 2 + i * ROWS] * V[3 * k + 2 + i * ROWS];

                            // set to a small number to avoid conflict with NaN==0
                            if (k == 0)
                            {
                                SLC[k].DISP_NS[i + j * PAR->WIDTH]=SMALL_NUM_FLOAT;
                                SLC[k].DISP_EW[i + j * PAR->WIDTH]=SMALL_NUM_FLOAT;
                                SLC[k].DISP_UD[i + j * PAR->WIDTH]=SMALL_NUM_FLOAT;
                            }
                        }
                    }
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
                                norm_tmp += TIME_MATRIX1[k + l * COLUMNS] * V[k + i * ROWS];
                            }
                            else if (DIM == 2)
                            {
                                norm_tmp += TIME_MATRIX1[2 * k + l * COLUMNS] * V[2 * k + i * ROWS] + TIME_MATRIX1[2 * k + 1 + l * COLUMNS] * V[2 * k + 1 + i * ROWS];
                            }
                            else if (DIM == 3)
                            {
                                norm_tmp += TIME_MATRIX1[3 * k + l * COLUMNS] * V[3 * k + i * ROWS] + TIME_MATRIX1[3 * k + 1 + l * COLUMNS] * V[3 * k + 1 + i * ROWS] + TIME_MATRIX1[3 * k + 2 + l * COLUMNS] * V[3 * k + 2 + i * ROWS];
                            }

                        norm_tmp -= pInSAR[l]->DATA[i + j * PAR->WIDTH];

                        NORM_AXY[i + j * PAR->WIDTH] += norm_tmp*norm_tmp;
                    }
                    NORM_AXY[i + j * PAR->WIDTH] = sqrt(NORM_AXY[i + j * PAR->WIDTH]);
                }
        }
        
        free(AT);
        free(V);
        
    
        #pragma omp critical
        {
            ++steps_completed;
            if ( (10*steps_completed % PAR->LENGTH) < 10)  WriteLog("completed "+ i2s(long(100.0 * steps_completed / PAR->LENGTH))+"%");
        }
    }

    free(TIME_MATRIX1);
    
    // compute average norms for the entire image
    float norm_x = 0, norm_axy = 0;
    long count = 0;

    for (int j = 0; j < PAR->LENGTH; j++)
        for (int i = 0; i < PAR->WIDTH; i++)
            if (ZSCORE_MASK[i + j * PAR->WIDTH] != 0)
            {
                norm_x += NORM_X[i + j * PAR->WIDTH];
                norm_axy += NORM_AXY[i + j * PAR->WIDTH];
                count++;
            }

    if (count > 0)
    {
        norm_x = norm_x / count;
        norm_axy = norm_axy / count;

        WriteLog("computed ||x|| and ||Ax-Y|| norms for lambda: " + f2s(norm_x) + " " + f2s(norm_axy) + " " + f2s(PAR->R_LAMBDA));
    }
    else
        throw std::invalid_argument(ERROR("no coherent in all interferograms pixels were found"));    
};

void CArea::InversionByPixel()
{
    WriteLog("\nstarting SVD inversion (by pixel) ... ");
    WriteLog("inversion parameters for problem Ax=Y:");
    WriteLog("A: " + i2s(COLUMNS) + " x " + i2s(ROWS));
    WriteLog("x: " + i2s(PAR->WIDTH) + " x " + i2s(COLUMNS));
    WriteLog("Y: " + i2s(PAR->WIDTH) + " x " + i2s(ROWS));
    WriteLog("in summary " + i2s(COLUMNS) + " unknowns & " + i2s(ROWS) + " equations");

    if (COLUMNS > ROWS) throw std::invalid_argument(ERROR("COLUMNS > ROWS, under-determined problem"));

    float_mem_alloc(NORM_X, PAR->WIDTH * PAR->LENGTH);
    float_mem_alloc(NORM_AXY, PAR->WIDTH * PAR->LENGTH);
    float_mem_alloc(RANK, PAR->WIDTH * PAR->LENGTH);
    float_mem_alloc(COND_NUM, PAR->WIDTH * PAR->LENGTH);

    long mem = 0;
    if (DIM == 1)
    {
        for (int k = 0; k < SLC.size(); k++)
        {
            float_mem_alloc(SLC[k].DISP_LOS, PAR->LENGTH * PAR->WIDTH);
        }
        mem = SLC.size() * PAR->WIDTH * PAR->LENGTH;
    }
    else if (DIM == 2)
    {
        for (int k = 0; k < SLC.size(); k++)
        {
            float_mem_alloc(SLC[k].DISP_EW, PAR->LENGTH * PAR->WIDTH);
            float_mem_alloc(SLC[k].DISP_UD, PAR->LENGTH * PAR->WIDTH);
        }
        mem = 2 * SLC.size() * PAR->WIDTH * PAR->LENGTH;
    }
    else if (DIM == 3)
    {
        for (int k = 0; k < SLC.size(); k++)
        {
            float_mem_alloc(SLC[k].DISP_NS, PAR->LENGTH * PAR->WIDTH);
            float_mem_alloc(SLC[k].DISP_EW, PAR->LENGTH * PAR->WIDTH);
            float_mem_alloc(SLC[k].DISP_UD, PAR->LENGTH * PAR->WIDTH);
        }
        mem = 3 * SLC.size() * PAR->WIDTH * PAR->LENGTH;
    }
    else if (DIM == 4)
    {
        for (int k = 0; k < SLC.size(); k++)
        {
            float_mem_alloc(SLC[k].DISP_NS, PAR->LENGTH * PAR->WIDTH);
            float_mem_alloc(SLC[k].DISP_EW, PAR->LENGTH * PAR->WIDTH);
            float_mem_alloc(SLC[k].DISP_UD, PAR->LENGTH * PAR->WIDTH);
            float_mem_alloc(SLC[k].DISP_UD1, PAR->LENGTH * PAR->WIDTH);
        }
        mem = 4 * SLC.size() * PAR->WIDTH * PAR->LENGTH;
    }

    WriteLog("successfully allocated " + i2s(mem * sizeof (float) / 1024000) + " MB of RAM for output data");

    // determine optimal size of LWORK (from LAPACK)
    sgelss_query_lwork(ROWS, COLUMNS, 1);
    
    omp_sched_t kind;
    int chunk_size;
    omp_get_schedule(&kind, &chunk_size);

    string schedule = "";
    switch(kind)
    {
        case omp_sched_static:
            schedule = "static";
            break;
        case omp_sched_dynamic:
            schedule = "dynamic";
            break;
        case omp_sched_guided:
            schedule = "guided";
            break;
        case omp_sched_auto:
            schedule = "auto";
            break;
        default:
            schedule = "implementatio nspecific";
            break;
    }

    WriteLog("default schedule type is " + schedule + " " + i2s(chunk_size) + ", changing to dynamic");

    int steps_completed = 0;
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int j = 0; j < PAR->LENGTH; j++)
        for (int i = 0; i < PAR->WIDTH; i++)
        {
            if (ZSCORE_MASK[i + j * PAR->WIDTH] != 0)
            {
                float *TIME_MATRIX1;
                float_mem_alloc(TIME_MATRIX1, COLUMNS * ROWS);
                memcpy(TIME_MATRIX1, TIME_MATRIX, COLUMNS * ROWS * sizeof (float));

                // allocating memory required for inversion
                float *V;
                float_mem_alloc(V, ROWS);

                // re-initialize V, column-major mode
                for (int k = 0; k < pInSAR.size(); k++) V[k] = pInSAR[k]->DATA[i + j * PAR->WIDTH];

                float *AT;
                float_mem_alloc(AT, COLUMNS * ROWS);

                int rows = 0;

                for (int k = 0; k < DSET.size(); k++)
                {
                    // constant geometry
                    float s_north = DSET[k].S[0], s_east = DSET[k].S[1], s_up = DSET[k].S[2];

                    // variable geometry, range set
                    if (DSET[k].LV_THETA && DSET[k].LV_PHI && DSET[k].TYPE == 0)
                    {

                        s_north = cos(DSET[k].LV_THETA[PAR->WIDTH * j + i]) * sin(DSET[k].LV_PHI[PAR->WIDTH * j + i]);
                        s_east = cos(DSET[k].LV_THETA[PAR->WIDTH * j + i]) * cos(DSET[k].LV_PHI[PAR->WIDTH * j + i]);
                        s_up = sin(DSET[k].LV_THETA[PAR->WIDTH * j + i]);
                    }
                    // variable geometry, azimuth set
                    else if (DSET[k].LV_THETA && DSET[k].LV_PHI && DSET[k].TYPE == 1)
                    {
                        s_north = -cos(DSET[k].LV_PHI[PAR->WIDTH * j + i]);
                        s_east = sin(DSET[k].LV_PHI[PAR->WIDTH * j + i]);
                        s_up = 0;
                    }

                    // populate time matrix
                    for (int jj = 0; jj < DSET[k].InSAR.size(); jj++)
                    {
                        for (int ii = 0; ii < DIM * (SLC.size() - 1); ii = ii + DIM)
                        {
                            if (DIM == 1)
                            {
                                TIME_MATRIX1[ii + 0 + rows * COLUMNS] = TIME_MATRIX[ii + 0 + rows * COLUMNS];

                            }
                            else if (DIM == 2)
                            {
                                TIME_MATRIX1[ii + 0 + rows * COLUMNS] = TIME_MATRIX[ii + 0 + rows * COLUMNS] * s_east;
                                TIME_MATRIX1[ii + 1 + rows * COLUMNS] = TIME_MATRIX[ii + 1 + rows * COLUMNS] * s_up;
                            }
                            else if (DIM == 3)
                            {
                                TIME_MATRIX1[ii + 0 + rows * COLUMNS] = TIME_MATRIX[ii + 0 + rows * COLUMNS] * s_north;
                                TIME_MATRIX1[ii + 1 + rows * COLUMNS] = TIME_MATRIX[ii + 1 + rows * COLUMNS] * s_east;
                                TIME_MATRIX1[ii + 2 + rows * COLUMNS] = TIME_MATRIX[ii + 2 + rows * COLUMNS] * s_up;
                            }
                            else if (DIM == 4)
                            {
                                TIME_MATRIX1[ii + 0 + rows * COLUMNS] = TIME_MATRIX[ii + 0 + rows * COLUMNS] * s_north;
                                TIME_MATRIX1[ii + 1 + rows * COLUMNS] = TIME_MATRIX[ii + 1 + rows * COLUMNS] * s_east;
                                TIME_MATRIX1[ii + 2 + rows * COLUMNS] = TIME_MATRIX[ii + 2 + rows * COLUMNS] * s_up;
                                TIME_MATRIX1[ii + 3 + rows * COLUMNS] = TIME_MATRIX[ii + 3 + rows * COLUMNS] * s_up;
                            }
                        }
                        rows++;
                    }
                }  
                
                if (PAR->DDNS && PAR->DDEW) // watch for a case with DIM=3 and DD files provided (i.e. still 3D but with the constraint)
                    for (int jj = rows; jj < (rows + SLC.size() - 1); jj++)
                    {
                        TIME_MATRIX1[jj * COLUMNS + DIM * (jj - rows) + 0] = PAR->DDNS[i + j * PAR->WIDTH];
                        TIME_MATRIX1[jj * COLUMNS + DIM * (jj - rows) + 1] = PAR->DDEW[i + j * PAR->WIDTH];
                        TIME_MATRIX1[jj * COLUMNS + DIM * (jj - rows) + 2] = -1;
                    } 
                
                    /*
                    if (i==0 && j==0)
                    {
                        FILE * pFile;
                        pFile = fopen("A.txt", "w");
                        if (pFile == 0) throw std::invalid_argument("Error: CArea::MakeTimeMatrix2D, pFile == 0. Cannot open file MSBAS_TIME_MATRIX.txt.");
                        for (int jj = 0; jj < ROWS; jj++)
                        {
                            for (int ii = 0; ii < COLUMNS; ii++)
                            {
                                fprintf(pFile, "%11.6f ", TIME_MATRIX1[ii + COLUMNS * jj]);
                            }
                            fprintf(pFile, "\n");
                        }
                        fclose(pFile);

                        pFile = fopen("Y.txt", "w");
                        if (pFile == 0) throw std::invalid_argument("Error: CArea::MakeTimeMatrix2D, pFile == 0. Cannot open file MSBAS_TIME_MATRIX.txt.");
                        for (int kk = 0; kk < ROWS; kk++) fprintf(pFile, "%11.6f \n", V[kk]);
                        fclose(pFile);
                        sleep(20);
                    }
                    */
                    
                
                // re-initialize AT, column-major mode
                for (int ii = 0; ii < COLUMNS; ii++) for (int jj = 0; jj < ROWS; jj++) AT[jj + ii * ROWS] = TIME_MATRIX1[ii + jj * COLUMNS];

                // perform inversion
                int rank_num=0; float cond_num=0;
                if (sgelss(AT, ROWS, COLUMNS, V, 1, rank_num, cond_num) != 0) throw std::invalid_argument(ERROR("LAPACK sgelss function signaled error"));

                // reconstruct time series and compute norm ||x||

                NORM_X[i + j * PAR->WIDTH] = 0;
                RANK[i + j * PAR->WIDTH] = float (rank_num);
                COND_NUM[i + j * PAR->WIDTH] = cond_num;

                for (int k = 0; k < (SLC.size() - 1); k++)
                {
                    double dt = (SLC[k + 1].DATE - SLC[k].DATE);
                    if (DIM == 1)
                    {
                        if(PAR->V_FLAG == 0)
                        {
                            SLC[k + 1].DISP_LOS[i + j * PAR->WIDTH] = SLC[k].DISP_LOS[i + j * PAR->WIDTH] + V[DIM * k + 0 + i] * dt;
                        }
                        else
                        {
                            SLC[k + 1].DISP_LOS[i + j * PAR->WIDTH] = V[DIM * k + 0 + i];
                        }
                         
                        NORM_X[i + j * PAR->WIDTH] += V[DIM * k + 0 + i] * V[DIM * k + 0 + i];

                        // set to a small number to avoid conflict with NaN==0    
                        if (k == 0)
                            SLC[k].DISP_LOS[i + j * PAR->WIDTH] = SMALL_NUM_FLOAT;
                    }
                    else if (DIM == 2)
                    {
                        if(PAR->V_FLAG == 0)
                        {
                            SLC[k + 1].DISP_EW[i + j * PAR->WIDTH] = SLC[k].DISP_EW[i + j * PAR->WIDTH] + V[DIM * k + 0 + i] * dt;
                            SLC[k + 1].DISP_UD[i + j * PAR->WIDTH] = SLC[k].DISP_UD[i + j * PAR->WIDTH] + V[DIM * k + 1 + i] * dt;
                        }
                        else
                        {
                            SLC[k + 1].DISP_EW[i + j * PAR->WIDTH] = V[DIM * k + 0 + i];
                            SLC[k + 1].DISP_UD[i + j * PAR->WIDTH] = V[DIM * k + 1 + i];
                        }                            
                        NORM_X[i + j * PAR->WIDTH] += V[DIM * k + 0 + i] * V[DIM * k + 0 + i] + V[DIM * k + 1 + i] * V[DIM * k + 1 + i];

                        // set to a small number to avoid conflict with NaN==0
                        if (k == 0)
                        {
                            SLC[k].DISP_EW[i + j * PAR->WIDTH] = SMALL_NUM_FLOAT;
                            SLC[k].DISP_UD[i + j * PAR->WIDTH] = SMALL_NUM_FLOAT;
                        }
                    }
                    else if (DIM == 3)
                    {
                        if(PAR->V_FLAG == 0)
                        {
                            SLC[k + 1].DISP_NS[i + j * PAR->WIDTH] = SLC[k].DISP_NS[i + j * PAR->WIDTH] + V[DIM * k + 0] * dt;
                            SLC[k + 1].DISP_EW[i + j * PAR->WIDTH] = SLC[k].DISP_EW[i + j * PAR->WIDTH] + V[DIM * k + 1] * dt;
                            SLC[k + 1].DISP_UD[i + j * PAR->WIDTH] = SLC[k].DISP_UD[i + j * PAR->WIDTH] + V[DIM * k + 2] * dt;
                        }
                        else
                        {
                            SLC[k + 1].DISP_NS[i + j * PAR->WIDTH] = V[DIM * k + 0];
                            SLC[k + 1].DISP_EW[i + j * PAR->WIDTH] = V[DIM * k + 1];
                            SLC[k + 1].DISP_UD[i + j * PAR->WIDTH] = V[DIM * k + 2];
                        }

                        NORM_X[i + j * PAR->WIDTH] += V[DIM * k + 0] * V[DIM * k + 0] + V[DIM * k + 1] * V[DIM * k + 1] + V[DIM * k + 2] * V[DIM * k + 2];

                        // set to a small number to avoid conflict with NaN==0
                        if (k == 0)
                        {
                            SLC[k].DISP_NS[i + j * PAR->WIDTH] = SMALL_NUM_FLOAT;
                            SLC[k].DISP_EW[i + j * PAR->WIDTH] = SMALL_NUM_FLOAT;
                            SLC[k].DISP_UD[i + j * PAR->WIDTH] = SMALL_NUM_FLOAT;
                        }
                    }
                    else if (DIM == 4)
                    {
                        if(PAR->V_FLAG == 0)
                        {
                            SLC[k + 1].DISP_NS[i + j * PAR->WIDTH] = SLC[k].DISP_NS[i + j * PAR->WIDTH] + V[DIM * k + 0] * dt;
                            SLC[k + 1].DISP_EW[i + j * PAR->WIDTH] = SLC[k].DISP_EW[i + j * PAR->WIDTH] + V[DIM * k + 1] * dt;
                            SLC[k + 1].DISP_UD[i + j * PAR->WIDTH] = SLC[k].DISP_UD[i + j * PAR->WIDTH] + V[DIM * k + 2] * dt;
                            SLC[k + 1].DISP_UD1[i + j * PAR->WIDTH] = SLC[k].DISP_UD1[i + j * PAR->WIDTH] + V[DIM * k + 3] * dt;
                        }
                        else
                        {
                            SLC[k + 1].DISP_NS[i + j * PAR->WIDTH] = V[DIM * k + 0];
                            SLC[k + 1].DISP_EW[i + j * PAR->WIDTH] = V[DIM * k + 1];
                            SLC[k + 1].DISP_UD[i + j * PAR->WIDTH] = V[DIM * k + 2];
                            SLC[k + 1].DISP_UD1[i + j * PAR->WIDTH] = V[DIM * k + 3];
                        }

                        NORM_X[i + j * PAR->WIDTH] += V[DIM * k + 0] * V[DIM * k + 0] + V[DIM * k + 1] * V[DIM * k + 1] + V[DIM * k + 2] * V[DIM * k + 2] + V[DIM * k + 3] * V[DIM * k + 3];

                        // set to a small number to avoid conflict with NaN==0
                        if (k == 0)
                        {
                            SLC[k].DISP_NS[i + j * PAR->WIDTH] = SMALL_NUM_FLOAT;
                            SLC[k].DISP_EW[i + j * PAR->WIDTH] = SMALL_NUM_FLOAT;
                            SLC[k].DISP_UD[i + j * PAR->WIDTH] = SMALL_NUM_FLOAT;
                            SLC[k].DISP_UD1[i + j * PAR->WIDTH] = SMALL_NUM_FLOAT;
                        }
                    }

                }
                NORM_X[i + j * PAR->WIDTH] = sqrt(NORM_X[i + j * PAR->WIDTH]);

                // compute norm ||Ax-Y||
                NORM_AXY[i + j * PAR->WIDTH] = 0;
                float norm_axy = 0;

                for (int l = 0; l < pInSAR.size(); l++)
                {
                    norm_axy = pInSAR[l]->DATA[i + j * PAR->WIDTH];
                    for (int k = 0; k < (SLC.size() - 1); k++)
                        if (DIM == 1)
                            norm_axy -= (TIME_MATRIX1[DIM * k + 0 + l * COLUMNS] * V[DIM * k + 0]);
                        else if (DIM == 2)
                            norm_axy -= (TIME_MATRIX1[DIM * k + 0 + l * COLUMNS] * V[DIM * k + 0] + TIME_MATRIX1[DIM * k + 1 + l * COLUMNS] * V[DIM * k + 1]);
                        else if (DIM == 3)
                            norm_axy -= (TIME_MATRIX1[DIM * k + 0 + l * COLUMNS] * V[DIM * k + 0] + TIME_MATRIX1[DIM * k + 1 + l * COLUMNS] * V[DIM * k + 1] + TIME_MATRIX1[DIM * k + 2 + l * COLUMNS] * V[DIM * k + 2]);
                        else if (DIM == 4)
                            norm_axy -= (TIME_MATRIX1[DIM * k + 0 + l * COLUMNS] * V[DIM * k + 0] + TIME_MATRIX1[DIM * k + 1 + l * COLUMNS] * V[DIM * k + 1] + TIME_MATRIX1[DIM * k + 2 + l * COLUMNS] * V[DIM * k + 2] + TIME_MATRIX1[DIM * k + 3 + l * COLUMNS] * V[DIM * k + 3]);

                    NORM_AXY[i + j * PAR->WIDTH] += norm_axy*norm_axy;
                }

                NORM_AXY[i + j * PAR->WIDTH] = sqrt(NORM_AXY[i + j * PAR->WIDTH]);

                free(TIME_MATRIX1);
                free(AT);
                free(V);
            }
            #pragma omp atomic
            steps_completed++;

            float local_steps_completed=steps_completed;
            float float_completed0 = 100 * (float) (local_steps_completed-1) / (float) (PAR->LENGTH*PAR->WIDTH);
            float float_completed = 100 * (float) (local_steps_completed) / (float) (PAR->LENGTH*PAR->WIDTH);
            if (floor(float_completed)>floor(float_completed0))
            {
                time_t now = time(0);
                tm* localtm = localtime(&now);
                WriteLog("completed "+ i2s(int(floor(float_completed)))+"% on " + asctime(localtm));
            }
    }

    WriteLog("SVD inversion completed successfully");

    // compute average norms for the entire image
    float norm_x = 0, norm_axy = 0;
    long count = 0;

    for (int j = 0; j < PAR->LENGTH; j++)
        for (int i = 0; i < PAR->WIDTH; i++)
            if (ZSCORE_MASK[i + j * PAR->WIDTH] != 0)
            {
                norm_x += NORM_X[i + j * PAR->WIDTH];
                norm_axy += NORM_AXY[i + j * PAR->WIDTH];
                count++;
            }

    if (count > 0)
    {
        norm_x = norm_x / count;
        norm_axy = norm_axy / count;

        WriteLog("computed ||x|| and ||Ax-Y|| norms for lambda: " + f2s(norm_x) + " " + f2s(norm_axy) + " " + f2s(PAR->R_LAMBDA));
    }
    else
        throw std::invalid_argument(ERROR("no coherent in all interferograms pixels were found"));
};

void CArea::PostProcessing()
{
    WriteLog("\nstarting post processing... ");

    ComputeLinearRate();

    WriteResultsToDisk();

    InteractiveMode();

};

void CArea::WriteResultsToDisk()
{
    WriteLog("writing results to a disk");

    // create file that will be used in post-processing for extracting time series from binary MSBAS output
    FILE * pFile0 = fopen("MSBAS_TSOUT.txt", "w");
    if (!pFile0) throw std::invalid_argument(ERROR("can not write to file MSBAS_TSOUT.txt"));
    if (PAR->FORMAT == 2)
        fprintf(pFile0, "%s=%i,%i\n", "FORMAT", PAR->FORMAT, PAR->FORMAT_INT_RADIUS);
    else
        fprintf(pFile0, "%s=%i\n", "FORMAT", PAR->FORMAT);

    fprintf(pFile0, "%s=%i,%i\n", "FILE_SIZE", PAR->FWIDTH, PAR->FLENGTH);
    fprintf(pFile0, "%s=%i", "C_FLAG", PAR->C_FLAG);
    if ((PAR->C_FLAG > 0)&&(PAR->C_FLAG < 10))
    {
        for (int i = 0; i < PAR->C_FLAG; i++)
            fprintf(pFile0, ", %i, %i", PAR->CPOS[i] + PAR->CSTART, PAR->LPOS[i] + PAR->LSTART);
        fprintf(pFile0, ", %i, %i\n", PAR->CSIZE, PAR->LSIZE);
    }
    else
        fprintf(pFile0, "\n");
    
    for (int k = 0; k < SLC.size(); k++)
    {
        if (DIM == 1)
        {
            string f1 = "MSBAS_" + SLC[k].NAME + "_LOS";
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
            string f1 = "MSBAS_" + SLC[k].NAME + "_EW";
            string f2 = "MSBAS_" + SLC[k].NAME + "_UD";

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
        else if (DIM == 3)
        {
            string f1 = "MSBAS_" + SLC[k].NAME + "_NS";
            string f2 = "MSBAS_" + SLC[k].NAME + "_EW";
            string f3 = "MSBAS_" + SLC[k].NAME + "_UD";

            if (PAR->FORMAT == 0 || PAR->FORMAT == 1)
            {
                f1 += ".bin";
                f2 += ".bin";
                f3 += ".bin";
            }
            else if (PAR->FORMAT == 2)
            {
                f1 += ".tif";
                f2 += ".tif";
                f3 += ".tif";
            }
            fprintf(pFile0, "%s %f %s %s %s\n", SLC[k].NAME.c_str(), SLC[k].DATE, f1.c_str(), f2.c_str(), f3.c_str());
        }
        else if (DIM == 4)
        {
            string f1 = "MSBAS_" + SLC[k].NAME + "_NS";
            string f2 = "MSBAS_" + SLC[k].NAME + "_EW";
            string f3 = "MSBAS_" + SLC[k].NAME + "_UD";
            string f4 = "MSBAS_" + SLC[k].NAME + "_UD1";

            if (PAR->FORMAT == 0 || PAR->FORMAT == 1)
            {
                f1 += ".bin";
                f2 += ".bin";
                f3 += ".bin";
                f4 += ".bin";
            }
            else if (PAR->FORMAT == 2)
            {
                f1 += ".tif";
                f2 += ".tif";
                f3 += ".tif";
                f4 += ".tif";
            }
            fprintf(pFile0, "%s %f %s %s %s %s\n", SLC[k].NAME.c_str(), SLC[k].DATE, f1.c_str(), f2.c_str(), f3.c_str(), f4.c_str());
        }
    }
    fclose(pFile0);

    //#pragma omp parallel for - disabled parallel writing, seems like a bug in GDAL
    for (int k = 0; k < SLC.size(); k++) SLC[k].WriteData();

    if (DIM == 1)
    {
        if (RATE_LOS) PAR->write_float_data("MSBAS_LINEAR_RATE_LOS", RATE_LOS);
        if (RATE_STD_LOS) PAR->write_float_data("MSBAS_LINEAR_RATE_STD_LOS", RATE_STD_LOS);
        if (RATE_R2_LOS) PAR->write_float_data("MSBAS_LINEAR_RATE_R2_LOS", RATE_R2_LOS);
    }
    else if (DIM == 2)
    {
        if (RATE_EW) PAR->write_float_data("MSBAS_LINEAR_RATE_EW", RATE_EW);
        if (RATE_STD_EW) PAR->write_float_data("MSBAS_LINEAR_RATE_STD_EW", RATE_STD_EW);
        if (RATE_R2_EW) PAR->write_float_data("MSBAS_LINEAR_RATE_R2_EW", RATE_R2_EW);
        if (RATE_UD) PAR->write_float_data("MSBAS_LINEAR_RATE_UD", RATE_UD);
        if (RATE_STD_UD) PAR->write_float_data("MSBAS_LINEAR_RATE_STD_UD", RATE_STD_UD);
        if (RATE_R2_UD) PAR->write_float_data("MSBAS_LINEAR_RATE_R2_UD", RATE_R2_UD);
    }
    else if (DIM == 3)
    {
        if (RATE_NS) PAR->write_float_data("MSBAS_LINEAR_RATE_NS", RATE_NS);
        if (RATE_STD_NS) PAR->write_float_data("MSBAS_LINEAR_RATE_STD_NS", RATE_STD_NS);
        if (RATE_R2_NS) PAR->write_float_data("MSBAS_LINEAR_RATE_R2_NS", RATE_R2_NS);
        if (RATE_EW) PAR->write_float_data("MSBAS_LINEAR_RATE_EW", RATE_EW);
        if (RATE_STD_EW) PAR->write_float_data("MSBAS_LINEAR_RATE_STD_EW", RATE_STD_EW);
        if (RATE_R2_EW) PAR->write_float_data("MSBAS_LINEAR_RATE_R2_EW", RATE_R2_EW);
        if (RATE_UD) PAR->write_float_data("MSBAS_LINEAR_RATE_UD", RATE_UD);
        if (RATE_STD_UD) PAR->write_float_data("MSBAS_LINEAR_RATE_STD_UD", RATE_STD_UD);
        if (RATE_R2_UD) PAR->write_float_data("MSBAS_LINEAR_RATE_R2_UD", RATE_R2_UD);
    }
    else if (DIM == 4)
    {
        if (RATE_NS) PAR->write_float_data("MSBAS_LINEAR_RATE_NS", RATE_NS);
        if (RATE_STD_NS) PAR->write_float_data("MSBAS_LINEAR_RATE_STD_NS", RATE_STD_NS);
        if (RATE_R2_NS) PAR->write_float_data("MSBAS_LINEAR_RATE_R2_NS", RATE_R2_NS);
        if (RATE_EW) PAR->write_float_data("MSBAS_LINEAR_RATE_EW", RATE_EW);
        if (RATE_STD_EW) PAR->write_float_data("MSBAS_LINEAR_RATE_STD_EW", RATE_STD_EW);
        if (RATE_R2_EW) PAR->write_float_data("MSBAS_LINEAR_RATE_R2_EW", RATE_R2_EW);
        if (RATE_UD) PAR->write_float_data("MSBAS_LINEAR_RATE_UD", RATE_UD);
        if (RATE_STD_UD) PAR->write_float_data("MSBAS_LINEAR_RATE_STD_UD", RATE_STD_UD);
        if (RATE_R2_UD) PAR->write_float_data("MSBAS_LINEAR_RATE_R2_UD", RATE_R2_UD);
        if (RATE_UD1) PAR->write_float_data("MSBAS_LINEAR_RATE_UD1", RATE_UD1);
        if (RATE_STD_UD1) PAR->write_float_data("MSBAS_LINEAR_RATE_STD_UD1", RATE_STD_UD1);
        if (RATE_R2_UD1) PAR->write_float_data("MSBAS_LINEAR_RATE_R2_UD1", RATE_R2_UD1);
    }

    PAR->write_float_data("MSBAS_NORM_X", NORM_X);
    PAR->write_float_data("MSBAS_NORM_AXY", NORM_AXY);
    PAR->write_float_data("MSBAS_RANK", RANK);
    PAR->write_float_data("MSBAS_COND_NUM", COND_NUM);
};

void CArea::ComputeLinearRate()
{
    WriteLog("computing linear rate(s)");

    if (DIM == 1)
    {
        float_mem_alloc(RATE_LOS, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_STD_LOS, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_R2_LOS, PAR->WIDTH * PAR->LENGTH);
    }
    else if (DIM == 2)
    {
        float_mem_alloc(RATE_EW, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_STD_EW, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_R2_EW, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_UD, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_STD_UD, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_R2_UD, PAR->WIDTH * PAR->LENGTH);
    }
    else if (DIM == 3)
    {
        float_mem_alloc(RATE_NS, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_STD_NS, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_R2_NS, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_EW, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_STD_EW, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_R2_EW, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_UD, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_STD_UD, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_R2_UD, PAR->WIDTH * PAR->LENGTH);
    }
    else if (DIM == 4)
    {
        float_mem_alloc(RATE_NS, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_STD_NS, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_R2_NS, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_EW, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_STD_EW, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_R2_EW, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_UD, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_STD_UD, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_R2_UD, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_UD1, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_STD_UD1, PAR->WIDTH * PAR->LENGTH);
        float_mem_alloc(RATE_R2_UD1, PAR->WIDTH * PAR->LENGTH);
    }

#pragma omp parallel for collapse(2)
    for (int j = 0; j < PAR->LENGTH; j++)
        for (int i = 0; i < PAR->WIDTH; i++)
        {
            float *x, *y;
            float_mem_alloc(x, SLC.size());
            float_mem_alloc(y, SLC.size());

            float Xa, Xb, Xsigmaa, Xsigmab, r2;

            if (DIM == 1 && ZSCORE_MASK[i + j * PAR->WIDTH] != 0)
            {
                for (int k = 0; k < SLC.size(); k++)
                {
                    x[k] = SLC[k].DATE;
                    y[k] = SLC[k].DISP_LOS[i + j * PAR->WIDTH];
                }

                linreg(x, y, SLC.size(), Xa, Xb, Xsigmaa, Xsigmab, r2);
                RATE_LOS[i + j * PAR->WIDTH] = Xb;
                RATE_STD_LOS[i + j * PAR->WIDTH] = Xsigmab;
                RATE_R2_LOS[i + j * PAR->WIDTH] = r2;
            }
            else if (DIM == 2 && ZSCORE_MASK[i + j * PAR->WIDTH] != 0)
            {
                for (int k = 0; k < SLC.size(); k++)
                {
                    x[k] = SLC[k].DATE;
                    y[k] = SLC[k].DISP_EW[i + j * PAR->WIDTH];
                }

                linreg(x, y, SLC.size(), Xa, Xb, Xsigmaa, Xsigmab, r2);
                RATE_EW[i + j * PAR->WIDTH] = Xb;
                RATE_STD_EW[i + j * PAR->WIDTH] = Xsigmab;
                RATE_R2_EW[i + j * PAR->WIDTH] = r2;

                for (int k = 0; k < SLC.size(); k++)
                {
                    x[k] = SLC[k].DATE;
                    y[k] = SLC[k].DISP_UD[i + j * PAR->WIDTH];
                }

                linreg(x, y, SLC.size(), Xa, Xb, Xsigmaa, Xsigmab, r2);
                RATE_UD[i + j * PAR->WIDTH] = Xb;
                RATE_STD_UD[i + j * PAR->WIDTH] = Xsigmab;
                RATE_R2_UD[i + j * PAR->WIDTH] = r2;
            }
            else if (DIM == 3 && ZSCORE_MASK[i + j * PAR->WIDTH] != 0)
            {
                for (int k = 0; k < SLC.size(); k++)
                {
                    x[k] = SLC[k].DATE;
                    y[k] = SLC[k].DISP_NS[i + j * PAR->WIDTH];
                }

                linreg(x, y, SLC.size(), Xa, Xb, Xsigmaa, Xsigmab, r2);
                RATE_NS[i + j * PAR->WIDTH] = Xb;
                RATE_STD_NS[i + j * PAR->WIDTH] = Xsigmab;
                RATE_R2_NS[i + j * PAR->WIDTH] = r2;

                for (int k = 0; k < SLC.size(); k++)
                {
                    x[k] = SLC[k].DATE;
                    y[k] = SLC[k].DISP_EW[i + j * PAR->WIDTH];
                }

                linreg(x, y, SLC.size(), Xa, Xb, Xsigmaa, Xsigmab, r2);
                RATE_EW[i + j * PAR->WIDTH] = Xb;
                RATE_STD_EW[i + j * PAR->WIDTH] = Xsigmab;
                RATE_R2_EW[i + j * PAR->WIDTH] = r2;

                for (int k = 0; k < SLC.size(); k++)
                {
                    x[k] = SLC[k].DATE;
                    y[k] = SLC[k].DISP_UD[i + j * PAR->WIDTH];
                }

                linreg(x, y, SLC.size(), Xa, Xb, Xsigmaa, Xsigmab, r2);
                RATE_UD[i + j * PAR->WIDTH] = Xb;
                RATE_STD_UD[i + j * PAR->WIDTH] = Xsigmab;
                RATE_R2_UD[i + j * PAR->WIDTH] = r2;
            }
            else if (DIM == 4 && ZSCORE_MASK[i + j * PAR->WIDTH] != 0)
            {
                for (int k = 0; k < SLC.size(); k++)
                {
                    x[k] = SLC[k].DATE;
                    y[k] = SLC[k].DISP_NS[i + j * PAR->WIDTH];
                }

                linreg(x, y, SLC.size(), Xa, Xb, Xsigmaa, Xsigmab, r2);
                RATE_NS[i + j * PAR->WIDTH] = Xb;
                RATE_STD_NS[i + j * PAR->WIDTH] = Xsigmab;
                RATE_R2_NS[i + j * PAR->WIDTH] = r2;

                for (int k = 0; k < SLC.size(); k++)
                {
                    x[k] = SLC[k].DATE;
                    y[k] = SLC[k].DISP_EW[i + j * PAR->WIDTH];
                }

                linreg(x, y, SLC.size(), Xa, Xb, Xsigmaa, Xsigmab, r2);
                RATE_EW[i + j * PAR->WIDTH] = Xb;
                RATE_STD_EW[i + j * PAR->WIDTH] = Xsigmab;
                RATE_R2_EW[i + j * PAR->WIDTH] = r2;

                for (int k = 0; k < SLC.size(); k++)
                {
                    x[k] = SLC[k].DATE;
                    y[k] = SLC[k].DISP_UD[i + j * PAR->WIDTH];
                }

                linreg(x, y, SLC.size(), Xa, Xb, Xsigmaa, Xsigmab, r2);
                RATE_UD[i + j * PAR->WIDTH] = Xb;
                RATE_STD_UD[i + j * PAR->WIDTH] = Xsigmab;
                RATE_R2_UD[i + j * PAR->WIDTH] = r2;

                for (int k = 0; k < SLC.size(); k++)
                {
                    x[k] = SLC[k].DATE;
                    y[k] = SLC[k].DISP_UD1[i + j * PAR->WIDTH];
                }

                linreg(x, y, SLC.size(), Xa, Xb, Xsigmaa, Xsigmab, r2);
                RATE_UD1[i + j * PAR->WIDTH] = Xb;
                RATE_STD_UD1[i + j * PAR->WIDTH] = Xsigmab;
                RATE_R2_UD1[i + j * PAR->WIDTH] = r2;

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
    float a3 = 0;
    float a4 = 0;
    float sd1 = 0;
    float sd2 = 0;
    float sd3 = 0;
    float sd4 = 0;

    if (PAR->I_FLAG == 1)
    {
        WriteLog("starting interactive mode");
        int flag = 1;
        while (flag == 1)
        {
            printf("\nenter five space delimited parameters - region location [column row], half-dimensions [column_radius row_radius] and repeat flag [0-no, 1-yes], e.g. for a region centered at [10 20] with dimensions of [2*5+1 2*8+1] enter '10 20 5 8 1': ");
            int retval = scanf("%i %i %i %i %i", &cref, &lref, &cradius, &lradius, &flag);
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
                    a3 = 0;
                    a4 = 0;
                    sd1 = 0;
                    sd2 = 0;
                    sd3 = 0;
                    sd4 = 0;

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
                                else if (DIM == 3)
                                {
                                    a1 = a1 + SLC[k].DISP_NS[ii + jj * PAR->WIDTH];
                                    a2 = a2 + SLC[k].DISP_EW[ii + jj * PAR->WIDTH];
                                    a3 = a3 + SLC[k].DISP_UD[ii + jj * PAR->WIDTH];
                                    num++;
                                }
                                else if (DIM == 4)
                                {
                                    a1 = a1 + SLC[k].DISP_NS[ii + jj * PAR->WIDTH];
                                    a2 = a2 + SLC[k].DISP_EW[ii + jj * PAR->WIDTH];
                                    a3 = a3 + SLC[k].DISP_UD[ii + jj * PAR->WIDTH];
                                    a4 = a4 + SLC[k].DISP_UD1[ii + jj * PAR->WIDTH];
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
                    else if (DIM == 3)
                    {
                        if (num > 0)
                        {
                            a1 = a1 / num;
                            a2 = a2 / num;
                            a3 = a3 / num;
                        }
                        else
                            a1 = a2 = a3 = 0;
                    }
                    else if (DIM == 4)
                    {
                        if (num > 0)
                        {
                            a1 = a1 / num;
                            a2 = a2 / num;
                            a3 = a3 / num;
                            a4 = a4 / num;
                        }
                        else
                            a1 = a2 = a3 = a4 = 0;
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
                                else if (DIM == 3)
                                {
                                    sd1 = sd1 + (SLC[k].DISP_NS[ii + jj * PAR->WIDTH] - a1)*(SLC[k].DISP_NS[ii + jj * PAR->WIDTH] - a1);
                                    sd2 = sd2 + (SLC[k].DISP_EW[ii + jj * PAR->WIDTH] - a2)*(SLC[k].DISP_EW[ii + jj * PAR->WIDTH] - a2);
                                    sd3 = sd3 + (SLC[k].DISP_UD[ii + jj * PAR->WIDTH] - a3)*(SLC[k].DISP_UD[ii + jj * PAR->WIDTH] - a3);
                                    num++;
                                }
                                else if (DIM == 4)
                                {
                                    sd1 = sd1 + (SLC[k].DISP_NS[ii + jj * PAR->WIDTH] - a1)*(SLC[k].DISP_NS[ii + jj * PAR->WIDTH] - a1);
                                    sd2 = sd2 + (SLC[k].DISP_EW[ii + jj * PAR->WIDTH] - a2)*(SLC[k].DISP_EW[ii + jj * PAR->WIDTH] - a2);
                                    sd3 = sd3 + (SLC[k].DISP_UD[ii + jj * PAR->WIDTH] - a3)*(SLC[k].DISP_UD[ii + jj * PAR->WIDTH] - a3);
                                    sd4 = sd4 + (SLC[k].DISP_UD1[ii + jj * PAR->WIDTH] - a4)*(SLC[k].DISP_UD1[ii + jj * PAR->WIDTH] - a4);
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
                    else if (DIM == 3)
                    {
                        if (num > 1)
                        {
                            sd1 = sqrt(sd1 / (num - 1));
                            sd2 = sqrt(sd2 / (num - 1));
                            sd3 = sqrt(sd3 / (num - 1));
                        }
                        else sd1 = sd2 = sd3 = 0;

                        fprintf(pFile, "%s %f %f %f %f %f %f %f\n", SLC[k].NAME.c_str(), SLC[k].DATE, a1, sd1, a2, sd2, a3, sd3);
                    }
                    else if (DIM == 4)
                    {
                        if (num > 1)
                        {
                            sd1 = sqrt(sd1 / (num - 1));
                            sd2 = sqrt(sd2 / (num - 1));
                            sd3 = sqrt(sd3 / (num - 1));
                            sd4 = sqrt(sd4 / (num - 1));
                        }
                        else sd1 = sd2 = sd3 = sd4 = 0;

                        fprintf(pFile, "%s %f %f %f %f %f %f %f %f %f\n", SLC[k].NAME.c_str(), SLC[k].DATE, a1, sd1, a2, sd2, a3, sd3, a4, sd4);
                    }
                }
                fclose(pFile);
            }
            else
            {
                printf("warning: this point exceeds boundaries of the region, try again\n");
            }
        }
    } // read from I_FLAG_FILE file
    else if (PAR->I_FLAG == 2)
    {
        FILE * pIFile = fopen(PAR->I_FLAG_FILE.c_str(), "r");
        if (!pIFile) throw std::invalid_argument(ERROR("can not open file" + PAR->I_FLAG_FILE));

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
                    a3 = 0;
                    a4 = 0;
                    sd1 = 0;
                    sd2 = 0;
                    sd3 = 0;
                    sd4 = 0;

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
                                else if (DIM == 3)
                                {
                                    a1 = a1 + SLC[k].DISP_NS[ii + jj * PAR->WIDTH];
                                    a2 = a2 + SLC[k].DISP_EW[ii + jj * PAR->WIDTH];
                                    a3 = a3 + SLC[k].DISP_UD[ii + jj * PAR->WIDTH];
                                    num++;
                                }
                                else if (DIM == 4)
                                {
                                    a1 = a1 + SLC[k].DISP_NS[ii + jj * PAR->WIDTH];
                                    a2 = a2 + SLC[k].DISP_EW[ii + jj * PAR->WIDTH];
                                    a3 = a3 + SLC[k].DISP_UD[ii + jj * PAR->WIDTH];
                                    a4 = a4 + SLC[k].DISP_UD1[ii + jj * PAR->WIDTH];
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
                    else if (DIM == 3)
                    {
                        if (num > 0)
                        {
                            a1 = a1 / num;
                            a2 = a2 / num;
                            a3 = a3 / num;
                        }
                        else
                            a1 = a2 = a3 = 0;
                    }
                    else if (DIM == 4)
                    {
                        if (num > 0)
                        {
                            a1 = a1 / num;
                            a2 = a2 / num;
                            a3 = a3 / num;
                            a4 = a4 / num;
                        }
                        else
                            a1 = a2 = a3 = a4 = 0;
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
                                else if (DIM == 3)
                                {
                                    sd1 = sd1 + (SLC[k].DISP_NS[ii + jj * PAR->WIDTH] - a1)*(SLC[k].DISP_NS[ii + jj * PAR->WIDTH] - a1);
                                    sd2 = sd2 + (SLC[k].DISP_EW[ii + jj * PAR->WIDTH] - a2)*(SLC[k].DISP_EW[ii + jj * PAR->WIDTH] - a2);
                                    sd3 = sd3 + (SLC[k].DISP_UD[ii + jj * PAR->WIDTH] - a3)*(SLC[k].DISP_UD[ii + jj * PAR->WIDTH] - a3);
                                    num++;
                                }
                                else if (DIM == 4)
                                {
                                    sd1 = sd1 + (SLC[k].DISP_NS[ii + jj * PAR->WIDTH] - a1)*(SLC[k].DISP_NS[ii + jj * PAR->WIDTH] - a1);
                                    sd2 = sd2 + (SLC[k].DISP_EW[ii + jj * PAR->WIDTH] - a2)*(SLC[k].DISP_EW[ii + jj * PAR->WIDTH] - a2);
                                    sd3 = sd3 + (SLC[k].DISP_UD[ii + jj * PAR->WIDTH] - a3)*(SLC[k].DISP_UD[ii + jj * PAR->WIDTH] - a3);
                                    sd4 = sd4 + (SLC[k].DISP_UD1[ii + jj * PAR->WIDTH] - a4)*(SLC[k].DISP_UD1[ii + jj * PAR->WIDTH] - a4);
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
                    else if (DIM == 3)
                    {
                        if (num > 1)
                        {
                            sd1 = sqrt(sd1 / (num - 1));
                            sd2 = sqrt(sd2 / (num - 1));
                            sd3 = sqrt(sd3 / (num - 1));
                        }
                        else sd1 = sd2 = sd3 = 0;

                        fprintf(pFile, "%s %f %f %f %f %f %f %f\n", SLC[k].NAME.c_str(), SLC[k].DATE, a1, sd1, a2, sd2, a3, sd3);
                    }
                    else if (DIM == 4)
                    {
                        if (num > 1)
                        {
                            sd1 = sqrt(sd1 / (num - 1));
                            sd2 = sqrt(sd2 / (num - 1));
                            sd3 = sqrt(sd3 / (num - 1));
                            sd4 = sqrt(sd4 / (num - 1));
                        }
                        else sd1 = sd2 = sd3 = sd4 = 0;

                        fprintf(pFile, "%s %f %f %f %f %f %f %f %f %f\n", SLC[k].NAME.c_str(), SLC[k].DATE, a1, sd1, a2, sd2, a3, sd3, a4, sd4);
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
    } // print all to a single file
    else if (PAR->I_FLAG == 3)
    {
        FILE * pFile = fopen("MSBAS_TIME_SERIES.txt", "w");
        if (!pFile) throw std::invalid_argument(ERROR("can not open MSBAS_TIME_SERIES.txt file."));

        if (DIM == 1)
        {
            fprintf(pFile, "%11s %11s %11s %11s ", "X", "Y", "RATE_LOS", "RATE_STD_LOS");
            for (int k = 0; k < SLC.size(); k++) fprintf(pFile, "%11.6f ", SLC[k].DATE);
            fprintf(pFile, "\n");
        }
        else if (DIM == 2)
        {
            fprintf(pFile, "%11s %11s %11s %11s %11s %11s ", "X", "Y", "RATE_EW", "RATE_STD_EW", "RATE_UD", "RATE_STD_UD");
            for (int k = 0; k < SLC.size(); k++) fprintf(pFile, "%11.6f %11.6f ", SLC[k].DATE, SLC[k].DATE);
            fprintf(pFile, "\n");
        }
        else if (DIM == 3)
        {
            fprintf(pFile, "%11s %11s %11s %11s %11s %11s %11s %11s ", "X", "Y", "RATE_NS", "RATE_STD_NS", "RATE_EW", "RATE_STD_EW", "RATE_UD", "RATE_STD_UD");
            for (int k = 0; k < SLC.size(); k++) fprintf(pFile, "%11.6f %11.6f %11.6f ", SLC[k].DATE, SLC[k].DATE, SLC[k].DATE);
            fprintf(pFile, "\n");
        }
        else if (DIM == 4)
        {
            fprintf(pFile, "%11s %11s %11s %11s %11s %11s %11s %11s %11s %11s ", "X", "Y", "RATE_NS", "RATE_STD_NS", "RATE_EW", "RATE_STD_EW", "RATE_UD", "RATE_STD_UD", "RATE_UD1", "RATE_STD_UD1");
            for (int k = 0; k < SLC.size(); k++) fprintf(pFile, "%11.6f %11.6f %11.6f %11.6f ", SLC[k].DATE, SLC[k].DATE, SLC[k].DATE, SLC[k].DATE);
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
                            fprintf(pFile, "%11i %11i %11.6f %11.6f ", i, j, RATE_LOS[i + j * PAR->WIDTH], RATE_STD_LOS[i + j * PAR->WIDTH]);
                            for (int k = 0; k < SLC.size(); k++) fprintf(pFile, "%11.6f ", SLC[k].DISP_LOS[i + j * PAR->WIDTH]);
                            fprintf(pFile, "\n");
                        }
                        else if (DIM == 2)
                        {
                            fprintf(pFile, "%11i %11i %11.6f %11.6f %11.6f %11.6f ", i, j, RATE_EW[i + j * PAR->WIDTH], RATE_STD_EW[i + j * PAR->WIDTH], RATE_UD[i + j * PAR->WIDTH], RATE_STD_UD[i + j * PAR->WIDTH]);

                            for (int k = 0; k < SLC.size(); k++) fprintf(pFile, "%11.6f %11.6f ", SLC[k].DISP_EW[i + j * PAR->WIDTH], SLC[k].DISP_UD[i + j * PAR->WIDTH]);
                            fprintf(pFile, "\n");
                        }
                        else if (DIM == 3)
                        {
                            fprintf(pFile, "%11i %11i %11f %11f %11.6f %11.6f %11.6f %11.6f ", i, j, RATE_NS[i + j * PAR->WIDTH], RATE_STD_NS[i + j * PAR->WIDTH], RATE_EW[i + j * PAR->WIDTH], RATE_STD_EW[i + j * PAR->WIDTH], RATE_UD[i + j * PAR->WIDTH], RATE_STD_UD[i + j * PAR->WIDTH]);

                            for (int k = 0; k < SLC.size(); k++) fprintf(pFile, "%11.6f %11.6f %11.6f ", SLC[k].DISP_NS[i + j * PAR->WIDTH], SLC[k].DISP_EW[i + j * PAR->WIDTH], SLC[k].DISP_UD[i + j * PAR->WIDTH]);
                            fprintf(pFile, "\n");
                        }
                        else if (DIM == 4)
                        {
                            fprintf(pFile, "%11i %11i %11f %11f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f ", i, j, RATE_NS[i + j * PAR->WIDTH], RATE_STD_NS[i + j * PAR->WIDTH], RATE_EW[i + j * PAR->WIDTH], RATE_STD_EW[i + j * PAR->WIDTH], RATE_UD[i + j * PAR->WIDTH], RATE_STD_UD[i + j * PAR->WIDTH], RATE_UD1[i + j * PAR->WIDTH], RATE_STD_UD1[i + j * PAR->WIDTH]);

                            for (int k = 0; k < SLC.size(); k++) fprintf(pFile, "%11.6f %11.6f %11.6f %11.6f ", SLC[k].DISP_NS[i + j * PAR->WIDTH], SLC[k].DISP_EW[i + j * PAR->WIDTH], SLC[k].DISP_UD[i + j * PAR->WIDTH], SLC[k].DISP_UD1[i + j * PAR->WIDTH]);
                            fprintf(pFile, "\n");
                        }
                    }
                }
        fclose(pFile);
    }
};


void CArea::sgelss_query_lwork(int m, int n, int nrhs)
{
    int info=0, rank=0, ldb=0;
    float rcond = 0.00001; // condition number ... choose a small value
    int lda = m; // leading dimension of A
    if (m > n) // leading dimension of B
        ldb = m;
    else
        ldb = n;
    
    float *work;
    float *s;
    
    LWORK = 10;
    float_mem_alloc(work, LWORK);
    float_mem_alloc(s, MINIMUM(m, n));
        
    LWORK = -1;
    sgelss_(&m, &n, &nrhs, 0, &lda, 0, &ldb, s, &rcond, &rank, work, &LWORK, &info);
    if (info != 0 || work[0] <= 0) throw std::invalid_argument(ERROR("cannot compute LWORK"));
    
    LWORK = work[0];
    free(work);
    free(s);

    WriteLog("optimal LWORK size: " + i2s(LWORK));
}        

inline int CArea::sgelss(float *a, int m, int n, float *b, int nrhs, int &rank_num, float &cond_num)
{
    rank_num=0; cond_num=0;
    int info=0, ldb=0, rank=0;
    float rcond = 0.00001; // condition number ... choose a small value
    int lda = m; // leading dimension of A
    if (m > n) // leading dimension of B
        ldb = m;
    else
        ldb = n;
    
    float *work;
    float *s;
    
    // use default LWORK when it is not defined above
    if (LWORK == 0) LWORK = 3*MINIMUM(m,n) + MAXIMUM(MAXIMUM(2*MINIMUM(m,n), MAXIMUM(m,n)), nrhs);
        
    float_mem_alloc(work, LWORK);
    float_mem_alloc(s, MINIMUM(m, n));

    sgelss_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &LWORK, &info);

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
    else if (info == 0)
    {
        int sdim = MINIMUM(m,n);    
        cond_num=s[0]/s[sdim-1];
        rank_num=rank;
        //printf("Rank of matrix A: %d\n", rank);
        //int sdim = MINIMUM(m,n);
        //printf("Singular values of A: ");
        //for(int i = 0; i < sdim; i++) printf("%f\t", s[i]);
        //printf("\n");
        //printf("Conditioning number of A: %f\n ", s[0]/s[sdim-1]);
        //printf("SVD solution is: ");
        //for(int i = 0; i < sdim; i++) printf("%f\t", b[i]);
        //printf("\n");

        //int sdim = MINIMUM(m,n);
        //RANK = rank;
        //COND_NUM = s[0]/s[sdim-1];
    }

    free(work);
    free(s);

    return info;
};

void CArea::linreg(float *x, float *y, int length, float &aa, float &bb, float &sigmaa, float &sigmab, float &r2)
{
    if (length < 2) throw std::invalid_argument(ERROR("at least two points required for computing linear trend"));

    aa = bb = sigmaa = sigmab = r2 = 0;

    float sx = 0, sy = 0, mx = 0, my = 0, ssxx = 0, ssyy = 0, ssxy = 0, s = 0;

    for (int i = 0; i < length; i++)
    {
        sx += x[i];
        sy += y[i];
    }

    mx = sx/length;
    my = sy/length;

    for (int i = 0; i < length; i++)
    {
        ssxx += (x[i]-mx)*(x[i]-mx);
        ssyy += (y[i]-my)*(y[i]-my);
        ssxy += (x[i]-mx)*(y[i]-my);
    }

    bb = ssxy/ssxx;
    aa = my - bb * mx;
    
    s = sqrt((ssyy-ssxy*ssxy/ssxx)/(length-2));
    sigmaa = s*sqrt(1/length + mx*mx/ssxx);
    sigmab = s/sqrt(ssxx);

    r2=ssxy*ssxy/ssxx/ssyy;
};

// from Numerical Recipes in C, Second Edition (1992) - do not use because it does not compute r2
// void CArea::linreg(float *x, float *y, int length, float &aa, float &bb, float &sigmaa, float &sigmab, float &chi2)
// {
// if (length < 2) throw std::invalid_argument(ERROR("at least two points required for computing linear trend"));

// aa = bb = sigmaa = sigmab = chi2 = 0;

// float ss = 0, sx = 0, sy = 0, st2 = 0, t = 0, sxoss = 0, sigdat = 0;

// for (int i = 0; i < length; i++)
// {
// sx += x[i];
// sy += y[i];
// }

// ss = length;
// sxoss = sx / ss;

// for (int i = 0; i < length; i++)
// {
// t = x[i] - sxoss;
// st2 += t*t;
// bb += t * y[i];
// }

// bb /= st2;
// aa = (sy - sx * bb) / ss;
// sigmaa = sqrt((1.0 + sx * sx / (ss * st2)) / ss);
// sigmab = sqrt(1.0 / st2);

// for (int i = 0; i < length; i++) chi2 += (y[i] - aa - bb * x[i])*(y[i] - aa - bb * x[i]);

// if (length > 2)
// {
// sigdat = sqrt(chi2 / (length - 2));
// sigmaa *= sigdat;
// sigmab *= sigdat;
// }
// };


