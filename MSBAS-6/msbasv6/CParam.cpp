#include "CParam.h"

CParam::CParam()
{
    NAME = "msbas_par.txt";
    
    FILE * pFile = fopen(NAME.c_str(), "w");
    if (!pFile) throw std::invalid_argument(ERROR("can not create " + NAME + " file"));
    
    fprintf(pFile, "# %s parameter file (commented values are optional)\n\n", NAME.c_str());
    fprintf(pFile, "FORMAT=0[,0]\n");
    fprintf(pFile, "FILE_SIZE=0,0\n");
    fprintf(pFile, "# WINDOW_SIZE=0,0,0,0\n");
    fprintf(pFile, "# C_FLAG=0[...]\n");
    fprintf(pFile, "# R_FLAG=0[...]\n");
    fprintf(pFile, "# I_FLAG=0[...]\n");
    fprintf(pFile, "# V_FLAG=0\n");
    fprintf(pFile, "# D_FLAG=0\n");
    fprintf(pFile, "# DD_NSEW_FILES=0,0\n");
    fprintf(pFile, "# DEM_FILE=0\n");
    fprintf(pFile, "# SET=0,0,0,0\n");
        
    fclose(pFile);
};

CParam::CParam(string name)
{
    NAME = name;
    FORMAT = 0;
    FORMAT_INT_RADIUS = 0;
    FWIDTH = FLENGTH = 0;
    CSTART = CSTOP = LSTART = LSTOP = WIDTH = LENGTH = 0;
    C_FLAG = CSIZE = LSIZE = 0;
    for (int i=0; i<10; i++) CPOS[i] = LPOS[i] = 0;
    
    R_FLAG = 0;
    R_LAMBDA = 0;
    I_FLAG = 0;
    I_FLAG_FILE.clear();
    V_FLAG = 0;
    D_FLAG = 0;
    
    poSrcDS = NULL;
    
    DEM = DDNS = DDEW = NULL;
    DEM_FILE.clear();
    DDNS_FILE.clear();
    DDEW_FILE.clear();
    
    ReadParFile();
    
};

CParam::~CParam()
{
    if (poSrcDS != NULL) GDALClose((GDALDatasetH) poSrcDS);
    if (DEM)  free (DEM);
    if (DDNS) free (DDNS);
    if (DDEW) free (DDEW);
    
};

void CParam::ReadParFile()
{
    string line, parameter, value;

    size_t n = 0;

    ifstream infile;
    infile.open(NAME.c_str());
    if (!infile.is_open()) throw std::invalid_argument(ERROR("cannot open parameter file " + NAME));

    WriteLog ("\nprovided parameters:");
    
    while (!infile.eof())
    {
        parameter.clear(); value.clear();
        getline(infile, line);

        // remove spaces and tabs
        for (unsigned int i = 0; i < line.length(); i++)
            if (line[i] == ' ' || line[i] == '\t')
            {
                line.erase(i, 1);
                i--;
            }
        
        // remove comments
        n = line.find_first_of('#');
        if (n != string::npos) line.erase(n, line.length()-n);

        // detect parameter
        n = line.find_first_of('=');
        if (n != string::npos)
        {
            parameter = line.substr(0, n);
            value = line.substr(n+1, line.length()-n);
        }
        if (!parameter.empty() && !value.empty())
        {
            WriteLog (parameter+"="+value);
            ProcessPar(parameter, value);
        }
    }
    ValidatePar();
};

void CParam::ProcessPar(string parameter, string value)
{
    string val;
    size_t n = 0;

    if (parameter == "FORMAT")
    {
        n = value.find_first_of(',');
        if (n == string::npos)
        {
            FORMAT = atoi(value.c_str());
            if (FORMAT != 0 && FORMAT != 1 && FORMAT != 2) throw std::invalid_argument(ERROR("incorrect value in FORMAT"));
            FORMAT_INT_RADIUS = 0;
        }
        else
        {
            val = value.substr(0, n);
            FORMAT = atoi(val.c_str());
            if (FORMAT != 2) throw std::invalid_argument(ERROR("incorrect value in FORMAT"));
            value.erase(0, n + 1);
            val = value.substr(0, value.length());
            FORMAT_INT_RADIUS = atoi(val.c_str());
            if (FORMAT_INT_RADIUS <0 || FORMAT_INT_RADIUS > 16) throw std::invalid_argument(ERROR("FORMAT_INT_RADIUS must be >0 && <=16"));
        }
    }
    else if (parameter == "FILE_SIZE")
    {
        n = value.find_first_of(',');
        if (n == string::npos) throw std::invalid_argument(ERROR("expecting ',' in FILE_SIZE"));
        val = value.substr(0, n);
        FWIDTH = atoi(val.c_str());
        value.erase(0, n + 1);
        
        val = value.substr(0, value.length());
        FLENGTH = atoi(value.c_str());
        
        if (FWIDTH <= 0 || FLENGTH <= 0) throw std::invalid_argument(ERROR("incorrect value(s) in FWIDTH, FLENGTH"));
    }
    else if (parameter == "WINDOW_SIZE")
    {
        n = value.find_first_of(',');
        if (n == string::npos) throw std::invalid_argument(ERROR("expecting ',' in WINDOW_SIZE"));
        val = value.substr(0, n);
        CSTART = atoi(val.c_str());
        value.erase(0, n + 1);

        n = value.find_first_of(',');
        if (n == string::npos) throw std::invalid_argument(ERROR("expecting ',' in WINDOW_SIZE"));
        val = value.substr(0, n);
        CSTOP = atoi(val.c_str());
        value.erase(0, n + 1);

        n = value.find_first_of(',');
        if (n == string::npos) throw std::invalid_argument(ERROR("expecting ',' in WINDOW_SIZE"));
        val = value.substr(0, n);
        LSTART = atoi(val.c_str());
        value.erase(0, n + 1);

        val = value.substr(0, value.length());
        LSTOP = atoi(value.c_str());
        
        if (CSTART < 0 || CSTOP < CSTART || LSTART < 0 || LSTOP < LSTART)
            throw std::invalid_argument(ERROR("incorrect value(s) in CSTART, CSTOP, LSTART, LSTOP"));
    }
    else if (parameter == "C_FLAG")
    {
        n = value.find_first_of(',');
        if (n == string::npos)
        {
            C_FLAG = atoi(value.c_str());
        }
        else
        {
            val = value.substr(0, n);
            C_FLAG = atoi(val.c_str());
            value.erase(0, n + 1);
            if ((C_FLAG > 0)&&(C_FLAG < 10))
            {
                for (int i = 0; i < C_FLAG; i++)
                {
                    n = value.find_first_of(',');
                    if (n == string::npos) throw std::invalid_argument(ERROR("expecting ',' in C_FLAG"));
                    val = value.substr(0, n);
                    CPOS[i] = atoi(val.c_str());
                    value.erase(0, n + 1);
                    
                    n = value.find_first_of(',');
                    if (n == string::npos) throw std::invalid_argument(ERROR("expecting ',' in C_FLAG"));
                    val = value.substr(0, n);
                    LPOS[i] = atoi(val.c_str());
                    value.erase(0, n + 1);
                }
                n = value.find_first_of(',');
                if (n == string::npos) throw std::invalid_argument(ERROR("expecting ',' in C_FLAG"));
                val = value.substr(0, n);
                CSIZE = atoi(val.c_str());
                value.erase(0, n + 1);
                
                val = value.substr(0, value.length());
                if (n == string::npos) throw std::invalid_argument(ERROR("expecting ',' in C_FLAG"));
                LSIZE = atoi(value.c_str());
            }
            else if (C_FLAG == 100)
            {
                n = value.find_first_of(',');
                if (n == string::npos) throw std::invalid_argument(ERROR("expecting ',' in C_FLAG"));
                val = value.substr(0, n);
                CSIZE = atoi(val.c_str());
                value.erase(0, n + 1);
                
                val = value.substr(0, value.length());
                if (n == string::npos) throw std::invalid_argument(ERROR("expecting ',' in C_FLAG"));
                LSIZE = atoi(value.c_str());
            }
        }
    }
    else if (parameter == "R_FLAG")
    {
        n = value.find_first_of(',');
        if (n == string::npos)
        {
            val = value.substr(0, n);
            R_FLAG = atoi(val.c_str());
            if (R_FLAG != 0) throw std::invalid_argument(ERROR("incorrect value in R_FLAG"));
        }
        else
        {
            val = value.substr(0, n);
            R_FLAG = atoi(val.c_str());
            if (R_FLAG < 1 || R_FLAG > 3) throw std::invalid_argument(ERROR("incorrect value in R_FLAG"));
            value.erase(0, n + 1);
            
            val = value.substr(0, value.length());
            if (n == string::npos) throw std::invalid_argument(ERROR("expecting ',' in R_FLAG"));
            R_LAMBDA = atof(val.c_str());
            if (R_LAMBDA <= 0) throw std::invalid_argument(ERROR("incorrect value in R_LAMBDA"));
        }
    }
    else if (parameter == "I_FLAG")
    {
        n = value.find_first_of(',');
        if (n == string::npos)
        {
            I_FLAG = atoi(value.c_str());
            if (I_FLAG != 0 && I_FLAG != 1 && I_FLAG != 3) throw std::invalid_argument(ERROR("incorrect value in I_FLAG"));
        }
        else
        {
            val = value.substr(0, n);
            I_FLAG = atoi(val.c_str());
            value.erase(0, n + 1);
            I_FLAG_FILE= value.substr(0, value.length());
            if (I_FLAG != 2 || I_FLAG_FILE.empty()) throw std::invalid_argument(ERROR("incorrect value in I_FLAG_FILE"));
        }
    }
    else if (parameter == "V_FLAG")
    {
        n = value.find_first_of(',');
        if (n == string::npos)
        {
            V_FLAG = atoi(value.c_str());
            if (V_FLAG < 0 || V_FLAG > 1) throw std::invalid_argument(ERROR("incorrect value in V_FLAG"));
        }
    }
    else if (parameter == "D_FLAG")
    {
        n = value.find_first_of(',');
        if (n == string::npos)
        {
            D_FLAG = atoi(value.c_str());
            if (D_FLAG < 0 || D_FLAG > 1) throw std::invalid_argument(ERROR("incorrect value in D_FLAG"));
        }
    }
    else if (parameter == "DEM_FILE")
    {
        n = value.find_first_of(',');
        if (n == string::npos) DEM_FILE = value;
    }
    else if (parameter == "DD_NSEW_FILES")
    {
        
        n = value.find_first_of(',');
        if (n == string::npos) throw std::invalid_argument(ERROR("expecting ',' in DD_NSEW_FILES"));
        DDNS_FILE = value.substr(0, n);
        value.erase(0, n + 1);
        DDEW_FILE = value.substr(0, value.length());
    }
    else if (parameter == "SET")
    {   
        string val0, val1, val2, val3, val4, val5, val6;
        n = value.find_first_of(',');
        if (n == string::npos) throw std::invalid_argument(ERROR("expecting ',' in SET"));
        val0 = value.substr(0, n);
        value.erase(0, n + 1);

        n = value.find_first_of(',');
        if (n == string::npos) throw std::invalid_argument(ERROR("expecting ',' in SET"));
        val1 = value.substr(0, n);
        value.erase(0, n + 1);

        n = value.find_first_of(',');
        if (n == string::npos) throw std::invalid_argument(ERROR("expecting ',' in SET"));
        val2 = value.substr(0, n);
        value.erase(0, n + 1);
        
        n = value.find_first_of(',');
        if (n == string::npos) throw std::invalid_argument(ERROR("expecting ',' in SET"));
        val3 = value.substr(0, n);
        value.erase(0, n + 1);
        
        n = value.find_first_of(',');
        // four values (constant geometry)
        if (n == string::npos)
            val4 = value.substr(0, value.length());
        // six values (varying geometry)
        else
        {
            val4 = value.substr(0, n);
            value.erase(0, n + 1);
            
            n = value.find_first_of(',');
            if (n == string::npos) throw std::invalid_argument(ERROR("expecting ',' in SET"));
            val5 = value.substr(0, n);
            value.erase(0, n + 1);
            
            val6 = value.substr(0, value.length());
            
        }
        if (val0.empty() || val1.empty() || val2.empty() || val3.empty() || val4.empty()) throw std::invalid_argument(ERROR("incorrect value(s) in SET"));
        DSET.push_back(SSet(val0, val1, val2, val3, val4, val5, val6));
    }
   
};

void CParam::ValidatePar()
{
    // initialize GDAL library
    if(FORMAT==2)
        GDALAllRegister();
    
    if (CSTOP > (FWIDTH-1) || LSTOP > (FLENGTH-1))
        throw std::invalid_argument(ERROR("incorrect value(s) in CSTOP, LSTOP"));
   
    // if WINDOW_SIZE is not defined then process entire image
    else if (CSTOP == 0 || LSTOP == 0)
    {
        CSTART = 0;
        CSTOP = FWIDTH-1;
        LSTART = 0;
        LSTOP = FLENGTH-1;
    }
    
    WIDTH = CSTOP - CSTART + 1;
    LENGTH = LSTOP - LSTART + 1;
    
    WriteLog ("\ncomputed parameters:");
    WriteLog ("WIDTH="+i2s(WIDTH)+ " LENGTH="+i2s(LENGTH));

    if ((C_FLAG > 0)&&(C_FLAG < 10))
    {
        for (int i = 0; i < C_FLAG; i++)
        {
            CPOS[i] = CPOS[i] - CSTART;
            LPOS[i] = LPOS[i] - LSTART;
            if (CPOS[i] < 0  || CPOS[i] > (WIDTH - 1) || LPOS[i] < 0  || LPOS[i] > (LENGTH - 1)) throw std::invalid_argument(ERROR("incorrect value(s) in CPOS[i], LPOS[i]"));
            
            WriteLog ("CPOS["+i2s(i)+"]="+i2s(CPOS[i]) + ", LPOS["+i2s(i)+"]="+i2s(LPOS[i]));
        }
        
    }
    
    if (DEM_FILE.size () != 0)
    {
        float_mem_alloc(DEM,WIDTH*LENGTH);
        read_float_data(DEM_FILE, DEM);
    }
    if (DDNS_FILE.size () != 0)
    {
        float_mem_alloc(DDNS,WIDTH*LENGTH);
        read_float_data(DDNS_FILE, DDNS);
    }
    if (DDEW_FILE.size () != 0)
    {
        float_mem_alloc(DDEW,WIDTH*LENGTH);
        read_float_data(DDEW_FILE, DDEW);
    }
    
};

void CParam::read_float_data(string name, float *data)
{
    if (FORMAT == 0 || FORMAT == 1)
    {
        ifstream in(name.c_str(), ifstream::in | ifstream::binary);
        if (in.fail()) throw std::invalid_argument(ERROR("can not open file " + name));

        char tmpstr[4];
        int i = 0, j = 0, num = 0;

        float *fdata;
        void *v;

        int floatwidth = 4 * FWIDTH;
        char* chdata = (char*) calloc(floatwidth, sizeof (char));

        if (!chdata) throw std::invalid_argument(ERROR("during reading file " + name + "can not allocate space for chdata"));

        j = 0;
        while (!in.eof())
        {
            in.read(chdata, floatwidth);
            num = in.gcount();

            if (num == floatwidth)
            {
                for (i = 0; i < FWIDTH; i++)
                {
                    if (FORMAT == 1)
                    {
                        tmpstr[3] = chdata[4 * i + 0];
                        tmpstr[2] = chdata[4 * i + 1];
                        tmpstr[1] = chdata[4 * i + 2];
                        tmpstr[0] = chdata[4 * i + 3];
                    }
                    else if (FORMAT == 0)
                    {
                        tmpstr[0] = chdata[4 * i + 0];
                        tmpstr[1] = chdata[4 * i + 1];
                        tmpstr[2] = chdata[4 * i + 2];
                        tmpstr[3] = chdata[4 * i + 3];
                    }
                    else
                        throw std::invalid_argument(ERROR("incorrect value in FORMAT, expecting 0 or 1"));

                    v = (void*) tmpstr;
                    fdata = reinterpret_cast<float*> (v);

                    //Set NaN values to 0
                    if (*fdata != *fdata) *fdata = 0;

                    if ((i >= CSTART)&&(i <= CSTOP)&&(j >= LSTART)&&(j <= LSTOP))
                        data[(i - CSTART)+(j - LSTART)*(CSTOP - CSTART + 1)] = *fdata;
                }
                j++;
            }
            else if (num != 0)
                throw std::invalid_argument(ERROR("incorrect number of columns detected in file " + name +  ", check FILE_SIZE"));
        }

        in.close();
        free(chdata);

    }
    else if (FORMAT == 2)
    {
        // check if GTiff format is supported
        GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
        if(poDriver == NULL) throw std::invalid_argument(ERROR("GTiff format is not supported"));

        // open dataset
        GDALDataset *poDataset = (GDALDataset *) GDALOpen(name.c_str(), GA_ReadOnly);
        if (poDataset == NULL) throw std::invalid_argument(ERROR("can not open file " + name));

        GDALDataset *poDataset1 = NULL;
        GDALRasterBand *poBand = NULL;
       
        // filling gaps by interpolation is requested
        if (FORMAT_INT_RADIUS>0)
        {    
            // create temp dataset to be interpolated
            poDataset1 = poDriver->CreateCopy((name+".temp").c_str(), poDataset, FALSE, NULL, NULL, NULL);
            if (poDataset1 == NULL) throw std::invalid_argument(ERROR("can not create temp file " + name + ".temp"));
            poBand = poDataset1->GetRasterBand(1);
        
            // fill gaps by interpolation
            GDALFillNodata(poBand, NULL, FORMAT_INT_RADIUS, 0, 1, NULL, NULL, NULL);
        }
        else
            poBand = poDataset->GetRasterBand(1);
        
        // extract data from dataset
        float *pafScanline = (float *) CPLMalloc(sizeof (float)* WIDTH * LENGTH);
        if (poBand->RasterIO(GF_Read, CSTART, LSTART, WIDTH, LENGTH, pafScanline, WIDTH, LENGTH, GDT_Float32, 0, 0) != 0)
            throw std::invalid_argument(ERROR("failed reading file " + name));

        // copy data to a permanent buffer
        memcpy(data, pafScanline, sizeof(float)*WIDTH*LENGTH);
        
        // remove temporary buffer
        CPLFree(pafScanline);

        // delete temp dataset
        if (poDataset1 != NULL)
        {
            GDALClose((GDALDatasetH) poDataset1);
            unlink((name+".temp").c_str());
            //poDriver->Delete((name+".temp").c_str()); # this command crushes when processing many *tif files, it is replaced by two commands above
        }

        // close dataset except one used in write_float_data
        if (poSrcDS == NULL)
            poSrcDS=poDataset;
        else 
            GDALClose((GDALDatasetH) poDataset);
    }
};

void CParam::write_float_data(string name, float *data)
{
    if (FORMAT == 0 || FORMAT == 1)
    {
        name=name+".bin";
        
        ofstream out;
        out.open (name.c_str(), ofstream::out | ofstream::binary);
        if (out.fail()) throw std::invalid_argument(ERROR("can not open file " + name));

        char tmpstr[4];

        int i = 0, j = 0;
        float a = 0, *fdata;
        char* pstr;

        while (j < FLENGTH)
        {
            for (i = 0; i < FWIDTH; i++)
            {
                if ((i >= CSTART)&&(i <= CSTOP)&&(j >= LSTART)&&(j <= LSTOP))
                    a = data[(i - CSTART)+(j - LSTART)*(CSTOP - CSTART + 1)];
                else
                    a = 0;

                fdata = &a;
                pstr = reinterpret_cast<char*> (fdata);
                if (FORMAT == 1)
                {
                    tmpstr[3] = *(pstr + 0);
                    tmpstr[2] = *(pstr + 1);
                    tmpstr[1] = *(pstr + 2);
                    tmpstr[0] = *(pstr + 3);
                }
                else if (FORMAT == 0)
                {
                    tmpstr[0] = *(pstr + 0);
                    tmpstr[1] = *(pstr + 1);
                    tmpstr[2] = *(pstr + 2);
                    tmpstr[3] = *(pstr + 3);
                }
                else
                    throw std::invalid_argument(ERROR("incorrect value in FORMAT, expecting 0 or 1"));
                out.write(tmpstr, 4);
            }
            j++;
        }
        out.close();
    }
    else if (FORMAT == 2)
    {
        name=name+".tif";

        // check if GTiff format is supported
        GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
        if(poDriver == NULL) throw std::invalid_argument(ERROR("GTiff format is not supported"));

        // pre-delete any dataset of the name soon to be created
        poDriver->QuietDelete(name.c_str()); 	
        
        if (poSrcDS == NULL) throw std::invalid_argument(ERROR("GTiff file template is unavailable"));
        GDALDataset *poDataset = poDriver->CreateCopy(name.c_str(), poSrcDS, FALSE,NULL, NULL, NULL);

        GDALRasterBand *poBand;
        poBand = poDataset->GetRasterBand(1);

        float *full_image;
        // if only part of image is used then artificially re-create full extent image so remaining buffer is wiped out
        if (WIDTH != FWIDTH || LENGTH != FLENGTH)
        {
            float_mem_alloc(full_image, FWIDTH * FLENGTH);
            for(int i=0; i<WIDTH; i++)
                for(int j=0; j<LENGTH; j++)
                    full_image[i+CSTART+(j+LSTART)*FWIDTH]=data[i+j*WIDTH];
        }
        else
            full_image=data;
                
        if(poBand->RasterIO(GF_Write, 0, 0, FWIDTH, FLENGTH, full_image, FWIDTH, FLENGTH, GDT_Float32, 0, 0) != 0) 
            throw std::invalid_argument(ERROR("cannot write to a GTiff file " + name));

        GDALClose((GDALDatasetH) poDataset);
    }
};

