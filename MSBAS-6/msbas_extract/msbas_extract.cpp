#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <sstream> 
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <stdexcept>
#include <algorithm> 
#include <gdal_priv.h>
#include <cpl_conv.h> // for CPLMalloc()

using namespace std;

#define ERROR(msg) ErrorMessage(__FILE__,__LINE__,__func__, msg)
#define float_mem_alloc(data,length) { data = (float*) calloc (length,sizeof(float)); if (!data) throw std::invalid_argument(ERROR("can not allocate memory for data")); memset(data, 0, sizeof (float)*length); }

string ErrorMessage (const char *fname, int lineno, const char *fxname, string msg);
void WriteLog (string message);
string i2s(long num);
string i2sf(long num); //formatted
string f2s(double num);

class CParam
{
public:
    string NAME;
    int FORMAT;
    int FWIDTH, FLENGTH;
    int X, Y, DX, DY;
    FILE * FILEOUT;

    CParam(string name, int x, int y, int dx, int dy);
    ~CParam(); 
    void ReadParFile();
    size_t split(const std::string &txt, std::vector<std::string> &strs, char ch);
    void ProcessPar(string parameter, string value);
    void read_float_data(string name, float &average, float &stdev);
    
private:

};

// input: file name [usually MSBAS_TSOUT.txt], region location [column row], half-dimensions [column_radius row_radius]
// output: file with time series for the region

int main(int argc, char** argv)
{
    try
    {
        WriteLog("hello, welcome to msbas_extract (v20200529)");

        if (argc < 6) throw std::invalid_argument(ERROR("missing (five) space delimited input parameters: file name [usually MSBAS_TSOUT.txt], region location [column row], half-dimensions [column_radius row_radius]"));
    
        char *name = argv[1];
        int x=atoi(argv[2]);
        int y=atoi(argv[3]);
        int dx=atoi(argv[4]);
        int dy=atoi(argv[5]);
        
        if (x<0 || y<0 || dx<0 || dy<0)
            throw std::invalid_argument(ERROR("parameters x, y, dx, dy cannot be negative"));
        
        CParam *par = new CParam (name, x,y, dx, dy);
        if (par == NULL) WriteLog("this should not have happened, exiting...");
        
        WriteLog("good bye");
    }
    catch(const std::invalid_argument& e)
    {
        WriteLog(e.what());
        return 1;
    }
    return 0;
}

string ErrorMessage (const char *fname, int lineno, const char *fxname, string msg)
{
    string fname1 = string(fname); 
    string lineno1 = i2s(lineno); 
    string fxname1 = string(fxname); 
    return "error in [file " + fname1 + ", line " + lineno1 + ", function " + fxname1 + "]: " + msg + ", exiting...";
}

void WriteLog (string message)
{
    printf("%s\n", message.c_str());
}

string i2s(long num)
{    
    ostringstream ostr;
    ostr << num;
    return ostr.str();
}

string i2sf(long num)
{    
    ostringstream ostr;
    ostr << setw(5) << num;
    return ostr.str();
}

string f2s(double num)
{    
    ostringstream ostr;
    ostr << showpoint << fixed << setprecision(8) << right << num;
    return ostr.str();
}

CParam::CParam(string name, int x, int y, int dx, int dy)
{
    NAME = name;
    FORMAT = 0;
    FWIDTH = FLENGTH = 0;
    X = x;
    Y = y;
    DX = dx;
    DY = dy;
    
    ReadParFile();
};

CParam::~CParam()
{
    fclose(FILEOUT);
};

void CParam::ReadParFile()
{
    string line, parameter, value;

    size_t n = 0;

    ifstream infile;
    infile.open(NAME.c_str());
    if (!infile.is_open()) throw std::invalid_argument(ERROR("cannot open parameter file " + NAME));

    WriteLog ("provided parameters:");
    
    while (getline(infile, line))
    {
        // avoid empty lines
        if (line.find_first_of('=') != string::npos || line.find_first_of(' ') != string::npos)
        {
            parameter.clear(); value.clear();

            // remove double spaces
            for (size_t i = 0; i < (line.length()-1); i++)
            {
                if (line[i] == ' ' && line[i+1] == ' ')
                {
                    line.erase(i+1, 1);
                    i--;
                }
            }

            // remove comments
            n = line.find_first_of('#');
            if (n != string::npos) line.erase(n, line.length()-n);

            // detect parameter 
            if (line.find_first_of('=') != string::npos)
            {
                n = line.find_first_of('=');
                parameter = line.substr(0, n);
                value = line.substr(n+1, line.length()-n);
                if (!parameter.empty() && !value.empty())
                {
                    WriteLog (parameter+"="+value);
                    ProcessPar(parameter, value);
                }
            }
            else if (line.find_first_of(' ') != string::npos)
            {
                // open file for writing if it is not already opened
                if (!FILEOUT)
                {
                    string name="MSBAS_EXTRACT_" + i2s(X)+"_" + i2s(Y) + "_" + i2s(DX) + "_" + i2s(DY) + ".txt";
                    FILEOUT = fopen(name.c_str(), "w");
                    WriteLog("writing results (YYMMDDTHHMMSS YYYY.YYYYYY [for each additional column] DISP DISP_ERROR) to a file " + name);
                }
                
                vector <string> v;
                int N=split(line, v, ' ');
                if (N>2)
                {
                    if (!FILEOUT) throw std::invalid_argument(ERROR("can not open FILEOUT "));
                    fprintf(FILEOUT, "%s %s", v[0].c_str(), v[1].c_str());    

                    int N=split(line, v, ' ');
                    for (int i=2; i<N; i++)
                    {
                        float average = 0, stdev=0;
                        read_float_data(v[i], average, stdev);
                        fprintf(FILEOUT, " %f %f", average, stdev);  
                    }

                    fprintf(FILEOUT, " \n");  
                }
            }
        }
    }
};

void CParam::ProcessPar(string parameter, string value)
{
    string val;
    size_t n = 0;

    if (parameter == "FORMAT")
    {
        FORMAT = atoi(value.c_str());
        if (FORMAT != 0 && FORMAT != 1 && FORMAT != 2 ) throw std::invalid_argument(ERROR("incorrect value in FORMAT"));
        
        // must be run once to initialize GDAL
        if (FORMAT == 2 ) GDALAllRegister();
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
        
        if (X>=FWIDTH || Y>=FLENGTH) throw std::invalid_argument(ERROR("incorrect value(s) in X, Y"));
        
        if ((X-DX)<0)
        {
            DX=X;
            WriteLog("adjusted DX to " + i2s(DX));
        }
        
        if ((Y-DY)<0)
        {
            DY=Y;
            WriteLog("adjusted DY to " + i2s(DY));
        }

        if ((X+DX)>=FWIDTH)
        {
            DX=FWIDTH-X-1;
            WriteLog("adjusted DX to " + i2s(DX));
        }
        
        if ((Y+DY)>=FLENGTH)
        {
            DY=FLENGTH-Y-1;
            WriteLog("adjusted DY to " + i2s(DY));
        }
    }
};

size_t CParam::split(const std::string &txt, vector<string> &strs, char ch)
{
    size_t pos = txt.find( ch );
    size_t initialPos = 0;
    strs.clear();
 
    // Decompose statement
    while(pos != string::npos)
    {
        strs.push_back( txt.substr( initialPos, pos - initialPos ));
        initialPos = pos + 1;
        pos = txt.find(ch, initialPos);
    }
 
    // Add the last one
    strs.push_back(txt.substr(initialPos, min(pos, txt.size()) - initialPos + 1));
 
    return strs.size();
};

void CParam::read_float_data(string name, float &average, float &stdev)
{
    float *data;
    
    if (FORMAT == 0 || FORMAT == 1)
    {
        float_mem_alloc(data, (2*DX+1)*(2*DY+1));
        
        ifstream in(name.c_str(), ifstream::in | ifstream::binary);
        if (in.fail()) throw std::invalid_argument(ERROR("can not open file " + name));

        char tmpstr[4];
        int i = 0, j = 0, num = 0;

        float *fdata;
        void *v;

        int floatwidth = 4 * FWIDTH;
        char* chdata = (char*) calloc(floatwidth, sizeof (char));

        if (!chdata) throw std::invalid_argument(ERROR("during reading file " + name + "can not allocate memory for chdata"));

        j = 0;
        while (!in.eof())
        {
            in.read(chdata, floatwidth);
            num = in.gcount();

            if (num == floatwidth)
            {
                if (j>=(Y-DY) && j<=(Y+DY))
                {
                    for (i = (X-DX); i <= (X+DX); i++)
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

                        int ii=i-(X-DX), jj=j-(Y-DY);

                        data[ii + jj*(2*DX+1)] = *fdata;

                    }
                }
            }
            else if (num != 0)
                throw std::invalid_argument(ERROR("incorrect number of columns detected in file " + name +  ", check FILE_SIZE"));
            
            j++;
        }
                        
        in.close();
        free(chdata);
    }
    else if (FORMAT == 2)
    {
        // open file
        GDALDataset *poDataset = (GDALDataset *) GDALOpen(name.c_str(), GA_ReadOnly);
        if (poDataset == NULL) throw std::invalid_argument(ERROR("can not open file " + name));

        GDALRasterBand *poBand;
        poBand = poDataset->GetRasterBand(1);

        /* extraction of useful parameters - presently not used
        double adfGeoTransform[6];
        printf("Driver: %s/%s\n", poDataset->GetDriver()->GetDescription(), poDataset->GetDriver()->GetMetadataItem(GDAL_DMD_LONGNAME));
        printf("Size is %dx%dx%d\n", poDataset->GetRasterXSize(), poDataset->GetRasterYSize(), poDataset->GetRasterCount());
        if (poDataset->GetProjectionRef() != NULL) printf("Projection is %s\n", poDataset->GetProjectionRef());

        if (poDataset->GetGeoTransform(adfGeoTransform) == CE_None)
        {
            printf("Origin = (%.6f,%.6f)\n", adfGeoTransform[0], adfGeoTransform[3]);
            printf("Pixel Size = (%.6f,%.6f)\n", adfGeoTransform[1], adfGeoTransform[5]);
        }

        int nBlockXSize, nBlockYSize;
        int bGotMin, bGotMax;
        double adfMinMax[2];
        poBand = poDataset->GetRasterBand(1);
        poBand->GetBlockSize(&nBlockXSize, &nBlockYSize);
        printf("Block=%dx%d Type=%s, ColorInterp=%s\n", nBlockXSize, nBlockYSize, GDALGetDataTypeName(poBand->GetRasterDataType()), GDALGetColorInterpretationName(poBand->GetColorInterpretation()));
        adfMinMax[0] = poBand->GetMinimum(&bGotMin);
        adfMinMax[1] = poBand->GetMaximum(&bGotMax);
        if (!(bGotMin && bGotMax)) GDALComputeRasterMinMax((GDALRasterBandH) poBand, TRUE, adfMinMax);
        printf("Min=%.3f, Max=%.3f\n", adfMinMax[0], adfMinMax[1]);

        if (poBand->GetOverviewCount() > 0) printf("Band has %d overviews.\n", poBand->GetOverviewCount());
        if (poBand->GetColorTable() != NULL) printf("Band has a color table with %d entries.\n", poBand->GetColorTable()->GetColorEntryCount());

        int nXSize = poBand->GetXSize();
        int nYSize = poBand->GetYSize();
         */

        data = (float *) CPLMalloc(sizeof (float)* (2*DX+1)* (2*DY+1));
        if (poBand->RasterIO(GF_Read, X-DX, Y-DY, 2*DX+1, 2*DY+1, data, 2*DX+1, 2*DY+1, GDT_Float32, 0, 0) != 0)
            throw std::invalid_argument(ERROR("failed reading file " + name));

        if( poDataset != NULL ) GDALClose( (GDALDatasetH) poDataset );
    }

    int count = 0;
    
    average=0;
    stdev=0;

    for (int i = 0; i < (2 * DX + 1); i++)
        for (int j = 0; j < (2 * DY + 1); j++)
            if (data[i + j * (2 * DX + 1)] != 0)
            {
                average += data[i + j * (2 * DX + 1)];
                count++;
            }

    if (count > 0)
        average /= count;
    else
        average = 0;

    if (count > 1)
    {
        for (int i = 0; i < (2 * DX + 1); i++)
            for (int j = 0; j < (2 * DY + 1); j++)
                if (data[i + j * (2 * DX + 1)] != 0)
                    stdev += (data[i + j * (2 * DX + 1)] - average)*(data[i + j * (2 * DX + 1)] - average);

        stdev = sqrt(stdev / (count - 1));
    }

    if (FORMAT == 0 || FORMAT == 1)
        free (data);
    else if (FORMAT == 2)
        CPLFree(data);
};
