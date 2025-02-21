#include "main.h"
#include "CParam.h"
#include "CArea.h"


int main(int argc, char** argv)
{
    try 
    {
        remove(LOG_FILE);
        WriteLog("hello, welcome to msbas (v20201008)");

        if (argc < 2)
        {
            CParam *par = new CParam ();
            throw std::invalid_argument(ERROR("missing parameter file, creating default parameter file " + par->NAME));
        }

        // initialize OpenMP
        WriteLog("number of threads available for parallel processing with OpenMP: " + i2s(omp_get_num_procs()));
        
        CParam *par = new CParam (argv[1]);

        CArea *area = new CArea (par);
        area->ComputeDIM();
        area->MakeTimeMatrix();
        area->ReadInterferograms();
        area->ComputeZScoreMask();
        area->CalibrateInterferograms();
        
        area->Inversion();
        area->PostProcessing();
        
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
    FILE * pFile;
    pFile = fopen(LOG_FILE, "a");
    if (pFile != NULL)
    {
        fprintf(pFile, "%s\n", message.c_str());
        fclose(pFile);
    }

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