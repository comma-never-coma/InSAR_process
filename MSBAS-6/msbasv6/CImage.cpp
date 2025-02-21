#include "CImage.h"

CImage::CImage(CParam *par,  const int setn, string date, string time)
{
    PAR = par;
    SETN = setn;
    MCOUNT = 0;
    SCOUNT = 0;
    
    DISP_LOS = NULL;
    DISP_NS = NULL;
    DISP_EW = NULL;
    DISP_UD = NULL;
    DISP_UD1 = NULL;
       
    InitDate(date, time);
}

CImage::CImage(const CImage &image)
{
    PAR = image.PAR;
    SETN = image.SETN;
    MCOUNT = image.MCOUNT;
    SCOUNT = image.SCOUNT;
    
    DISP_LOS = image.DISP_LOS;
    DISP_NS = image.DISP_NS;
    DISP_EW = image.DISP_EW;
    DISP_UD = image.DISP_UD;
    DISP_UD1 = image.DISP_UD1;
   
    NAME = image.NAME;
    DATE = image.DATE;
}

CImage::~CImage()
{ 
    if (DISP_LOS) free(DISP_LOS);
    if (DISP_NS) free(DISP_NS);
    if (DISP_EW) free(DISP_EW);
    if (DISP_UD) free(DISP_UD);
    if (DISP_UD1) free(DISP_UD1);
} 

void CImage::InitDate(string date, string time)
{
    if (date.size() != 8 || time.size() != 6) throw std::invalid_argument(ERROR("date.size() != 8 || time.size() != 6, incorrect format in date or time"));

    int y = atoi(date.substr(0, 4).c_str());
    int m = atoi(date.substr(4, 2).c_str());
    int d = atoi(date.substr(6, 2).c_str());

    int H = atoi(time.substr(0, 2).c_str());
    int M = atoi(time.substr(2, 2).c_str());
    int S = atoi(time.substr(4, 2).c_str());

    if (y<1990 || y>2050 ) throw std::invalid_argument(ERROR("incorrect format - year"));
    if (m<1 || m>12 ) throw std::invalid_argument(ERROR("incorrect format - month"));
    if (d<1 || d>31 ) throw std::invalid_argument(ERROR("incorrect format - day"));
    if (H<0 || H>23 ) throw std::invalid_argument(ERROR("incorrect format - hour"));
    if (M<0 || M>59 ) throw std::invalid_argument(ERROR("incorrect format - minute"));
    if (S<0 || S>59 ) throw std::invalid_argument(ERROR("incorrect format - second"));


    double t = (H * 3600.0 + M * 60.0 + S * 1.0) / (3600.0 * 24);

    bool isLeapYear = ((y % 4 == 0 && y % 100 != 0) || (y % 400 == 0));
    int daysInYear = isLeapYear ? 366 : 365;

    int mday = 0;

    if (m >= 1) mday = d;
    if (m >= 2) mday += 31;
    if (m >= 3 && !isLeapYear) mday += 28;
    if (m >= 3 && isLeapYear) mday += 29;
    if (m >= 4) mday += 31;
    if (m >= 5) mday += 30;
    if (m >= 6) mday += 31;
    if (m >= 7) mday += 30;
    if (m >= 8) mday += 31;
    if (m >= 9) mday += 31;
    if (m >= 10) mday += 30;
    if (m >= 11) mday += 31;
    if (m == 12) mday += 30;

    DATE = y + (mday - 1 + t) / daysInYear;
    NAME = date + "T" + time;
};

void CImage::WriteData()
{
    if (DISP_LOS) PAR->write_float_data("MSBAS_" + NAME + "_LOS", DISP_LOS);
    if (DISP_NS)  PAR->write_float_data("MSBAS_" + NAME + "_NS", DISP_NS);
    if (DISP_EW)  PAR->write_float_data("MSBAS_" + NAME + "_EW", DISP_EW);
    if (DISP_UD)  PAR->write_float_data("MSBAS_" + NAME + "_UD", DISP_UD);
    if (DISP_UD1)  PAR->write_float_data("MSBAS_" + NAME + "_UD1", DISP_UD1);
};