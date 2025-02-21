#include "CImage.h"

CImage::CImage(CParam *par,  const int snum, string date, string time)
{
    PAR = par;
    SNUM = snum;
    MCOUNT = 0;
    SCOUNT = 0;
    
    DISP_LOS = NULL;
    DISP_EW = NULL;
    DISP_UD = NULL;
       
    Init(date, time);
}

CImage::CImage(const CImage &image)
{
    PAR = image.PAR;
    SNUM = image.SNUM;
    MCOUNT = image.MCOUNT;
    SCOUNT = image.SCOUNT;
    
    DISP_LOS = image.DISP_LOS;
    DISP_EW = image.DISP_EW;
    DISP_UD = image.DISP_UD;
   
    NAME = image.NAME;
    DATE = image.DATE;
}

void CImage::Init(string date, string time)
{
    if (date.size() != 8 || time.size() != 6) throw std::invalid_argument(ERROR("date.size() != 8 || time.size() != 6, incorrect format in date or time"));

    int y = atoi(date.substr(0, 4).c_str());
    int m = atoi(date.substr(4, 2).c_str());
    int d = atoi(date.substr(6, 2).c_str());

    int H = atoi(time.substr(0, 2).c_str());
    int M = atoi(time.substr(2, 2).c_str());
    int S = atoi(time.substr(4, 2).c_str());

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
    if (DISP_EW)  PAR->write_float_data("MSBAS_" + NAME + "_EW", DISP_EW);
    if (DISP_UD)  PAR->write_float_data("MSBAS_" + NAME + "_UD", DISP_UD);
};