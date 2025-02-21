#ifndef CIMAGE_H
#define	CIMAGE_H

#include "main.h"
#include "CParam.h"

class CImage
{
public:
    string NAME;
    double DATE;
    int MCOUNT;
    int SCOUNT;

    float *DISP_LOS;
    float *DISP_EW;
    float *DISP_UD;
    
    CImage(CParam *par,  const int snum, string date, string time);
    CImage(const CImage &image);
    void Init(string date, string time);
    void WriteData();
    
    bool operator < (const CImage &image) const {return (DATE < image.DATE);}
    bool operator == (const CImage &image) const {return (DATE == image.DATE);}

    int SNUM;
   
private:
    CParam *PAR;
};

#endif	/* CIMAGE_H */

