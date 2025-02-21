#ifndef CIMAGE_H
#define	CIMAGE_H

#include "main.h"
#include "CParam.h"

class CImage
{
public:
    int SETN;
    int MCOUNT;
    int SCOUNT;

    string NAME;
    double DATE;

    float *DISP_LOS;
    float *DISP_NS;
    float *DISP_EW;
    float *DISP_UD;
    float *DISP_UD1;
    
    CImage(CParam *par,  const int setn, string date, string time);
    CImage(const CImage &image);
    ~CImage();
    
    void InitDate(string date, string time);
    void WriteData();
    
    bool operator < (const CImage &image) const {return (DATE < image.DATE);}
    bool operator == (const CImage &image) const {return (DATE == image.DATE);}
   
private:
    CParam *PAR;
};

#endif	/* CIMAGE_H */

