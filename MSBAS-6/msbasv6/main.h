#ifndef MAIN_H
#define	MAIN_H

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
#include "gdal_alg.h"
#include <omp.h>

using namespace std;

#define ERROR(msg) ErrorMessage(__FILE__,__LINE__,__func__, msg)

#define float_mem_alloc(data,columns) { data = (float*) calloc (columns,sizeof(float)); if (!data) throw std::invalid_argument(ERROR("can not allocate memory for data")); memset(data, 0, sizeof (float)*columns); }
#define float_2d_mem_alloc(data,rows,columns) { data = (float**) calloc (rows,sizeof(float*)); for(int i = 0; i < rows; i++) { data[i] = (float*) calloc(columns, sizeof(float));  if (!data[i]) throw std::invalid_argument(ERROR("can not allocate memory for data")); memset(data[i], 0, sizeof (float)*columns); }}
 
#define MAXIMUM(m, n) (m) < (n) ? (n) : (m)
#define MINIMUM(m, n) (m) < (n) ? (m) : (n)
//#define ABS(x) (x > 0) ? (x) : (-x) //now defined in cpl_conv.h
#define LARGE_NUM_FLOAT 100000000.0
#define SMALL_NUM_FLOAT 0.000001
#define LARGE_NUM_INT 100000000
#define PI 3.1415926
#define MAX_LENGTH_FILENAME 256 // max length of file name
#define LOG_FILE "MSBAS_LOG.txt"

string ErrorMessage (const char *fname, int lineno, const char *fxname, string msg);
void WriteLog (string message);
string i2s(long num);
string i2sf(long num); //formatted
string f2s(double num);

extern "C"
{
    // minimum norm solution
    void sgelss_(int *m, int *n, int *nrhs, float *a, int *lda,
                        float *b, int *ldb, float *s, float *RCOND, int *RANK, float *work,
                        int *lwork, int *info);
    
    // LU decomoposition of a general matrix
    void sgetrf_(int* m, int *n, float* a, int* lda, int* ipiv, int* info);

    // generate inverse of a matrix given its LU decomposition
    void sgetri_(int* n, float* a, int* lda, int* ipiv, float* work, int* lwork, int* info);
}

#endif	/* MAIN_H */

