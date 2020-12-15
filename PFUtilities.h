#ifndef GUARD_PFUtilities_h
#define GUARD_PFUtilities_h

#include<iostream>
#include<fstream>
#include<fftw3.h>
#include<vector>
#include<cmath>
#include "stdlib.h"
#include "time.h"
#include "Parameters.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>

#define REAL 0
#define IMAG 1

const double PI = 3.141592653589793238463;

//Functions related to mathematical manipulation and processing in 2D

void apply2DPeriodicBC(int *, int *, int *, int *, const computeParam *);
void apply2DZeroFluxBC(int *, int *, int *, int *, const computeParam *);
double **get2DLaplacian(double **, const computeParam *);
double get2DLaplacian(double **, int, int, \
                      int, int, int, int,
                      const computeParam *);

void prepare2DFFTBasis(double *, double*, const computeParam *);
void prepare2DFFTSquaredBasis(double *, double *, const computeParam *);
double *getK2Power2DFFTBasis(double , const computeParam *);

fftw_complex *get2DFFT(double **, long, long);
fftw_complex *get2DiFFT(fftw_complex *, long, long);
fftw_complex *get2DiFFTAfterREALManipulation(double **, double **, long, long);
void set2DiFFTArray(fftw_complex *, double **, long, long);
void display2DFFTdata(fftw_complex *, long, long);

double **get2DREAL(fftw_complex *, long, long);
double **get2DCMPLX(fftw_complex *, long, long);

// Fuctions related to Array memory allocation
double **create2DField(const computeParam *);
double ***create3DField(const computeParam *);
void setZero2DField(double **, const computeParam *);

void delete2DArray(double**, long);
void displayArray(double **, const computeParam *);



// Functions related to file I/O
std::map<std::string, double> readParameterFile(const std::string& );
double **read2DVTKFile(std::string);
void write2DData(double **, const computeParam *,std::string);
void write2DVTKfile(double **, const computeParam *, long, std::string);
void write2DVTKfile(double **, const computeParam *, long);
void write1DTXTfile(double *, double *, long, std::string);
void write2DVTKfile(double **, double **, const computeParam *, long);
void write2DVTKfile(double **, double **, const computeParam *,\
                    long, std::string);

void write2DVTKfile(double **, std::string,\
                    double **, std::string, \
                    const computeParam *,\
                    long , std::string);

void write2DVTKfile(double **, std::string,\
                    double **, std::string,\
                    double **, std::string,\
                    const computeParam *,\
                    long , std::string);

void write2DVTKfile(double **, std::string,\
                    double **, std::string,\
                    double **, std::string,\
                    double **, std::string,\
                    const computeParam *,\
                    long , std::string);

//Functions related to mathematical manipulation and processing in 3D
void setZero3DField(double **, const computeParam *);

void display3DArray(double ***, const computeParam *);

void write3DVTKfile(double ***, std::string,\
                    double ***, std::string, \
                    const computeParam *,\
                    long , std::string);

#endif
