#ifndef GUARD_InitialStructures_h
#define GUARD_InitialStructures_h

#include<iostream>
#include "Parameters.h"
#include "PFUtilities.h"


struct circle
{
  double xc, yc; //Center coordinates
  double rad;
  //step sizes
};

/**

IMPORTANT FUNCTIONS

**/

circle *createCircleStruct();

void placeCircles(std::vector<circle *>,int, int,\
                  double **, std::vector<double **>, int);

void prepareFromLastRunVTKFiles(double **, std::vector<double **>,\
                                const computeParam *, \
                                const materialParam *);

void prepareUECircleTriangle(double **, std::vector<double **>,\
                             const computeParam *, \
                             const materialParam *);

void prepareTwoCircle(double **, std::vector<double **>,\
                      const computeParam *, \
                      const materialParam *);

void prepareECircleTriangle(double **, std::vector<double **>,\
                          const computeParam *, \
                          const materialParam *);


void prepareTwoSurfaceOneCircle(double **, std::vector<double **>,\
                                const computeParam *, \
                                const materialParam *);


void prepareTwoSurfaces(double **, std::vector<double **>,\
                        const computeParam *, \
                        const materialParam *);

void prepareTwoSurfaceCircleTriangle(double **, std::vector<double **>,\
                                    const computeParam *, \
                                    const materialParam *);

void prepareSurfaceCircleTriangle(double **, std::vector<double **>,\
                                 const computeParam *, \
                                 const materialParam *);


void prepareSurfaceCircle(double **, std::vector<double **>,\
                          const computeParam *, \
                          const materialParam *);



void prepareTwoSurfaceTwoCircle(double **, std::vector<double **>,\
                                const computeParam *, \
                                const materialParam *);

void prepareSurfaceTwoCircle(double **, std::vector<double **>,\
                             const computeParam *, \
                             const materialParam *);

void prepareSurfaceThreeCircle(double **, std::vector<double **>,\
                               const computeParam *, \
                               const materialParam *);

void prepareTwoSurfaceThreeCircle(double **, std::vector<double **>,\
                                  const computeParam *, \
                                  const materialParam *);


void prepareSurfaceOneBigTwoSmallCircles(double **, std::vector<double **>,\
                                  const computeParam *, \
                                  const materialParam *);

/**

TEST functions for Nx=100 and Ny=100

**/



void makeCircles(int, int, \
                 long,double *,double *, double,\
                 double **, std::vector<double **>, int);


//Funtions related to microstrucutre creation
void insertSurface(double **, std::vector<double **>,\
                   const computeParam *, \
                   const materialParam *);

//Funtions related to microstrucutre creation
void insertCircle(double **, std::vector<double **>,\
                  const computeParam *, \
                  const materialParam *);

//Funtions related to microstrucutre creation
void prepareInitialTestMicrostructures(double **, std::vector<double **>,\
                                       const computeParam *, \
                                       const materialParam *);


#endif
