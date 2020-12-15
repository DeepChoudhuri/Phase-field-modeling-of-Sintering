#ifndef GUARD_CoupledEqFDSinter_h
#define GUARD_CoupledEqFDSinter_h

#include <vector>
#include "Parameters.h"


//Funtions related to microstrucutre creation
void prepareInitialMicrostructure(double **, std::vector<double **>,\
                                  const computeParam *,const materialParam *);

//This function does the following:
// 1. Takes all the Eta data from the Vector and compresses them
//    to a single writable matrix, and
// 2. Writes the concetration file also to a the same file as the Etas.
void writeAllDataToVTKFile(double **, std::vector<double **>,\
                           const computeParam *, const materialParam *,\
                           long);

//Derivative of the Free energy with respect to the concentration.
//Takes in one concentration value,
//and multiple non-conservative parameter fields
double dFdCon(double, std::vector<double **>,\
              const materialParam *,long, long);

//Derivative of the Free energy with respect
//to the non-conserved parameter(s) eta
////Takes in one concentration value,
//and MULTIPLE non-conservative (Eta) parameter fields
double dFdEtai(double, double, std::vector<double **>,\
               const materialParam *,long, long);

//Total  bulk free energy at a given time step
//Takes in ONE concentration field,
//and MULTIPLE non-conservative (Eta) parameter fields
double BulkFreeEnergy(double **, std::vector<double **>,\
                      const computeParam *,const materialParam *);


// Concentration dependent mobility
double mobility(double, std::vector<double **>,\
                const materialParam *,long, long);


// Solve differential equations to evolve concentration and the Eta fiels
void evolveMicrostrucutre(double **, std::vector<double **>,\
                          const computeParam *,const materialParam *);


//Get the area of the smallest cicle
double getSmallCircleArea(std::vector<double **>,
                          const computeParam *,const materialParam *);

//Get the area of the largest cicle
double getLargestCircleArea(std::vector<double **>,
                          const computeParam *,const materialParam *);

//Get neck length of the small particle/circle
double getSmallCircleNeckLength(std::vector<double **>,
                                const computeParam *,const materialParam *);


#endif
