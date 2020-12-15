#ifndef GUARD_Parameters_h
#define GUARD_Parameters_h

/* Parameters for numerical solutions*/
struct computeParam
{
  long Nx, Ny, Nz; //defines the 3D grid
  long nstep, nprint;
  double dtime;
  double dx, dy, dz; //step sizes
  int dataWRITE;
  int rerunVTKWRITE;
  long lastIterNum;
  long writeAFTERIter;

};

/*Material parameters*/
struct materialParam
{
  int initStructureType;
  // Parameters for coupled equation
  int etaNum; // Number of non-conserved eta variables. IN case of sintering number of particles
  double Kc; //Gradient coefficient for Cahn-Hilliard (CH) equation
  //double M; //Mobility coefficient for CH will be defined latter using sintering formulations
  double Ke; //Gradient coefficient for Allen-Cahn (AC) equation
  double L; // Mobility coefficient for AC equation

  //Parameters for sintering
  double dvol; //Volume diffusivity
  double dvap; //Diffusivity of the vapor phase
  double dsur; //Surface Diffusivity
  double dgb;  //Grain boundary Diffusivity

  // Parameters in the Free energy expression
  double A, B;

  double R1Factor, R2Factor, R3Factor, LEVELFactor, NumFLATSurf;
  int rerunETANum,rerunSTRUCTURE;

};

// Define parameters
void defineParameters(computeParam *, materialParam *);

#endif
