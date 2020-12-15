#include "Parameters.h"
#include "PFUtilities.h"

// Initializing the compute parameters - MODIFY these variables!!
void defineParameters(computeParam *cp, materialParam *mp)
{
  std::map<std::string, double> list = readParameterFile("parameters.in");
  //Computational parameters
  cp->Nx    = list["NX"];
  cp->Ny    = list["NY"];
  cp->Nz    = list["NZ"];
  cp->nprint = list["NPRINT"];
  cp->nstep = list["NSTEPS"] + 1;
  cp->dx    = list["DX"]; // anything below 0.4 is junk!!
  cp->dy    = cp->dx;
  cp->dz    = 0.0;
  cp->dtime = list["DT"];
  cp->dataWRITE = list["dataWRITE"];

  cp->rerunVTKWRITE = list["rerunVTKWRITE"];
  cp->lastIterNum = list["lastIterNum"];
  cp->writeAFTERIter = list["writeAFTERIter"];


  /* Parameters for coupled equation
  */


  mp->R1Factor = list["R1Factor"];
  mp->R2Factor = list["R2Factor"];
  mp->R3Factor = list["R3Factor"];
  mp->LEVELFactor = list["LEVELFactor"];
  mp->NumFLATSurf = list["NumFLATSurf"];

  mp->rerunETANum=list["rerunETANum"];
  mp->rerunSTRUCTURE=list["rerunSTRUCTURE"];

  mp->initStructureType = list["STRUCUTRE"];

  if(mp->initStructureType == 1){
     mp->etaNum = 3;
    }else if (mp->initStructureType == 2) {
      mp->etaNum = 2;
    }else if (mp->initStructureType == 3) {
      mp->etaNum = 4;
    }else if (mp->initStructureType == 4) {
      mp->etaNum = 5;
    }else if (mp->initStructureType == 5) {
      mp->etaNum = 3;
    }else if (mp->initStructureType == 6) {
      mp->etaNum = 2;
    }else if (mp->initStructureType == 7) {
      mp->etaNum = 3;
    }else if (mp->initStructureType == 8) {
      mp->etaNum = 3;
    }else if (mp->initStructureType == 9) {
      mp->etaNum = 4;
    }else if (mp->initStructureType == 10) {
      mp->etaNum = 2;
    }else if (mp->initStructureType == 11) {
      mp->etaNum = 4;
    }else if (mp->initStructureType == 12) {
      mp->etaNum = 5;
    }else if (mp->initStructureType == 13) {
      mp->etaNum = 4;
    }else if (mp->initStructureType == 100) {
      mp->etaNum = mp->rerunETANum;
    }


  mp->Kc = list["KC"]; //Concentration gradient
  mp->Ke = list["KE"];//non-conserved orderparameter gradient
  mp->L  = list["L"];//Allen-Cahn constant

  //Parameters for sintering. dvol, dval, dsur and dgb
  //are used for defining "M" in CH equation
  mp->dvol = list["DVOL"];
  mp->dvap = list["DVAP"];
  mp->dsur = list["DSUR"];
  mp->dgb  = list["DGB"];

  // Parameters in the Free energy expression
  mp->A = list["A"];
  mp->B = list["B"]; //Works for 1,2,5


  /*
  //mp->etaNum = 5;
  mp->Kc   = 5.0; //Concentration gradient
  mp->Ke   = 2.0;//non-conserved orderparameter gradient
  mp->L    = 5.0;//Allen-Cahn constant

  //Parameters for sintering. dvol, dval, dsur and dgb
  //are used for defining "M" in CH equation
  mp->dvol  = 0.04;
  mp->dvap  = 0.002;
  mp->dsur  = 16.0;
  mp->dgb   = 1.6;

  // Parameters in the Free energy expression
  mp->A = 16;
  mp->B = 1; //Works for 1,2,5
  */
  //mp->etaNum = 5;

}//END function
