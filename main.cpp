#include<iostream>
#include "CoupledEqFDSinter.h" //CHANGE THIS header file per your NEED*/
#include "PFUtilities.h"
#include "InitialStructures.h" //CHANGE THIS header file per your NEED*/



int main(int argc, char const *argv[])
{
  /*
  Initializing microstrucural parameters
  */
  computeParam *cp = new computeParam(); //Define strucutre for computational parameters
  materialParam *mp = new materialParam(); // Define strucutre for material parameters
  defineParameters(cp,mp);

  // Define 2D concetration field
  double **Cxy = create2DField(cp);

  // Vector containing 2D slices of all eta xy fields
  std::vector<double **> Evec;
  //Create 2D slices inside the vector

  for(int i = 0; i < mp->etaNum; ++i){
      Evec.push_back(create2DField(cp));//Fill with ZERO valued matrix
  }
  //std::cout << "Not reading VTK file " << std::endl;
  std::cout << "Number of Etas: " << Evec.size() << std::endl;


    if (mp->initStructureType == 100) //from an input file
    {
      //std::cout << "Reading data from VTK file.... " << std::endl;
      //prepareFromLastRunVTKFiles(Cxy, Evec, cp, mp);//Type-10

      std::string concFname = "Conc" + std::string("_") + \
                              std::to_string(cp->lastIterNum) + ".vtk";
      Cxy = read2DVTKFile(concFname);
      std::string etaFname;
      for(int i = 0; i < mp->etaNum; i++)
      {
        etaFname = "Eta" + std::string("-") + std::to_string(i)+\
                    std::string("_") + std::to_string(cp->lastIterNum) + ".vtk";
        Evec[i] = read2DVTKFile(etaFname);
        mp->initStructureType = mp->rerunSTRUCTURE;
      }
    }
    else if(mp->initStructureType == 1){
     prepareECircleTriangle(Cxy, Evec, cp, mp); //Type-1
    }else if (mp->initStructureType == 2) {
      prepareSurfaceCircle(Cxy, Evec, cp, mp); //Type-2
    }else if (mp->initStructureType == 3) {
      prepareSurfaceCircleTriangle(Cxy, Evec, cp, mp); //Type-3
    }else if (mp->initStructureType == 4) {
      prepareTwoSurfaceCircleTriangle(Cxy, Evec, cp, mp);//Type-4
    }else if (mp->initStructureType == 5) {
      prepareSurfaceTwoCircle(Cxy, Evec, cp, mp);//Type-5
    }else if (mp->initStructureType == 6) {
      prepareTwoCircle(Cxy, Evec, cp, mp);//Type-6
    }else if (mp->initStructureType == 7) {
      prepareTwoSurfaceOneCircle(Cxy, Evec, cp, mp);//Type-7
    }else if (mp->initStructureType == 8) {
      prepareUECircleTriangle(Cxy, Evec, cp, mp);//Type-8
    }else if (mp->initStructureType == 9) {
      prepareTwoSurfaceTwoCircle(Cxy, Evec, cp, mp);//Type-9
    }else if (mp->initStructureType == 10) {
      prepareTwoSurfaces(Cxy, Evec, cp, mp);//Type-10
    }else if (mp->initStructureType == 11) {
      prepareSurfaceThreeCircle(Cxy, Evec, cp, mp); //Type-11
    }else if (mp->initStructureType == 12) {
      prepareTwoSurfaceThreeCircle(Cxy, Evec, cp, mp); //Type-11
    }else if (mp->initStructureType == 13) {
      prepareSurfaceOneBigTwoSmallCircles(Cxy, Evec, cp, mp); //Type-11
    }else {
      std::cout<< "WRONG input structure!!" << std::endl;
      exit(1);
    }

  //writeAllDataToVTKFile(Cxy, Evec, cp, mp, 0); // For testing VTK output

  /*
  END intial microstrucutre preparation
  */

  evolveMicrostrucutre(Cxy, Evec, cp, mp);

  //Garbage cleaning - delete dangling pointers
  delete cp;
  delete mp;
  //delete2DArray(Cxy,cp->Nx);
  //Evec.clear();
  return 0;
} //END function
