#include "InitialStructures.h"

/**

IMPORTANT FUNCTIONS

**/

circle *createCircleStruct()
{
  circle *cr = new circle();
  cr->xc = 0.0;
  cr->yc = 0.0;
  cr->rad = 0.0;
  return cr;
}//END function


/*
This is the main function for placing circles within initial microstrucutre
*/
void placeCircles(std::vector<circle *> cVec, int i, int j, \
                  double **Cxy, std::vector<double **> Evec, int ns)
{
  long nCircles = cVec.size();
  double dist;
  for(int p = 0; p < nCircles; p++){
    dist = std::sqrt(std::pow((i-cVec[p]->xc),2.0) + std::pow((j-cVec[p]->yc),2.0));
    if(dist <= cVec[p]->rad){ //If the distance is less than the Radius
      Cxy[i][j] = 0.5*(1.0+std::tanh(cVec[p]->rad-dist)); // Assign max. values for concentration
      Evec[p+ns][i][j] = 0.5*(1.0+std::tanh(cVec[p]->rad-dist)); // Assign max. values for etas for each particle
    }
  }
}// END FUNCTION

/*
Use old VTK files to fillup Eta and Conc field
*/
void prepareFromLastRunVTKFiles(double **Cxy, std::vector<double **> Evec,\
                                const computeParam *cp, \
                                const materialParam *mp)
{
  std::string concFname = "Conc" + std::string("_") + \
                          std::to_string(cp->lastIterNum) + ".vtk";
  Cxy = read2DVTKFile(concFname);
  //std::cout << read2DVTKFile(concFname) << std::endl;

  std::string etaFname;
  for(int i = 0; i < mp->etaNum; i++)
  {
    etaFname = "Eta" + std::string("-") + std::to_string(i)+\
                std::string("_") + std::to_string(cp->lastIterNum) + ".vtk";
    Evec[i] = read2DVTKFile(etaFname);
  }
  //std::cout << "Number of Etas: " << Evec.size() << std::endl;
}//END function




/*
THREE circles of differnt of VARIABLE radius
*/
void prepareUECircleTriangle(double **Cxy, std::vector<double **> Evec,\
                             const computeParam *cp, \
                             const materialParam *mp)
{

  double r1 = mp->R1Factor*cp->Nx;
  double r2 = mp->R2Factor*cp->Nx;
  double r3 = mp->R3Factor*cp->Nx;

  double ys = mp->LEVELFactor*cp->Ny;

  long nCircles = mp->etaNum;
  std::vector<circle *> circleVec;
  for(int i = 0; i < nCircles; ++i){
    circleVec.push_back(createCircleStruct());//Fill with ZERO valued matrix
  }

  circleVec[0]->rad= r1 ;
  circleVec[0]->xc = 0.5*cp->Nx - 0.98*r1 ;
  circleVec[0]->yc = ys;

  circleVec[1]->rad= r2;
  circleVec[1]->xc = 0.5*cp->Nx + r2;
  circleVec[1]->yc = ys;

  circleVec[2]->rad= r3;

  double d = (r1*(r3+r1)- r2*(r3-r1))/(r1+r2);
  double h = std::sqrt(std::pow(r3+r1,2.0) - std::pow(d,2.0));

  double x3 = circleVec[0]->xc + \
              (d/(r1+r2))*(circleVec[1]->xc - circleVec[0]->xc) + \
              (h/(r1+r2))*(-circleVec[1]->yc + circleVec[0]->yc);

  double y3 = circleVec[0]->yc + \
              (d/(r1+r2))*(circleVec[1]->yc - circleVec[0]->yc) + \
              (h/(r1+r2))*(circleVec[1]->xc - circleVec[0]->xc);


  circleVec[2]->xc = x3;
  circleVec[2]->yc = y3;

  for(int i = 0; i < cp->Nx; ++i)
  {
    for(int j = 0; j < cp->Nx; ++j)
    {
      //makeCircles(i,j,nCircles,xc,yc,Rad,Cxy,Evec,0);
      placeCircles(circleVec,i,j,Cxy,Evec,0);
    }
  }

}//END Function



/*
TWO circles of differnt of equal radius
*/
void prepareTwoCircle(double **Cxy, std::vector<double **> Evec,\
                          const computeParam *cp, \
                          const materialParam *mp)
{
  long nCircles = mp->etaNum;

  std::vector<circle *> circleVec;
  for(int i = 0; i < nCircles; ++i){
    circleVec.push_back(createCircleStruct());//Fill with ZERO valued matrix
  }

  circleVec[0]->rad= mp->R1Factor*cp->Nx;
  circleVec[0]->xc = 0.5*cp->Nx - 0.98*circleVec[0]->rad;
  circleVec[0]->yc = 0.5*cp->Ny;


  circleVec[1]->rad= mp->R3Factor*cp->Nx;
  circleVec[1]->xc = 0.5*cp->Nx + circleVec[1]->rad;
  circleVec[1]->yc = 0.5*cp->Ny;

  //placeCircles(circleVec);
  for(int i = 0; i < cp->Nx; ++i)
  {
    for(int j = 0; j < cp->Nx; ++j)
    {
      placeCircles(circleVec,i,j,Cxy,Evec,0);
    }
  }

}//END Function


void prepareECircleTriangle(double **Cxy, std::vector<double **> Evec,\
                          const computeParam *cp, \
                          const materialParam *mp)
{
  long nCircles = mp->etaNum;
  std::vector<circle *> circleVec;
  for(int i = 0; i < nCircles; ++i){
    circleVec.push_back(createCircleStruct());//Fill with ZERO valued matrix
  }

  double Rad = mp->R1Factor*cp->Nx;
  circleVec[0]->rad= Rad;
  circleVec[0]->xc = 0.5*cp->Nx - Rad;
  circleVec[0]->yc = 0.5*cp->Ny - Rad/std::sqrt(3);

  circleVec[1]->rad= Rad;
  circleVec[1]->xc = 0.5*cp->Nx + Rad;
  circleVec[1]->yc = 0.5*cp->Ny - Rad/std::sqrt(3);

  circleVec[2]->rad= Rad;
  circleVec[2]->xc = 0.5*cp->Nx;
  circleVec[2]->yc = 0.5*cp->Ny + (2.0/std::sqrt(3))*Rad;

  for(int i = 0; i < cp->Nx; ++i)
  {
    for(int j = 0; j < cp->Nx; ++j)
    {
      //makeCircles(i,j,nCircles,xc,yc,Rad,Cxy,Evec,0);
      placeCircles(circleVec,i,j,Cxy,Evec,0);
    }
  }

}//END Function



void prepareTwoSurfaceCircleTriangle(double **Cxy, std::vector<double **> Evec,\
                                    const computeParam *cp, \
                                    const materialParam *mp)
{
  int numFlatSurfaces = 2;

  double ys = mp->LEVELFactor*cp->Ny;
  double xs = 0.5*cp->Nx;
  double slice = 0.5*cp->Ny - ys;
  double Rad = 0.5*slice;

  long nCircles = mp->etaNum-numFlatSurfaces;
  std::vector<circle *> circleVec;
  for(int i = 0; i < nCircles; ++i){
    circleVec.push_back(createCircleStruct());//Fill with ZERO valued matrix
  }

  circleVec[0]->rad= Rad;
  circleVec[0]->xc = 0.5*cp->Nx - 0.49*slice;
  circleVec[0]->yc = ys + 0.5*slice;

  circleVec[1]->rad= Rad;
  circleVec[1]->xc = 0.5*cp->Nx + 0.5*slice;
  circleVec[1]->yc = ys + 0.5*slice;

  circleVec[2]->rad= Rad;
  circleVec[2]->xc = 0.5*cp->Nx;
  circleVec[2]->yc = (0.5*std::sqrt(3.0)+0.5)*slice + ys;

  /*double *xc = new double[nCircles];
  double *yc = new double[nCircles];

  xc[0] = 0.5*cp->Nx - 0.49*slice;
  yc[0] = ys + 0.5*slice;

  xc[1] = 0.5*cp->Nx + 0.5*slice;
  yc[1] = ys + 0.5*slice;

  xc[2] = 0.5*cp->Nx;
  yc[2] = (0.5*std::sqrt(3.0)+0.5)*slice + ys;*/

  //double dist;

  for(int i = 0; i < cp->Nx; ++i)
  {
    for(int j = 0; j < cp->Nx; ++j)
    {
      if(j <= ys && i < xs)
      {
        Cxy[i][j] = 0.9999;
        Evec[0][i][j] = 0.999;
      }
      else if(j <= ys && i > xs+cp->dx){
        Cxy[i][j] = 0.9999;
        Evec[1][i][j] = 0.9999;
      }else{
        Cxy[i][j] = 0.0;
        Evec[0][i][j] = 0.0;
        Evec[1][i][j] = 0.0;
      }

      //makeCircles(i,j,nCircles,xc,yc,Rad,Cxy,Evec,2);
      placeCircles(circleVec,i,j,Cxy,Evec,2);
    }
  }
} //END Function


void prepareSurfaceCircleTriangle(double **Cxy, std::vector<double **> Evec,\
                                 const computeParam *cp, \
                                 const materialParam *mp)
{
  int numFlatSurfaces = 1;
  double ys = mp->LEVELFactor*cp->Ny;
  double slice = 0.5*cp->Ny - ys;
  double Rad = 0.5*slice;

  int nCircles = mp->etaNum - numFlatSurfaces;
  std::vector<circle *> circleVec;
  for(int i = 0; i < nCircles; ++i){
    circleVec.push_back(createCircleStruct());//Fill with ZERO valued matrix
  }

  circleVec[0]->rad= Rad;
  circleVec[0]->xc = 0.5*cp->Nx - 0.49*slice;
  circleVec[0]->yc = ys + 0.5*slice;

  circleVec[1]->rad= Rad;
  circleVec[1]->xc = 0.5*cp->Nx + 0.5*slice;
  circleVec[1]->yc = ys + 0.5*slice;

  circleVec[2]->rad= Rad;
  circleVec[2]->xc = 0.5*cp->Nx;
  circleVec[2]->yc = (0.5*std::sqrt(3.0)+0.5)*slice + ys;

  /*double *xc = new double[nCircles];
  double *yc = new double[nCircles];

  xc[0] = 0.5*cp->Nx - 0.49*slice;
  yc[0] = ys + 0.5*slice;

  xc[1] = 0.5*cp->Nx + 0.5*slice;
  yc[1] = ys + 0.5*slice;

  xc[2] = 0.5*cp->Nx;
  yc[2] = (0.5*std::sqrt(3.0)+0.5)*slice + ys;*/

  for(int i = 0; i < cp->Nx; ++i)
  {
    for(int j = 0; j < cp->Nx; ++j)
    {
      if(j <= ys)
      {
        Cxy[i][j] = 0.9999;
        Evec[0][i][j] = 0.999;
      }
      else{
        Cxy[i][j] = 0.0;
        Evec[0][i][j] = 0.0;
      }

      //makeCircles(i,j,nCircles,xc,yc,Rad,Cxy,Evec,1);
      placeCircles(circleVec,i,j,Cxy,Evec,1);
    }
  }
  //delete[] xc;
  //delete[] yc;
} // END function

/*
Type-13
*/
void prepareSurfaceOneBigTwoSmallCircles(double **Cxy, std::vector<double **> Evec,\
                                         const computeParam *cp, \
                                         const materialParam *mp)
{
  int numFlatSurfaces = 1;
  double ys = mp->LEVELFactor*cp->Ny;

  double r1 = mp->R1Factor*cp->Ny;
  double r3 = mp->R3Factor*cp->Ny;

  int nCircles = mp->etaNum - numFlatSurfaces;
  std::vector<circle *> circleVec;
  for(int i = 0; i < nCircles; ++i){
    circleVec.push_back(createCircleStruct());//Fill with ZERO valued matrix
  }

  circleVec[0]->rad= r1;
  circleVec[0]->xc = 0.5*cp->Nx;
  circleVec[0]->yc = ys + r1;

  circleVec[1]->rad= r3;
  circleVec[1]->yc = ys + r3;
  circleVec[1]->xc = circleVec[0]->xc + std::sqrt(std::pow(r1+r3,2.0) - \
                                                  std::pow(circleVec[0]->yc - \
                                                  circleVec[1]->yc,2.0));
  circleVec[2]->rad= r3;
  circleVec[2]->yc = ys + r3;
  circleVec[2]->xc = circleVec[0]->xc - std::sqrt(std::pow(r1+r3,2.0) - \
                                                  std::pow(circleVec[0]->yc - \
                                                  circleVec[1]->yc,2.0));
  /*double *xc = new double[nCircles];
  double *yc = new double[nCircles];

  xc[0] = 0.5*cp->Nx - 0.49*slice;
  yc[0] = ys + 0.5*slice;

  xc[1] = 0.5*cp->Nx + 0.5*slice;
  yc[1] = ys + 0.5*slice;*/

  for(int i = 0; i < cp->Nx; ++i)
  {
    for(int j = 0; j < cp->Nx; ++j)
    {
      if(j <= ys)
      {
        Cxy[i][j] = 0.9999;
        Evec[0][i][j] = 0.999;
      }
      else{
        Cxy[i][j] = 0.0;
        Evec[0][i][j] = 0.0;
      }

      //makeCircles(i,j,nCircles,xc,yc,Rad,Cxy,Evec,1);
      placeCircles(circleVec,i,j,Cxy,Evec,1);
    }
  }
} // END function


/*
Type-5
*/
void prepareSurfaceTwoCircle(double **Cxy, std::vector<double **> Evec,\
                             const computeParam *cp, \
                             const materialParam *mp)
{
  int numFlatSurfaces = 1;
  double ys = mp->LEVELFactor*cp->Ny;

  double r1 = mp->R1Factor*cp->Ny;
  double r3 = mp->R3Factor*cp->Ny;

  int nCircles = mp->etaNum - numFlatSurfaces;
  std::vector<circle *> circleVec;
  for(int i = 0; i < nCircles; ++i){
    circleVec.push_back(createCircleStruct());//Fill with ZERO valued matrix
  }

  circleVec[0]->rad= r1;
  circleVec[0]->xc = 0.5*cp->Nx - 0.49*r1;
  circleVec[0]->yc = ys + r1;

  circleVec[1]->rad= r3;
  circleVec[1]->yc = ys + r3;
  circleVec[1]->xc = circleVec[0]->xc + std::sqrt(std::pow(r1+r3,2.0) - \
                                                  std::pow(circleVec[0]->yc - \
                                                  circleVec[1]->yc,2.0));

  /*double *xc = new double[nCircles];
  double *yc = new double[nCircles];

  xc[0] = 0.5*cp->Nx - 0.49*slice;
  yc[0] = ys + 0.5*slice;

  xc[1] = 0.5*cp->Nx + 0.5*slice;
  yc[1] = ys + 0.5*slice;*/

  for(int i = 0; i < cp->Nx; ++i)
  {
    for(int j = 0; j < cp->Nx; ++j)
    {
      if(j <= ys)
      {
        Cxy[i][j] = 0.9999;
        Evec[0][i][j] = 0.999;
      }
      else{
        Cxy[i][j] = 0.0;
        Evec[0][i][j] = 0.0;
      }

      //makeCircles(i,j,nCircles,xc,yc,Rad,Cxy,Evec,1);
      placeCircles(circleVec,i,j,Cxy,Evec,1);
    }
  }
} // END function


/*
Type-11
*/
void prepareTwoSurfaceThreeCircle(double **Cxy, std::vector<double **> Evec,\
                                  const computeParam *cp, \
                                  const materialParam *mp)
{
  int numFlatSurfaces = 2;
  double xs = 0.5*cp->Nx;
  double ys = mp->LEVELFactor*cp->Ny;

  double r1 = mp->R1Factor*cp->Ny;
  double r2 = mp->R2Factor*cp->Ny;
  double r3 = mp->R3Factor*cp->Ny; //smallest circle

  int nCircles = mp->etaNum - numFlatSurfaces;
  std::vector<circle *> circleVec;
  for(int i = 0; i < nCircles; ++i){
    circleVec.push_back(createCircleStruct());//Fill with ZERO valued matrix
  }

  circleVec[0]->rad= r1;
  circleVec[0]->xc = 0.5*cp->Nx - 0.49*r1*2;
  circleVec[0]->yc = ys + r1;


  circleVec[1]->rad= r2;
  circleVec[1]->xc = 0.5*cp->Nx + 0.49*r2*2;
  circleVec[1]->yc = ys + r2;


  // Smallest circle
  circleVec[2]->rad= r3;
  circleVec[2]->yc = ys + r3;
  circleVec[2]->xc = circleVec[0]->xc + 0.98*std::sqrt(std::pow(r1+r3,2.0) - \
                                                  std::pow(circleVec[2]->yc - \
                                                           circleVec[0]->yc ,2.0));
  //circleVec[2]->xc = circleVec[0]->xc + std::sqrt(std::pow(r1+r3,2.0) - \
  //                                                std::pow(circleVec[0]->yc - \
  //                                                circleVec[1]->yc,2.0));

  /*double *xc = new double[nCircles];
  double *yc = new double[nCircles];

  xc[0] = 0.5*cp->Nx - 0.49*slice;
  yc[0] = ys + 0.5*slice;

  xc[1] = 0.5*cp->Nx + 0.5*slice;
  yc[1] = ys + 0.5*slice;*/

  for(int i = 0; i < cp->Nx; ++i)
  {
    for(int j = 0; j < cp->Nx; ++j)
    {
      if(j <= ys && i < xs)
      {
        Cxy[i][j] = 0.9999;
        Evec[0][i][j] = 0.999;
      }
      else if(j <= ys && i > xs/*+cp->dx*/){
        Cxy[i][j] = 0.9999;
        Evec[1][i][j] = 0.9999;
      }else{
        Cxy[i][j] = 0.0;
        Evec[0][i][j] = 0.0;
        Evec[1][i][j] = 0.0;
      }

      //makeCircles(i,j,nCircles,xc,yc,Rad,Cxy,Evec,2);
      placeCircles(circleVec,i,j,Cxy,Evec,numFlatSurfaces);
    }
  }
} // END function



/*
Type-11
*/
void prepareSurfaceThreeCircle(double **Cxy, std::vector<double **> Evec,\
                               const computeParam *cp, \
                               const materialParam *mp)
{
  int numFlatSurfaces = 1;
  double ys = mp->LEVELFactor*cp->Ny;

  double r1 = mp->R1Factor*cp->Ny;
  double r2 = mp->R2Factor*cp->Ny;
  double r3 = mp->R3Factor*cp->Ny; //smallest circle

  int nCircles = mp->etaNum - numFlatSurfaces;
  std::vector<circle *> circleVec;
  for(int i = 0; i < nCircles; ++i){
    circleVec.push_back(createCircleStruct());//Fill with ZERO valued matrix
  }

  circleVec[0]->rad= r1;
  circleVec[0]->xc = 0.5*cp->Nx - 0.49*r1*2;
  circleVec[0]->yc = ys + r1;


  circleVec[1]->rad= r2;
  circleVec[1]->xc = 0.5*cp->Nx + 0.49*r2*2;
  circleVec[1]->yc = ys + r2;


  // Smallest circle
  circleVec[2]->rad= r3;
  circleVec[2]->yc = ys + r3;
  circleVec[2]->xc = circleVec[0]->xc + 0.98*std::sqrt(std::pow(r1+r3,2.0) - \
                                                  std::pow(circleVec[2]->yc - \
                                                           circleVec[0]->yc ,2.0));
  //circleVec[2]->xc = circleVec[0]->xc + std::sqrt(std::pow(r1+r3,2.0) - \
  //                                                std::pow(circleVec[0]->yc - \
  //                                                circleVec[1]->yc,2.0));

  /*double *xc = new double[nCircles];
  double *yc = new double[nCircles];

  xc[0] = 0.5*cp->Nx - 0.49*slice;
  yc[0] = ys + 0.5*slice;

  xc[1] = 0.5*cp->Nx + 0.5*slice;
  yc[1] = ys + 0.5*slice;*/

  for(int i = 0; i < cp->Nx; ++i)
  {
    for(int j = 0; j < cp->Nx; ++j)
    {
      if(j <= ys)
      {
        Cxy[i][j] = 0.9999;
        Evec[0][i][j] = 0.999;
      }
      else{
        Cxy[i][j] = 0.0;
        Evec[0][i][j] = 0.0;
      }

      //makeCircles(i,j,nCircles,xc,yc,Rad,Cxy,Evec,1);
      placeCircles(circleVec,i,j,Cxy,Evec,1);
    }
  }
} // END function





void prepareSurfaceCircle(double **Cxy, std::vector<double **> Evec,\
                         const computeParam *cp, \
                         const materialParam *mp)
{
  //double slice = 0.333333;
  int numFlatSurfaces = 1;

  double ys = mp->LEVELFactor*cp->Ny;

  int nCircles = mp->etaNum - numFlatSurfaces;
  std::vector<circle *> circleVec;
  for(int i = 0; i < nCircles; ++i){
    circleVec.push_back(createCircleStruct());//Fill with ZERO valued matrix
  }

  circleVec[0]->rad= mp->R3Factor*cp->Ny;
  circleVec[0]->xc = 0.5*cp->Nx;
  circleVec[0]->yc = ys + circleVec[0]->rad;

  /*double *xc = new double[1];
  double *yc = new double[1];

  xc[0] = 0.5*cp->Nx;
  yc[0] = ys + 0.5*slice;*/

  for(int i = 0; i < cp->Nx; ++i)
  {
    for(int j = 0; j < cp->Nx; ++j)
    {
      if(j <= ys)
      {
        Cxy[i][j] = 0.9999;
        Evec[0][i][j] = 0.999;
      }
      else{
        Cxy[i][j] = 0.0;
        Evec[0][i][j] = 0.0;
      }
      /*dist = std::sqrt(std::pow((i-xc),2.0) + std::pow((j-yc),2.0));

      if(dist <= Rad){ //If the distance rx1 is less than the Radius
        Cxy[i][j] = 0.5*(1.0+std::tanh(Rad-dist)); // Assign max. values for concentration
        Evec[1][i][j] = 0.5*(1.0+std::tanh(Rad-dist)); // Assign max. values for etas for each particle*/
      //makeCircles(i,j,1,xc,yc,Rad,Cxy,Evec,1);
        placeCircles(circleVec,i,j,Cxy,Evec,1);
      }
    }
}//END Function


//Type-9
void prepareTwoSurfaceTwoCircle(double **Cxy, std::vector<double **> Evec,\
                                const computeParam *cp, \
                                const materialParam *mp)
{
  int numFlatSurfaces = 2;
  double factor = 0.4;
  double ys = factor*cp->Ny;
  double xs = 0.5*cp->Nx;

  double r1 = mp->R1Factor*cp->Ny;
  double r3 = mp->R3Factor*cp->Ny;

  int nCircles = mp->etaNum - numFlatSurfaces;
  std::vector<circle *> circleVec;
  for(int i = 0; i < nCircles; ++i){
    circleVec.push_back(createCircleStruct());//Fill with ZERO valued matrix
  }

  circleVec[0]->rad= r1;
  circleVec[0]->xc = 0.5*cp->Nx - 0.49*r1;
  circleVec[0]->yc = ys + r1;

  circleVec[1]->rad= r3;
  circleVec[1]->yc = ys + r3;
  circleVec[1]->xc = circleVec[0]->xc + std::sqrt(std::pow(r1+r3,2.0) - \
                                                  std::pow(circleVec[0]->yc - \
                                                  circleVec[1]->yc,2.0));


  for(int i = 0; i < cp->Nx; ++i)
  {
    for(int j = 0; j < cp->Nx; ++j)
    {
      if(j <= ys && i < xs)
      {
        Cxy[i][j] = 0.9999;
        Evec[0][i][j] = 0.999;
      }
      else if(j <= ys && i > xs/*+cp->dx*/){
        Cxy[i][j] = 0.9999;
        Evec[1][i][j] = 0.9999;
      }else{
        Cxy[i][j] = 0.0;
        Evec[0][i][j] = 0.0;
        Evec[1][i][j] = 0.0;
      }

      //makeCircles(i,j,nCircles,xc,yc,Rad,Cxy,Evec,2);
      placeCircles(circleVec,i,j,Cxy,Evec,2);
    }
  }
} //END Function


void prepareTwoSurfaces(double **Cxy, std::vector<double **> Evec,\
                        const computeParam *cp, \
                        const materialParam *mp)
{

  double ys = mp->LEVELFactor*cp->Ny;
  double xs = 0.5*cp->Nx;


  for(int i = 0; i < cp->Nx; ++i)
  {
    for(int j = 0; j < cp->Nx; ++j)
    {
      if(j <= ys && i < xs)
      {
        Cxy[i][j] = 0.9999;
        Evec[0][i][j] = 0.9999;
      }
      else if(j <= ys && i > xs/*+cp->dx*/){
        Cxy[i][j] = 0.9999;
        Evec[1][i][j] = 0.9999;
      }else{
        Cxy[i][j] = 0.0;
        Evec[0][i][j] = 0.0;
        Evec[1][i][j] = 0.0;
      }
    }
  }
} //END Function


//Type 7
void prepareTwoSurfaceOneCircle(double **Cxy, std::vector<double **> Evec,\
                                const computeParam *cp, \
                                const materialParam *mp)
{
  int numFlatSurfaces = 2;
  double ys = mp->LEVELFactor*cp->Ny;
  double xs = 0.5*cp->Nx;

  int nCircles = mp->etaNum - numFlatSurfaces;
  std::vector<circle *> circleVec;
  for(int i = 0; i < nCircles; ++i){
    circleVec.push_back(createCircleStruct());//Fill with ZERO valued matrix
  }

  circleVec[0]->rad= mp->R3Factor*cp->Ny;
  circleVec[0]->xc = 0.5*cp->Nx;
  circleVec[0]->yc = ys + circleVec[0]->rad;

  for(int i = 0; i < cp->Nx; ++i)
  {
    for(int j = 0; j < cp->Nx; ++j)
    {
      if(j <= ys && i < xs)
      {
        Cxy[i][j] = 0.9999;
        Evec[0][i][j] = 0.999;
      }
      else if(j <= ys && i > xs/*+cp->dx*/){
        Cxy[i][j] = 0.9999;
        Evec[1][i][j] = 0.9999;
      }else{
        Cxy[i][j] = 0.0;
        Evec[0][i][j] = 0.0;
        Evec[1][i][j] = 0.0;
      }

      //makeCircles(i,j,nCircles,xc,yc,Rad,Cxy,Evec,2);
      placeCircles(circleVec,i,j,Cxy,Evec,2);
    }
  }
} //END Function


/**

TEST functions for Nx=100 and Ny=100

**/


void makeCircles(int i, int j, \
                 long nCircles,double *xc,double *yc, double Rad,\
                 double **Cxy, std::vector<double **> Evec, int ns)
{
  double dist;
  for(int p = 0; p < nCircles; p++){
    dist = std::sqrt(std::pow((i-xc[p]),2.0) + std::pow((j-yc[p]),2.0));
    if(dist <= Rad){ //If the distance is less than the Radius
      Cxy[i][j] = 0.5*(1.0+std::tanh(Rad-dist)); // Assign max. values for concentration
      Evec[p+ns][i][j] = 0.5*(1.0+std::tanh(Rad-dist)); // Assign max. values for etas for each particle
    }
  }
} // END FUNCTION


void insertSurface(double **Cxy, std::vector<double **> Evec,\
                                  const computeParam *cp, \
                                  const materialParam *mp)
{
  double ys = 0.333*cp->Ny;

  for(int i = 0; i < cp->Nx; ++i)
  {
    for(int j = 0; j < cp->Nx; ++j)
    {
      if(j <= ys)
      {
        Cxy[i][j] = 0.9999;
        Evec[0][i][j] = 0.999;
      }
      else{
        Cxy[i][j] = 0.0;
        Evec[0][i][j] = 0.0;
      }
    }
  }
}


void insertCircle(double **Cxy, std::vector<double **> Evec,\
                                  const computeParam *cp, \
                                  const materialParam *mp)
{
  double xc = 0.5*cp->Nx;
  double yc = 0.5*cp->Ny;
  double dist = 0.0;
  double Rad = 20.0;

  for(int i = 0; i < cp->Nx; ++i)
  {
    for(int j = 0; j < cp->Nx; ++j)
    {
      double distX = std::pow(i - xc, 2.0);
      double distY = std::pow(j - yc, 2.0);
      dist = std::sqrt(distX+distY);
      //std::cout << dist << " ";
      //std::cout << 0.5*(1.0+std::tanh(Rad-dist)) << " ";

      Cxy[i][j] = 0.5*(1.0+std::tanh(Rad-dist));
      Evec[0][i][j] = 0.5*(1.0+std::tanh(Rad-dist));
    }
  }

} //END function

//Funtions related to microstrucutre creation - for (100x100 matrix ONLY!!)
void prepareInitialTestMicrostructures(double **Cxy, std::vector<double **> Evec,\
                                      const computeParam *cp, \
                                      const materialParam *mp)
{
  std::cout << "Preparing initial microstructutre ..." << std::endl;

  //Coordinates of particle - dynamically assigned
  double *xc;
  double *yc;

  // used when more than two etas or particles
  double RN = 10.0;
  double RNnew = RN;

  // used when more than two etas or particles
  double R1 = 20.0;
  double R2 = 0.5*R1;

  // For evaluating is the point falls within a particle
  double rxy1 = 0;
  double rxy2 = 0; // This one is used for bi-crystal configuration


  if(mp->etaNum !=2){ // Here etaNum is treated as number of particles
    xc = new double[mp->etaNum];
    yc = new double[mp->etaNum];

    //Equilaterial Triangle - 3 particles of equal size
    /*RNnew = 20.0;
    xc[0] = 30.0; yc[0] = 30.0;
    xc[1] = 70.0; yc[1] = 30.0;
    xc[2] = 50.0; yc[2] = 64.64;*/

    //Square - 9 particles
    RNnew = 10.0;
    xc[0] = 29.0; yc[0] = 50.0;
    xc[1] = 50.0; yc[1] = 50.0;
    xc[2] = 71.0; yc[2] = 50.0;
    xc[3] = 50.0; yc[3] = 29.0;
    xc[4] = 50.0; yc[4] = 71.0;
    xc[5] = 39.0; yc[5] = 39.0;
    xc[6] = 61.0; yc[6] = 39.0;
    xc[7] = 39.0; yc[7] = 61.0;
    xc[8] = 61.0; yc[8] = 61.0;

    for(int ei = 0; ei < mp->etaNum; ++ei){
      if(ei > 4){
        RNnew = 0.5*RN;
      }

      srand((unsigned)time( NULL));
      for(int i = 0; i < cp->Nx; ++i){
        for(int j = 0; j < cp->Ny; ++j){
          rxy1 = std::sqrt(std::pow((i-xc[ei]),2.0) + std::pow((j-yc[ei]),2.0));
          if(rxy1 <= RNnew){ //If the distance rx1 is less than the Radius
            //Cxy[i][j] = 0.999; // Assign max. values for concentration
            Cxy[i][j] = 0.5*(1.0+std::tanh(RNnew-rxy1));
            Evec[ei][i][j] = 0.999; // Assign max. values for etas for each particle
          }
        }
      }
    }
  } // This ends microstrucutre with more than 2 particles

  if(mp->etaNum == 2){ // When we have a Bi-crystal, i.e. only TWO particles
    xc = new double[mp->etaNum];
    yc = new double[mp->etaNum];
    //Coordinates assignment
    xc[0] = cp->Nx*0.5;
    yc[0] = 40.0;
    yc[1] = 70.0;

    for(int i = 0; i < cp->Nx; ++i){
      for (int j = 0; j < cp->Ny; ++j) {
        rxy1 = std::sqrt(std::pow((i-xc[0]),2.0) + std::pow((j-yc[0]),2.0));
        rxy2 = std::sqrt(std::pow((i-xc[0]),2.0) + std::pow((j-yc[1]),2.0));

        if(rxy1 <= R1){ //If the distance rx1 is less than the Radius
          Cxy[i][j] = 0.5*(1.0+std::tanh(R1-rxy1)); // Assign max. values for concentration
          Evec[0][i][j] = 0.5*(1.0+std::tanh(R1-rxy1)); // Assign max. values for etas for each particle
        }

        if(rxy2 <= R2){ //If the distance rx1 is less than the Radius
          Cxy[i][j] = 0.5*(1.0+std::tanh(R2-rxy2)); // Assign max. values for concentration
          Evec[1][i][j] = 0.5*(1.0+std::tanh(R2-rxy2)); // Assign max. values for etas for each particle
        }
      }
    }
    //xc = NULL; yc = NULL;
  } //End bi-crystal configuration
  std::cout << "Intial microstructure created!" << std::endl;
  //Garbage cleaning
  delete[] xc;
  delete[] yc;
}//END function
