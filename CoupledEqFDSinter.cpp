#include "CoupledEqFDSinter.h"
#include "PFUtilities.h"


//This function does the following:
// 1. Takes all the Eta data from the Vector and compresses them
//    to a single writable matrix, and
// 2. Writes the concetration file also to a the same file as the Etas.
void writeAllDataToVTKFile(double ** cxy, std::vector<double **> evec,\
                           const computeParam *cp, const materialParam *mp,\
                           long istep)
{
  double **exy2 = create2DField(cp);
  double **muxy = create2DField(cp);

  double sumei2 = 0.0;
  double sumei3 = 0.0;

  for(int i = 0; i < cp->Nx; ++i){
    for(int j = 0; j < cp->Ny; ++j){
      for(int ne = 0; ne < mp->etaNum; ++ne){
        exy2[i][j] = exy2[i][j] + std::pow(evec[ne][i][j], 2.0);

        sumei2 = sumei2 + std::pow(evec[ne][i][j],2.0);
        sumei3 = sumei3 + std::pow(evec[ne][i][j],3.0);
      }
      muxy[i][j] = 2.0*mp->A*(cxy[i][j]*(1.0-cxy[i][j])*(1.0-cxy[i][j]) - \
                   cxy[i][j]*cxy[i][j]*(1.0-cxy[i][j])) + \
                   mp->B*(2.0*cxy[i][j]-6.0*sumei2+4.0*sumei3);
      sumei2 = 0.0;
      sumei3 = 0.0;
    }
  }
  write2DVTKfile(cxy, "DEN", \
                 exy2, "ALL-ETA", \
                 muxy, "MU", \
                 cp, istep, "COMBINED");
}//END function



//Derivative of the Free energy with respect to the concentration.
////Takes in one concentration value,
//and MULTIPLE non-conservative (Eta) parameter fields
double dFdCon(double cxy, std::vector<double **> evec,\
              const materialParam *mp,long p, long q)
{
  double sumei2 = 0.0;
  double sumei3 = 0.0;
  for(int ei = 0; ei < mp->etaNum; ++ei){
    sumei2 = sumei2 + std::pow(evec[ei][p][q],2.0);
    sumei3 = sumei3 + std::pow(evec[ei][p][q],3.0);
  }
  return 2.0*mp->A*(cxy*(1.0-cxy)*(1.0-cxy) - cxy*cxy*(1.0-cxy)) + \
         mp->B*(2.0*cxy-6.0*sumei2+4.0*sumei3);
//         2.0*mp->A*(1.0 - cxy)*(1.0 - 2.0*cxy);
//         2.0*mp->A*(cxy*std::pow(1.0-cxy,2.0) - std::pow(cxy,2.0)*(1.0-cxy));

} //END function



//Derivative of the Free energy with respect
//to the non-conserved parameter(s) eta
////Takes in one concentration value,
//and MULTIPLE non-conservative (Eta) parameter fields
double dFdEtai(double cxy, double eixy, std::vector<double **> evec,\
               const materialParam *mp, long p, long q)
{
  double sumei2 = 0.0;
  for(int i=0; i < mp->etaNum; ++i){
    sumei2 = sumei2 + std::pow(evec[i][p][q],2.0);
  }
  return 12.0*mp->B*((1.0-cxy)*eixy -\
         (2.0-cxy)*std::pow(eixy,2.0) +\
         eixy*sumei2);//Expression was analytically derived
}//END function



//Total  bulk free energy at a given time step
//Takes in ONE concentration field,
//and MULTIPLE non-conservative (Eta) parameter fields
double BulkFreeEnergy(double **cxy, std::vector<double **> evec,\
                      const computeParam *cp,const materialParam *mp)
{
  double totFreeEng = 0.0;
  double sumei2, sumei3;

  for(int i = 0; i < cp->Nx; ++i){
    for(int j = 0; j < cp->Ny; ++j){

      sumei2 = 0.0;sumei3 = 0.0;
      for(int ei = 0; ei < mp->etaNum; ++ei){
        sumei2 = sumei2 + pow(evec[ei][i][j],2.0);
        sumei3 = sumei3 + pow(evec[ei][i][j],3.0);
      }

      totFreeEng = totFreeEng + mp->A*std::pow((cxy[i][j]*(1.0-cxy[i][j])),2.0)\
                  + mp->B*(std::pow(cxy[i][j], 2.0) + 6.0*(1.0-cxy[i][j])*sumei2\
                  - 4.0*(2.0-cxy[i][j])*sumei3 + 3.0*std::pow(sumei2,2.0));
    }
  }
  return totFreeEng;
}//END function



// Concentration dependent mobility
double mobility(double cxy, std::vector<double **> evec,\
              const materialParam *mp,long p, long q)
{
  double mob = 0.0;
  double phi = std::pow(cxy,3.0)*(10.0-15.0*cxy+6.0*std::pow(cxy,2.0));
  double sumProdij = 0.0;
  for(int i = 0; i < mp->etaNum; ++i){
    for(int j = 0; j < mp->etaNum; ++j){
      if(i!=j){
        sumProdij = sumProdij + evec[i][p][q]*evec[j][p][q];
      }
    }
  }
  mob = mp->dvol*phi + \
        mp->dvap*(1.0-phi) + \
        mp->dsur*std::pow(cxy,2.0)*std::pow((1.0-cxy),2.0) +\
        //mp->dsur*cxy*(1.0-cxy) +
        mp->dgb*sumProdij;

  return mob;
}



// Solve differential equations to evolve concentration and the Eta fiels
void evolveMicrostrucutre(double **Cxy, std::vector<double **>Evec,\
                          const computeParam *cp,const materialParam *mp)
{
  std::cout << "Solving CH and AC coupled differential equations ... "<< std::endl;

  // Stores bulk Free Energy at each time step
  double *FBulk = new double[cp->nstep];
  // For storing time data
  double *step = new double[cp->nstep];
  double *time = new double[cp->nstep]; // To be used for plotting F vs. t
  double tt = 0; // time variable

  //For concentration fields for CH equations
  double **lapConc; // Laplacian of concetration
  double **idFdc = create2DField(cp); // Matrix evaluating inner derivative
  double **lapIdFdc; // Laplacian of the inner derivative

  //For non-conservative ith eta field for AC equations
  double **lapEtai; // Laplacian of ith eta
  double idFdetai; // Evaluation of inner derivative
  double **etaxyi  = create2DField(cp); // ith eta matrix for temporary storage

  std::ofstream outData;
  std::string fname = "small_circle_data_structureType" + std::string("-") + \
                      std::to_string(mp->initStructureType)+ ".txt";
  outData.open(fname);

  /*out << "TIME" << " " << "STEP" << " " <<  "ENERGY" << " " \
      << "DIAMETER" << "\n";*/

  // START time integration
  for(int t = 0; t < cp->nstep; ++t){
    step[t] = t;
    time[t] = tt;

  /*-------------------Evolve concentration fields-----------------------------*/
    lapConc = get2DLaplacian(Cxy, cp);// Laplacian of the concentration field

    //START loop for evaluating DF/Dcon - k.Laplacian(con)
    for(int i = 0; i < cp->Nx; i++){
      for(int j = 0; j < cp->Ny; j++){
        idFdc[i][j] =  dFdCon(Cxy[i][j],Evec,mp,i,j) - 0.5*mp->Kc*lapConc[i][j];
      }
     }//END loop for evaluating DF/Dcon - k.Laplacian(con)

     lapIdFdc = get2DLaplacian(idFdc,cp); //Laplacian of th einner derivative
     //START time integration loop for concentration
     for(int i = 0; i < cp->Nx; ++i){
       for(int j = 0; j < cp->Ny; ++j){

         Cxy[i][j] = Cxy[i][j] + mobility(Cxy[i][j],Evec,mp,i,j)\
                                 *cp->dtime*lapIdFdc[i][j];

         if(Cxy[i][j] >= 0.9999){ //Set MAX limits for small variations
           Cxy[i][j] = 0.9999;
         }

         if(Cxy[i][j] < 0.00001){ //Set MIN limits for small variations
           Cxy[i][j] = 0.00001;
         }
       }
      }//END time integration loop for concentration
  /*-------------------END evolving concentration fields-------------------*/


  /*-------------------Evolve non-consevatinve Eta fields-------------------*/
    for(int ei = 0; ei < mp->etaNum; ei++){//Run though number of Etas

      for(int i = 0; i < cp->Nx; i++){
        for(int j = 0; j < cp->Ny; j++){
          etaxyi[i][j] = Evec[ei][i][j];
        }
      }

      lapEtai = get2DLaplacian(etaxyi, cp);
      //START time integration loop for etai
      for(int i = 0; i < cp->Nx; i++){
        for(int j = 0; j < cp->Ny; j++){
          idFdetai = dFdEtai(Cxy[i][j],etaxyi[i][j],Evec,mp,i,j) - \
                           0.5*mp->Ke*lapEtai[i][j]; // DF/Detai - k.Laplacian(etai)
          etaxyi[i][j] = etaxyi[i][j] - cp->dtime*mp->L*idFdetai; // time step

          if(etaxyi[i][j] >= 0.9999){ //Set MAX limits for small variations
             etaxyi[i][j] = 0.9999;
          }
          if(etaxyi[i][j] < 0.0001){ //Set MIN limits for small variations
            etaxyi[i][j] = 0.0001;
          }
        }
      } //END time integration loop

      //HARD copying of temporary eta values to the Parent eta data strucutre
      for(int i = 0; i < cp->Nx; i++){
        for(int j = 0; j < cp->Ny; j++){
          Evec[ei][i][j] = etaxyi[i][j];
        }
      }
      //displayArray(Evec[0], cp);
    } // END etas run
  /*-------------------END evolving non-consevatinve Eta fields--------------*/

    if((t%cp->nprint == 0) && (t < cp->nstep)){
        writeAllDataToVTKFile(Cxy,Evec,cp,mp,t);
    }

    FBulk[t] = BulkFreeEnergy(Cxy,Evec,cp,mp);

    if(t%100 == 0 && cp->dataWRITE==1){
      if(mp->initStructureType==2 || mp->initStructureType==7){
         outData << time[t] << " " << step[t] << " " <<  FBulk[t] << " " \
         << getSmallCircleArea(Evec,cp,mp) << " " \
         << getSmallCircleNeckLength(Evec,cp,mp) << "\n";
      }else if(mp->initStructureType==5 || mp->initStructureType==6 || \
               mp->initStructureType==8 || mp->initStructureType==9 || \
               mp->initStructureType==11 || mp->initStructureType==12)
      {
        outData << time[t] << " " << step[t] << " " <<  FBulk[t] << " " \
        << getSmallCircleArea(Evec,cp,mp) << " " \
        << getSmallCircleNeckLength(Evec,cp,mp) << " " \
        << getLargestCircleArea(Evec,cp,mp) <<"\n";
      }

      outData.flush();

    if(cp->rerunVTKWRITE==1 && t!=0 && t==cp->writeAFTERIter)
    {
      int nVec = Evec.size();
      for(int i =0; i < nVec; i++)
      {
        std::string etaFname = "Eta" + std::string("-") + std::to_string(i);
        write2DVTKfile(Evec[i],cp,t,etaFname);
      }
      write2DVTKfile(Cxy,cp,t,"Conc");

    }
      /*std::cout << time[t] << " " << step[t] << " " <<  FBulk[t] << " " \
               << getSmallCircleArea(Evec,cp,mp) << " " \
               << getSmallCircleNeckLength(Evec,cp,mp) << std::endl;*/
    }

    tt = tt + cp->dtime;
  }//END MAIN time integration loop


//write1DTXTfile(time, FBulk, cp->nstep, "energy-time.txt");
//write1DTXTfile(step, FBulk, cp->nstep, "energy-step.txt");
//writeAllDataToVTKFile(Cxy,Evec,cp,mp,cp->nstep);

std::cout << "Finished solving CH and AC coupled differential equations "<< std::endl;

//Garbage cleaning
delete[] FBulk;
delete[] time;
delete2DArray(lapConc,cp->Nx);
delete2DArray(lapIdFdc,cp->Nx);
delete2DArray(etaxyi,cp->Nx);

outData.close();

}//END function


/*

Get neck length of the small particle/circle

*/
double getSmallCircleNeckLength(std::vector<double **> eVec,
                                const computeParam *cp,const materialParam *mp)
{
  double neckLength = 0.0;
  int istrucType = mp->initStructureType;

  for(int i = 0; i < cp->Nx ; ++i)
  {
    for(int j = 0; j < cp->Ny ; ++j)
    {
      // eVec[2] is the small circle in 8, 5, 2, 6, 7
      // eVec[3] is the small circle in 9

      if(istrucType==8)
        neckLength = neckLength + eVec[2][i][j]*eVec[0][i][j] +\
                     eVec[2][i][j]*eVec[1][i][j];
        //neckLength = neckLength + eVec[1][i][j]*eVec[0][i][j];
      else if(istrucType==9)
        neckLength = neckLength + eVec[3][i][j]*eVec[1][i][j] +\
                     eVec[3][i][j]*eVec[2][i][j];
      else if(istrucType==5)
        neckLength = neckLength + eVec[2][i][j]*eVec[0][i][j] +\
                     eVec[2][i][j]*eVec[1][i][j];
      else if(istrucType==6)
        neckLength = neckLength + eVec[1][i][j]*eVec[0][i][j];
      else if(istrucType==2)
        neckLength = neckLength + eVec[1][i][j]*eVec[0][i][j];
      else if(istrucType==7)
        neckLength = neckLength + eVec[2][i][j]*eVec[0][i][j] + \
                     eVec[2][i][j]*eVec[1][i][j];
      else if(istrucType==11)
              neckLength = neckLength +\
                           eVec[3][i][j]*eVec[0][i][j] +\
                           eVec[3][i][j]*eVec[1][i][j] +\
                           eVec[3][i][j]*eVec[2][i][j];
      else if(istrucType==12)
              neckLength = neckLength +\
                           eVec[4][i][j]*eVec[0][i][j] +\
                           eVec[4][i][j]*eVec[1][i][j] +\
                           eVec[4][i][j]*eVec[2][i][j] +\
                           eVec[4][i][j]*eVec[3][i][j];
    }
  }
  return neckLength;
}

/*
Get the area of the smallest cicle
*/
double getSmallCircleArea(std::vector<double **> eVec,
                          const computeParam *cp,const materialParam *mp)
{
  int istrucType = mp->initStructureType;
  double area = 0.0;

  int circleNo = 0;

  if(istrucType==5)
  {
    circleNo = 2;
  }else if(istrucType==8)
  {
    circleNo = 2;
  }else if(istrucType==6)
  {
    circleNo = 1;
  }else if(istrucType==9  || istrucType==11)
  {
    circleNo = 3;
  }else if(istrucType==2)
  {
    circleNo = 1;
  }else if(istrucType==7)
  {
    circleNo = 2;
  }else if(istrucType==12)
  {
    circleNo = 4;
  }


  for(int i = 0; i < cp->Nx ; ++i)
  {
    for(int j = 0; j < cp->Ny ; ++j)
    {
      if(eVec[circleNo][i][j] > 0.01 && eVec[circleNo][i][j] < 1.0)
        area = area + cp->dx*cp->dy;
    }
  }

  return area;
}


/*
Get the area of the Largest cicle
*/
double getLargestCircleArea(std::vector<double **> eVec,
                            const computeParam *cp,const materialParam *mp)
{
  int istrucType = mp->initStructureType;
  double area = 0.0;

  int circleNo = 0;

  if(istrucType==5)
  {
    circleNo = 1;
  }else if(istrucType==8)
  {
    circleNo = 0;
  }else if(istrucType==6)
  {
    circleNo = 0;
  }else if(istrucType==9)
  {
    circleNo = 2;
  }

  for(int i = 0; i < cp->Nx ; ++i)
  {
    for(int j = 0; j < cp->Ny ; ++j)
    {
      if(eVec[circleNo][i][j] > 0.01 && eVec[circleNo][i][j] < 1.0)
        area = area + cp->dx*cp->dy;
    }
  }

  return area;
}










































//
