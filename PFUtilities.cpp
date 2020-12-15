#include "PFUtilities.h"

/**--------------------------------------------------------------------------

For processing 2D datasets

--------------------------------------------------------------------------**/


/*
Apply 2D periodic boundary condition
*/
void apply2DPeriodicBC(int *im, int *ip, int *jm, int *jp, const computeParam *cp)
{
  if(*im <= 0){
    *im = cp->Nx - 1;
  }
  if(*ip >= cp->Nx){
    *ip = 0;
  }

  if(*jm <= 0){
    *jm = cp->Ny - 1;
  }
  if(*jp >= cp->Ny){
    *jp = 0;
  }
} //END Function


/*
Apply 2D zero flux boundary condition
*/
void apply2DZeroFluxBC(int *im, int *ip, int *jm, int *jp, const computeParam *cp)
{
  if(*im <= 0){
    *im = 2;
  }
  if(*ip >= cp->Nx){
    *ip = cp->Nx-2;
  }

  if(*jm <= 0){
    *jm = 1;
  }
  if(*jp >= cp->Ny){
    *jp = cp->Nx-2;
  }
} //END Function


/*
Get centered differnece for a partial derivative
*/
double getCenteredDifference(double **Pxy, int i, int j, \
                             int im, int ip, double delta)
{
  double val = (Pxy[ip][j] - Pxy[im][j])/(2.0*delta);
  return val;
}//END FUNCTION


/*
Get Laplacian using the 5-point stencil method.Returns one number
*/
double get2DLaplacian(double **Pxy, int i, int j, \
                      int im, int ip, int jm, int jp,
                      const computeParam *cp)
{
  double PxyaE = Pxy[ip][j];
  double PxyaW = Pxy[im][j];
  double PxyaN = Pxy[i][jp];
  double PxyaS = Pxy[i][jm];
  double PxyaC = Pxy[i][j];

  double val = (PxyaN + PxyaS + PxyaE + PxyaW - 4.0*PxyaC)/(cp->dx*cp->dy);

  return val;

}//END FUNCTION



/* Get Lapnacian using the 5-point stencil method.
Also check for periodic boundary conditions*/
double **get2DLaplacian(double **Axy, const computeParam *cp)
{
  int im, ip, jm, jp;

  double **dummy = create2DField(cp);
  /*dummy = new double*[cp->Nx];
  for(int i = 0; i < cp->Nx; i++){
    dummy[i] = new double[cp->Ny];
  }*/

  for(int i = 0; i< cp->Nx; i++){
    for(int j = 0; j< cp->Ny; j++){
      jp = j + 1;
      jm = j - 1;
      ip = i + 1;
      im = i - 1;

      // Check for periodic boundary conditions
      apply2DPeriodicBC(&im, &ip, &jm, &jp, cp);
      /*if(im < 0){
        im = cp->Nx - 1;
      }
      if(ip == cp->Nx){
        ip = 0;
      }

      if(jm < 0){
        jm = cp->Ny - 1;
      }
      if(jp == cp->Ny){
        jp = 0;
      }*/
      //Creating five point stencil points
      double aE = Axy[ip][j];
      double aW = Axy[im][j];
      double aN = Axy[i][jp];
      double aS = Axy[i][jm];
      double aC = Axy[i][j];

      //evaluating values usinf all the five stencil points
      dummy[i][j] = (aN + aS + aE + aW - 4.0*aC)/(cp->dx*cp->dy);
    }
  }
  return dummy;
}//END function



/*Prepare the X and Y -axis basis set for FFT computation
NOTE that the double arrays Kx and Ky are modifued within the function*/
void prepare2DFFTBasis(double *Kx, double*Ky, const computeParam *cp)
{
  long nx21 = cp->Nx/2 + 1;
  long nx2 = cp->Nx + 2;

  long ny21 = cp->Ny/2 + 1;
  long ny2 = cp->Ny + 2;

  double delKx = (2.0*PI)/(cp->Nx*cp->dx);
  double delKy = (2.0*PI)/(cp->Ny*cp->dy);

  double tempx = 0.0;
  for(int i = 0; i < nx21; i++){
    tempx = i*delKx;
    Kx[i] = tempx;
    Kx[nx2-i-2] = (-1.0)*tempx;
  }

  double tempy = 0.0;
  for(int i = 0; i < ny21; i++){
    tempy = i*delKy;
    Ky[i] = tempy;
    Ky[ny2-i-2] = (-1.0)*tempy;
  }
} //END Function


/*Returns the the nth power of K2 = Kx^2 + Ky^2*/
double *getK2Power2DFFTBasis(double nthPower, const computeParam *cp)
{
  double *Kx = new double[cp->Nx];
  double *Ky = new double[cp->Ny];
  double *K2 = new double[cp->Nx * cp->Ny];

  double *Kn = new double[cp->Nx * cp->Ny];

  prepare2DFFTBasis(Kx, Ky, cp); //Prepares the X and Y basis

  int k = 0;
  for(int i = 0; i< cp->Nx; i++){
    for(int j = 0; j< cp->Ny; j++){
      K2[k] = std::pow(Kx[i],2.0) + std::pow(Ky[j],2.0);
      Kn[k] = std::pow(K2[k],nthPower);
      k++;
    }
  }

  delete []Kx;
  delete []Ky;

  return Kn;

}//END function



/*Prepare the squared version of the X and Y -axis basis set for FFT computation
i.e. K2 = Kx^2 + Ky^2 and K4 = K2^2
NOTE: K2 and K4 are 1D arrays*/
void prepare2DFFTSquaredBasis(double *K2, double *K4, const computeParam *cp)
{
  double *Kx = new double[cp->Nx];
  double *Ky = new double[cp->Ny];

  prepare2DFFTBasis(Kx, Ky, cp); //Prepares the X and Y basis

  int k = 0;
  for(int i = 0; i< cp->Nx; i++){
    for(int j = 0; j< cp->Ny; j++){
      K2[k] = std::pow(Kx[i],2.0) + std::pow(Ky[j],2.0);
      K4[k] = std::pow(K2[k],2.0);
      k++;
    }
  }

  delete []Kx;
  delete []Ky;

}//END function




/* Function returns FFT of a 2D feild contanind real values of the FFT.
Use the get2DREAL function to get a 2D array if real values*/
fftw_complex *get2DFFT(double **Axy, long Nx, long Ny)
{
  //Create the 2D array for FFT
  double **fftData = new double*[3];

  //Allocate memory to the array
  for (int i = 0; i < 3; i++) {
    fftData[i] = new double[3];
  }
  //Create input and output FFTs
  fftw_complex *in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
  fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);

  //Create plan
  fftw_plan plan = fftw_plan_dft_2d(Nx, Ny, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Fill input array
    int k = 0;
    for(int i = 0; i < Nx; i++){
      for(int j = 0; j < Ny; j++){
        in[k][REAL] = Axy[i][j];
        in[k++][IMAG] = 0.0;
      }
    }

  //Execute FFT plan
  fftw_execute(plan);
  //Garbage cleaning
  fftw_destroy_plan(plan);
  fftw_free(in);
  //fftw_free(out);

  return out;
} //End function



/* Function returns the complete inverse FFT of a 2D feild.
Use the get2DREAL function to get a 2D array if real values*/
fftw_complex *get2DiFFT(fftw_complex *inFFT, long Nx, long Ny)
{
  double normalization = Nx*Ny; //Normalization factor

  //Create output FFTs
  fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
  //Create plan
  fftw_plan plan = fftw_plan_dft_2d(Nx, Ny, inFFT, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  //Execute FFT plan
  fftw_execute(plan);

  //scale the output to the obtain the exact result
  for(int i = 0; i < Nx*Ny; i++){
    out[i][REAL] /= normalization;
    out[i][IMAG] /= normalization;
  }
  //Garbage cleaning
  fftw_destroy_plan(plan);
  return out;
} //End Function



/* Function returns INVERSE FFT of a 2D feild by taking two arrays as inputs
First ayyaryu contains real values while the second contains complex values*/
fftw_complex *get2DiFFTAfterREALManipulation(double **REALxy, double **CMPLXxxy, long Nx, long Ny)
{
  double normalization = Nx*Ny; //Normalization factor
  //Create input and output FFTs
  fftw_complex *in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
  fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);

  //Create plan
  fftw_plan plan = fftw_plan_dft_2d(Nx, Ny, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

  // Fill input array
  int k = 0;
  for(int i = 0; i < Nx; i++){
    for(int j = 0; j < Ny; j++){
      in[k][REAL]   = REALxy[i][j];
      in[k++][IMAG] = CMPLXxxy[i][j];
    }
  }
  //Execute FFT plan
  fftw_execute(plan);
  //scale the output to the obtain the exact result
  for(int i = 0; i < Nx*Ny; i++){
    out[i][REAL] /= normalization;
    out[i][IMAG] /= normalization;
  }
  //Garbage cleaning
  fftw_destroy_plan(plan);
  fftw_free(in);
  //fftw_free(out);

  return out;
} //End function



/* Function sets up a 2D array of double
NOTE: Be sure to create/allocate memory to the 2D array before passing to this function*/
void set2DiFFTArray(fftw_complex *inFFT, double **ifftArray, long Nx, long Ny)
{
  double normalization = Nx*Ny; //Normalization factor

  //Create output FFTs
  fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
  //Create plan
  fftw_plan plan = fftw_plan_dft_2d(Nx, Ny, inFFT, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  //Execute FFT plan
  fftw_execute(plan);
  //scale the output to the obtain the exact result
  for(int i = 0; i < Nx*Ny; i++){
    out[i][REAL] /= normalization;
    out[i][IMAG] /= normalization;
  }
  // Copy real portion of the FFT to the new array
  int k = 0;
  for(int i = 0; i < Nx; i++){
    for(int j = 0; j < Ny; j++){
        ifftArray[i][j] =  out[k++][REAL];
    }
  }
  //Garbage cleaning
  fftw_destroy_plan(plan);
  fftw_free(out);
} //End Function



/*This function extracts the REAL values from the FFT, fills up those values in
a 2D array and returns to the calling point*/
double **get2DREAL(fftw_complex *inFFT, long Nx, long Ny)
{
  //Allocate memory to the 2D array
  double **realArray = new double*[Nx];
  for (int i = 0; i < Nx; i++) {
    realArray[i] = new double[Ny];
  }
  // Copy real portion of the FFT to the new array
  int k = 0;
  for(int i = 0; i < Nx; i++){
    for(int j = 0; j < Ny; j++){
        realArray[i][j] =  inFFT[k++][REAL];
    }
  }
  return realArray;
} //END function



/*This function extracts the  COMPLEX values from the FFT, fills up those values in
a 2D array and returns to the calling point*/
double **get2DCMPLX(fftw_complex *inFFT, long Nx, long Ny)
{
  double **cmplxArray = new double*[Nx];
  for (int i = 0; i < Nx; i++) {
    cmplxArray[i] = new double[Ny];
  }
  // Copy real portion of the FFT to the new array
  int k = 0;
  for(int i = 0; i < Nx; i++){
    for(int j = 0; j < Ny; j++){
        cmplxArray[i][j] =  inFFT[k++][IMAG];
    }
  }
  return cmplxArray;
}//END function




void display2DFFTdata(fftw_complex *fftData, long Nx, long Ny)
{
  std::cout << std::endl;
  int k = 0;
  for(int i = 0; i < Nx; i++){
    for(int j = 0; j < Ny; j++){
        std::cout << fftData[k][REAL] << " + i*" << fftData[k][IMAG] <<"   ";
        k++;
    }
    std::cout << std::endl;
  }
} //END Function



/* Set a double 2D array to zero*/
void setZero2DField(double **Fxy, const computeParam *cp){
    for(int i = 0; i < cp->Nx; ++i){
      for(int j = 0; j < cp->Nx; ++j){
        Fxy[i][j] = 0.0;
      }
    }
}//END function



/*Funciton for deleting 2D SQUARE array of typedouble*/
void delete2DArray(double** Arr, long arraySize)
{
  for(int i = 0; i < arraySize; ++i){ //For 2D array of pointers
    delete[] Arr[i];
  }
}//END function


//Display the 2D array oin the screen. MUST be a square matrix
void displayArray(double **Arr, const computeParam *cp)
{
  std::cout << std::endl;
  for (int i = 0; i < cp->Nx; i++) {
    for (int j = 0; j < cp->Ny; j++) {
      std::cout << Arr[i][j] << " ";
    }
    std::cout << std::endl;
  }
}//END function


/* Creates a 2D array and sets the value of individual elements to ZERO*/
double **create2DField(const computeParam * cp){
  //Allocate memory for the concentration field
  double **Fxy = new double*[cp->Nx];
  for(int i = 0; i < cp->Nx; i++){
    Fxy[i] = new double[cp->Ny];
  }
  setZero2DField(Fxy, cp);
  return Fxy;
}//END function



/*
Print file for two arrays contnaing double - basaically x y for plotting
*/
void write1DTXTfile(double *x, double *y, long Nx, std::string filename)
{
  std::ofstream out;
  out.open(filename);
  //std::cout << "Number of items: " << Nx << std::endl;
  for(int i = 0; i < Nx; i++){
    out << x[i] << "  " << y[i] << "\n";
  }
  out.close();
}//END function



//Print  2D array to a file
void write2DData(double **Arr, const computeParam *cp,std::string filename)
{
  std::ofstream out;
  out.open(filename);
  for (int i = 0; i < cp->Nx; i++) {
    for (int j = 0; j < cp->Ny; j++) {
      out << Arr[i][j] << " ";
      //std::cout <<  Arr[i][j] << " ";
    }
    out << std::endl;
  }
  out.close();
}//END function


/***
 *
 * This program reads a 2D VTK file and outputs a pointer to an array of 2D data
 * Author: Logan Blake
 *
 **/

double **read2DVTKFile(std::string filePath) {

    //Read file and set global variables
    std::string line;
    std::ifstream inFile (filePath);
    int nx = 0;
    int ny = 0;
    int nz = 0;
    std::vector<std::vector<float> > data;
    int lines_read = 0;
    int data_line = 0;

    //Read DTK file and process lines
    if (inFile.is_open()) {
        while (std::getline (inFile,line)) {
            lines_read++;
            std::stringstream lineStream(line);
            //If line contains "dimensions" saves them to nx,ny,nz variables
            if (line.find("DIMENSIONS") != std::string::npos){
                std::vector<int> dimensions;
                std::string dimensions_temp;
                int dimensions_found;

                while(!lineStream.eof()){
                    lineStream >> dimensions_temp;
                    if (std::stringstream(dimensions_temp)>>dimensions_found)
                        dimensions.push_back(dimensions_found);
                }

                nx = dimensions[0];
                ny = dimensions[1];
                nz = dimensions[2];
            }
            //Find where the point data is
            if (line.find("POINT_DATA") != std::string::npos){
                data_line = lines_read+2;
            }

            //If line is > specified skipped lines (header information)
            //Push data to vector then push vector to data vector
            if (data_line != 0 and lines_read>data_line and lines_read<(nx*ny*nz)+data_line+1) {
                std::vector<float> line_data;
                std::string temp;
                float found;
                while(!lineStream.eof()) {
                    lineStream >> temp;
                    if (std::stringstream(temp)>>found) {
                        line_data.push_back(found);
                    }
                }
                data.push_back(line_data);
            }
        }
    }
    inFile.close();

  //Create a 2D array of pointers -- Deep added
   double **vtkData = new double*[nx];
   for(int i = 0; i < ny; i++){
	  vtkData[i] = new double[ny];
   }

    //Process vectors from data vector and fill array -- Deep added
    int rowCount = 0;
    for (int i = 0; i<nx; i++) {
        for (int j = 0; j<ny; j++) {
            vtkData[i][j] = data[rowCount][0];
	    rowCount++;
        }
    }
    return vtkData;
}//END function


/***
 *
 * This program reads a variable-value pairs from an input files
 *
 **/

std::map<std::string, double>readParameterFile(const std::string& inputFile)
{
  typedef std::map<std::string, double> ConfigInfo;
  ConfigInfo keyValue;

  //std::cout << "Inside the function" <<std::endl;
  std::ifstream inFile(inputFile);
  std::string line;

   while (std::getline(inFile, line)) {
     //std::cout << "what" << std::endl;
     std::istringstream is_line(line);
     std::string key;
     if (std::getline(is_line, key,'=')) {
       std::string value;
       if(key[0] == '#')
        continue;
       if(std::getline(is_line, value))
       {
         std::string::size_type sz;
         keyValue[key] = std::stod(value, &sz);
       }
     }
   }

  inFile.close();

  return keyValue;
}//END function


/**--------------------------------------------------------------------------

Write 2D VTK files

--------------------------------------------------------------------------**/


/* Write data to VTK file with on ONE type of dataset*/
void write2DVTKfile(double **Cxy, const computeParam *cp, \
                    long tstep, std::string filename)
{
  std::ofstream out;
  std::string fname = filename + std::string("_") + std::to_string(tstep) + ".vtk";
  out.open(fname);

  long nx, ny, nz;
  nx = cp->Nx;
  ny = cp->Ny;
  nz = cp->Nz;

  long nTot = nx*ny*nz; //total number of data or grid points

  //Start writting to file

  //Header information
  out << "# vtk DataFile Version 2.0\n";
  out << "time_10.vtk\n";
  out << "ASCII\n";
  out << "DATASET STRUCTURED_GRID\n";

  //Coordinates of grid points
  out << "DIMENSIONS " << nx << " " << ny << " " << nz <<"\n";
  out << "POINTS " << nTot << " float"<<"\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << (i-1)*cp->dx << " " << (j-1)*cp->dy << " " << 0.0 <<"\n";
    }
  }//END function

  //Actual grid point values
  out << "POINT_DATA " << nTot <<"\n";
  out << "SCALARS PHI float 1\n";
  out << "LOOKUP_TABLE default\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << Cxy[i][j] << "\n";
    }
  }

  out.close();
}//END function



/* Write data to VTK file with on ONE type of dataset*/
void write2DVTKfile(double **Cxy, const computeParam *cp, long tstep)
{
  std::ofstream out;
  std::string fname = std::string("time_") + std::to_string(tstep) + ".vtk";
  out.open(fname);

  long nx, ny, nz;
  nx = cp->Nx;
  ny = cp->Ny;
  nz = cp->Nz;

  long nTot = nx*ny*nz; //total number of data or grid points

  //Start writting to file

  //Header information
  out << "# vtk DataFile Version 2.0\n";
  out << "time_10.vtk\n";
  out << "ASCII\n";
  out << "DATASET STRUCTURED_GRID\n";

  //Coordinates of grid points
  out << "DIMENSIONS " << nx << " " << ny << " " << nz <<"\n";
  out << "POINTS " << nTot << " float"<<"\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << (i-1)*cp->dx << " " << (j-1)*cp->dy << " " << 0.0 <<"\n";
    }
  }//END function

  //Actual grid point values
  out << "POINT_DATA " << nTot <<"\n";
  out << "SCALARS CON float 1\n";
  out << "LOOKUP_TABLE default\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << Cxy[i][j] << "\n";
    }
  }

  out.close();
}//END function


/* Write data to VTK file with on TWO types of dataset*/
void write2DVTKfile(double **Cxy, double **Exy, const computeParam *cp, long tstep)
{
  std::ofstream out;
  std::string fname = std::string("time_") + std::to_string(tstep) + ".vtk";
  out.open(fname);

  long nx, ny, nz;
  nx = cp->Nx;
  ny = cp->Ny;
  nz = cp->Nz;

  long nTot = nx*ny*nz; //total number of data or grid points

  //Start writting to file

  //Header information
  out << "# vtk DataFile Version 2.0\n";
  out << "time_10.vtk\n";
  out << "ASCII\n";
  out << "DATASET STRUCTURED_GRID\n";

  //Coordinates of grid points
  out << "DIMENSIONS " << nx << " " << ny << " " << nz <<"\n";
  out << "POINTS " << nTot << " float"<<"\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << (i-1)*cp->dx << " " << (j-1)*cp->dy << " " << 0.0 <<"\n";
    }
  }//END function

  //Actual grid point concentration values
  out << "POINT_DATA " << nTot <<"\n";
  out << "SCALARS CON float 1\n";
  out << "LOOKUP_TABLE default\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << Cxy[i][j] << "\n";
    }
  }

  //Actual grid point non-conservative Eta values
  //out << "POINT_DATA " << nTot <<"\n";
  out << "SCALARS ETAs float 1\n";
  out << "LOOKUP_TABLE default\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << Exy[i][j] << "\n";
    }
  }

  out.close();
}//END function



/* Write data to VTK file with on TWO types of dataset*/
void write2DVTKfile(double **Cxy, double **Exy, const computeParam *cp,\
                    long tstep, std::string filename)
{
  std::ofstream out;
  std::string fname = filename + std::string("_") + std::to_string(tstep) + ".vtk";
  out.open(fname);

  long nx, ny, nz;
  nx = cp->Nx;
  ny = cp->Ny;
  nz = cp->Nz;

  long nTot = nx*ny*nz; //total number of data or grid points

  //Start writting to file

  //Header information
  out << "# vtk DataFile Version 2.0\n";
  out << "time_10.vtk\n";
  out << "ASCII\n";
  out << "DATASET STRUCTURED_GRID\n";

  //Coordinates of grid points
  out << "DIMENSIONS " << nx << " " << ny << " " << nz <<"\n";
  out << "POINTS " << nTot << " float"<<"\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << (i-1)*cp->dx << " " << (j-1)*cp->dy << " " << 0.0 <<"\n";
    }
  }//END function

  //Actual grid point concentration values
  out << "POINT_DATA " << nTot <<"\n";
  out << "SCALARS CON float 1\n";
  out << "LOOKUP_TABLE default\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << Cxy[i][j] << "\n";
    }
  }

  //Actual grid point non-conservative Eta values
  //out << "POINT_DATA " << nTot <<"\n";
  out << "SCALARS ETAs float 1\n";
  out << "LOOKUP_TABLE default\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << Exy[i][j] << "\n";
    }
  }

  out.close();
}//END function


/* Write data to VTK file with on TWO types of dataset*/
void write2DVTKfile(double **P1xy, std::string paramName1,\
                    double **P2xy, std::string paramName2, \
                    const computeParam *cp,\
                    long tstep, std::string filename)
{
  std::ofstream out;
  std::string fname = filename + std::string("_") + std::to_string(tstep) + ".vtk";
  out.open(fname);

  long nx, ny, nz;
  nx = cp->Nx;
  ny = cp->Ny;
  nz = cp->Nz;

  long nTot = nx*ny*nz; //total number of data or grid points

  //Start writting to file

  //Header information
  out << "# vtk DataFile Version 2.0\n";
  out << "DATA\n";
  out << "ASCII\n";
  out << "DATASET STRUCTURED_GRID\n";

  //Coordinates of grid points
  out << "DIMENSIONS " << nx << " " << ny << " " << nz <<"\n";
  out << "POINTS " << nTot << " float"<<"\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << (i-1)*cp->dx << " " << (j-1)*cp->dy << " " << 0.0 <<"\n";
    }
  }//END function

  //Actual grid point concentration values
  out << "POINT_DATA " << nTot <<"\n";
  out << "SCALARS " << paramName1 << " float 1\n";
  out << "LOOKUP_TABLE default\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << P1xy[i][j] << "\n";
    }
  }

  //Actual grid point non-conservative Eta values
  //out << "POINT_DATA " << nTot <<"\n";
  //out << "SCALARS ETAs float 1\n";
  out << "SCALARS " << paramName2 << " float 1\n";
  out << "LOOKUP_TABLE default\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << P2xy[i][j] << "\n";
    }
  }

  out.close();
}//END function



/* Write data to VTK file with on THREE types of dataset*/
void write2DVTKfile(double **P1xy, std::string paramName1,\
                    double **P2xy, std::string paramName2,\
                    double **P3xy, std::string paramName3,\
                    const computeParam *cp,\
                    long tstep, std::string filename)
{
  std::ofstream out;
  std::string fname = filename + std::string("_") + std::to_string(tstep) + ".vtk";
  out.open(fname);

  long nx, ny, nz;
  nx = cp->Nx;
  ny = cp->Ny;
  nz = cp->Nz;

  long nTot = nx*ny*nz; //total number of data or grid points

  //Start writting to file

  //Header information
  out << "# vtk DataFile Version 2.0\n";
  out << "time_10.vtk\n";
  out << "ASCII\n";
  out << "DATASET STRUCTURED_GRID\n";

  //Coordinates of grid points
  out << "DIMENSIONS " << nx << " " << ny << " " << nz <<"\n";
  out << "POINTS " << nTot << " float"<<"\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << (i-1)*cp->dx << " " << (j-1)*cp->dy << " " << 0.0 <<"\n";
    }
  }//END function

  //Actual grid point concentration values
  out << "POINT_DATA " << nTot <<"\n";
  out << "SCALARS " << paramName1 << " float 1\n";
  out << "LOOKUP_TABLE default\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << P1xy[i][j] << "\n";
    }
  }

  //Actual grid point non-conservative Eta values
  //out << "POINT_DATA " << nTot <<"\n";
  //out << "SCALARS ETAs float 1\n";
  out << "SCALARS " << paramName2 << " float 1\n";
  out << "LOOKUP_TABLE default\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << P2xy[i][j] << "\n";
    }
  }

  out << "SCALARS " << paramName3 << " float 1\n";
  out << "LOOKUP_TABLE default\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << P3xy[i][j] << "\n";
    }
  }

  out.close();
}//END function



/* Write data to VTK file with on THREE types of dataset*/
void write2DVTKfile(double **P1xy, std::string paramName1,\
                    double **P2xy, std::string paramName2,\
                    double **P3xy, std::string paramName3,\
                    double **P4xy, std::string paramName4,\
                    const computeParam *cp,\
                    long tstep, std::string filename)
{
  std::ofstream out;
  std::string fname = filename + std::string("_") + std::to_string(tstep) + ".vtk";
  out.open(fname);

  long nx, ny, nz;
  nx = cp->Nx;
  ny = cp->Ny;
  nz = cp->Nz;

  long nTot = nx*ny*nz; //total number of data or grid points

  //Start writting to file

  //Header information
  out << "# vtk DataFile Version 2.0\n";
  out << "time_10.vtk\n";
  out << "ASCII\n";
  out << "DATASET STRUCTURED_GRID\n";

  //Coordinates of grid points
  out << "DIMENSIONS " << nx << " " << ny << " " << nz <<"\n";
  out << "POINTS " << nTot << " float"<<"\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << (i-1)*cp->dx << " " << (j-1)*cp->dy << " " << 0.0 <<"\n";
    }
  }//END function

  out << "POINT_DATA " << nTot <<"\n";
  out << "SCALARS " << paramName1 << " float 1\n";
  out << "LOOKUP_TABLE default\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << P1xy[i][j] << "\n";
    }
  }

  out << "SCALARS " << paramName2 << " float 1\n";
  out << "LOOKUP_TABLE default\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << P2xy[i][j] << "\n";
    }
  }

  out << "SCALARS " << paramName3 << " float 1\n";
  out << "LOOKUP_TABLE default\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << P3xy[i][j] << "\n";
    }
  }

  out << "SCALARS " << paramName4 << " float 1\n";
  out << "LOOKUP_TABLE default\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << P4xy[i][j] << "\n";
    }
  }

  out.close();
}//END function



/**--------------------------------------------------------------------------

For processing 3D datasets

--------------------------------------------------------------------------**/


/* Prints the 3D array on screen*/
void display3DArray(double ***Fxyz, const computeParam *cp){
    for(int i = 0; i < cp->Nx; ++i){
      for(int j = 0; j < cp->Nx; ++j){
        for(int k = 0; k < cp->Nz; ++k)
          std::cout<< Fxyz[i][j][k] << std::endl;
      }
    }
}//END function



/* Set a double 2D array to zero*/
void setZero3DField(double ***Fxyz, const computeParam *cp){
    for(int i = 0; i < cp->Nx; ++i){
      for(int j = 0; j < cp->Nx; ++j){
        for(int k = 0; k < cp->Nz; ++k)
        Fxyz[i][j][k] = 0.0;
      }
    }
}//END function



/* Creates a 3D array and sets the value of individual elements to ZERO*/
double ***create3DField(const computeParam * cp){
  //Allocate memory for the concentration field
  double ***Fxyz = new double**[cp->Nx];
  for(int i = 0; i < cp->Nx; i++){
    Fxyz[i] = new double*[cp->Ny];
    for(int j = 0; j < cp->Ny; j++){
      Fxyz[i][j] = new double[cp->Nz];
    }
  }
  setZero3DField(Fxyz, cp);
  return Fxyz;
}//END function


/* Write data to VTK file with on TWO types of dataset*/
void write3DVTKfile(double ***P1xyz, std::string paramName1,\
                    double ***P2xyz, std::string paramName2, \
                    const computeParam *cp,\
                    long tstep, std::string filename)
{
  std::ofstream out;
  std::string fname = filename + std::string("_") + std::to_string(tstep) + ".vtk";
  out.open(fname);

  long nx, ny, nz;
  nx = cp->Nx;
  ny = cp->Ny;
  nz = cp->Nz;

  long nTot = nx*ny*nz; //total number of data or grid points

  //Start writting to file

  //Header information
  out << "# vtk DataFile Version 2.0\n";
  out << "time_10.vtk\n";
  out << "ASCII\n";
  out << "DATASET STRUCTURED_GRID\n";

  //Coordinates of grid points
  out << "DIMENSIONS " << nx << " " << ny << " " << nz <<"\n";
  out << "POINTS " << nTot << " float"<<"\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      for(int k = 0; k < nz; k++){
        out << (i-1)*cp->dx << " " << (j-1)*cp->dy << " " << (k-1)*cp->dz <<"\n";
      }
    }
  }//END function

  //Actual grid point concentration values
  out << "POINT_DATA " << nTot <<"\n";
  out << "SCALARS " << paramName1 << " float 1\n";
  out << "LOOKUP_TABLE default\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      for(int k = 0; k < nz; k++){
        out << P1xyz[i][j][k] << "\n";
      }
    }
  }

  //Actual grid point non-conservative Eta values
  //out << "POINT_DATA " << nTot <<"\n";
  //out << "SCALARS ETAs float 1\n";
  out << "SCALARS " << paramName2 << " float 1\n";
  out << "LOOKUP_TABLE default\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      for(int k = 0; k < nz; k++){
        out << P2xyz[i][j][k] << "\n";
      }
    }
  }

  out.close();
}//END function













//
