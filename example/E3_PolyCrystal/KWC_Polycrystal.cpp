#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <vector>

//Relevant Library
#include "DataOut.h"
#include "InitCrystal.h"
#include "PostProcessor.h"
#include "PrimalDual.h"
#include "KWCThresholding.h"

using namespace std;
using namespace DataOut;


int main(int argc, char *argv[]){

  int n1=atoi(argv[1]);
  int n2=atoi(argv[2]);
  int n3=atoi(argv[3]);
  double epsilon=atof(argv[4]);

  
  /* Construct Classes */
  int const DIM=3;
  
  /* number of grain in the domain*/
  const unsigned int lcount=50;
    
  int pcount = n1*n2*n3;
  double dt= epsilon * epsilon;
  double *Xangles = new double[lcount]();
  double *Yangles = new double[lcount]();
  double *Zangles = new double[lcount]();
  double *eta = new double[pcount]();
  int *labels = new int[pcount]();
  double *JField = new double[pcount]();

  int Nthread = 1; //Number total thread to be used for FFTW
  
  /* Create 2D voronoi crystal, Z orientaion lies in [0,maxZangle]*/
  double maxZangle = 70.0* M_PI/180.0;
  InitializeCrystal::RandomCrystalConfiguration2D(n3,n2,n1,lcount,labels, maxZangle, Xangles,
                                                       Yangles,Zangles);
 
  //Set material type
  char materialType='S';
  Material material(n3,n2,n1,materialType);
  material.simpleKWC_s=1.0; //set material constant value s
  material.setUpClass(Xangles,Yangles,Zangles,labels,JField);
 
  //Construct PD Algorithm class
  
  double PDerror=1e-6; // tolerance of Primal-dual algorithm
  int PDmaxIters=10000; // allowable iteration number of the algorithm
  PrimalDual<DIM> EtaSubProblem(n3, n2, n1, PDerror, PDmaxIters,lcount,epsilon, Nthread);
  
  //Criteria for identifying interior regions of grain
  double initThresCriteria = 2*epsilon ;
  KWCThreshold<DIM> FastMarching(n3,n2,n1,lcount,initThresCriteria);

  EtaSubProblem.setUpClass(eta, Xangles, Yangles, Zangles, labels, JField);
  /* 'P' stands for Dirichlet boundary condition will be enfored during thresholding */
  FastMarching.setUpClass(eta, Xangles, Yangles, Zangles, labels, JField,'P');
  
  
  /* FFMEPG movie data */
  unsigned char *pixels=new unsigned char[pcount];
  unsigned char *colors=new unsigned char[lcount];
  for(int l=0;l<lcount;l++){
      // Distribute colors to angles
        colors[l]=255*Zangles[l]/maxZangle;
  }
  
  //Start Dynamics
  char *string;
  asprintf(&string,"ffmpeg -y -f rawvideo -vcodec rawvideo -pix_fmt gray -s %dx%d -r 30 -i - -f mp4 -q:v 5 -an -vcodec mpeg4 out.mp4", n1,n2);
  
  //open an output pipe
  FILE *pipeout = popen(string, "w");
  
  clock_t t;
  t = clock();
 
  for (int i =0 ; i<20 ; i++)
  {
    std::cout << "  " << i << "time-step begins...." << std::endl;
    /* add movie frame to file */
    PrepareFFMPEG2DPixels(n1,n2,0,pixels, labels, colors);
    fwrite(pixels, 1, pcount, pipeout);
    
    material.calculateFieldJ('P');
    EtaSubProblem.run(epsilon);
    FastMarching.run(epsilon);
    
  }

  t = clock()-t;
  std::cout << "Total time = " << ((float)t)/CLOCKS_PER_SEC <<std::endl;

  EtaSubProblem.freeMemory();
  FastMarching.freeMemory();
	
  
  fflush(pipeout);
  pclose(pipeout);
 
  delete[] eta; eta=NULL;
  delete[] labels; labels=NULL;
  delete[] Xangles; Xangles=NULL;
  delete[] Yangles; Yangles=NULL;
  delete[] Zangles; Zangles=NULL;
    
}

