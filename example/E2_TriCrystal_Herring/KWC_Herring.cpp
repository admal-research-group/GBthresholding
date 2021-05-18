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
#include "Metrics.h"
#include "InitCrystal.h"
#include "PostProcessor.h"
#include "PrimalDual_Neumann.h"
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
  /* total number of grain */
  const unsigned int lcount=3;
  /* total number of grid point */
  int pcount = n1*n2*n3;
  double dt= epsilon * epsilon; //initial choice of dt
	
  double *Xangles = new double[lcount]();
  double *Yangles = new double[lcount]();
  double *Zangles = new double[lcount]();
  double *eta = new double[pcount]();
  int *labels = new int[pcount]();
  double *JField = new double[pcount]();
  
  
  int Nthread = 1; //Number total thread to be used for FFTW
  
  /* 2D Crystal has only Z angles */
  Xangles[0]=0; Xangles[1]=0; Xangles[2]=0;
  Yangles[0]=0; Yangles[1]=0; Yangles[2]=0;
  Zangles[0]=0.0 * (M_PI)/6 ; Zangles[1]=(M_PI)/6; Zangles[2]=2*(M_PI)/6;
  
  InitializeCrystal::tripleConfiguration(n3,n2,n1,labels);
  
  /* running with original KWC model*/
  char materialType='S';
  Material material(n3,n2,n1,materialType);
  material.simpleKWC_s=1.0; //set material constant value s
  material.setUpClass(Xangles,Yangles,Zangles,labels,JField);
  
  //Construct PD Algorithm class
  double PDerror=1e-6; // tolerance of Primal-dual algorithm
  int PDmaxIters=10000; // allowable iteration number of the algorithm
  PrimalDual<DIM> EtaSubProblem(n3, n2, n1, PDerror, PDmaxIters,lcount,epsilon, Nthread);
  
  //Criteria for identifying interior regions of grain
  double initThresCriteria = 0.9;
  
  KWCThreshold<DIM> FastMarching(n3,n2,n1,lcount,initThresCriteria);

  EtaSubProblem.setUpClass(eta, Xangles, Yangles, Zangles, labels, JField);
  /* 'D' stands for Dirichlet boundary condition will be enfored during thresholding */
  FastMarching.setUpClass(eta, Xangles, Yangles, Zangles, labels, JField,'D');
    
  
  /* FFMPEG movie data */
  unsigned char *pixels=new unsigned char[pcount];
  unsigned char *colors=new unsigned char[lcount];
  for(int l=0;l<lcount;l++){
      // Distribute colors to angles
        colors[l]=255*Zangles[l]/(M_PI) * 3;
  }
  
  //Start Dynamics
  char *string;
  asprintf(&string,"ffmpeg -y -f rawvideo -vcodec rawvideo -pix_fmt gray -s %dx%d -r 30 -i - -f mp4 -q:v 5 -an -vcodec mpeg4 out.mp4", n1,n2);
  //open an output pipe
  FILE *pipeout = popen(string, "w");
  
  //Start Dynamics
  for (int i =0 ; i<200 ; i++)
  {
    PrepareFFMPEG2DPixels(n1,n2,0,pixels, labels, colors);
    fwrite(pixels, 1, pcount, pipeout);
    
    material.calculateFieldJ('N');
    EtaSubProblem.run(epsilon);
    FastMarching.run(epsilon);
  }

  fflush(pipeout);
  pclose(pipeout);

  EtaSubProblem.freeMemory();
	FastMarching.freeMemory();
  
  delete[] eta; eta=NULL;
  delete[] labels; labels=NULL;
  delete[] Xangles; Xangles=NULL;
  delete[] Yangles; Yangles=NULL;
  delete[] Zangles; Zangles=NULL;
    
}

