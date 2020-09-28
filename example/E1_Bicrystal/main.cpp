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
#include "PrimalDual_Neumann.h"


using namespace std;
using namespace DataOut;

/*
 
 Input : Two orientation angles of two crystal
 Output : crystal order phase (eta) solution and KWC Energy of the bicrystal

 This code uses a Primal-dual algorithm
 to solve $\eta(x)$ for a given $\theta(x)$.
 
 In case of Bicrystal, both the analytic solution of $\eta$
 and KWC energy are known.
 Therefore, we can use this configuration to veryfy the Primal-dual method
 
*/

int main(int argc, char *argv[]){

  //Read global variables from bash
  //Define grid size and set model parameter epsilon
    
  int n1=atoi(argv[1]);
  int n2=atoi(argv[2]);
  int n3=atoi(argv[3]);
  double epsilon=atof(argv[4]);
 
	/* Construct Classes */
  int const DIM=3;
  /* total number of grain */
  const unsigned int lcount=2;
  /* total number of grid point */
	int pcount = n1*n2*n3; 
	double dt= epsilon * epsilon; //initial choice of dt 
	  
  double *Xangles = new double[lcount]();
  double *Yangles = new double[lcount]();
  double *Zangles = new double[lcount]();
  double *eta = new double[pcount]();
  int *labels = new int[pcount]();
  double *JField = new double[pcount]();
  
  int Nthread = 1; //The number of threads to be used for FFTW
  
  //Input
  double orientationAngleLeftCrystal = 0.0;
  double orientationAngleRightCrystal = 30.0* M_PI/180;
  
  //Initialize Polycrystal angles
  Zangles[0]=orientationAngleLeftCrystal;
  Zangles[1]=orientationAngleRightCrystal;

  //Initialize crystal: set initial condition of $\theta(x)$
  InitializeCrystal::oneD_Bicrystal_configuration(n3,n2,n1,labels);
  
  char materialType='S';
  Material material(n3,n2,n1,materialType);
  material.simpleKWC_s=1.0; //set material constant value s
  material.setUpClass(Xangles,Yangles,Zangles,labels,JField);
  
  //Run Primal-dual algorithm
  material.calculateFieldJ('N');
 
  double PDerror=1e-6; // tolerance of Primal-dual algorithm
  int PDmaxIters=10000; // allowable iteration number of the algorithm
  //Construct PD Algorithm class
  PrimalDual<DIM> EtaSubProblem(n3, n2, n1, PDerror, PDmaxIters,lcount,epsilon, Nthread);
  
  //Link vaiables' pointers to the PD algortihm class
  EtaSubProblem.setUpClass(eta, Xangles, Yangles, Zangles, labels, JField);
  EtaSubProblem.run(epsilon);
   
  //Export $\eta$ solution in the vtk format
  //DataOut::Output2DvtuScalar(n1, n2, n3, eta, "eta", "etaSolution_", 0);

  double energy = computeKWCEnergy(material.simpleKWC_s,n1,n2,epsilon,eta,Zangles,labels);
  
  std::cout <<"  GB Energy:  " << energy << std::endl;
    
  EtaSubProblem.freeMemory();
	
  delete[] eta; eta=NULL;
  delete[] labels; labels=NULL;
  delete[] Xangles; Xangles=NULL;
  delete[] Yangles; Yangles=NULL;
  delete[] Zangles; Zangles=NULL;
}

