//
// KWCJumpFunctoin.h
// myproject
//
// Class KWC Optimizer
// This code construct J(\theta) of extended KWC model from
// given GB data set W(\theta).
// W(\theta) is distributed to
// 1)  f(\eta)
// 2)  \laplace \eta
// 3)   g(\eta)|\grad \theta|
// using Newton's Method
//
// Created by Jaekwang Kim on 3/11/20.
// Copyright © 2020 Jaekwang Kim. All rights reserved.

#ifndef KWCJumpFunction_h
#define KWCJumpFunction_h
#include <iomanip>   //std::setw
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <fstream>      // std::ifstream
#include <vector>


using namespace std;

double computeEtaCuspWithJ(double J,int max_iterator);
double computeTotalEnergy(double eta_cusp, double J);

class KWCDesignJ  {
public:
   
  int nData;
  string inputFileName;
  string outputFileName;
    
  KWCDesignJ(int n, string sInput, string sOutput)  {
    std::cout <<"  KWC Optimizer Class constructed" << std::endl;
    nData=n;
    inputFileName = sInput;
    outputFileName = sOutput;
    
    std::cout <<"  Input data is.... "<< inputFileName <<std::endl;
    std::cout <<"  Output data will be.... "<< outputFileName <<std::endl;
  }
  void run();
    
  private:
};

void KWCDesignJ::run()
{
  double max_iterator=10000;
    
  //input
  double W_data[nData]; // Given data
    
  //output
  double J[nData]; //Objective Function
  double tilt_angle[nData];; // unit [radian]
  double W_kwc[nData];
  double eta_cusp[nData];
    
  //read data
  std::ifstream readFile;
  readFile.open(inputFileName);
  
  std::string varstr;
 
  if (readFile.is_open())
  {
    getline(readFile,varstr); //pass 2 lines
    getline(readFile,varstr);
    for(int i=0; i< nData; i ++)  {
    readFile >> tilt_angle[i] >> W_data[i];
    }
  }else  {
   cout << "file cannot be open" <<endl; exit(1);
  }
    
  readFile.close();
    
  double F0; // error between W_exp and W_kwc at step 0
  double F1; // error at step 1
  
  
  double J_step0; //Guess J
  double J_step1;
  
  double W_kwc_step0;
  double W_kwc_step1;
  
  double eta_cusp_t=0.0;
  
  //Loop over data point
  for(int i=0; i<nData; i++)
  {
        //[1] First two points guess, choices are arbitrary
        J_step0=W_data[i];
        J_step1=0.5 * W_data[i];
        
        eta_cusp_t = computeEtaCuspWithJ(J_step0 , max_iterator);
        W_kwc_step0=computeTotalEnergy(eta_cusp_t,J_step0);
        F0= W_data[i]-W_kwc_step0;
        
        eta_cusp_t = computeEtaCuspWithJ(J_step1 , max_iterator);
        W_kwc_step1=computeTotalEnergy(eta_cusp_t,J_step1);
        F1= W_data[i]-W_kwc_step1;
        
        //[2] Execute Newton iteration
        double small_tol =1e-9; //Newtonain solver tolerance
        
        do
        {
            //update J
            double J_temp = J_step1;
            //Compute next step
            J_step1 = J_step1-(J_step0-J_step1)/(F0-F1) * F1;
            J_step0 = J_temp;
            
            //update W_kwc
            W_kwc_step0 = W_kwc_step1;
            eta_cusp_t = computeEtaCuspWithJ(J_step1 , max_iterator);
            W_kwc_step1=computeTotalEnergy(eta_cusp_t,J_step1);
            
            F0= W_data[i]-W_kwc_step0;
            F1= W_data[i]-W_kwc_step1;
            
        }while(fabs(F1)> small_tol );
        
        //save solutions
        eta_cusp[i]=eta_cusp_t;
        W_kwc[i]=W_kwc_step1;
        J[i]=J_step1;
    
        if(J[i] >= DBL_MIN && J[i] <= DBL_MAX )  { // you do not have a coverence solution
         //do nothing
        }else{
           J[i]=0;
        }
    
        if(W_data[i]==0)
        {
          eta_cusp[i]=1.0;
          W_kwc[i]=0;
          J[i]=0;
        }
    }//End:Loop over data point
    
    
    //Check
  
    for(int i=0; i<nData; i++)
    {
        //std::cout << "( tilt_angle(degree), W_data , W_kwc, J[theta], etacusp )=" <<
        std::cout <<
        "("<< std::setw(5) << tilt_angle[i]<<", "
        << std::setw(10)<< W_data[i] <<", "
        << std::setw(10)<<W_kwc[i]<<", "
        << std::setw(10)<< J[i]<<","
        << std::setw(10)<< eta_cusp[i] <<")"
        <<std::endl;
    }
  
    
    
    std::ofstream pipe;
    pipe.open(outputFileName);
    
    pipe    << std::setw(15)  << " # tiltAngle"
    << std::setw(15)  << "# Jtheta"
    << std::setw(15)  << "# Total Energy"
    << std::endl;
    pipe << std::setw(15)  << "# ------------"
    << std::setw(15)  << " -------------" <<" "
    << std::setw(15)  << " -------------" << std::endl;
  
    for(int i=0; i<nData; i++)
    {
        pipe << std::setw(15) << tilt_angle[i]
        << std::setw(15) << J[i]
        << std::setw(15) << W_data[i]
        << std::endl;
      
    }
    pipe.close(); 
}

//Class declaration ends here

//Simple functions
double computeEtaCuspWithJ(double J,int max_iterator)
{
  //find eta from J(\theta) using jump condition
  // 2*(1-eta_cusp) = s * g,eta (1-eta_cusp) * J
  // It uses a Newton's method to iteratively solve the above equation
  
  double s=1.0;
  double f=FLT_MAX; //Error function: this should be minimized.
  double grad_f=FLT_MAX;
  
  double tol=1e-10;  // Newtonian solver criteria
  
  double eta_t=0.9;   // Initial guess
  double eta_cusp=0.0;
    
  for(int l=0 ; l<max_iterator; l++)  {
    
    f= 2.0*(1.0-eta_t) - s/(1-eta_t) * J;
    grad_f = -2.0 -  s*J/(1-eta_t);
    eta_t = eta_t - f/grad_f;
        
    if(fabs(f)<tol)  {
      eta_cusp=eta_t;
      break;
    }
  }
    
  if(fabs(f) > tol)
    cout << "Warning Eta_cusp value may not be found correctly " <<endl;
  
  return eta_cusp;
}

//KWC Toal Energy
double computeTotalEnergy(double eta_cusp, double J)
{
  //Given J and eta_cusp, compute the KWC total energy (Energy in bicrystal)
  double s=1.0;
  double W_kwc= (1.0-eta_cusp)*(1.0-eta_cusp)-s*log(1-eta_cusp)*J;
  return W_kwc;
}

#endif /* kwc_jumpfunctoin_h */
