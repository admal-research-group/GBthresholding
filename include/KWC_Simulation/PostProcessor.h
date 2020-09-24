// post_processor.h
// Created by Jaekwang Kim on 11/10/19.
// Header files include PostProcessors of KWC model simulation
// It read KWC (phase, orientation) data from outer .txt file and
// does postprocessing (e.g. computing KWC energy)

/*
This code should be prepared to run on the fly while the main code is running
or separately from result.txt file
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <math.h>


#include "DataOut.h"
#include "Material.h"

#ifndef Post_processor_h
#define Post_processor_h

using namespace std;
using namespace DataOut;


//2d energy first, from a pointer
double computeKWCEnergy(double sParam, int n1,int n2, double epsilon,
                        double *eta, double *angles, int *labels)
{
  cout<< "  Computes KWC Energy...\n";
  int pcount = n1*n2;
  
  double dA=1.0/n1/n1;
  double dl=1.0/n1;
  
  double *energyField;
  energyField = (double*) calloc(pcount,sizeof(double));
 
  double jumpEng=0.0;
  double gradientEng=0.0;
  
  int k=0;
  
  for(int i=0;i<n2;i++){
    for(int j=0;j<n1;j++){
     
      double x_val =  ((j*1.0)/n1); double y_val =  ((i*1.0)/n2);
      int xm=j-1; int ym=i-1; int xp=j+1; int yp=i+1;
      
      if(xm<0){
        xm=0;
      }
      
      if(xp>n1-1){
        xp=n1-1;
      }
      
      if(yp>n2-1){
        yp=n2-1;
      }
      
      if(ym<0)
      {
        ym=0;
      }
      
      //Gradient of eta //second order
      double etax=n1*(eta[i*n1+xp]-eta[i*n1+xm])/2;
      double etay=n2*(eta[yp*n1+j]-eta[ym*n1+j])/2;
    
      //Ideally, Jump of $\theta$, Exists as a delta function
      //Finite Difference Central Scheme diffuses it to a half magnitude on 2 grid point
      double localJump;
      
      localJump = fabs(angles[labels[i*n1+xp]]-angles[labels[i*n1+xm]]);
      double Jumpx=fabs(localJump);
      
      localJump = fabs(angles[labels[yp*n1+j]]-angles[labels[ym*n1+j]]);
      double Jumpy=fabs(localJump);
      
      double meanJump = sqrt(Jumpx*Jumpx+Jumpy*Jumpy);
      
      // [a] (epsilon/2) * (\nabla eta)^2
      double engGradEta =0.5 * epsilon *( etax*etax + etay*etay ) * dA;
      // [b] ( 1/(2epsilon) ) * f(\eta)
      double engHmgEta = (1.0/(2*epsilon)) * (1.0-eta[i*n1+j])*(1.0-eta[i*n1+j]) * dA;
      // [c] coupling term // This does not exists in the integrand (integrated quantitiy)
      double gIntegrated = (-1.0) * log(1.00-eta[i*n1+j]); // g(eta)= -log(1-eta).
     
      if(gIntegrated < 0) //melted
      {
        gIntegrated = 0.0;
        engHmgEta = (1.0/(2*epsilon)) * dA;
      }
     
      double engJmpTheta;
      double tol = pow(10,-10);
      
      // as g(\eta) diverges to infty when eta=1, we need some trick.
      
      if(meanJump>tol){
         engJmpTheta = 0.5* sParam * gIntegrated * meanJump * dl;}
      else  {
         engJmpTheta = 0;
      }
      
       energyField[k]= engGradEta + engHmgEta + engJmpTheta;
       jumpEng += engJmpTheta;
       gradientEng += engGradEta+engHmgEta;
       k=k+1;
      
        }
    }
  
    double totalEnergy=0;
  
    for(int i=0; i<pcount ; i++)
    {
        totalEnergy+=energyField[i];
    }
  
    cout << "   Total energy: " << totalEnergy << endl;
    cout << "   Jump eng " << jumpEng << endl;
    cout << "   Gradient eng " << gradientEng << endl;
  
    return totalEnergy;
}


//2d energy first, from a pointer
double computeKWCEnergyCovarianceModel(Material &material, int n1,int n2, double epsilon,
                        double *eta, double *angles, int *labels)
{
  cout<< "  Computes KWC Energy...\n";
  int pcount = n1*n2;
  
  double dA=1.0/n1/n1;
  double dl=1.0/n1;
  
  double *energyField;
  energyField = (double*) calloc(pcount,sizeof(double));
 
  double jumpEng=0.0;
  double gradientEng=0.0;
  
  int k=0;
  
  for(int i=0;i<n2;i++){
    for(int j=0;j<n1;j++){
     
      double x_val =  ((j*1.0)/n1); double y_val =  ((i*1.0)/n2);
      int xm=j-1; int ym=i-1; int xp=j+1; int yp=i+1;
      
      if(xm<0){
        xm=0;
      }
      
      if(xp>n1-1){
        xp=n1-1;
      }
      
      if(yp>n2-1){
        yp=n2-1;
      }
      
      if(ym<0)
      {
        ym=0;
      }
      
      //Gradient of eta
      double etax=n1*(eta[i*n1+xp]-eta[i*n1+xm])/2;
      double etay=n2*(eta[yp*n1+j]-eta[ym*n1+j])/2;
    
      //Ideally, Jump of $\theta$, Exists as a delta function
      //Finite Difference Central Scheme diffuses it to a half magnitude on 2 grid point
      double localJump;
      
      localJump = fabs(angles[labels[i*n1+xp]]-angles[labels[i*n1+xm]]);
      double Jumpx=localJump;
      
      localJump = fabs(angles[labels[yp*n1+j]]-angles[labels[ym*n1+j]]);
      double Jumpy=localJump;
      
      double meanJump = sqrt(Jumpx*Jumpx+Jumpy*Jumpy);
      
      double Jtheta = material.computeJthetaCovarianceModel(meanJump);
     
      
      // [a] (epsilon/2) * (\nabla eta)^2
      double engGradEta =(epsilon/2.0) *( etax*etax + etay*etay ) * dA;
      // [b] ( 1/(2epsilon) ) * f(\eta)
      double engHmgEta = (1.0/(2*epsilon)) * (1.0-eta[i*n1+j])*(1.0-eta[i*n1+j]) * dA;
      // [c] coupling term // This does not exists in the integrand (integrated quantitiy)
      double gIntegrated = (-1.0) * log(1.02-eta[i*n1+j]) * Jtheta; // g(eta)= -log(1-eta).
     
      double tol = pow(10,-10);
      double engJmpTheta;
      // as g(\eta) diverges to infty when eta=1, we need some trick.
      if(meanJump>tol){
         engJmpTheta = 0.5 * gIntegrated * dl;}
      else  {
         engJmpTheta = 0.5 * gIntegrated * dl;;
      }
      
       energyField[k]= engGradEta + engHmgEta + engJmpTheta;
       jumpEng += engJmpTheta;
       gradientEng += engGradEta+engHmgEta;
       k=k+1;
      
        }
    }
  
    double totalEnergy=0;
  
    for(int i=0; i<pcount ; i++)
    {
        totalEnergy+=energyField[i];
    }
  
    cout << "   Total energy: " << totalEnergy << endl;
    cout << "   Jump eng " << jumpEng << endl;
    cout << "   Gradient eng " << gradientEng << endl;
  
    return totalEnergy;
}


#endif /* post_processor_h */
