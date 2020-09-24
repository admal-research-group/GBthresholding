// algo_primal_dual.h
// Created by Jaekwang Kim on 11/10/19.
// Header files include EtaSubproblem solvers for KWC Grainboundary Models simulation

/*
It  reads \theta (x) and determines minimum solutin of phi
	Input :: n1,n2,n3, labels 
	Xangles(x) Yangles(x) Zangles(x), Label(x) 

	At time t=0,  Xangles(x) and Label(x)
	can be generated from functions such as init_1d_bicrystal
 
Output :: eta(x)

(Note) :: This code applies Neumann boundary condition on $\eta$.
To do this, it uses a discrete cosine fourier transform
*/
#ifndef algo_primal_dual_Neumann_h
#define algo_primal_dual_Neumann_h

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <math.h>
#include <fftw3.h>

#include "DataOut.h"
#include "Metrics.h"
#include "Material.h"

using namespace std;

// [Goal] Solves eta_sub problem
// :Minimizing \eta KWC energy for given discrete \theta configuration.
// using Primal Dual algorithm
// Dualize the minimizing problem in terms of eta_star, psi_star


//Poisson equation Real DFT solver struct
typedef struct{
    fftw_plan Forward;
    fftw_plan Backward;
    //fftw_complex *workspace;
    double *workspace;
    double *kernel;
}poissonSolver;

template<int DIM>
class PrimalDual
{
    
  public:
  
  PrimalDual(const unsigned int n3,
    const unsigned int n2,
    const unsigned int n1,
    double error,
    int maxIters,
    const unsigned int lcount,
    double epsilon,
    int Nthread
   ); //Consructor
    
  poissonSolver fftps; /* The"workspace"-stucture */
    
  void setUpClass(double *eta_pointer,
    double *Xangle_pointer,
    double *Yangle_pointer,
    double *Zangle_pointer,
    int *label_pointer,
    double *energyField_pointer,
    char materialType);
    
  void run(Material &material, double epsilon); // Code driver
  void freeMemory ();
		
  //The Primal-dual algorithms begins with taking these variables from external program
  unsigned int n1,n2,n3;

  int *labels;
  double *Xangles;
  double *Yangles;
  double *Zangles;
  double *eta;
  double *energyField;
  double epsilon;
	
  char materialType;
  
  private:
    
  int Nthread; // denotes the number of thread being used for FFT and iFFT
  double error;
  int maxIters;
  int pcount;
  int lcount;
  double *psi;
  double *etaOld;

  
  void calculateEnergy (Material &material);
  void calculateEnergySmooth (double epsilon);
  double solveEta(double *eta, double *psi, double tau);
  double solvePsi (double *eta, double *etaOld,double *psi, double sigma, double alpha);
        

};

//Constructor
template<int DIM>
PrimalDual<DIM>::PrimalDual(
  const unsigned int n3,
  const unsigned int n2,
  const unsigned int n1,
  double error,
  int maxIters,
  const unsigned int lcount,
  double epsilon,
  int Nthread)
:n3(n3), n2(n2), n1(n1), error(error), maxIters(maxIters),lcount(lcount),epsilon(epsilon),
Nthread(Nthread)
{
  std::cout << std::endl;
  std::cout << " [Constructor] C++ Primal Dual Algorithm Class is being created..." << std::endl;
  std::cout << "  n1:   " << n1 <<std::endl;
  std::cout << "  n2:   " << n2 <<std::endl;
  std::cout << "  n3:   " << n3 <<std::endl;
	std::cout << "  epsilon:   " << epsilon <<std::endl;
  std::cout << "  initial dt:   " << epsilon *epsilon  <<std::endl;
  std::cout << "  Target Error:   " << error <<std::endl;
  std::cout << "  MaxIters:   " << maxIters <<std::endl;
	std::cout << "  lcount:   " << lcount <<std::endl;
  std::cout << "  Nthread(FFTW):   " << Nthread <<std::endl;
  std::cout << " [Constructor] Class construction ended!" << std::endl;
}


template<int DIM>
void PrimalDual<DIM>::setUpClass(double *eta_pointer,
						   double *Xangle_pointer,
						   double *Yangle_pointer,
						   double *Zangle_pointer,
						   int *label_pointer,
						   double *energyField_pointer,
               char materialType_input)
{
  // Pointer variable memory spaces will be linked from outer pointers
  // ::Global variable within PrimalDual Class
    
  std::cout << std::endl;
  std::cout <<"   Setting PD Algorithm workspace..." << std::endl;
  pcount = n1 * n2 * n3;
	
  materialType = materialType_input;
  
	//link pointers
	eta= eta_pointer;
	Xangles= Xangle_pointer; 
	Yangles= Yangle_pointer; 
	Zangles= Zangle_pointer; 
	labels= label_pointer; 
	energyField= energyField_pointer;
	
	psi = new double[pcount];
	etaOld = new double[pcount];
	fftps.workspace= new double[pcount];
	fftps.kernel = new double[pcount];

	int fftw_init_threads(void); 
	fftw_plan_with_nthreads (Nthread);

  fftps.Forward=fftw_plan_r2r_2d(n2, n1, fftps.workspace, fftps.workspace,
                                 FFTW_REDFT10, FFTW_REDFT10,
                                 FFTW_MEASURE);
  
  fftps.Backward=fftw_plan_r2r_2d(n2, n1, fftps.workspace, fftps.workspace,
                                  FFTW_REDFT01, FFTW_REDFT01,
                                  FFTW_MEASURE);

  for(int i=0; i<n3; i++)  { // z loop
    for(int j=0; j<n2; j++)  { //y loop
      for(int k=0; k<n1; k++)  { //x loop
        
        double x = M_PI*k/(n1*1.0);
        double y = M_PI*j/(n2*1.0);
        double z = M_PI*i/(n3*1.0);
        double Eigenvalues=2*n1*n1*(1-cos(x))+2*n2*n2*(1-cos(y))
                           +2*n3*n3*(1-cos(z));

        fftps.kernel[i*n2*n1+j*n1+k]=Eigenvalues;
      }
     }
  }

}

	

template<int DIM>
void PrimalDual<DIM>::calculateEnergy(Material &material)
{
  for(int i=0; i<n3; i++)  {
    for(int j=0;j<n2;j++)  {
      for(int k=0;k<n1;k++){
                
        int xp=k+1; int xm=k-1;
        int yp=j+1; int ym=j-1;
        int zp=i+1; int zm=i-1;
                
        if(xp>n1-1)
          xp=n1-1;
                
        if(yp>n2-1)
          yp=n2-1;
                
        if(zp>n3-1)
          zp=n3-1;
                
        if(xm<0)
          xm=0;
                
        if(ym<0)
          ym=0;
                
        if(zm<0)
          zm=0;
                
        double local_jump_x; double local_jump_y; double local_jump_z;

        local_jump_x= 1.0 * fabs(Zangles[labels[i*n1*n2+j*n1+xp]]-Zangles[labels[i*n1*n2+j*n1+xm]]);
        local_jump_y= 1.0 * fabs(Zangles[labels[i*n1*n2+yp*n1+k]]-Zangles[labels[i*n1*n2+ym*n1+k]]);
        local_jump_z= 1.0 * fabs(Zangles[labels[zp*n1*n2+j*n1+k]]-Zangles[labels[zm*n1*n2+j*n1+k]]);

        double Misorientation =sqrt(local_jump_x * local_jump_x +
                                    local_jump_y * local_jump_y +
                                    local_jump_z * local_jump_z);
        
        double misorientationENG=0.0;
        
        //Unavoidably, the local jump is diffused by central Finite difference into two grid point
        //It must have a same energy in the sense of weak form, so we take half of Jump energy
        //on each grid point
        
        if(materialType=='c')  {
         misorientationENG = 0.5 * material.computeJthetaCovarianceModel(Misorientation);
        }else  {
        misorientationENG = 0.5* material.computeJthetaSimpleMaterial(Misorientation);
        }
        
        energyField[i*n1*n2+j*n1+k]= misorientationENG;
				
      }
    }
  }
}

template<int DIM>
void PrimalDual<DIM>::calculateEnergySmooth(double epsilon)
{
  double dt = epsilon * epsilon;
  //calculate diffused \(theta) for small times from discrete values of theta
  double e=pow(epsilon,4.0);
    
  for(int i=0;i<pcount;i++)  {
    fftps.workspace[i]=energyField[i];
  }
    
  fftw_execute(fftps.Forward);
    
  for(int i=0;i<pcount;i++)  {
    fftps.workspace[i]*=exp(-fftps.kernel[i]*e)/(4*pcount);
  }
    
  fftw_execute(fftps.Backward);
    
  for(int i=0;i<pcount;i++){
    energyField[i]=fftps.workspace[i];
  }

}


template<int DIM>
void PrimalDual<DIM>::run(Material &material, double epsilon)
{
  cout << "  Execute Primal-dual algorithm for eta sub problem."<<endl;

  double tau=epsilon;
  double sigma= 1/tau;
  double etaChange;
  double psiChange;

  calculateEnergy(material);
  calculateEnergySmooth (epsilon);
    
  // Control Iteration
  double strongConvexity=2/epsilon;
  double alpha=1/sqrt(1+strongConvexity*tau);
	
	//eta initialization 
	for(int i=0;i< pcount ;i++){
		eta[i]=0.0; 
		psi[i]=0.0;
	}
    
  for(int i=0;i<  maxIters ;i++)  {
    
    memcpy(etaOld,eta,pcount*sizeof(double));
    etaChange= solveEta(eta, psi, tau);
    psiChange= solvePsi(eta, etaOld, psi, sigma, alpha);
        
    alpha=1/sqrt(1+strongConvexity*tau);
    tau*=alpha;
    sigma/=alpha;
     
    if(i%30==0)  {
      cout<< "  (etaChange, psiChange)=(" <<etaChange <<","<<psiChange<<")"<<endl;
		}
        
    if(etaChange<error)  {
      cout << "    PD Algorithm Tolerance is achieved\n"<<endl;
      break;
    }
    
  }
}


template<int DIM>
double PrimalDual<DIM>::solveEta(double *eta, double *psi, double tau)
{
  // We take form for g(x) = -log(1-x), f(x)= (0.5/epsilon) * (1-x)^2
  double maxEtaChange = 0;
    
  for(int i=0;i<n3;i++)  {
    for(int j=0;j<n2;j++)  {
      for(int k=0;k<n1;k++)  {
    
        double norm = n1 * energyField[i*n1*n2+j*n1+k];
        double qa= 1.0/epsilon+1/tau;
        double qb= psi[i*n1*n2+j*n1+k]-2.0/epsilon-1/tau-eta[i*n1*n2+j*n1+k]/tau;
        double qc= 1.0/epsilon-psi[i*n1*n2+j*n1+k] + eta[i*n1*n2+j*n1+k]/tau - norm;
                
        double discrim= qb*qb-4*qa*qc;
        double qz= sqrt(discrim);
                
        double old= eta[i*n2*n1+j*n1+k];
                
        eta[i*n1*n2+j*n1+k] =-(qb+qz)/(2*qa);
        maxEtaChange = fmax( fabs(eta[i*n1*n2+j*n1+k]-old), maxEtaChange);
                
      }
    }
  }
    return maxEtaChange;
}


template<int DIM>
double PrimalDual<DIM>::solvePsi (double *eta, double *etaOld, double *psi, double sigma,
                                  double alpha)
{
	//double *second_derivative_values;
  //second_derivative_values= (double*) calloc(3,sizeof(double));
    
  for(int i=0;i<n3;i++)  {
    for(int j=0;j<n2;j++)  {
      for(int k=0;k<n1;k++)  {
        double etaBar = eta[i*n2*n1+j*n1+k]+ alpha * (eta[i*n2*n1+j*n1+k]-etaOld[i*n2*n1+j*n1+k]);
        fftps.workspace[i*n2*n1+j*n1+k]=epsilon*( psi[i*n2*n1+j*n1+k]/sigma + etaBar);
      }
    }
  }
    
  fftw_execute(fftps.Forward);
  fftps.workspace[0]=0.0;  //zeroing

  for(int i=0;i<pcount;i++){
    fftps.workspace[i]*= fftps.kernel[i]/((1+epsilon*fftps.kernel[i]/sigma)*4*pcount);
  }
    
  fftw_execute(fftps.Backward);
    
  double max_psiChange=0;
    
  for(int i=0;i<pcount;i++){
    double old=psi[i];
    psi[i]=fftps.workspace[i];
    max_psiChange=fmax(fabs(psi[i]-old),max_psiChange);
  }
    
  return max_psiChange;
}



template<int DIM>
void PrimalDual<DIM>::freeMemory()
{
  delete[] psi; psi=NULL;
  delete[] fftps.workspace; fftps.workspace=NULL;
  delete[] fftps.kernel; fftps.kernel=NULL;
  delete[] etaOld; etaOld=NULL;
  
  fftw_cleanup_threads();
  fftw_destroy_plan(fftps.Forward);
  fftw_destroy_plan(fftps.Backward);

  std::cout<<"  Memories of FFTW has freed successfully"<< "\n"<< std::endl;
}




#endif /* algo__primal_dual_h */
