// Primal_dual.h
// Created by Jaekwang Kim on 11/10/19.
// Header files include EtaSubproblem solvers for KWC Grainboundary Models simulation

/*
It  reads \theta (x) and determines minimum solutin of phi
	Input :: n1,n2,n3, labels 
	Xangles(x) Yangles(x) Zangles(x), Label(x) 

	At time t=0,  Xangles(x) and Label(x)
	can be generated from functions such as init_1d_bicrystal
 
Output :: eta(x)
*/

#ifndef PrimalDual_h
#define PrimalDual_h

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

//Poisson equation FFT solver struct
typedef struct  {
  fftw_plan Forward;
  fftw_plan Backward;
  double *kernel;

  fftw_complex *workspace;
  
}poissonSolver;


template<int DIM>
class PrimalDual{
    
  public:
    
  PrimalDual  (const unsigned int n3,
    const unsigned int n2,
    const unsigned int n1,
    double error,
    int maxIters,
    const unsigned int lcount,
    double epsilon,
    int Nthread
  ); //Consructor
    
  poissonSolver fftps; /* The"workspace"-stucture */
    
  void setUpClass  (double *eta_pointer,
    double *Xangle_pointer,
    double *Yangle_pointer,
    double *Zangle_pointer,
    int *label_pointer,
    double *JField_pointer
    );
    
  void run(double epsilon); // Class Driver
  void freeMemory ();
		
  //The Primal-dual algorithms begins with taking these variables from external program
  unsigned int n1,n2,n3;
    
  int *labels;
  double *Xangles;
  double *Yangles;
  double *Zangles;
  double *eta;
  double *JField;
  double epsilon;
  
  private:
    
  int Nthread;  // denotes the number of thread being used for FFT and iFFT
  double error;
  int maxIters;
  int pcount;
  int lcount;
  double *psi;
  double *etaOld;
  
  void calculateEnergySmooth (double epsilon);
  double solveEta (double *eta, double *psi, double tau);
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
:n3(n3), n2(n2), n1(n1),error(error), maxIters(maxIters),lcount(lcount),epsilon(epsilon),
 Nthread(Nthread)
{
  std::cout << std::endl;
  std::cout << "  [Constructor] C++ Primal Dual Algorithm Class is being created..." << std::endl;
  std::cout << "  n1:   " << n1 <<std::endl;
  std::cout << "  n2:   " << n2 <<std::endl;
  std::cout << "  n3:   " << n3 <<std::endl;
  std::cout << "  epsilon:   " << epsilon <<std::endl;
  std::cout << "  initial dt:   " << epsilon * epsilon <<std::endl;
  std::cout << "  Target Error:   " << error <<std::endl;
  std::cout << "  MaxIters:   " << maxIters <<std::endl;
  std::cout << "  lcount:   " << lcount <<std::endl;
  std::cout << "  Nthread(FFTW):   " << Nthread <<std::endl;
  std::cout << "  [Constructor] Class construction ended!" << std::endl;
}


template<int DIM>
void PrimalDual<DIM>::setUpClass(double *eta_pointer,
  double *Xangle_pointer,
  double *Yangle_pointer,
  double *Zangle_pointer,
  int *label_pointer,
  double *JField_pointer
  )
{
  // [Goal] memory spaces will be prepared
  // for variables used throughout the program
  // ::Global variable within PrimalDual Class
    
  std::cout << std::endl;
  std::cout <<"   Setting PD Algorithm workspace..." << std::endl;

  pcount = n1 * n2 * n3;
	
  //Link pointers
  eta= eta_pointer;
  Xangles= Xangle_pointer;
  Yangles= Yangle_pointer;
  Zangles= Zangle_pointer;
  labels= label_pointer;
  JField= JField_pointer;
	
  psi = new double[pcount]();
  etaOld = new double[pcount]();
  fftps.kernel = new double[pcount]();

  int fftw_init_threads(void);
  fftw_plan_with_nthreads (Nthread);
  fftps.workspace= new fftw_complex[pcount]();

  fftps.Forward=fftw_plan_dft_3d(n3,n2,n1,fftps.workspace,fftps.workspace,
	                                       -1,FFTW_MEASURE);

  fftps.Backward=fftw_plan_dft_3d(n3,n2,n1,fftps.workspace,fftps.workspace,
	                                        1, FFTW_MEASURE);
	
  for(int i=0; i<n3; i++)  { // z loop
    for(int j=0; j<n2; j++) { //y loop
      for(int k=0; k<n1; k++) { //x loop
        
         //They are negative Laplacian operator...
          double x = 2*M_PI*k/(n1*1.0);
          double y = 2*M_PI*j/(n2*1.0);
          double z = 2*M_PI*i/(n3*1.0);
          double Eigenvalues=2*n1*n1*(1-cos(x))+2*n2*n2*(1-cos(y))
                                   +2*n3*n3*(1-cos(z));

                fftps.kernel[i*n2*n1+ j*n1+k]=Eigenvalues;
      }
    }
  }

}



template<int DIM>
void PrimalDual<DIM>::calculateEnergySmooth(double epsilon)
{
  double dt = epsilon * epsilon;
  double t=pow(dt,2.0);
    
  for(int i=0;i<pcount;i++){
    fftps.workspace[i][0]=JField[i];
    fftps.workspace[i][1]=0.0;
  }
    
  fftw_execute(fftps.Forward);
    
  for(int i=0;i<pcount;i++){
    fftps.workspace[i][0]*=exp(-fftps.kernel[i]*t)/pcount;
    fftps.workspace[i][1]*=exp(-fftps.kernel[i]*t)/pcount;
  }
    
  fftw_execute(fftps.Backward);
    
  for(int i=0;i<pcount;i++){
    	JField[i]=fftps.workspace[i][0];
      if(JField[i]<0){
        JField[i]=0.0;}
    
  }
	
}


template<int DIM>
void PrimalDual<DIM>::run(double epsilon)
{
  cout << "  Execute Primal-dual algorithm for eta sub problem."<<endl;
  double tau=epsilon;
  double sigma= 1/tau;
  double eta_change;
  double psi_change;

  calculateEnergySmooth (epsilon);
    
  // Control Speed of Iteration
  double strongConvexity=2/epsilon;
  double alpha=1/sqrt(1+strongConvexity*tau);
	
	//eta initialization 
  for(int i=0;i< pcount ;i++){
    eta[i]=0.0; psi[i]=0.0;
  }
    
  for(int i=0;i<  maxIters ;i++){
        
    memcpy(etaOld,eta,pcount*sizeof(double));
        
    eta_change= solveEta(eta, psi, tau);
    psi_change= solvePsi(eta, etaOld, psi, sigma, alpha);
        
    alpha=1/sqrt(1+strongConvexity*tau);
    tau*=alpha;
    sigma/=alpha;
     
    if(i%20==0)
      cout<< "Iteration :" << i << " step.... " <<
       "  (Eta_change, Psi_change)=(" <<eta_change <<","<<psi_change<<")"<<endl;
		
        
    if(eta_change<error)  {
      cout << "  PD Algorithm Tolerance is achieved\n"<<endl;
            break;
    }
  }
}


template<int DIM>
double PrimalDual<DIM>::solveEta (double *eta, double *psi, double tau)
{
  // We take form for g(x) = -log(1-x), f(x)= (1/epsilon) * (1-x)^2
  double maxEtaChange = 0;
    
  for(int i=0; i<n3; i++)  {
    for(int j=0;j<n2;j++)  {
      for(int k=0;k<n1;k++){
                
        //evaluate nabla (\theta)
        //The front factor (n1*n2*n3)^(1/3) will be valid only when "n1=n2=n3"
        double norm = n1 * JField[i*n1*n2+j*n1+k];
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

  double *second_derivative_values;
  second_derivative_values= (double*) calloc(3,sizeof(double));
    
  for(int i=0; i<n3; i++)  {
    for(int j=0;j<n2;j++)  {
      for(int k=0;k<n1;k++)  {

        //Neumann Condition
        //double etaBar = eta[i*n2*n1+j*n1+k]+ alpha * (eta[i*n2*n1+j*n1+k]-etaOld[i*n2*n1+j*n1+k]);
        //fftps.workspace[i*n2*n1+j*n1+k]=epsilon*( psi[i*n2*n1+j*n1+k]/sigma + etaBar);

        //Periodic condition
        calc_second_derivative (psi,i,j,k,n3,n2,n1,second_derivative_values);
        double psi_xx=second_derivative_values[0];
        double psi_yy=second_derivative_values[1];
        double psi_zz=second_derivative_values[2];

        calc_second_derivative (eta,i,j,k,n3,n2,n1,second_derivative_values);

        double eta_xx=second_derivative_values[0];
        double eta_yy=second_derivative_values[1];
        double eta_zz=second_derivative_values[2];

        calc_second_derivative (etaOld,i,j,k,n3,n2,n1,second_derivative_values);
        double etaOld_xx=second_derivative_values[0];
        double etaOld_yy=second_derivative_values[1];
        double etaOld_zz=second_derivative_values[2];

        // etaxx bar is linear interpolation between etaxx_old and etaxx
        double etaxxBar=(1+alpha)*eta_xx-alpha*etaOld_xx;
        double etayyBar=(1+alpha)*eta_yy-alpha*etaOld_yy;
        double etazzBar=(1+alpha)*eta_zz-alpha*etaOld_zz;

        fftps.workspace[i*n2*n1+j*n1+k][0]=-epsilon*
            (
             psi_xx+psi_yy+psi_zz +
             sigma*(etaxxBar+etayyBar+etazzBar)
            ) /sigma;
        
        fftps.workspace[i*n2*n1+j*n1+k][1]=0;
                
      }
    }
  }
    
  fftw_execute(fftps.Forward);
  
  //zeroing
  fftps.workspace[0][0]=0.0;
  fftps.workspace[0][1]=0.0;
    
  for(int i=0;i<pcount;i++){
    fftps.workspace[i][0]/=(1+epsilon*fftps.kernel[i]/sigma)*pcount;
    fftps.workspace[i][1]/=(1+epsilon*fftps.kernel[i]/sigma)*pcount;
   }

  fftw_execute(fftps.Backward);
    
  double maxPsiChange=0;
    
  for(int i=0;i<pcount;i++){
    double old=psi[i];
    psi[i]=fftps.workspace[i][0];
    maxPsiChange=fmax(fabs(psi[i]-old),maxPsiChange);
  }
    
  return maxPsiChange;
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

  std::cout<<"  Memories of FFTW has freed successfully."<< std::endl;
}



#endif /* PrimalDual_h */
