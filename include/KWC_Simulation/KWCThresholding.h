// algo_primal_dual.h
// Created by Jaekwang Kim on 11/10/19.
// Header file includes KWC Thresholding dynamics using Fast Marching Method 

// Being updated:: Outsourcing thresholidng value

/*
It  reads \eta (x) and the grain boundary 
	Input :: n1,n2,n3, 
	Xangles(x) Yangles(x) Zangles(x), Label(x) 

	At time t=0,  Xangles(x) and Label(x)
	can be generated from functions such as init_1d_bicrystal
 
Output :: eta(x)
*/


#ifndef algo_KWCThreshold_h
#define algo_KWCThreshold_h
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <math.h>
#include <fftw3.h>
#include "DataOut.h"
#include "Heap.h"
#include "watch_fastmarching.h"

using namespace std;
using namespace DataOut;


template<int dim>
class KWCThreshold
{
public:

  KWCThreshold(const unsigned int n3,
               const unsigned int n2,
               const unsigned int n1,
               const unsigned int lcount,
               double initThresCriteria);
    	
  double *eta;
  double *Xangles;
  double *Yangles;
  double *Zangles;
  int *labels;
  double *energy;
  
  
  void run (double epsilon);
  void setUpClass  (
         double *eta_pointer,
         double *Xangle_pointer,
         double *Yangle_pointer,
         double *Zangle_pointer,
         int *label_pointer,
         double *Energy_pointer,
         char boundary);
						   
  void freeMemory ();
  double thresCriteria;
  int watch_fastmarching=0;
  bool STUCK;

private:
  
  unsigned int n1,n2,n3;
  int pcount;
  int lcount;
    
  double *grid; //temporary grid where Fast Marching Method proceeds.
  double *indicator;
  i_heap heap;
		
  vector<vector<int> > threeDNeighbors; //26 neighbor points in 3d
  char thetaBoundary;
  
  void identifyInitBoundary (double epsilon);
  bool thresholding ();
  int executeFastMarchingIteration();
  double eikonalUpdate (int k, int j, int i);
  double eikonalUpdate_2d (double h1, double h2, double U1, double U2, double f);
  	
};

//Constructor
template<int dim>
KWCThreshold<dim>::KWCThreshold(const unsigned int n3,
 const unsigned int n2,
 const unsigned int n1,
 const unsigned int lcount,
 double initThresCriteria)
:n3(n3), n2(n2), n1(n1), lcount(lcount), thresCriteria(initThresCriteria)
{
  std::cout << std::endl;
  std::cout << "  [Constructor] C++ Thresholding Class is being created..." << std::endl;
  std::cout << "  n1:   " << n1 <<std::endl;
  std::cout << "  n2:   " << n2 <<std::endl;
  std::cout << "  n3:   " << n3 <<std::endl;
	std::cout << "  lcount:   " << lcount <<std::endl;
	std::cout << "  init_thresCriteria :   " << thresCriteria <<std::endl;
  std::cout << "  [Constructor] Class construction ended!" << std::endl;
    
}


template<int dim>
void KWCThreshold<dim>::setUpClass(double *eta_pointer,
						   double *Xangle_pointer,
						   double *Yangle_pointer,
						   double *Zangle_pointer,
						   int *label_pointer,
						   double *Energy_pointer,
						   char boundary)
{
  std::cout << "   Setting Fast Marching workspace..." << std::endl;
  thetaBoundary = boundary;
  pcount = n1 * n2 * n3;
    
	//link pointers
	eta= eta_pointer;
	Xangles= Xangle_pointer; 
	Yangles= Yangle_pointer; 
	Zangles= Zangle_pointer; 
	labels= label_pointer; 
	energy= Energy_pointer;

	grid = new double[pcount]();
	indicator = new double[pcount]();

  heap = i_heap_create_empty_heap(pcount, pcount);
    
	//Create neighbor vector indexes
  std::vector<int>  vector3d;
  for(int i=-1;i<2;i++)  {
    for(int j=-1;j<2;j++)  {
      for(int k=-1;k<2;k++)  {
      
        vector3d.push_back(k);
        vector3d.push_back(j);
        vector3d.push_back(i);
                
         if(i==0 && j==0 && k==0){
           vector3d.clear();
         }else{
           threeDNeighbors.push_back(vector3d);
           vector3d.clear();
          }
       }
     }
    }
  
}


template<int dim>
void KWCThreshold<dim>::identifyInitBoundary(double epsilon)
{
  cout << "  Fast Marching Begins"<<endl;

   // When we conduct poly crystal simulation
   // as simulation goes, the J(||\theta||) becomes smaller as GBs of large misorientation
   // disappears. This results a dynamics become stuck
   // as we compare "J(||\theta||)" to some thresholding criteria
   // to identify Initial grain boundary.
   // Thus, we have to release thresholding criteria to some degree as simulation goes.
    
  bool stuckTest = true;
    
	if(thetaBoundary=='D')
	{ 
		cout<< "  Dirichlet boundary condition is implemented" <<endl;
	}else
  {
    cout<< "  Periodic boundary condition is implemented" <<endl;
  }
	
  //double *testField;
  //testField=new double[pcount];

  //Start search GBs
  do  {
    for(int i=0; i<n3; i++)  {
      for(int j=0;j<n2;j++){
        for(int k=0;k<n1;k++)  {

          int xm=(k-1+n1)%n1; int xp=(k+1+n1)%n1;
			    int ym=(j-1+n2)%n2; int yp=(j+1+n2)%n2;
			    int zm=(i-1+n3)%n3; int zp=(i+1+n3)%n3;
    
          double txp = n1  * 1.0 * (energy[i*n2*n1+j*n1+xp]) ;
          double txm = n1  * 1.0 * (energy[i*n2*n1+j*n1+xm]) ;
          double typ = n2  * 1.0 * (energy[i*n2*n1+yp*n1+k]) ;
          double tym = n2  * 1.0 * (energy[i*n2*n1+ym*n1+k]) ;
                
          double norm= 0.5*sqrt(txp*txp+txm*txm) + 0.5*sqrt(typ*typ+tym*tym);
              
          //double metric = norm/ sqrt(1+norm*norm);
			    //testField[i*n2*n1+j*n1+k]=norm;
			   
          if( norm < thresCriteria ){
          //if( metric < 2*epsilon ){
          // IF this is too narrow, It cannot go , // This result stucking motion
          //Initial condition of FM, Basins of grain
            indicator[i*n2*n1+j*n1+k]=0.0;
            grid[i*n2*n1+j*n1+k]=0.0;
           }else  {
            stuckTest = false;
					  indicator[i*n2*n1+j*n1+k]=FLT_MAX;
            grid[i*n2*n1+j*n1+k]=FLT_MAX;
           }
				
        
				  if(thetaBoundary=='D')  {
  
					  if(k<1 || k>n1-2 || j < 1 || j>n2-2)  {
					  // if \theta is Dirichlet, you don't need to consider them
					    indicator[i*n2*n1+j*n1+k]=0.0;
					    grid[i*n2*n1+j*n1+k]=0.0;
					  }
				  }
      }
    }
  }
  
        
  }while(stuckTest==true);
  
  //Make the below line alive if you want to watch how interior region is recognized
  //Output2DvtuScalar(n1, n2, n3, indicator, "indicator", "vtu/indicator_",0);
 
}


template<int dim>
bool KWCThreshold<dim>::thresholding()
{
  bool stuckTest = true; //Assume Stuck first
  int changeCount =0;
    
  i_heap_clear_heap(&heap);
  
  for(int i=0; i<n3; i++){
    for(int j=0; j<n2; j++){
      for(int k=0; k<n1; k++){
        
        int index=i*n2*n1+j*n1+k;
        if(grid[index]==0)  {
          
            int interior=1;
            for(int m=0; m<26; m++)  {
            
              int x=k+threeDNeighbors[m][0]; if(x>n1-1){x=0;}else if(x<0){x=n1-1;}
              int y=j+threeDNeighbors[m][1]; if(y>n2-1){y=0;}else if(y<0){y=n2-1;}
              int z=i+threeDNeighbors[m][2]; if(z>n3-1){z=0;}else if(z<0){z=n3-1;}
						
              if(grid[z*n2*n1 + y*n1 + x]==FLT_MAX)  {
                interior=0;
              }
                    
            }
            if(!interior) {//if 'interior==0' , they are 'the FRONT'
              i_heap_insert_node_with_location(&heap, index, grid[index], index);
            }
          }
	    }//k_loop
	  }//j_loop
  }//i_loop
  
  WatchFastMarching FMwatch(n1,n2,lcount,labels,grid);
  
  if(watch_fastmarching==1)  {
    FMwatch.videoInitialize();
  }

  //Short neighbor
  int xgrid[6]={1, 0, -1,  0, 0, 0};
  int ygrid[6]={0, 1, 0 , -1, 0, 0};
  int zgrid[6]={0, 0, 0 ,  0, 1,-1};
    
  while(!i_heap_empty(&heap))  {
       
    int root_index= executeFastMarchingIteration();

    if(grid[root_index]!=0) {
    
      int x = root_index % n1;
      int y = ( root_index % (n1*n2)) /n1;
      int z = root_index / (n1*n2);
    
      int xp=(x+1)%n1;
      int xm=(x-1+n1)%n1;
      int yp=(y+1)%n2;
      int ym=(y-1+n2)%n2;
      double dx1=grid[y*n1+x]-grid[y*n1+xm];
      double dx2=grid[y*n1+xp]-grid[y*n1+x];
      double dy1=grid[y*n1+x]-grid[ym*n1+x];
      double dy2=grid[yp*n1+x]-grid[y*n1+x];
      double dx, dy;
	  
	  // Root does not have minimum value in this loop,
      // since you are considering visited neighbor cell together in this step.
      // First, we judge where has closer approximation (where is upwind)
      // Here, "upwind" only means which is the relevent direction on grid
      // it does not consider the direction of the flow
      
      if(dx1+dx2>0){
            dx=dx1;
      }else{
            dx=dx2;
      }
      if(dy1+dy2>0){
            dy=dy1;
      }else{
            dy=dy2;
      }

	  // After we decided which the upwind,
      // you consider whether or not accepting the update from each direction
      // This may consider update from other direction
	  
      double norm=sqrt(dx*dx+dy*dy);
      int a,b;
      
      if(dx/norm>.5){
              a=1;
      }else if(dx/norm>-.5){
              a=0;
      }else{
              a=-1;
      }
      if(dy/norm>.5){
              b=1;
      }else if(dy/norm>-.5){
              b=0;
      }else{
              b=-1;
      }

      //calculate the direction of the characteristic
      int cx=(x-a+n1)%n1;
      int cy=(y-b+n2)%n2;

			//Threshold
      double previousLabel = labels[root_index];
      labels[root_index]=labels[cy*n1+cx];

      //We will count how many labels are change
      if(previousLabel != labels[root_index])  {
        changeCount +=1;
      }
           
      int updatedIndex = cy*n1+cx;
      
      if(watch_fastmarching==1)   {
        FMwatch.videoUpdate_pixel(updatedIndex);
      }

    }
  }//do until heap becomes empty

  if(watch_fastmarching==1)  {
    FMwatch.videoClose();
  }
  
  //(IF not enough points is updated...it speakouts warning
  // and suggests increase epsilon
   if(changeCount > n2*n3 * 0.01)  {
     stuckTest =false;
   }
    
  std::cout << "  FMM ended."<<endl;
  return stuckTest;
   
}

template<int dim>
double KWCThreshold<dim>::eikonalUpdate(int x, int y, int z)
{
  double result =0.0;
  //eikonal update in 3d with non-equal spacing for each direction
  double hx=1/(n1*1.0);
  double hy=1/(n2*1.0);
  double hz=1/(n3*1.0);
    
  int index = z*n2*n1 + y*n1 + x;
    
  int x_west=(x-1+n1)%n1; int x_east=(x+1+n1)%n1;
  int y_south=(y-1+n2)%n2; int y_north=(y+1+n2)%n2;
  int z_down = (z-1+n3)%n3; int z_up=(z+1+n3)%n3;
  
  if(thetaBoundary=='D')
  {
    x_west=x-1; if(x_west<0){x_west=0;}
    x_east=x+1; if(x_east>n1-1){x_east=n1-1;}
    y_south=y-1; if(y_south<0){y_south=0;}
    y_north=y+1; if(y_north>n2-1){y_north=n2-1;}
    z_down=z-1; if(z_down<0){z_down=0;}
    z_up=z+1;  if(z_up>n3-1){z_up=n3-1;}
  }
	
	
  int west= z*n2*n1+ y*n1 + x_west;
  int east= z*n2*n1+ y*n1 + x_east;
    
  int south = z*n2*n1 + y_south * n1 + x;
  int north = z*n2*n1 + y_north * n1 + x;
    
  int down = z_down*n2*n1 + y*n1 + x;
  int up = z_up*n2*n1 + y*n1 +x;

  //f inverse in slowness
  double finverse = (1-eta[index]) * (1-eta[index]);
    
  //testgrid is grid, true grid is indicator -- this is equivalent to upwind
  double U_H=fmin(grid[west],grid[east]);
  double U_V=fmin(grid[north],grid[south]);
  double U_Z=fmin(grid[up], grid[down]);
    
  //We solve the quadratic equation for U_ij, current cell, using up wind
    
  double hxs = hx*hx;
  double hys = hy*hy;
  double hzs = hz*hz;

  double a = hxs*hys + hys*hzs + hxs*hzs;
  double b = (-2.0) * (hys*hzs * U_H + hxs*hzs * U_V + hxs*hys* U_Z);
  double c = hys*hzs*pow(U_H,2)+ hxs*hzs*pow(U_V,2) +
			       hzs*hys*pow(U_Z,2)- hxs*hys*hzs*pow(finverse,2);
    
  double discrim = b*b - 4.0 *a *c;
    
  if(discrim>0) {
    result=(sqrt(b*b-4.0*a*c)-b)/(2.0*a);
  }else  {   //if discrim < 0,  then we take 2 dimensional updates in each plane
		// and choose the smallest.
    double result_xy;
    double result_yz;
    double result_zx;
    
    //in x-y plane, y-z, z-x plane separatedly
    result_xy = eikonalUpdate_2d(hxs,hys,U_H,U_V,finverse);
    result_yz = eikonalUpdate_2d(hys,hzs,U_V,U_Z,finverse);
    result_zx = eikonalUpdate_2d(hzs,hxs,U_Z,U_H,finverse);
    
    result=fmin(result_xy,result_yz);
		result=fmin(result, result_zx); 
        
		result=result_xy;
  }
	 
	 return result;
}


 
template<int dim>
double KWCThreshold<dim>::eikonalUpdate_2d
(double h1s, double h2s, double U1, double U2, double finverse)
{
	//conduct 2d eikonalUpdate, when discrim <0
  //h1s,h2s: h1^2 where h1= 1.0/n1;
  //U1,U2 : minimum value in 1 and 2 direction
    
  double result=0.0;
  double alpha = h1s/h2s;
  
	double a = h1s + h2s; 
	double b = -2.0 * (h1s*U1 + h2s*U2); 
	double c = h1s*U1*U1 + h2s*U2*U2 - h1s*h2s*pow(finverse,2);
		
  double discrim = b*b - 4*a*c;
    
  if(discrim>0)  {
    result = (sqrt(b*b-4*a*c)-b)/(2*a);
  }else if (U2 < U1)  {
    result = U2 + sqrt(h2s) * finverse;
  }else  {
    result = U1 + sqrt(h1s) * finverse;
  }
  
  return result;
}


template<int dim>
int KWCThreshold<dim>::executeFastMarchingIteration()
{
	//Purpose: return the index of root 
	//and eikonal update around the neighbros of root
	//The root will be deleted after this function 
	
  int *currentIndex=&heap.root[0].originalIndex;
  double *min=&heap.root[0].key;
    
  //index(i,j,k) is descrbied as "i*n2*n1+j*n1+k"
  // current point
  int cx = *currentIndex % n1;
  int cy = (*currentIndex % (n1*n2)) / n1;
	int cz = (*currentIndex)/ (n1*n2);
    
  int xgrid[6]={1, 0, -1,  0, 0, 0};
  int ygrid[6]={0, 1, 0 , -1, 0, 0};
  int zgrid[6]={0, 0, 0 ,  0, 1,-1};
  
  for(int m=0; m<6; m++)  {
    //periodicity applied.
    int x=(cx+xgrid[m]+n1)%n1;
    int y=(cy+ygrid[m]+n2)%n2;
    int z=(cz+zgrid[m]+n3)%n3;
		int index= z*n2*n1 + y*n1 + x;
		double possible = fmax(eikonalUpdate(x,y,z),*min);
		
		if(thetaBoundary=='D')  {
    
      int xn=cx+xgrid[m]; if(xn<0){xn=0;}else if(xn>n1-1){xn=n1-1;}
      int yn=cy+ygrid[m]; if(yn<0){yn=0;}else if(yn>n2-1){yn=n2-1;}
      int zn=cz+zgrid[m]; if(zn<0){zn=0;}else if(zn>n3-1){zn=n3-1;}
      index= zn*n2*n1 + yn*n1 + xn;
      possible = fmax(eikonalUpdate(xn,yn,zn),*min);
    }
    
    if(possible<grid[index])  {
      if(grid[index]==FLT_MAX)  {
        grid[index]=possible;
          if(possible<=1.1*grid[index])  {
            i_heap_insert_node_with_location(&heap, index, grid[index], index);
          }
       }else{
          int heapLocation=heap.locations[index];
          i_heap_decrease_key(&heap,heapLocation,grid[index]);
            }
    }
    
  }
  
  i_heap_delete_min(&heap);
  return *currentIndex;
	
}

template<int dim>
void KWCThreshold<dim>::run(double epsilon)
{
  identifyInitBoundary (epsilon);
  STUCK=thresholding ();
  std::cout <<"  Eaxmine Stuck Test :" << STUCK << std::endl;
}

template<int dim>
void KWCThreshold<dim>::freeMemory()
{
  i_heap_destroy_heap(&heap);
  delete [] grid;
  delete [] indicator;
}



#endif 
