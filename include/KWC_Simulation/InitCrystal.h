//  init_crystal.h
//  Created by Jaekwang Kim on 7/10/19.
//  outputting data for visualization

#include <iostream>
#include <fstream>
#include <ios>
#include "Metrics.h"
#include <time.h>


#ifndef initCrystal_h
#define initCrystal_h

using namespace std;

namespace InitializeCrystal
{


void RandomCrystalConfiguration3D(int n3, int n2, int n1, int lcount, int *labels,
		double *Xangles,double *Yangles, double *Zangles)
{
  std::cout<< "   Initialize with a 3D total " << lcount<< "-number of Polycrystal...\n";
  double seeds[3*lcount];
  double max=0;

  	  // Generate random seeds and angles
      for(int l=0;l<lcount;l++){

          seeds[2*l]=rand()/(RAND_MAX*1.0);
          seeds[2*l+1]=rand()/(RAND_MAX*1.0);
          seeds[2*l+2]=rand()/(RAND_MAX*1.0);

          Xangles[l]=rand()/(RAND_MAX*1.0);
          Yangles[l]=rand()/(RAND_MAX*1.0);
          Zangles[l]=rand()/(RAND_MAX*1.0);

          double misorientation =sqrt( Xangles[l]*Xangles[l]
									  +Yangles[l]*Yangles[l]
								      +Zangles[l]*Zangles[l]);
          if(misorientation>max){
              max=misorientation;}
      }

      //find the nearest seed and assign the label
      for(int i=0; i<n3; i++){
              for(int j=0;j<n2;j++){
                  for(int k=0;k<n1;k++){

              //present coordinate
              double x=k/(n1*1.0); double y=j/(n2*1.0); double z=i/(n3*1.0);

              double min=FLT_MAX;
              int minIndex=0;

              for(int l=0;l<lcount;l++){

                  double xdist=periodic_distance(x,seeds[2*l]);
                  double ydist=periodic_distance(y,seeds[2*l+1]);
                  double zdist=periodic_distance(z,seeds[2*l+2]);

                  double dist=xdist*xdist+ydist*ydist+zdist*zdist;

                  if(dist<min){
                  min=dist;
                  minIndex=l;}
              }
              labels[i*n2*n1+j*n1+k]=minIndex;
          }}} //close i,j,k
}



    
void RandomCrystalConfiguration2D(int n3, int n2, int n1, int lcount, int *labels, double maxZangle,
                                    double *Xangles,double *Yangles, double *Zangles)
{
    srand((unsigned int)time(NULL));
  
    std::cout<< "   Initialize with a 3D total " << lcount<< "-number of Polycrystal...\n";
    double seeds[3*lcount];
    double max=0;
        
    // Generate random seeds and angles
        for(int l=0;l<lcount;l++){
            
            seeds[2*l]=rand()/(RAND_MAX*1.0);
            seeds[2*l+1]=rand()/(RAND_MAX*1.0);
            seeds[2*l+2]=rand()/(RAND_MAX*1.0);
            
            Xangles[l]=0.0;
            Yangles[l]=0.0;
            Zangles[l]= maxZangle * rand()/(RAND_MAX*1.0);
            
        }
    
    /* Seed out*/
  
    std::ofstream pipe;
    pipe.open("crystal.txt");
    
    for(int l=0;l<lcount;l++) {
        
        pipe << std::setw(5) << seeds[2*l] << " "
        << std::setw(5) << seeds[2*l+1] << " "
        << std::setw(5) << seeds[2*l+2] << " "
        << std::setw(5) << Zangles[l] << std::endl;
    }
    
    pipe.close();
    
        //find the nearest seed and assign the label
        for(int i=0; i<n3; i++){
            for(int j=0;j<n2;j++){
                for(int k=0;k<n1;k++){
                    
                    //present coordinate
                    double x=k/(n1*1.0); double y=j/(n2*1.0); double z=i/(n3*1.0);
                    
                    double min=FLT_MAX;
                    int minIndex=0;
                    
                    for(int l=0;l<lcount;l++){
                        
                        double xdist=periodic_distance(x,seeds[2*l]);
                        double ydist=periodic_distance(y,seeds[2*l+1]);
                        double zdist=periodic_distance(z,seeds[2*l+2]);
                        
                        double dist=xdist*xdist+ydist*ydist; //2D distance
                        //double dist=xdist*xdist+ydist*ydist+zdist*zdist;
                        
                        if(dist<min){
                            min=dist;
                            minIndex=l;}
                    }
                    labels[i*n2*n1+j*n1+k]=minIndex;
    }}} //close i,j,k
}


  
void tempRandomCrystalConfiguration2D(int n3, int n2, int n1, int lcount, int *labels, double maxZangle,
                                    double *Xangles,double *Yangles, double *Zangles)
{
    // srand((unsigned int)time(NULL));
    srand(14);
  
    std::cout<< "   Initialize with a 3D total " << lcount<< "-number of Polycrystal...\n";
    double seeds[3*lcount];
    double max=0;
  
    // Generate random seeds and angles
        for(int l=0;l<lcount;l++){
            seeds[2*l]=rand()/(RAND_MAX*1.0);
            seeds[2*l+1]=rand()/(RAND_MAX*1.0);
            seeds[2*l+2]=rand()/(RAND_MAX*1.0);
            Xangles[l]=0.0;
            Yangles[l]=0.0;
            Zangles[l]=rand()/(RAND_MAX*1.0);
        }
  
        //find the nearest seed and assign the label
        for(int i=0; i<n3; i++){
            for(int j=0;j<n2;j++){
                for(int k=0;k<n1;k++){
                  
                    //present coordinate
                    double x=k/(n1*1.0); double y=j/(n2*1.0); double z=i/(n3*1.0);
                  
                    double min=FLT_MAX;
                    int minIndex=0;
                  
                    for(int l=0;l<lcount;l++){
                      
                        double xdist=periodic_distance(x,seeds[2*l]);
                        double ydist=periodic_distance(y,seeds[2*l+1]);
                        double zdist=periodic_distance(z,seeds[2*l+2]);
                      
                        double dist=xdist*xdist+ydist*ydist; //2D distance
                        //double dist=xdist*xdist+ydist*ydist+zdist*zdist;
                      
                        if(dist<min){
                            min=dist;
                            minIndex=l;}
                    }
                    labels[i*n2*n1+j*n1+k]=minIndex;
    }}} //close i,j,k
}

    
void Shirinking_configuration_3D(int n3, int n2, int n1, int *labels)
{
  std::cout<< "   Initialize with a shirinking sphere Grain Boundary...\n";

    for(int i=0; i<n3; i++) // z loop
    {
        for(int j=0; j<n2; j++) //y loop
        {
            for(int k=0; k<n1; k++) //x loop
            {
                double x = k/(n1*1.0); double y = j/(n2*1.0); double z = i/(n3*1.0);

                double min=FLT_MAX; int minIndex=0;

                if( pow(x-0.5,2)+pow(y-0.5,2)+pow(z-0.5,2) < 0.25*0.25)
                {
                    minIndex=1;
                }
                else{
                    minIndex=0;
                }
                labels[i*n1*n2+j*n1+k]=minIndex; //labels denotes
            }
        }
    }
}


void ShirinkingConfiguration2D(int n3, int n2, int n1, int *labels)
{
  std::cout<< "  Initialize with a shirinking sphere Grain Boundary...\n";

  for(int i=0; i<n3; i++)  { // z loop
    for(int j=0; j<n2; j++)  { //y loop
     for(int k=0; k<n1; k++)  { //x loop
           
       double x = k/(n1*1.0); double y = j/(n2*1.0);
       double min=FLT_MAX; int minIndex=0;

       if( pow(x-0.5,2)+pow(y-0.5,2) < 0.25*0.25) {
         minIndex=1;
       }else{
         minIndex=0;
       }
      
       labels[i*n1*n2+j*n1+k]=minIndex;
      }
    }
  }
}


void oneD_Bicrystal_configuration(int n3, int n2, int n1, int *labels)
{
    std::cout<< "   Initialize with a 1D bicrystal \n";
    
    for(int i=0; i<n3; i++) // z loop
    {
        for(int j=0; j<n2; j++) //y loop
        {
            for(int k=0; k<n1; k++) //x loop
            {
                double x = k/(n1*1.0); double y = j/(n2*1.0); double z = i/(n3*1.0);
                
                double min=FLT_MAX; int minIndex=0;
                
                if( x > 0.5)
                {
                    minIndex=1;
                }
                else{
                    minIndex=0;
                }
                labels[i*n1*n2+j*n1+k]=minIndex; //labels denotes
            }
        }
    }
}



void Temp_configuration(int n3, int n2, int n1, int *labels)
{
    std::cout<< "   Initialize with a 1D bicrystal \n";
    
    for(int i=0; i<n3; i++) // z loop
    {
        for(int j=0; j<n2; j++) //y loop
        {
            for(int k=0; k<n1; k++) //x loop
            {
                double x = k/(n1*1.0); double y = j/(n2*1.0); double z = i/(n3*1.0);
                
                double min=FLT_MAX; int minIndex=0;
                
                if( x >0.25 & x< 0.75 & y>0.25 & y<0.75)
                {
                    minIndex=1;
                }
                else{
                    minIndex=0;
                }
                labels[i*n1*n2+j*n1+k]=minIndex; //labels denotes
            }
        }
    }
}

void Triple_configuration(int n3, int n2, int n1, int *labels)
{
    std::cout<< "   Initialize with a triple configuration \n";
    for(int i=0; i<n3; i++) // z loop
    {
        for(int j=0; j<n2; j++) //y loop
        {
            for(int k=0; k<n1; k++) //x loop
            {	
                double x = k/(n1*1.0); double y = j/(n2*1.0); double z = i/(n3*1.0);
                double min=FLT_MAX; int minIndex=0;
				
            if(x<0.8){
                minIndex=1;
            }else
            {
                if(y>=0.5){minIndex=0;}
				else{minIndex=2;}    
            }
				labels[i*n1*n2+j*n1+k]=minIndex; //labels denotes
           } 
        }
    }
}


void HerringConfiguration(int n3, int n2, int n1, int *labels) {

  std::cout<< "   Initialize with a triple configuration \n";
    
  for(int i=0; i<n3; i++) { // z loop
    for(int j=0; j<n2; j++) { //y loop
      for(int k=0; k<n1; k++) { //x loop
            
        double x = k/(n1*1.0); double y = j/(n2*1.0); double z = i/(n3*1.0);
        double min=FLT_MAX; int minIndex=0;
				
        if(x<1.0){
          double value1 = (1.0-0.5)/(0.25-1.0) * (x-1.0) + 0.5;
          double value2 = -(1.0-0.5)/(0.25-1.0) * (x-1.0) + 0.5;
				
          if(y>value1)
            minIndex=2;
          else if(y<=value1 && y>value2)
            minIndex=1;
          else
            minIndex=0;
        }else {
            if(y>=0.5){minIndex=2;}
				else{minIndex=0;}    
        }
        labels[i*n1*n2+j*n1+k]=minIndex; //labels denotes
      }
    }
  }
}
    
void HerringConfigurationSecondType(int n3, int n2, int n1, int *labels) {
  std::cout<< "   Initialize with a triple configuration \n";
  
  for(int i=0; i<n3; i++) { // z loop
    for(int j=0; j<n2; j++) { //y loop
      for(int k=0; k<n1; k++) { //x loop
        
        double x = k/(n1*1.0); double y = j/(n2*1.0); double z = i/(n3*1.0);
        double min=FLT_MAX; int minIndex=0;
        
        if(x<0.5){
           minIndex=0;
        }else {
          if(y>0.5)
            minIndex=1;
          else
            minIndex=2;
        }
        
        labels[i*n1*n2+j*n1+k]=minIndex; //labels denotes
      }
    }
  }
}

    
    


void Shirinking_Compare(int n3, int n2, int n1, int *labels){

    std::cout<< "   Initialize with a triple configuration \n";

    for(int i=0; i<n3; i++) // z loop
    {
        for(int j=0; j<n2; j++) //y loop
        {
            for(int k=0; k<n1; k++) //x loop
            {
                double x = k/(n1*1.0); double y = j/(n2*1.0); double z = i/(n3*1.0);
                double min=FLT_MAX; int minIndex=0;

            if( (x-0.5)*(x-0.5) + (y-0.25)*(y-0.25) <0.15*0.15){
				minIndex=1;

            }else if ((x-0.5)*(x-0.5) + (y-0.75)*(y-0.75) <0.15*0.15)
            {
            	minIndex=2;
            }else
            {
            	minIndex=0;
            }

			labels[i*n1*n2+j*n1+k]=minIndex; //labels denotes

           }
        }
    }
}

}




#endif /* Init_crystal_h */
