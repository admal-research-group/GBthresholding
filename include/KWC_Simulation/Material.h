//
//  Material.h
//  myproject
//
//  Created by Jaekwang Kim on 3/4/20.
//  Copyright © 2020 Jaekwang Kim. All rights reserved.
//

#ifndef Material_h
#define Material_h

#include "Metrics.h"

//Module calculates Jtheta, which should be dependent on Material

class Material {
  
public :
  
  char materialType;
  double simpleKWC_s; // you should put this in the constructor
  
  //Constructor
  Material(const unsigned int n3,
    const unsigned int n2,
    const unsigned int n1,
    char materialType);

  
  void setUpClass(
    double *Xangle_pointer,
    double *Yangle_pointer,
    double *Zangle_pointer,
    int *label_pointer,
    double *fieldJ_pointer);
  
  void calculateFieldJ (char boundary);
  
  
  //Original KWC model
  double computeJthetaSimpleMaterial (double misorientation)  {
  
        return simpleKWC_s * misorientation;
  };
    
  //Covariance model: it will use linear interpolation to the available discrete data set
  double computeJthetaCovarianceModel (double misorientation)  {
    if(boolCovarianceModelSet==false)  {
      std::cout << "covariance model is not set properly "<< std::endl;
      exit(1);
    }
      return linear_interpolate(misorientation, nData, misorientationDatum, JfunDatum);
  }
    
  //Covariacne model: read discrete data of $J(theta)$ (inputs format)
  void setCovarianceModel(int dataNum, string filename)  {
    nData= dataNum;
    
    JfunDatum = new double[nData];
    energyDatum = new double[nData];
    misorientationDatum = new double[nData];
        
    ifstream readFile(filename);
        
    double trashValue;
    std::string varstr;
        
    if (readFile.is_open())  {
    
      getline(readFile,varstr); //pass 2 lines
      getline(readFile,varstr);
      
      for(int i=0; i< nData; i ++)  {
                readFile >> misorientationDatum[i]  >> JfunDatum[i] >> energyDatum[i] ;
                misorientationDatum[i] *= M_PI/180;
      }
      // Check whether it read correctly.
      for(int i=0; i< nData; i ++)  {
        std::cout << "i=" <<i <<": "<<  std::setw(8) <<  misorientationDatum[i]*180/M_PI
        << " " << std::setw(8) << JfunDatum[i] << std::setw(8) <<  energyDatum[i] << std::endl;
       }

      //exit(1); 

     }else  {
       cout << "Covariance Data File cannot be open" <<endl;
       exit(1);
    }
       cout << "Covariance Data files are read correctly...."<<endl;
       boolCovarianceModelSet=true;
    }
    
private:
  
    //"Function J of misorienation"
    bool boolCovarianceModelSet=false;
    int nData;
    double *misorientationDatum;
    double *energyDatum;
    double *JfunDatum;
  
    //"Field Quantitiy J"
    unsigned int n1,n2,n3; 
    int *labels;
    double *Xangles;
    double *Yangles;
    double *Zangles;
    double *fieldJ;
  
};


//Constructor
Material::Material(
  const unsigned int n3,
  const unsigned int n2,
  const unsigned int n1,
  char materialType)
:n3(n3), n2(n2), n1(n1),materialType(materialType)
{
  //Initialy set the parameter s=1. It is public variable, so can be modified
  //at any time
  simpleKWC_s=1.0;
}


void Material::setUpClass(
  double *Xangle_pointer,
  double *Yangle_pointer,
  double *Zangle_pointer,
  int *label_pointer,
  double *fieldJ_pointer)
{
 
  std::cout <<"....Setting Material Class pointers..." << std::endl;
  
  Xangles= Xangle_pointer;
  Yangles= Yangle_pointer;
  Zangles= Zangle_pointer;
  labels= label_pointer;
  fieldJ= fieldJ_pointer;

}

void Material::calculateFieldJ(char boundary)
{
  //set n1,n2,n3

  for(int i=0; i<n3; i++){
    for(int j=0;j<n2;j++){
      for(int k=0;k<n1;k++){
        
        //Periodic Neighbors
        /*
        int xm=(k-1+n1)%n1; int xp=(k+1+n1)%n1;
        int ym=(j-1+n2)%n2; int yp=(j+1+n2)%n2;
        int zm=(i-1+n3)%n3; int zp=(i+1+n3)%n3;
        */
        //Neighboring points
        int xp,xm,yp,ym,zp,zm;
        
        /*Neumann Condition*/
        if(boundary=='N')
        {
          xp=k+1; xm=k-1;
          yp=j+1; ym=j-1;
          zp=i+1; zm=i-1;
        
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
            
        }else if(boundary == 'P')
        {
          xm=(k-1+n1)%n1; xp=(k+1+n1)%n1;
          ym=(j-1+n2)%n2; yp=(j+1+n2)%n2;
          zm=(i-1+n3)%n3; zp=(i+1+n3)%n3;
        }else
        {
          std::cout << "Boundary condition for 'calculateFieldJ' is not set properly"<< std::endl;
          exit(1);
        }
        
        double local_jump_x,local_jump_y,local_jump_z;

        //Xangles
        local_jump_x= 1.0 * fabs(Xangles[labels[i*n1*n2+j*n1+xp]]-Xangles[labels[i*n1*n2+j*n1+xm]]);
        local_jump_y= 1.0 * fabs(Xangles[labels[i*n1*n2+yp*n1+k]]-Xangles[labels[i*n1*n2+ym*n1+k]]);
        local_jump_z= 1.0 * fabs(Xangles[labels[zp*n1*n2+j*n1+k]]-Xangles[labels[zm*n1*n2+j*n1+k]]);

        double misorientationX=0.0;
        misorientationX=  sqrt(local_jump_x * local_jump_x +
                               local_jump_y * local_jump_y +
                               local_jump_z * local_jump_z);
        
        //Yangles
        local_jump_x= 1.0 * fabs(Yangles[labels[i*n1*n2+j*n1+xp]]-Yangles[labels[i*n1*n2+j*n1+xm]]);
        local_jump_y= 1.0 * fabs(Yangles[labels[i*n1*n2+yp*n1+k]]-Yangles[labels[i*n1*n2+ym*n1+k]]);
        local_jump_z= 1.0 * fabs(Yangles[labels[zp*n1*n2+j*n1+k]]-Yangles[labels[zm*n1*n2+j*n1+k]]);

        double misorientationY=0.0;
        misorientationY=  sqrt(local_jump_x * local_jump_x +
                               local_jump_y * local_jump_y +
                               local_jump_z * local_jump_z);
        
        //Zangles
        local_jump_x= 1.0 * fabs(Zangles[labels[i*n1*n2+j*n1+xp]]-Zangles[labels[i*n1*n2+j*n1+xm]]);
        local_jump_y= 1.0 * fabs(Zangles[labels[i*n1*n2+yp*n1+k]]-Zangles[labels[i*n1*n2+ym*n1+k]]);
        local_jump_z= 1.0 * fabs(Zangles[labels[zp*n1*n2+j*n1+k]]-Zangles[labels[zm*n1*n2+j*n1+k]]);
        
        double misorientationZ=0.0;
        misorientationZ=  sqrt(local_jump_x * local_jump_x +
                               local_jump_y * local_jump_y +
                               local_jump_z * local_jump_z);

        double Misorientation =sqrt(misorientationX * misorientationX +
                                    misorientationY * misorientationY +
                                    misorientationZ * misorientationZ);
        
        //Pass one single misorientation value, regardless of DIMensioality of the system
        double misorientationENG=0.0;
        
        //Unavoidably, the local jump is diffused by central Finite difference into two grid point
        //It must have a same energy in the sense of weak form, so we take half of Jump energy
        //on each grid point
        
        if(materialType=='S'){
          misorientationENG = 0.5 * computeJthetaSimpleMaterial(Misorientation);
          
        }else if(materialType=='C')
        {
          misorientationENG = 0.5 * computeJthetaCovarianceModel(Misorientation);
        }else
        {
          cout << "Material Type is not properly set. Exit !!" <<endl;
          exit(1);
        }
        //For now, it is simply jump of \theta
        
        fieldJ[i*n1*n2+j*n1+k]= misorientationENG;
        if(misorientationENG >0.01)
        {
          //cout << "Showing non trivial energy" <<endl;
        }
        
      }
    }
  }

}





#endif /* Material_h */
