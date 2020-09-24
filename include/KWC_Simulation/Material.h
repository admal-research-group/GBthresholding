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
  
  double s=1.0;
  //Original KWC model
  double computeJthetaSimpleMaterial (double misorientation)  {
  
        return s * misorientation;
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
    
    bool boolCovarianceModelSet=false;
    int nData;
    double *misorientationDatum;
    double *energyDatum;
    double *JfunDatum;
  
};





#endif /* Material_h */
