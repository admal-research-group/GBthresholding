//  my_metrics.h
//  Created by Jaekwang Kim on 7/10/19.
//  Header files include simple function metrics..

#ifndef my_metrics_h
#define my_metrics_h

double periodic_distance(double x1, double x2);
double linear_interpolate(double x, unsigned int data_number, double *x_k, double *y_k );
double calc_h1_diff(double *data1, double *data2, int n1, int n2);
void calc_second_derivative (double *function_values, int i, int j, int k, int n3, int n2, int n1, double *second_derivative_values);
void examineNeighborLabel (int n3, int n2, int n1, int cx, int cy, int cz, int &xm, int &xp, int &ym, int &yp, int &zm, int &zp, char labelBoundary);

void examineNeighborLabel (int n3, int n2, int n1, int cx, int cy, int cz, int &xm, int &xp, int &ym, int &yp, int &zm, int &zp, char labelBoundary)
{
  if(labelBoundary == 'N' || 'D')
  {
   xm=cx-1; if(xm<0){xm=0;}
   xp=cx+1; if(xp>n1-1){xp=n1-1;}
   
   ym=cy-1; if(ym<0){ym=0;}
   yp=cy+1; if(yp>n2-1){yp=n2-1;}
   
   zm=cz-1; if(zm<0){zm=0;}
   zp=cz+1; if(zp>n3-1){zp=n3-1;}
  }
  else if(labelBoundary == 'P')
  {
   xp=(cx+1+n1)%n1;
   xm=(cx-1+n1)%n1;
   yp=(cy+1+n2)%n2;
   ym=(cy-1+n2)%n2;
   zp=(cz+1+n3)%n3;
   zm=(cz-1+n3)%n3;
  }
  else
  {std::cout << "Label boundary is not well defined"<< std::endl; exit(1);}
  
}




double linear_interpolate(double x, unsigned int dataNumber, double *x_k, double *y_k )
{
    //I assume the x_k are sorted in a increaseing order
    //Find n, such that x_k[n], x_k[n+1] that includes x
    
    unsigned int n=WINT_MAX;
    
    for(int i=0; i<dataNumber; i++)
    {
        if(x_k[i] <= x && x_k[i+1]> x)
        n=i;
    }
    
    if(n==WINT_MAX){  // if not found, apply periodic data
       x = x - x_k[dataNumber-1];
     
      // search again
      for(int i=0; i<dataNumber; i++)
     {
        if(x_k[i] <= x && x_k[i+1]> x)
        n=i;
     }
     
    }
  
    if(n==WINT_MAX){ //still not found, terminate program
      std::cout<<"  Significant Error Detected in linear interpolation function"; exit(1);
    }

    double slope=(y_k[n+1]-y_k[n])/(x_k[n+1]-x_k[n]);
    double value = y_k[n] + slope * (x - x_k[n]);
    return value;
}



double calc_h1_diff(double *data1, double *data2, int n1, int n2){
    
    double h1=0;
    
    for(int i=0;i<n2;i++){
        for(int j=0;j<n1;j++){
            
            int xm=(j-1+n1)%n1;
            int ym=(i-1+n2)%n2;
            
            double data1x=n1*(data1[i*n1+j]-data1[i*n1+xm]);
            double data1y=n2*(data1[i*n1+j]-data1[ym*n1+j]);
            
            double data2x=n1*(data2[i*n1+j]-data2[i*n1+xm]);
            double data2y=n2*(data2[i*n1+j]-data2[ym*n1+j]);
            
            h1+=pow(data2x-data1x,2)+pow(data1y-data2y,2);
            
        }
    }
    
    h1/=n1*n2;
    
    return h1;
}


double periodic_distance(double x1, double x2){
    
    //assume Lx=1;
    
    double diff=fabs(x1-x2);
    double dist=.5-fabs(.5-diff);
    
    return dist;
}


void calc_first_derivative (double *function_values, int k, int j, int i, int n1, int n2, int n3, double *first_derivative_values,char boundary)
{
    int xp,xm,yp,ym,zm,zp;
  
    if(boundary=='P')
    {
      xm=(k-1+n1)%n1; ym=(j-1+n2)%n2; zm=(i-1+n3)%n3;
      xp=(k+1+n1)%n1; yp=(j+1+n2)%n2; zp=(i+1+n3)%n3;
    }else if(boundary=='N')
    {
    
      xp=(k+1); xm=(k-1); yp=(j+1); ym=(j-1); zp=(i+1); zm=(i-1);
    
      if(xp>n1-1)
      {xp=n1-1;}
    
      if(xm<0)
      {xm=0;}
    
      if(yp>n2-1)
      {yp=n2-1;}
    
      if(ym<0)
      {ym=0;}
    
      if(zp>n3-1)
      {zp=n3-1;}
    
      if(zm<0)
      {zm=0;}
    
    }else
    {
      std::cout << "Boundary condition for 'computeGradientEta' is not set properly"<< std::endl;
      exit(1);
    }
  
  
    first_derivative_values[0]= 0.5 * n1 * (function_values[i*n2*n1 + j*n1 + xp]
                                          - function_values[i*n2*n1 + j*n1 + xm]);
  
 
    first_derivative_values[1]= 0.5 * n2 * (function_values[i*n2*n1 + yp*n1 + k]
                                          -function_values[i*n2*n1 + ym*n1 + k]);
  
    first_derivative_values[2]= 0.5 * n3 * (function_values[zp*n2*n1 + j*n1 + k]
                                          -function_values[zm*n2*n1 + j*n1 + k]);
}



void calc_first_derivative_2 (double *function_values, int k, int j, int i, int n1, int n2, int n3, double *first_derivative_values,char boundary)
{
    int xp,xm,yp,ym,zm,zp;
  
  
    xm=(k-2+n1)%n1; ym=(j-2+n2)%n2; zm=(i-2+n3)%n3;
    xp=(k+2+n1)%n1; yp=(j+2+n2)%n2; zp=(i+2+n3)%n3;
   
  
    first_derivative_values[0]= 0.25 * n1 * (function_values[i*n2*n1 + j*n1 + xp]
                                          - function_values[i*n2*n1 + j*n1 + xm]);
  
 
    first_derivative_values[1]= 0.25 * n2 * (function_values[i*n2*n1 + yp*n1 + k]
                                          -function_values[i*n2*n1 + ym*n1 + k]);
  
    first_derivative_values[2]= 0.25 * n3 * (function_values[zp*n2*n1 + j*n1 + k]
                                          -function_values[zm*n2*n1 + j*n1 + k]);
}


//Very accurate estimate of first graident 
void calcFirstDerivative (double *function_values, int k, int j, int i, int n1, int n2, int n3,
                          double *first_derivative_values)
{  //periodic condition only
  int xp1,xp2,xp3,xp4,xm1,xm2,xm3,xm4;
  int yp1,yp2,yp3,yp4,ym1,ym2,ym3,ym4;
  
  xp1=(k+1+n1)%n1; xp2=(k+2+n1)%n1; xp3=(k+3+n1)%n1; xp4=(k+4+n1)%n1;
  xm1=(k-1+n1)%n1; xm2=(k-2+n1)%n1; xm3=(k-3+n1)%n1; xm4=(k-4+n1)%n1;
  
  yp1=(j+1+n2)%n2; yp2=(j+2+n2)%n2; yp3=(j+3+n2)%n2; yp4=(j+4+n2)%n2;
  ym1=(j-1+n2)%n2; ym2=(j-2+n2)%n2; ym3=(j-3+n2)%n2; ym4=(j-4+n2)%n2;
  
  //accuracy is eighth order
  first_derivative_values[0]= n1 * (
  (4.0/5.0) * function_values[i*n2*n1 + j*n1 + xp1]
  -(1.0/5.0) * function_values[i*n2*n1 + j*n1 + xp2]
  +(4.0/105) * function_values[i*n2*n1 + j*n1 + xp3]
  -(1.0/280) * function_values[i*n2*n1 + j*n1 + xp4]
  
  -(4.0/5.0) * function_values[i*n2*n1 + j*n1 + xm1]
  +(1.0/5.0) * function_values[i*n2*n1 + j*n1 + xm2]
  -(4.0/105) * function_values[i*n2*n1 + j*n1 + xm3]
  +(1.0/280) * function_values[i*n2*n1 + j*n1 + xm4]
  );
  
  
  first_derivative_values[1]= n1 * (
  (4.0/5.0) * function_values[i*n2*n1 + yp1*n1 + k]
  -(1.0/5.0) * function_values[i*n2*n1 + yp2*n1 + k]
  +(4.0/105) * function_values[i*n2*n1 + yp3*n1 + k]
  -(1.0/280) * function_values[i*n2*n1 + yp4*n1 + k]
  
  -(4.0/5.0) * function_values[i*n2*n1 + ym1*n1 + k]
  +(1.0/5.0) * function_values[i*n2*n1 + ym2*n1 + k]
  -(4.0/105) * function_values[i*n2*n1 + ym3*n1 + k]
  +(1.0/280) * function_values[i*n2*n1 + ym4*n1 + k]
  );
  
  
}


 
void calc_second_derivative (double *function_values, int i, int j, int k, int n3, int n2, int n1, double *second_derivative_values)
{
    
    //it will be good if you add some...assert on the size of
    //second_derivative_valaues;
    
    //periodic
  
    int xm=(k-1+n1)%n1; int ym=(j-1+n2)%n2; int zm=(i-1+n3)%n3;
    int xp=(k+1+n1)%n1; int yp=(j+1+n2)%n2; int zp=(i+1+n3)%n3;
    
    
    //finite length
    //theta- finite length
    /*
    int xp=(k+1); int xm=(k-1);
    int yp=(j+1); int ym=(j-1);
    int zp=(i+1); int zm=(i-1);
    
    if(xp>n1-1)
    {xp=n1-1;}
    
    if(xm<0)
    {xm=0;}
    
    if(yp>n2-1)
    {yp=n2-1;}
    
    if(ym<0)
    {ym=0;}
    
    if(zp>n3-1)
    {zp=n3-1;}
    
    if(zm<0)
    {zm=0;}
    */
    
    second_derivative_values[0]= n1*n1 * ( function_values[i*n2*n1 + j*n1 + xm]
                                 + function_values[i*n2*n1 + j*n1 + xp]
                                 - 2*function_values[i*n2*n1 + j*n1 + k]);
    
    second_derivative_values[1]= n2*n2 * ( function_values[i*n2*n1 + ym*n1 + k]
                                 + function_values[i*n2*n1 + yp*n1 + k]
                                 - 2*function_values[i*n2*n1 + j*n1 + k]);
    
    second_derivative_values[2]= n3*n3 * ( function_values[zm*n2*n1 + j*n1 + k]
                                + function_values[zp*n2*n1 + j*n1 + k]
                                - 2*function_values[i*n2*n1 + j*n1 + k]);
    
    
}


#endif /* my_metrics_h */
