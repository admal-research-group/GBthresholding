//  my_metrics.h
//  Created by Jaekwang Kim on 7/10/19.
//  Header files include simple function metrics..

#ifndef my_metrics_h
#define my_metrics_h

double periodic_distance(double x1, double x2);
double linear_interpolate(double x, unsigned int data_number, double *x_k, double *y_k );
double calc_h1_diff(double *data1, double *data2, int n1, int n2);
void calc_second_derivative (double *function_values, int i, int j, int k, int n3, int n2, int n1, double *second_derivative_values);


double linear_interpolate(double x, unsigned int data_number, double *x_k, double *y_k )
{
    //I assume the x_k are sorted in a increaseing order
    //Find n, such that x_k[n], x_k[n+1] that includes x
    
    unsigned int n=WINT_MAX;
    
    for(int i=0; i<data_number; i++)
    {
        if(x_k[i] <= x && x_k[i+1]> x)
        n=i;
    }
    
    if(n==WINT_MAX)
    {
        std::cout<<"Significant Error Detected in linear interpolation function"; exit(1);
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

 
void calc_second_derivative (double *function_values, int i, int j, int k, int n3, int n2, int n1, double *second_derivative_values)
{
    
    //it will be good if you add some...assert on the size of
    //second_derivative_valaues;
    
    //periodic
    /*
    int xm=(k-1+n1)%n1; int ym=(j-1+n2)%n2; int zm=(i-1+n3)%n3;
    int xp=(k+1+n1)%n1; int yp=(j+1+n2)%n2; int zp=(i+1+n2)%n2;
    */
    
    //finite length
    //theta- finite length
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
