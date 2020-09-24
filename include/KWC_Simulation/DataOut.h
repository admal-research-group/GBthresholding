
//  Created by Jaekwang Kim on 7/10/19.
//  Outputting data for visualization

#include <iostream>
#include <fstream>
#include <sstream>
#include <ios>
#include <iomanip>
#include "Metrics.h"

#ifndef DataOut_h
#define DataOut_h

using namespace std;

namespace DataOut{

void exportOrientationInfo(int ny, int nx, int lcount, int *labels ,double *Zangles, string file_name)
{
  cout <<"Outputs Orientatoin_Info"<<endl;
  fstream fout;
  fout.open(file_name,ios::out);

  //Angle info
  fout << lcount << std::endl;
  for(int i=0; i<lcount; i++) {
    fout << Zangles[i] <<std::endl; }
    
  //label info
  for(int j=0;j<ny;j++) {
    for(int k=0;k<nx;k++) {
      fout << labels[j*nx+k] << " " ;
    }
    fout <<std::endl;
  }

  fout.close();
}
    
    
void Output3DvtuBinary(int nx, int ny, int nz, double *Xangles, double *Yangles, double *Zangles,
			   int *labels, string variable_name, string s,int data_num)
{
  char file_number[4];
  sprintf(file_number,"%03d",data_num);
  s+=file_number;
  s+=".vtu";

  fstream fout;
  char file_name[s.size()+1];
  strcpy(file_name,s.c_str());

  fout.open(file_name,ios::out| ios::binary);

  //Provide Data
  {
  fout << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << endl;
  //fout << "<Binary>" << endl;
  fout << "<UnstructuredGrid>" << endl;
  fout << "<Piece NumberOfPoints=\" "<< nx*ny*nz<<" \" " <<
  "NumberOfCells=\" " << (nx-1)*(ny-1)*(nz-1) <<" \"> " << endl;
  fout << "<PointData Scalars=\"scalars\">" << endl;
  fout << "<DataArray type=\"Float32\" Name=\" " << variable_name << "\" format=\"ascii\">" << endl;


  for(int i=0; i<nz; i++)  {
    for(int j=0; j<ny; j++)  {
      for(int k=0; k<nx; k++)  {
        //theta is magnitude of orientation respect to global rotation system,
        //which is chosen as the Identity tensor
        //You will lose information of the direction of misorientation
        
        double theta=sqrt(Xangles[labels[i*ny*nx+j*nx+k]]*Xangles[labels[i*ny*nx+j*nx+k]]
                        +Yangles[labels[i*ny*nx+j*nx+k]]*Yangles[labels[i*ny*nx+j*nx+k]]
                        +Zangles[labels[i*ny*nx+j*nx+k]]*Zangles[labels[i*ny*nx+j*nx+k]]);

        fout << std::setprecision(2) << theta << " " ;
        //fout.write((unsigned char *) &theta, sizeof(double) );
        //fout.write(reinterpret_cast<unsigned char *>(&theta),sizeof(double));
      }
      fout << endl;
    }
  }
  fout << "</DataArray>" << endl;
  fout << "</PointData>" << endl;
  }


    //Provide Coordinate Points
    {
        fout << "<Points>" <<endl;
        fout << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" <<endl;
        //Increase x->y>-z
        for(int i=0; i<nz; i++)
        {for(int j=0; j<ny; j++)
        {for(int k=0; k<nx; k++)
        {
        	/*
        	double x,y,z;
        	x=(1.0*k)/nx;
        	y=(1.0*j)/ny;
        	z=(1.0*i)/nx;

        	fout.write((char *) &x, sizeof(double) );
        	fout.write((char *) &y, sizeof(double) );
        	fout.write((char *) &z, sizeof(double) );
			    */
          
        	fout << (1.0*k)/nx << " " << (1.0*j)/ny << " " << (1.0*i)/nz << endl;

        }}}

        fout << "</DataArray>" << endl;
        fout << "</Points>" <<endl;
    }

    // Cell Connectivity
    {
        fout << "<Cells>" <<endl;
        fout << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" <<endl;

        // If it works, change to cell type 12
        for(int i=0; i<nz-1; i++)
        {
            for(int j=0; j<ny-1; j++)
            {
                for(int k=0; k<nx-1; k++)
                {

                    int org= i*nx*ny + j*nx +k;

                    fout << org << " " << org+1 << " " << org+1+nx << " " << org+nx << " "
                    << org+nx*ny << " " << org+nx*ny+1 << " " << org+nx*ny+nx+1 << " " << org+nx*ny+nx << endl;

                    /*
                	int org= i*nx*ny + j*nx +k;
                	int vertex2= org+1;
                	int vertex3= org+1+nx;
                	int vertex4= org+nx ;
                	int vertex5= org+nx*ny	;
                	int vertex6= org+nx*ny+1 ;
                	int vertex7= org+nx*ny+nx+1;
					        int vertex8= org+nx*ny+nx;

					        fout.write((char *) &org, sizeof(int) );
					        fout.write((char *) &vertex2, sizeof(int) );
					        fout.write((char *) &vertex3, sizeof(int) );
					        fout.write((char *) &vertex4, sizeof(int) );
					        fout.write((char *) &vertex5, sizeof(int) );
					        fout.write((char *) &vertex6, sizeof(int) );
					        fout.write((char *) &vertex7, sizeof(int) );
					        fout.write((char *) &vertex8, sizeof(int) );
					        */
        }}}
        fout << "</DataArray>" << endl;
    }


    // Cell Type Offset
    {
        fout << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" <<endl;
        int count=1;
        for(int i=0; i<nz-1; i++)
        {   for(int j=0; j<ny-1; j++)
            {       for(int k=0; k<nx-1; k++)
                    {
                      //int out_number = count * 8;
                      //fout.write((char *) &out_number, sizeof(int) );
            		  fout << count*8 << " ";
                      count+=1;
                    }
                    fout << std::endl;
            }
        }
        fout << "</DataArray>" << endl;
    }

    // Cell Type information
    {
      fout << "<DataArray type=\"UInt8\" Name=\"types\" format=\"binary\">" <<endl;

        for(int i=0; i<nz-1; i++)
        {   for(int j=0; j<ny-1; j++)
            {for(int k=0; k<nx-1; k++)
                {
            		int out_number = 12;
            		fout.write((char *) &out_number, sizeof(int) );
            		//fout << 12 << " ";
                }
                fout << std::endl;}
        }

        fout << std::endl;
        fout << "</DataArray>" << endl;
        fout << "</Cells>" << endl;
        fout << "</Piece>" << endl;
        fout << "</UnstructuredGrid>" << endl;
        fout << "</VTKFile>" << endl;
    }

    //VTK file format end
	fout.close();

}

//output double type data
void Output3Dvtu(int nx, int ny, int nz, double *Xangles, double *Yangles, double *Zangles,
			   int *labels, string variable_name, string s,int data_num)
{
	char file_number[4];
	sprintf(file_number,"%03d",data_num);
	s+=file_number;
	s+=".vtu";

	fstream fout;
	char file_name[s.size()+1];
	strcpy(file_name,s.c_str());
	fout.open(file_name,ios::out);

	//Provide Data
	{
		fout << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << endl;
		fout << "<UnstructuredGrid>" << endl;
		fout << "<Piece NumberOfPoints=\" "<< nx*ny*nz<<" \" " <<
				"NumberOfCells=\" " << (nx-1)*(ny-1)*(nz-1) <<" \"> " << endl;
		fout << "<PointData Scalars=\"scalars\">" << endl;
		fout << "<DataArray type=\"Float32\" Name=\" " << variable_name << "\" format=\"ascii\">" << endl;


	    for(int i=0; i<nz; i++)
	    	{    for(int j=0; j<ny; j++)
	            {
	                for(int k=0; k<nx; k++)
	                {
	                	//theta is magnitude of orientation respect to global rotation system,
	                	//which is chosen as the Identity tensor
	                	//You will lose information of the direction of misorientation
	                	double theta=sqrt(Xangles[labels[i*ny*nx+j*nx+k]]*Xangles[labels[i*ny*nx+j*nx+k]]
										 +Yangles[labels[i*ny*nx+j*nx+k]]*Yangles[labels[i*ny*nx+j*nx+k]]
										 +Zangles[labels[i*ny*nx+j*nx+k]]*Zangles[labels[i*ny*nx+j*nx+k]]);

	                    fout << theta << " " ;
	                }
	                fout << endl;
	            }

	    	}
	    fout << "</DataArray>" << endl;
	    fout << "</PointData>" << endl;
	}


    //Provide Coordinate Points
    {
        fout << "<Points>" <<endl;
        fout << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" <<endl;
        //Increase x->y>-z
        for(int i=0; i<nz; i++)
        {for(int j=0; j<ny; j++)
        {for(int k=0; k<nx; k++)
        {fout << (1.0*k)/nx << " " << (1.0*j)/ny << " " << (1.0*i)/nz << endl;
        }}}

        fout << "</DataArray>" << endl;
        fout << "</Points>" <<endl;
    }

    // Cell Connectivity
    {
        fout << "<Cells>" <<endl;
        fout << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" <<endl;

        // If it works, change to cell type 12
        for(int i=0; i<nz-1; i++)
        {
            for(int j=0; j<ny-1; j++)
            {
                for(int k=0; k<nx-1; k++)
                {
                    int org= i*nx*ny + j*nx +k;

                    fout << org << " " << org+1 << " " << org+1+nx << " " << org+nx << " "
                    << org+nx*ny << " " << org+nx*ny+1 << " " << org+nx*ny+nx+1 << " " << org+nx*ny+nx << endl;
        }}}
        fout << "</DataArray>" << endl;
    }


    // Cell Type Offset
    {
        fout << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" <<endl;
        int count=1;
        for(int i=0; i<nz-1; i++)
        {   for(int j=0; j<ny-1; j++)
            {       for(int k=0; k<nx-1; k++)
                    {
                      fout << count*8 << " ";
                      count+=1;
                    }
                    fout << std::endl;
            }
        }
        fout << "</DataArray>" << endl;
    }

    // Cell Type information
    {
      fout << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" <<endl;

        for(int i=0; i<nz-1; i++)
        {   for(int j=0; j<ny-1; j++)
            {for(int k=0; k<nx-1; k++)
                { fout << 12 << " ";}
                fout << std::endl;}
        }

        fout << std::endl;
        fout << "</DataArray>" << endl;
        fout << "</Cells>" << endl;
        fout << "</Piece>" << endl;
        fout << "</UnstructuredGrid>" << endl;
        fout << "</VTKFile>" << endl;
    }

  //VTK file format end
	fout.close();

}


void Output3Dvtu(int nx, int ny, int nz, double *data, string variable_name, string s,int data_num)
{
	char file_number[4];
	sprintf(file_number,"%03d",data_num);
	s+=file_number;
	s+=".vtu";

	fstream fout;
	char file_name[s.size()+1];
	strcpy(file_name,s.c_str());
	fout.open(file_name,ios::out);

	//Provide Data
	{
		fout << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << endl;
		fout << "<UnstructuredGrid>" << endl;
		fout << "<Piece NumberOfPoints=\" "<< nx*ny*nz<<" \" " <<
				"NumberOfCells=\" " << (nx-1)*(ny-1)*(nz-1) <<" \"> " << endl;
		fout << "<PointData Scalars=\"scalars\">" << endl;
		fout << "<DataArray type=\"Float32\" Name=\" " << variable_name << "\" format=\"ascii\">" << endl;


	    for(int i=0; i<nz; i++)
	    	{    for(int j=0; j<ny; j++)
	            {
	                for(int k=0; k<nx; k++)
	                {
	                    fout << data[i*ny*nx+j*nx+k] << " " ;
	                }
	                fout << endl;
	            }

	    	}
	    fout << "</DataArray>" << endl;
	    fout << "</PointData>" << endl;
	}


    //Provide Coordinate Points
    {
        fout << "<Points>" <<endl;
        fout << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" <<endl;
        //Increase x->y>-z
        for(int i=0; i<nz; i++)
        {for(int j=0; j<ny; j++)
        {for(int k=0; k<nx; k++)
        {fout << (1.0*k)/nx << " " << (1.0*j)/ny << " " << (1.0*i)/nz << endl;
        }}}

        fout << "</DataArray>" << endl;
        fout << "</Points>" <<endl;
    }

    // Cell Connectivity
    {
        fout << "<Cells>" <<endl;
        fout << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" <<endl;

        // If it works, change to cell type 12
        for(int i=0; i<nz-1; i++)
        {
            for(int j=0; j<ny-1; j++)
            {
                for(int k=0; k<nx-1; k++)
                {
                    int org= i*nx*ny + j*nx +k;

                    fout << org << " " << org+1 << " " << org+1+nx << " " << org+nx << " "
                    << org+nx*ny << " " << org+nx*ny+1 << " " << org+nx*ny+nx+1 << " " << org+nx*ny+nx << endl;
        }}}
        fout << "</DataArray>" << endl;
    }


    // Cell Type Offset
    {
        fout << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" <<endl;
        int count=1;
        for(int i=0; i<nz-1; i++)
        {   for(int j=0; j<ny-1; j++)
            {       for(int k=0; k<nx-1; k++)
                    {
                      fout << count*8 << " ";
                      count+=1;
                    }
                    fout << std::endl;
            }
        }
        fout << "</DataArray>" << endl;
    }

    // Cell Type information
    {
      fout << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" <<endl;

        for(int i=0; i<nz-1; i++)
        {   for(int j=0; j<ny-1; j++)
            {for(int k=0; k<nx-1; k++)
                { fout << 12 << " ";}
                fout << std::endl;}
        }

        fout << std::endl;
        fout << "</DataArray>" << endl;
        fout << "</Cells>" << endl;
        fout << "</Piece>" << endl;
        fout << "</UnstructuredGrid>" << endl;
        fout << "</VTKFile>" << endl;
    }

    //VTK file format end
	fout.close();

}


// select z=a and output (x,y,a)
void Output2D(int n1, int n2, int zvalue, double *data, string s)
{
  char file_name[s.size()+1];
  strcpy(file_name,s.c_str());
  FILE *fptr = fopen(file_name,"w");
    
  int i = zvalue;
  for(int j=0;j<n2;j++){
    for(int k=0;k<n1;k++){
      fprintf(fptr, "%f ", data[i*n1*n2+j*n1+k]);
    }
    fprintf(fptr, "\n");
  }
  fclose(fptr);
}

//output int type data
void Output2D(int n1, int n2, int zvalue, int *data, string s)
{
    char file_name[s.size()+1];
    strcpy(file_name,s.c_str());
    
    FILE *fptr = fopen(file_name,"w");
    
    int i = zvalue;
    
    for(int j=0;j<n2;j++){
        for(int k=0;k<n1;k++){
            
            fprintf(fptr, "%d ", data[i*n1*n2+j*n1+k]);
        }
        fprintf(fptr, "\n");
    }
    fclose(fptr);
}

void Output2DvtuAngle(int nx, int ny, int nz, double *Zangles,
                     int *labels, string variable_name, string s,int data_num)
{
  char file_number[4];
  sprintf(file_number,"%03d",data_num);
  s+=file_number;
  s+=".vtu";
        
  fstream fout;
  char file_name[s.size()+1];
  strcpy(file_name,s.c_str());
  fout.open(file_name,ios::out);
        
  //Provide Data
{
  fout << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << endl;
  fout << "<UnstructuredGrid>" << endl;
  fout << "<Piece NumberOfPoints=\" "<< nx*ny*nz<<" \" "
         <<  "NumberOfCells=\" " << (nx-1)*(ny-1)<<" \"> " << endl;
  fout << "<PointData Scalars=\"scalars\">" << endl;
  fout << "<DataArray type=\"Float32\" Name=\" " << variable_name << "\" format=\"ascii\">" << endl;
        
    for(int j=0; j<ny; j++){
      for(int k=0; k<nx; k++){
        //theta is magnitude of orientation respect to global rotation system,
        //which is chosen as the Identity tensor
        //You will lose information of the direction of misorientation
        double theta=sqrt(Zangles[labels[j*nx+k]]*Zangles[labels[j*nx+k]]);
        fout << theta << " " ;
       }
      fout << endl;
    }
      
  fout << "</DataArray>" << endl;
  fout << "</PointData>" << endl;
}
        
        
//Provide Coordinate Points
{
  fout << "<Points>" <<endl;
  fout << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" <<endl;
  
  //Increase x->y>-z
  for(int i=0; i<nz; i++){
    for(int j=0; j<ny; j++){
      for(int k=0; k<nx; k++){
        fout << (1.0*k)/nx << " " << (1.0*j)/ny << " " << (1.0*i)/nz << endl;
      }
    }
  }
            
  fout << "</DataArray>" << endl;
  fout << "</Points>" <<endl;
}
        
// Cell Connectivity
{
  fout << "<Cells>" <<endl;
  fout << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" <<endl;
            
  for(int j=0; j<ny-1; j++) {
    for(int k=0; k<nx-1; k++) {
      int org= j*nx +k;
      fout << org << " " << org+1 << " " << org+1+nx << " " << org+nx << endl;
    }
  }
  fout << "</DataArray>" << endl;
}
        
// Cell Type Offset
{
  fout << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" <<endl;
  int count=1;
            
  for(int j=0; j<ny-1; j++){
    for(int k=0; k<nx-1; k++){
      fout << count*4 << " ";
      count+=1;
    }
    fout << std::endl;
  }
    
  fout << "</DataArray>" << endl;
}
        
// Cell Type information
{
  fout << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" <<endl;
 
  for(int j=0; j<ny-1; j++) {
    for(int k=0; k<nx-1; k++){
      fout << 9 << " ";
    }
    fout << std::endl;
  }
            
  fout << std::endl;
  fout << "</DataArray>" << endl;
  fout << "</Cells>" << endl;
  fout << "</Piece>" << endl;
  fout << "</UnstructuredGrid>" << endl;
  fout << "</VTKFile>" << endl;
}
  //VTK file format end
  fout.close();
}


void Output2DvtuScalar(int nx, int ny, int nz, double *data, string variable_name, string s, int data_num)
{
    char file_number[4];
    sprintf(file_number,"%03d",data_num);
    s+=file_number;
    s+=".vtu";
        
    fstream fout;
    char file_name[s.size()+1];
    strcpy(file_name,s.c_str());
    fout.open(file_name,ios::out);
        
    //Provide Data
    {
        fout << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << endl;
        fout << "<UnstructuredGrid>" << endl;
        fout << "<Piece NumberOfPoints=\" "<< nx*ny*nz<<" \" " <<
        "NumberOfCells=\" " << (nx-1)*(ny-1)<<" \"> " << endl;
        fout << "<PointData Scalars=\"scalars\">" << endl;
        fout << "<DataArray type=\"Float32\" Name=\" " << variable_name << "\" format=\"ascii\">" << endl;
            
            for(int j=0; j<ny; j++)
            {
                for(int k=0; k<nx; k++)
                {
                    fout << data[j*nx+k] << " " ;
                }
                fout << endl;
            }
            fout << "</DataArray>" << endl;
            fout << "</PointData>" << endl;
        }
    
        //Provide Coordinate Points
        {
            fout << "<Points>" <<endl;
            fout << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" <<endl;
            //Increase x->y>-z
            for(int i=0; i<nz; i++)
            {for(int j=0; j<ny; j++)
            {for(int k=0; k<nx; k++)
            {fout << (1.0*k)/nx << " " << (1.0*j)/ny << " " << (1.0*i)/nz << endl;
            }}}
            
            fout << "</DataArray>" << endl;
            fout << "</Points>" <<endl;
        }
        
        // Cell Connectivity
        {
            fout << "<Cells>" <<endl;
            fout << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" <<endl;
            
            for(int j=0; j<ny-1; j++)
            {
                for(int k=0; k<nx-1; k++)
                {
                    int org= j*nx +k;
                    fout << org << " " << org+1 << " " << org+1+nx << " " << org+nx << endl;
                }
                
            }
            fout << "</DataArray>" << endl;
        }
        
        
        // Cell Type Offset
        {
            fout << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" <<endl;
            int count=1;
            
            for(int j=0; j<ny-1; j++)
            {
                for(int k=0; k<nx-1; k++)
                {
                    fout << count*4 << " ";
                    count+=1;
                }
                fout << std::endl;
            }
            fout << "</DataArray>" << endl;
        }
        
        // Cell Type information
        {
            fout << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" <<endl;
            for(int j=0; j<ny-1; j++)
            {
                for(int k=0; k<nx-1; k++)
                {
                    fout << 9 << " ";
                }
                fout << std::endl;
            }
            
            fout << std::endl;
            fout << "</DataArray>" << endl;
            fout << "</Cells>" << endl;
            fout << "</Piece>" << endl;
            fout << "</UnstructuredGrid>" << endl;
            fout << "</VTKFile>" << endl;
        }
        
        //VTK file format end
        fout.close();
        
}
    
    
    
//output int type data
void Output2dAngles(int n1, int n2, int zvalue,
                     int *labels, double *angles, string s, int data_num)
{
  char file_number[4];
  sprintf(file_number,"%03d",data_num);
  s+=file_number; s+=".txt";
    
  char file_name[s.size()+1];
  strcpy(file_name,s.c_str());
    
  FILE *fptr = fopen(file_name,"w");
    
  int i = zvalue;
    
  for(int j=0;j<n2;j++){
    for(int k=0;k<n1;k++){
      fprintf(fptr, "%f ", angles[labels[i*n1*n2+j*n1+k]]);
    }
      
    fprintf(fptr, "\n");
  }
    
  fclose(fptr);
}


//output int type data
void Output1Dsolution(int n1, int n2, int zvalue,
                     double *values, string s, int data_num)
{
  char file_number[4];
  sprintf(file_number,"%03d",data_num);
  s+=file_number; s+=".txt";
  
  char file_name[s.size()+1];
  strcpy(file_name,s.c_str());
  
  FILE *fptr = fopen(file_name,"w");
  
  int i = 0;
  int j = n2/2;
  
  for(int k=0;k<n1;k++){
      double x = (1.0)*k/n1;
      fprintf(fptr, "%f %f ", x, values[i*n1*n2+j*n1+k]);
      fprintf(fptr, "\n");
  }
  
  
  fclose(fptr);
}

void PrepareFFMPEG2DPixels(int n1, int n2,int zvalue,
                           unsigned char *pixels, int *labels,
                           unsigned char *colors)
{
  int i = zvalue;
  for(int j=0;j<n2;j++) {
    for(int k=0;k<n1;k++){
      pixels[j*n1+k]=colors[labels[i*n2*n1+j*n1+k]];
    }
  }
    
}
    
}; //End of Name Space


#endif /* output_data_h */
