/*
	Author: Crispin Contreras
	Class: Physics 905
	Purpose: Solves the 3D radial Schrodinger equation in a harmonic potetial
	with two electrons. One way is using the Jacobi method taking into account
	the interaction between electrons and the other using the Armadillo library. 

*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include "armadillo"

using namespace arma;
using namespace std;
ofstream ofile;


void rotate(mat& , mat& , int , int , int);

int main()
{	
	int n; //number of steps
	int iterations=0;
	int inter; 
	char outfile[30]; //Rho and eigenfunctions
	char Eigen[30]; // Eigenvalues, elapsed time, iterations
	double epsilon = 1.0e-8;
	double rho_max, rho_min, w_r; 
	double *rho_i, *V_Pot, *V_Pin; 
	double h; 

	//Read Input
	cout<<"Enter the number of steps n: ";
	cin>>n;
	cout<<"Enter the value of rho_max: ";
	cin>>rho_max; 
	cout<<"Enter the value of rho_min: ";
	cin>>rho_min;
	cout<<"Do you want non-interacting(0) or interacting (1): ";
	cin>>inter;
	cout<<"Enter the value of the frequency w_r: ";
	cin>>w_r; 
	cout<<"Enter the name of the outputfile for Rho and Eigenfunctions: ";
	cin>>outfile;
	cout<<"Enter the name of file for Eigenvalues, elapse time, and total number of iterations: ";
	cin>>Eigen; 
		

	//calculate the step size
	h = (rho_max - rho_min)/n; 

	//Allocate Memory 
	rho_i = new double [n+1];
	V_Pot = new double [n-1];
	V_Pin = new double [n-1]; 
	
	//Fillout value of potential 	
	for(int i=0; i<=n; i++)
	{
	   rho_i[i] = rho_min + i*h; 

	}
	
	for (int i=0; i< (n-1); i++)
	{
	    //Non-Interacting Potential 
	    V_Pot[i]=w_r*w_r*rho_i[i+1]*rho_i[i+1];
	
	    //Interacting Potential 
	    V_Pin[i] = V_Pot[i] + (1/rho_i[i+1]);	
		
	}


	//Declare variables
	mat A(n-1,n-1);
	mat R(n-1,n-1);  
	
	//Filling the matrices
	A.zeros(n-1, n-1); 
	R.eye(n-1,n-1);

	for (int i=0 ; i<(n-1); i++)
	{	

	     if(inter == 0)
	     {
	        A(i, i) = 2/(h*h) + V_Pot[i];
	     }
	     if(inter == 1)
	     {
		A(i,i) = 2/(h*h) + V_Pin[i]; 	

	     }

	     if( i<(n-2))
	     {
		A(i,i+1) = A(i+1, i)= -1/(h*h);		
					
	     }
		
	}

	
	//Copy for armadillo
	mat C = A;
		

	//Implementing Jacobi Method
	double  max = fabs(A(0,1)); 
	int l; //indices
	int k; 

	clock_t start, finish;	
	start = clock();
		
	while(max > epsilon)
	{
		max =0.0;
	        for (int i = 0 ; i< (n-2) ; i++)
       		{
                    for(int j =i+1 ; j<(n-1); j++)
                    {
                       if (fabs(A(i, j)) > max)
                       {                  
 			 max =fabs(A(i,j)); 		 
	  
			//Find the indices
                          k=i;
                          l =j;	
                       }
                     }
                }//End of i loop
		
		//Call function to rotate
		rotate(A,R,k,l,n);  

		iterations++;
	}

	finish=clock();
	double time_j = (double (finish-start)/CLOCKS_PER_SEC); 

	//Assorting the eigevalues and eigenvectors
	vec Max(n-1); 
	int Loc[3]={0,1,2};
	double temp;
 
	for(int i=0; i<(n-1);i++)
	{
	   Max(i)=A(i,i); 
	}

	for(int i=0; i<(n-1);i++)
	{
	   for(int j=i+1; j<(n-1);j++)
	   {
	      if(Max(i)>Max(j))
	      {
		 temp=Max(i); 
		 Max(i)=Max(j); 
		 Max(j)=temp;
		 if(i<3)
		 {
		   Loc[i]=j;  
		 }
	      }

	   }

	}
	
	//Solving using Armadillo
        vec eigval(n-1);
        mat eigvec(n-1,n-1);

        start = clock();
        eig_sym(eigval, eigvec, C);
        finish = clock();

        double time_ar =(double(finish-start)/CLOCKS_PER_SEC);

	//three lowest states Armadillo
	vec V0 = eigvec.col(0);
	vec V1 = eigvec.col(1); 
	vec V2 = eigvec.col(2); 
	//From Jacobi
	vec R0 = R.col(Loc[0]);
	vec R1 = R.col(Loc[1]); 
	vec R2 = R.col(Loc[2]);

	//Unit Test, the norm should be equal to one 
	double V0Sum=0.0, V1Sum=0.0, V2Sum=0.0;
	double R0Sum=0.0, R1Sum=0.0, R2Sum=0.0;

	for(int i=0; i<(n-1); i++)
	{
  	   V0Sum+=V0(i)*V0(i); 
	   V1Sum+=V1(i)*V1(i); 
	   V2Sum+=V2(i)*V2(i);
	   R0Sum+=R0(i)*R0(i); 
	   R1Sum+=R1(i)*R1(i); 
	   R2Sum+=R2(i)*R2(i); 	

	} 
	
	cout<<"Norm of Ground (Ar) "<<V0Sum<<endl;
        cout<<"Norm of 1st (AR) "<<V1Sum<<endl;
        cout<<"Norm of 2nd (AR) "<<V2Sum<<endl;
        cout<<"Norm of Ground (Ja) "<<R0Sum<<endl;
        cout<<"Norm of 1st (Ja) "<<R1Sum<<endl;
        cout<<"Norm of 2nd (Ja) "<<R2Sum<<endl;


	//Writing to file, Rho, EigenFunctions
	ofile.open(outfile);
	ofile<<setiosflags(ios::showpoint | ios::uppercase);
	ofile<<"#Rho"<<setw(23)<<"Ground(Ar)";
	ofile<<setw(20)<<"1st(Ar)"<<setw(20)<<"2nd(Ar)";
	ofile<<setw(20)<<"Ground(Ja)"<<setw(20)<<"1st(Ja)";
	ofile<<setw(20)<<"2nd(Ja)"<<endl; 

	for(int i=0; i<(n-1); i++)
	{
   	   ofile<<rho_i[i+1]; 
	   ofile<<setw(20)<<(1/h)*V0(i)*V0(i);
	   ofile<<setw(20)<<(1/h)*V1(i)*V1(i);
	   ofile<<setw(20)<<(1/h)*V2(i)*V2(i);
	   ofile<<setw(20)<<(1/h)*R0(i)*R0(i);
	   ofile<<setw(20)<<(1/h)*R1(i)*R1(i);
	   ofile<<setw(20)<<(1/h)*R2(i)*R2(i)<<endl;
	}

	ofile.close();
	
	//Wrting to file with Eigenvalues  

	ofile.open(Eigen);
	ofile<<setiosflags(ios::showpoint | ios::uppercase); 
	ofile<<"#Eig(Ja)"<<setw(25)<<"Time(Ja)";
	ofile<<setw(25)<<"iterations"<<setw(25)<<"Eig(Ar)";
	ofile<<setw(25)<<"Time(Ar)"<<endl; 

	for (int i=0; i<5; i++)
	{
		ofile<<setprecision(6)<<Max(i);
		ofile<<setw(25)<<setprecision(6)<<time_j;
		ofile<<setw(25)<<iterations;
		ofile<<setw(25)<<eigval(i); 
		ofile<<setw(25)<<time_ar<<endl; 
	}


	ofile.close();

	delete [] V_Pot;
	delete [] V_Pin; 
	delete [] rho_i;

}//End Main Program 


void rotate(mat &A, mat &R, int k, int l, int n)
{

	//Rotation of the Matrix
        double tau;
        double t;
        double s;
        double c;

        if(A(k,l) != 0.0)
        {
          tau = (A(l, l)-A(k,k))/(2*A(k,l));

          if(tau > 0)
          {
            t = -tau + sqrt(1 + tau*tau);
          }
          else
          {
            t = -tau - sqrt(1+ tau*tau);
          }
          c = 1/sqrt(1 + t*t);
          s = t*c;
        }
        else
        {
          c =1.0;
          s= 0.0;
	}
	
	double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
        a_kk = A(k,k);
        a_ll = A(l,l);

        //Changing the matrix elemensts with indeces k an l 
        A(k, k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
        A(l, l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
        A(k,l) = 0.0;
        A(l,k) = 0.0;

        //Changing the remaining elements 
        for(int i= 0; i< (n-1); i++)
        {
           if( i != k && i != l)
           {
              a_ik = A(i,k);
              a_il = A(i,l);
              A(i,k) = a_ik*c - a_il*s;
              A(k,i) = A(i,k);
              A(i,l) = a_il*c + a_ik*s;
              A(l,i) = A(i,l);
            }

           //Compute the new eigenvectors
           r_ik = R(i,k);
           r_il = R(i,l);
           R(i,k) = c*r_ik - s*r_il;
           R(i,l) = c*r_il + s*r_ik;
	}
	return; 
}
