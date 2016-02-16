/***************************************************************
***** Author: Crispin Contreras	
***** Class:  Physics 905
***** Purpose: solve tridiagonal matrix with LU decomposition (Armadillo) and my own method.
****************************************************************/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include "time.h"
#include "armadillo"

using namespace arma; 
using namespace std;
ofstream ofile;

int main()
{
        //DECLARE VARIABLES     
        int size_matrix;
	double size_h; 
        double *diagonal_element, *known_fun, *numerical_f, *x;
	double *analytic_f, *epsilon;
	char outfilename[30];

        //READ INPUT 
        cout<<"Enter the size of the matrix: ";
        cin>>size_matrix;
	cout<<"Enter the name of the output file: ";
	cin>>outfilename;


        //ALLOCATE MEMORY 
        size_h = 1.0/(1.0+ size_matrix);
        diagonal_element = new double[size_matrix];
        known_fun = new double[size_matrix+2];
	numerical_f = new double[size_matrix + 2];
	x = new double [size_matrix + 2];
	analytic_f = new double [size_matrix + 2];



	//CALCULATE h^2*f(X) AND x[i]
	for (int i=0; i<= (size_matrix+1); i++)
	{
		//Define x
		x[i]=i*size_h;
		
		//Calculate analytic function		 
                analytic_f[i] = 1.0 -(1-exp(-10))*x[i] - exp(-10*x[i]);

		//Calculating source of charge h^2*f
		if( (i>0) && (i< (size_matrix+1)))
		{
			known_fun [i-1] = size_h*size_h*(100.0*exp(-10.0*x[i]));
		}
	}


	//Clock
	clock_t start, finish;
	start = clock();


	//FORWARD SUBSTITUTION
	diagonal_element[0] = 2.0; // Defines First Diagonal element	
	for(int i=0; i <= ( size_matrix+1); i++)
        {
                //Calculate diagonal terms
		if(i < size_matrix)
		{
                	diagonal_element[i+1] = 2.0-(1.0/diagonal_element[i]);
	
		}//end if statement

	
	        //Calculate known functions
                if( (i>0) && (i<size_matrix))
		{
			known_fun[i] += known_fun[i-1]/diagonal_element[i-1] ;
	
		}
		
	
	}//End of Loop
	

	//BACKWARD SUBSTITUTION 
        numerical_f[size_matrix]= known_fun[size_matrix-1]/diagonal_element[size_matrix -1];
	numerical_f[0]=0.0;
	numerical_f[size_matrix+1]=0;
	
	
 	for(int m = (size_matrix-1); m > 0; m-- )
	{	
			numerical_f[m] = (known_fun[m-1] + numerical_f[m+1])/diagonal_element[m-1];

	}//End of Loop

	finish = clock(); // End clock 
	double elapsed_t = (double(finish-start)/CLOCKS_PER_SEC); 



	//OPENING, WRITING TO FILE, AND FINDIND THE RELATIVE ERROR
        ofile.open(outfilename);
        ofile<<setiosflags(ios::showpoint | ios::uppercase);
        ofile<<"#    x:			Analytic:		Numerical:		Epsilon Largest:		log(h)"<<endl;
	
	//Define Variable to find relative error
        epsilon = new double [size_matrix+2];
	epsilon[0] = epsilon[size_matrix+2]=0;//set error here to 0 since it's undefined	
	
	for(int i=0; i<=(size_matrix+1); i++)
	{
		if(i>0)
		{
                epsilon[i]=log10(fabs((numerical_f[i]-analytic_f[i])/analytic_f[i]));
		}

		//writing to file
		ofile<<setw(15)<<setprecision(8)<< x[i];
                ofile<<setw(20)<<setprecision(8)<<analytic_f[i]; 
		ofile<<setw(25)<<setprecision(8)<<numerical_f[i];
                ofile<<setw(28)<<setprecision(8)<<epsilon[i];
		ofile<<setw(29)<<setprecision(8)<<log10(size_h)<<endl;
	}//End Of Loop
		

	//USING LU DECOMPOSITION WITH ARMADILLO	
	mat A(size_matrix, size_matrix);
	mat L(size_matrix, size_matrix);
	mat U(size_matrix, size_matrix);
	colvec V(size_matrix);	
	colvec  W(size_matrix);

	
	A.zeros(size_matrix, size_matrix); //Filling the matrix with zeros
	A(0, 0)=2.0;
	W(0) =  size_h*size_h*(100.0*exp(-10.0*size_h));


	//filling the rest of the matrix	
	for(int i=1; i<size_matrix; i++)
	{
		A(i, i) = 2.0;
		A(i,i-1)= A(i-1, i) = -1;
		W(i) = size_h*size_h*(100.0*exp(-10.0*size_h*(i+1)));
	}	
	
	
	//Clock
        start = clock();

	lu(L,U,A);
	colvec Y = solve(L, W);
	V = solve(U, Y);
	
	finish = clock(); // End clock 
        double elap_tLU = (double(finish-start)/CLOCKS_PER_SEC);
	
	cout << setiosflags(ios::showpoint | ios::uppercase);
	cout << setprecision(50) << "Elapsed time mine: "<<elapsed_t<<" LU time: "<<elap_tLU<<endl
	ofile << "# LU decomp (used to verify if the results are the same)" << endl;
	for(int i=0; i<size_matrix; i++)
	{
                epsilon[i]=log10(fabs((V[i]-analytic_f[i+1])/analytic_f[i+1]));
		//writing to file
		ofile<<setw(7)<<setprecision(8)<<"# "<<x[i+1];
                ofile<<setw(10)<<setprecision(8)<<"# "<<analytic_f[i+1]; 
		ofile<<setw(15)<<setprecision(8)<<"# "<<V[i];
                ofile<<setw(20)<<setprecision(8)<<"# "<<epsilon[i]<< endl;
	}//End Of Loop


	ofile.close();

        //free memory
        delete [] diagonal_element;
        delete [] known_fun;
	delete [] x; 
        delete [] numerical_f;
	delete [] analytic_f;
	delete [] epsilon;
        return 0;

}//End main program 



