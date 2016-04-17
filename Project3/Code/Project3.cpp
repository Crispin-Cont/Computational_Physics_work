/*
	Author:Crispin Contreras
	Class: Physics 905 Computational Physics
	Purpose: Code for small scale solar system 

*/

#include <iostream>
#include <iomanip>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;
ofstream ofile; 

void initialize(double*& , double*& ,double*& , double*& , double&, double&, int&);
void Verlet(double*&, double*&, double*&, double*&, double*&, double&, double&, double&, int&);

int main ()
{
	double *x, *y, *r,*Vx, *Vy;
	double t_f, t_i;
	double h; 
	int n; //Step size

	//Allocate memory
	x = new double [n+1];
	y = new double [n+1]; 
	r = new double [n+1]; 
	Vx = new double [n+1];
	Vy = new double [n+1]; 

	
	//Call function to get the values 
	initialize(x, y, Vx, Vy,t_f, t_i,n); 
	
	//Call Function Verlet
	Verlet(x, y,Vx,Vy,r,h,t_f,t_i,n);
	
	
	for(int i=0; i<n; i++)
	{
	   cout<<"Radius :"<<r[i]<<endl;	
	}

	//Delete arrays 
	delete [] x; 
	delete [] y;
	delete [] r;  
	delete [] Vx; 
	delete [] Vy; 


} //END OF MAIN FUNCTION

//This Function prints out the prompts to get the initial velocities and positions
void initialize(double *&x , double *&y, double *&Vx, double *&Vy, double &t_f, double &t_i, int &n)
{
	cout<<"Enter the number of steps: ";
	cin>>n; 
	cout<<"Enter the Value of x_0: ";
	cin>>x[0];
	cout<<"Enter the Value of y_0: "; 
	cin>>y[0];
	cout<<"Enter the Value of Vx_0: ";
	cin>>Vx[0]; 
	cout<<"Enter the Value of Vy_0: ";
	cin>>Vy[0]; 
	cout<<"Enter the initial time t_i: ";
	cin>>t_i;
	cout<<"Enter the Final time t_f: ";
	cin>>t_f;

	return;
} 

//Function that uses the Verlet Method
void Verlet(double *&x, double *&y,double *&Vx, double *&Vy, double *&r, double &h, double &t_f, double &t_i, int &n)
{
	h = (t_f - t_i)/n; 
	r[0]=sqrt(pow(x[0],2) + pow(y[0],2));

	for(int i=0; i<n; i++)
	{
	   x[i+1] = x[i] + h*Vx[i] -x[i]*((h*h)/2.0)*((4.0*pow(M_PI,2))/(pow(r[i],3)));
	   y[i+1] = y[i] + h*Vy[i] -y[i]*((h*h)/2.0)*((4.0*pow(M_PI,2))/(pow(r[i],3)));
	   r[i+1] = sqrt(pow(x[i+1],2) + pow(y[i+1],2));
	   Vx[i+1] = Vx[i] -(h/2.0)*((x[i+1]*4.0*pow(M_PI,2))/(pow(r[i+1],3)) + x[i]*(4.0*pow(M_PI,2))/(pow(r[i],3)));
	   Vy[i+1] = Vy[i] -(h/2.0)*((y[i+1]*4.0*pow(M_PI,2))/(pow(r[i+1],3)) + y[i]*(4.0*pow(M_PI,2))/(pow(r[i],3)));
	}

}
