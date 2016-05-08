#ifndef PLANET_H
#define PLANET_H
#include <cmath>
#define _USE_MATH_DEFINES

class planet 
{
  	
	
 public:
 	//Properties
	int n; 
 	double mass,t_i,t_f;  
 	double pos[3];
	double vel[3];  
	double L[3];
	
 	//Initializers
 	planet();
 	planet(double M,int step);

	//Functions
 	double Kin_E();
	double Ang_M(); 
};
#endif
