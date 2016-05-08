#ifndef SOLVER_H
#define SOLVER_H
#include <cmath>
#include "/home/quetzalcoatl/Computational_Physics_work/Project3/Code/planet.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

using std::vector; 

class solver
{

 public:
 	friend class planet; 

	//Properties 
	int tot_p; //Total planets
	vector<planet> all_planets; 
	double G; 
	
	//Initializers 
	solver();

	//Functions
	void add(planet newplanet); //adds new planet
	void verlet(planet &N,int type);  
	void RK4(planet &N,int type);
	void Force(double &Fx, double &Fy,double &Fz, double x, double y,double z,double m);
	void Header_Pos(int type,std::ofstream& ofile); 
	void Header_Energy(int type, std::ofstream& ofile); 
	void Write_Pos(std::ofstream& ofile, double time); 
	void Write_Energy(std::ofstream& ofile, double time, int type); 
	void Potential(double &pot, double x, double y,double z, double m1, double m2); 

}; 
#endif
