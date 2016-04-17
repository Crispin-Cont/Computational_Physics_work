#ifndef SOLVER_H
#define SOLVER_H
#include <cmath>
#include "/home/quetzalcoatl/Scratch/Project3/Code/planet.h"
#include <vector>
using std::vector; 

class solver
{

 public:
 	friend class planet; 

	//Properties 
	int tot_p; //Total planets
	double total_Kinetic; 
	double total_Potential;
	vector<planet> all_planets; 
	double G; 

	//Initializers 
	solver();

	//Functions
	void add(planet newplanet); //adds new planet
	void verlet(planet &N);  
	void RK4(planet &N);
	void ForceRK4(double &Fx, double &Fy,double&Fz, double x, double y,double z);
	void writeRK4(planet &N); 
	void writeVe(planet &N); 
}; 
#endif
