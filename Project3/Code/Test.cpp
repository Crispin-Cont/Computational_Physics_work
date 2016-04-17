#include <iostream>
#include "/home/quetzalcoatl/Scratch/Project3/Code/planet.h"
#include "/home/quetzalcoatl/Scratch/Project3/Code/solver.h"
#include <cmath>

using namespace std; 

int main()
{
  solver system;
  int type; 

  /*Chose what system will be solved. 
    type = 0 solves the Earth-Sun system with the sun stationary
    type =1 solves the Earth-Sun-Jupiter system with the sun stationary
    type = 2 solves the Earth-Sun-Jupiter system with the sun, all with repect to the center of mass
    type =3 solves the whole solar sytem with the sun stationary			
  */ 			

  cout<<"Which system will you like to be solve, type =: "
  cin>>type; 
  

  if(type == 0)
  {			
  	planet Earth(3e-6, 100);  
	system.add(Earth);
  } 
  else if(type ==1)
  {
	planet Earth(3.0034e-6, 100); 
	system.add(Earth);
	planet Jupiter(9.5e-4,100); 
	system.add(Jupiter);
  }
  else if(type == 2)
  {
	planet Sun(1,100);
	system.add(Sun);
 	planet Earth(3.0034e-6, 100);
        system.add(Earth);
        planet Jupiter(9.5e-4,100);
        system.add(Jupiter);
  }
  else if(type == 3)
  {
	planet Mercury(1.6605e-7,100);
	system.add(Mercury);
	planet Venus(2.4483e-6,100);
	system.add(Venus);
        planet Earth(3.0034e-6, 100);
        system.add(Earth);
	planet Mars(3.2278e-7,100);
	system.add(Mars); 
        planet Jupiter(9.5449e-4,100);
        system.add(Jupiter);
	planet Saturn(2.8580e-4,100); 
	system.add(Saturn);
	planet Uranus(4.3656e-5,100);
	system.add(Uranus); 
	planet Neptune(5.1506e-5,100);
	system.add(Neptune);
	planet Pluto(6.5728e-8,100 );
	system.add(Pluto);
  }

	
  system.RK4(Earth); 
  system.verlet(Earth);
  system.writeRK4(Earth);
  system.writeVe(Earth); 	

}
