#include <iostream>
#include "/home/quetzalcoatl/Computational_Physics_work/Project3/Code/planet.h"
#include "/home/quetzalcoatl/Computational_Physics_work/Project3/Code/solver.h"
#include <cmath>

using namespace std; 

int main()
{
  solver systemRK;
  solver systemVV; 
  int type; 
  int n; 

  /*Chose what system will be solved. 
    type = 0 solves the Earth-Sun system with the sun stationary
    type =1 solves the Earth-Sun-Jupiter system with the sun stationary
    type = 2 solves the Earth-Sun-Jupiter system with the sun, all with repect to the center of mass
    type =3 solves the whole solar sytem with the sun stationary			
  */ 			

  cout<<"Which system will you like to be solve, type = ";
  cin>>type; 
  cout<<"Enter the step size n: "; 
  cin>>n; 

  
  planet Earth(3.0034e-6, n);

  if(type == 0)
  {			
	systemRK.add(Earth);
	systemVV.add(Earth); 
  } 
  else if(type ==1)
  {
	systemRK.add(Earth);
	systemVV.add(Earth);
	int typeJ; 
 	cout<<"Mass of Jupiter; 0:normal, 1:10times, 2:E3times: ";
        cin>>typeJ; 
	if(typeJ == 0)
        {
	   planet Jupiter(9.5449e-4,100); 
	   systemRK.add(Jupiter);
	   systemVV.add(Jupiter); 
        }
	if(typeJ == 1)
        {
	   planet Jupiter(9.5449e-3,100);
           systemRK.add(Jupiter);
           systemVV.add(Jupiter);
	} 
	if(typeJ == 2)
  	{
 	   planet Jupiter(9.5449e-1,100);
           systemRK.add(Jupiter);
           systemVV.add(Jupiter);
	}
  
  }
  else if(type == 2)
  {
	planet Sun(1,100);
	systemRK.add(Sun);
	systemVV.add(Sun); 
        systemRK.add(Earth);
	systemVV.add(Earth); 
        planet Jupiter(9.5449e-4,100);
        systemRK.add(Jupiter);
	systemVV.add(Jupiter); 
  }
  else if(type == 3)
  {
	planet Mercury(1.6605e-7,100);
	systemRK.add(Mercury);
	systemVV.add(Mercury); 
	planet Venus(2.4483e-6,100);
	systemRK.add(Venus);
	systemVV.add(Venus); 
        systemRK.add(Earth);
	systemVV.add(Earth); 
	planet Mars(3.2278e-7,100);
	systemRK.add(Mars); 
	systemVV.add(Mars); 
        planet Jupiter(9.5449e-4,100);
        systemRK.add(Jupiter);
	systemVV.add(Jupiter); 
	planet Saturn(2.8580e-4,100); 
	systemRK.add(Saturn);
	systemVV.add(Saturn); 
	planet Uranus(4.3656e-5,100);
	systemRK.add(Uranus); 
	systemVV.add(Uranus); 
	planet Neptune(5.1506e-5,100);
	systemRK.add(Neptune);
	systemVV.add(Neptune); 
	planet Pluto(6.5728e-8,100 );
	systemRK.add(Pluto);
	systemVV.add(Pluto); 
  }

  systemRK.RK4(Earth,type); 
  systemVV.verlet(Earth,type);
}
