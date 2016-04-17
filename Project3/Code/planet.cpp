#include "/home/quetzalcoatl/Scratch/Project3/Code/planet.h"
#include <iostream>

using namespace std; 

 planet::planet()	
 {}
 

 planet::planet(double M, int step)
 {
	n = step;
	mass = M; 

	if( M ==1.0)
	{
	   cout<<"Enter information for the SUN"<<endl;
	}
	else if( M== 1.2e-7)
	{
	   cout<<"Enter infomation for MERCURY"<<endl; 
	}
        else if( M== 2.4e-6)
        {
           cout<<"Enter infomation for VENUS"<<endl;
        }
	else if( M ==1.5e-6)
        {
           cout<<"Enter information for the EARTH"<<endl;
        }
        else if( M== 3.3e-7)
        {
           cout<<"Enter infomation for MARS"<<endl;
        }
        else if( M==9.5e-4)
        {
           cout<<"Enter infomation for JUPITER"<<endl;
        }
	else if( M==2.75e-4)
        {
           cout<<"Enter infomation for SATURN"<<endl;
        }
        else if( M==4.4e-5)
        {
           cout<<"Enter infomation for URANUS"<<endl;
        }
        else if( M ==5.1e-5)
        {
           cout<<"Enter information for the NEPTUNE"<<endl;
        }
        else if( M==5.6e-9)
        {
           cout<<"Enter infomation for PLUTO"<<endl;
        }
       
	cout<<"Intial positon x_i: ";
	cin>>pos[0];
        cout<<"Intial position y_i: ";
        cin>>pos[1];
	cout<<"Intial position z_i: ";
        cin>>pos[2];
        cout<<"Intial velocity Vx_i: ";
        cin>>vel[0];
        cout<<"Inital velocity Vy_i: ";
        cin>>vel[1];
	cout<<"Inital velocity Vz_i: ";
        cin>>vel[2];
	cout<<"Intial time t_i: "; 
	cin>>t_i; 
	cout<<"Final time t_f: "; 
	cin>>t_f; 	

 }
 
 planet::~planet()
 {}

