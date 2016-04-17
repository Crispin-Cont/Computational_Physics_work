#include "/home/quetzalcoatl/Scratch/Project3/Code/planet.h"
#include "/home/quetzalcoatl/Scratch/Project3/Code/solver.h"
#include "time.h" 
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>


using namespace std; 

 solver::solver()
 {
 	tot_p = 0; //Total planets
	total_Kinetic = 0; 
	total_Potential = 0; 
	G = 4*M_PI*M_PI; 
 }

 void solver::add(planet newplanet)
 {
 	tot_p += 1; 
	all_planets.push_back(newplanet); 
 }

 //This calculates the position of the planets using the Verlet Method
 void solver::verlet(planet &N)
 {
 	double h = (N.t_f - N.t_i)/N.n;
	double Fx,Fy,Fz,Fx1,Fy1,Fz1;
	double acc[tot_p][3];	
	double Nacc[tot_p][3];
	double rel_pos[3]

        for(int i=0; i< (N.n-1); i++)
        {
	   for(int j=0; j<tot_p; j++)
	   {	
              planet &este = all_planets[j];
	      Fx=Fy=Fz=Fx1=Fy1=Fz1=0; 
	      ForceRK4(Fx,Fy,Fz,este.pos[0],este.pos[1],este.pos[2]); 	
	      if(type>0)
	      {
		for(int l=1+j; l<tot_p;l++)
		{
		    planet &otro = all_planets[l];  	
		    for(int d =0; d<3; d++)
		    {
			rel_pos[d] = este.pos[d] - otro.pos[d]; 
		    }
		    ForceRK4(Fx,Fy,Fz,rel_pos[0],rel_pos[1],rel_pos[2]); 
		}
	
	      }
	      acc[j][0] = Fx;
	      acc[j][1] = Fy;
	      acc[j][2] = Fz; 
	      for(int l=0; l<3; l++)
	      {
		este.pos[l] += h*este.vel[l] +0.5*h*h*acc[j][l];
	      }
	   N.x_v[i+1] = N.x_v[i] + h*N.Vx_v[i] -N.x_v[i]*((h*h)/2.0)*(G/(pow(N.r_v[i],3)));
           N.r_v[i+1] = sqrt(pow(N.x_v[i+1],2) + pow(N.y_v[i+1],2));
           N.Vx_v[i+1] = N.Vx_v[i] -(h/2.0)*((N.x_v[i+1]*G)/(pow(N.r_v[i+1],3)) + N.x_v[i]*(G)/(pow(N.r_v[i],3)));
           N.Vy_v[i+1] = N.Vy_v[i] -(h/2.0)*((N.y_v[i+1]*G)/(pow(N.r_v[i+1],3)) + N.y_v[i]*(G)/(pow(N.r_v[i],3)));
	   }	
        }

 }

 //Calculates the position using the RK4 method. 
 void solver::RK4(planet &N)
 {	

        double h = (N.t_f - N.t_i)/N.n;	
	double Fx,Fy,Fz; 
	double rel_pos[3];
	double r_R[tot_p]; 
	

	//Setting up ks, the first [] is the total planets to be solved for and 
	//the second is the number of dimenstions 
	double k1_v[tot_p][3], k2_v[tot_p][3],k3_v[tot_p][3],k4_v[tot_p][3]; 
        double k1_x[tot_p][3], k2_x[tot_p][3],k3_x[tot_p][3],k4_x[tot_p][3];
	
		
	//Initializes Writting 	
	writeRK4(); 

	//Setting up the k valuesas 
	for(int i=0; i<N.n;i++)
	{
	   time=(i+1)*h; 
 	   ofile<<time; 
  	
	   //Seting up K1	
	   for(int j=0; j<tot_p; j++)
	   {
		planet &este=all_planets[j];
		Fx=Fy=Fz=0.0;
		if(type == 0)
		{	
                   ForceRK4(Fx,Fy,Fz,este.pos[0],este.pos[1],este.pos[2]);	
		}
		else
		{
		   for(int l=1+j; l<tot_p; l++)
		   {
		      planet &otro=all_planets[l];
		      for(int a=0; a<3;a++){rel_pos[a]=-(otro.pos[a]-este.pos[a]);}
	              ForceRK4(Fx,Fy,Fz,rel_pos[0],rel_pos[1],rel_pos[2]);
		   }
		}
	

	      	k1_v[j][0] = h*Fx;
	        k1_v[j][1] = h*Fy;
		k1_v[j][2] = h*Fz; 
	        for(int l= 0; l<3;l++)
		{
		   k1_x[j][l] = h*este.vel[l];
		}
	        
		 
	   }//End of loop	

	   //Setting up K2
	   for(int j=0; j<tot_p; j++)
           {
                planet &este=all_planets[j];
                Fx=Fy=Fz=0.0;
		
		if(type == 0)
		{
		   for(int a=0; a<3;a++)
                   {
                      rel_pos[a]=este.pos[a]+k1_x[j][a]/2.0;
                   }
                   ForceRK4(Fx,Fy,Fz,rel_pos[0],rel_pos[1],rel_pos[2]);

		}	
		else
		{
                   for(int l=1+j; l<tot_p; l++)
                   {   
                      planet &otro=all_planets[l];
                      for(int a=0; a<3;a++)
                      {
                         rel_pos[a]=-((otro.pos[a]+k1_x[l][a]/2.0)-(este.pos[a]+k1_x[j][a]/2.0));
                      }
                      ForceRK4(Fx,Fy,Fz,rel_pos[0],rel_pos[1],rel_pos[2]);
                   }
		}

                k2_v[j][0] = h*Fx;
                k2_v[j][1] = h*Fy;
                k2_v[j][2] = h*Fz;
                for(int l= 0; l<3;l++)
                {
                   k2_x[j][l] = h*(este.vel[l]+ k1_v[j][l]/2.0);
                }


            }//End of loop

	    //Setting up K3
	    for(int j=0; j<tot_p; j++)
            {
                planet &este=all_planets[j];
                Fx=Fy=Fz=0.0;
		
		if(type == 0)
                {
                   for(int a=0; a<3;a++)
                   {
                      rel_pos[a]=este.pos[a]+k2_x[j][a]/2.0;
                   }
                   ForceRK4(Fx,Fy,Fz,rel_pos[0],rel_pos[1],rel_pos[2]);

                }
		else
		{
                    for(int l=1+j; l<tot_p; l++)
                    {
                       planet &otro=all_planets[l];
                       for(int a=0; a<3;a++)
                       {
                          rel_pos[a]=-((otro.pos[a]+k2_x[l][a]/2.0)-(este.pos[a]+k2_x[j][a]/2.0));
                       }
                       ForceRK4(Fx,Fy,Fz,rel_pos[0],rel_pos[1],rel_pos[2]);
                    }
		}


                k3_v[j][0] = h*Fx;
                k3_v[j][1] = h*Fy;
                k3_v[j][2] = h*Fz;
                for(int l= 0; l<3;l++)
                {
                   k3_x[j][l] = h*(este.vel[l]+ k2_v[j][l]/2.0);
                }


            }//End of loop

	    //Setting up K4
	    for(int j=0; j<tot_p; j++)
            {
                planet &este=all_planets[j];
                Fx=Fy=Fz=0.0;
		if(type == 0)
                {
                   for(int a=0; a<3;a++)
                   {
                      rel_pos[a]=este.pos[a]+k3_x[j][a];
                   }
                   ForceRK4(Fx,Fy,Fz,rel_pos[0],rel_pos[1],rel_pos[2]);

                }
		else
		{
                   for(int l=1+j; l<tot_p; l++)
                   {
                      planet &otro=all_planets[l];
                      for(int a=0; a<3;a++)
                      {
                         rel_pos[a]=-((otro.pos[a]+k3_x[l][a])-(este.pos[a]+k3_x[j][a]));
                      }
                      ForceRK4(Fx,Fy,Fz,rel_pos[0],rel_pos[1],rel_pos[2]);
                   }
		}


                k4_v[j][0] = h*Fx;
                k4_v[j][1] = h*Fy;
                k4_v[j][2] = h*Fz;
                for(int l= 0; l<3;l++)
                {
                   k4_x[j][l] = h*(este.vel[l]+ k3_v[j][l]);
                }


             }//End of loop	

	     //Update the functions
	   
             //This calculates the radius of the planets
             for(int j=0;j<tot_p;j++)
             {
                planet &este =all_planets[j];
                r_R[j] = sqrt(este.pos[0]*current.pos[0] + este.pos[1]*este.pos[1] + este.pos[2]*este.pos[2]);
		
		ofile<<setw(20)<<r_R[j];
			
		for(int l=0; l<3; l++)
		{
		   este.pos[l] += (k1_x[j][l] +2*(k2_x[j][l]+k3_x[j][l]) +k4_x[j][l])/6.0;
		   ofile<<setw(20)<<este.pos[l]<<endl;
		   este.vel[l] += (k1_v[j][l] +2*(k2_v[j][l]+k3_v[j][l]) +k4_v[j][l])/6.0;
		}
             }
   	}

	ofile.close();
 }

 
 void solver::ForceRK4(double &Fx, double &Fy,double &Fz, double x, double y, double z)
 {
	double r= sqrt(x*x+y*y+z*z); 

	Fx -= (G/(r*r*r))*x; 
	Fy -= (G/(r*r*r))*y; 
	Fz -= (G/(r*r*r))*z;
 }	
 
/*
 void solver::ForceVe(double &Fx,double &Fy, double &Fz)
 {
	Fx -= (G/(r*r*r))*x; 
	Fy -= (G/(r*r*r))*y; 
	Fz -= (G/(r*r*r))*z; 

 }
 */
 void solver::writeRK4(planet &N)
 {
	ofstream ofile;
        char RK4[30];

        cout<<"Enter File Name: ";
        cin>>RK4;
	ofile.open(RK4);
	ofile<<setiosflags(ios::showpoint | ios::uppercase);
	ofile<<"#Time(Yr)";
	
	if(type == 0)
	{	
	   ofile<<setw(20)<<"X_Ea"<<setw(20)<<"Y_Ea"<<setw(20)<<"Z_Ea"<<setw(20)<<"R_earth"<<endl;
	}
	else if(type == 1)
	{
	   ofile<<setw(20)<<"X_Ea"<<setw(20)<<"Y_Ea"<<setw(20)<<"Z_Ea"<<setw(20)<<"R_Earth";
	   ofile<<setw(20)<<"X_Ju"<<setw(20)<<"Y_Ju"<<setw(20)<<"Z_Ju"<<setw(20)<<"R_Jupiter"<<endl;	
	}
	else if(type == 2)
	{
	   ofile<<setw(20)<<"X_Su"<<setw(20)<<"Y_Su"<<setw(20)<<"Z_Su"<<setw(20)<<"R_Sun";
	   ofile<<setw(20)<<"X_Ea"<<setw(20)<<"Y_Ea"<<setw(20)<<"Z_Ea"<<setw(20)<<"R_Earth";
           ofile<<setw(20)<<"X_Ju"<<setw(20)<<"Y_Ju"<<setw(20)<<"Z_Ju"<<setw(20)<<"R_Jupiter"<<endl;
	}
	else if(type == 3)
	{
	   ofile<<setw(20)<<"X_Me"<<setw(20)<<"Y_Me"<<setw(20)<<"Z_Me"<<setw(20)<<"R_Mercury";
           ofile<<setw(20)<<"X_Ve"<<setw(20)<<"Y_Ve"<<setw(20)<<"Z_Ve"<<setw(20)<<"R_Venus";
           ofile<<setw(20)<<"X_Ea"<<setw(20)<<"Y_Ea"<<setw(20)<<"Z_Ea"<<setw(20)<<"R_Earth";
	   ofile<<seww(20)<<"X_Ma"<<setw(20)<<"Y_Ma"<<setw(20)<<"Z_Ma"<<setw(20)<<"R_Mars";
           ofile<<setw(20)<<"X_Ju"<<setw(20)<<"Y_Ju"<<setw(20)<<"Z_Ju"<<setw(20)<<"R_Jupiter";
	   ofile<<setw(20)<<"X_Sa"<<setw(20)<<"Y_Sa"<<setw(20)<<"Z_Sa"<<setw(20)<<"R_Saturn";
           ofile<<setw(20)<<"X_Ur"<<setw(20)<<"Y_Ur"<<setw(20)<<"Z_Ur"<<setw(20)<<"R_Uranus";
           ofile<<setw(20)<<"X_Ne"<<setw(20)<<"Y_Ne"<<setw(20)<<"Z_Ne"<<setw(20)<<"R_Neptune";
           ofile<<setw(20)<<"X_Pu"<<setw(20)<<"Y_Pu"<<setw(20)<<"Z_Pu"<<setw(20)<<"R_Pluto"<<endl;
	
	}	
 }

 void solver::writeVe(planet &N)
 {

        ofstream ofile;
        double time;
        double h = (N.t_f - N.t_i)/N.n;
        char Ve[30];


        cout<<"Enter File Name: ";
        cin>>Ve;

        ofile.open(Ve);
        ofile<<setiosflags(ios::showpoint | ios::uppercase);
        ofile<<"#Time(Yr)"<<setw(20)<<"X(au)"<<setw(20)<<"Y(au)"<<setw(20)<<"r_earth"<<endl;
        for(int i =0; i<= N.n; i++)
        {
           time=i*h;
           ofile<<time;
           if(i< N.n)
           {
              ofile<<setw(20)<<N.x_v[i];
              ofile<<setw(20)<<N.y_v[i];
	      ofile<<setw(20)<<N.r_v[i]<<endl;
           }
        }
        ofile.close();


 }

