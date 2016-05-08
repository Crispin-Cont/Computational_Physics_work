#include "/home/quetzalcoatl/Computational_Physics_work/Project3/Code/planet.h"
#include "/home/quetzalcoatl/Computational_Physics_work/Project3/Code/solver.h"
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
	G = 4*M_PI*M_PI; 
 }

 void solver::add(planet newplanet)
 {
 	tot_p += 1; 
	all_planets.push_back(newplanet); 
 }

 //This calculates the position of the planets using the Verlet Method
 void solver::verlet(planet &N,int type)
 {
 	double h = (N.t_f - N.t_i)/N.n;
	double Fx,Fy,Fz,Fx1,Fy1,Fz1;
	double acc[tot_p][3];	
	double Nacc[tot_p][3];
	double rel_pos[3];
	double time; 
	
	
	char Ve[30];
	char Energy2[30]; //contains the Kinetic energy, potential energy, and angular mom. 
        cout<<"Enter the name of the Verlet file: ";
        cin>>Ve;
	cout <<"Enter the name of the Energy file: ";
	cin>>Energy2; 	

	std::ofstream ofile(Ve); 
	std::ofstream E_output(Energy2); 
	Header_Pos(type,ofile); 
	Header_Energy(type,E_output); 

	//Write initial values to the file 
	time= 0.0; 
	Write_Pos(ofile,time);
	Write_Energy(E_output,time,type); 

	//Start the clock 
	clock_t start_VV,finish_VV;  
	start_VV = clock(); 

        for(int i=0; i< N.n; i++)
        {  
	   time=(i+1)*h;
	   for(int j=0; j<tot_p; j++)
	   {	
              planet &este = all_planets[j];
	      Fx=Fy=Fz=Fx1=Fy1=Fz1=0; 
	      Force(Fx,Fy,Fz,este.pos[0],este.pos[1],este.pos[2],1); 	
	      if(type>0)
	      {
		for(int l=0; l<tot_p;l++)
		{  
		    if(l == j)
		    {
			Fx+=0;
			Fy+=0; 
			Fz+=0; 
	            }
		    else
		    {
		    	planet &otro = all_planets[l];  	
		    	for(int d =0; d<3; d++)
		    	{
				rel_pos[d] = este.pos[d] - otro.pos[d]; 
		    	}
		    	Force(Fx,Fy,Fz,rel_pos[0],rel_pos[1],rel_pos[2],otro.mass);
		    } 
		}
	
	      }
	      acc[j][0] = Fx;
	      acc[j][1] = Fy;
	      acc[j][2] = Fz; 
	      
	      //Update the position . 
	      for(int l=0; l<3; l++)
	      {
		este.pos[l] += h*este.vel[l] +0.5*h*h*acc[j][l];
	      }

	      //Values for the velocity 
	      Force(Fx1,Fy1,Fz1,este.pos[0],este.pos[1],este.pos[2],1);
	      if(type>0)
	      {	
	      	for(int l=0; l<tot_p;l++)
                {
                    if(l == j)
                    {
                        Fx1+=0;
                        Fy1+=0;
			Fz1+=0; 
                    }
                    else
                    {
                        planet &otro = all_planets[l];
                        for(int d =0; d<3; d++)
                        {
                                rel_pos[d] = este.pos[d] - otro.pos[d];
                        }
                        Force(Fx1,Fy1,Fz1,rel_pos[0],rel_pos[1],rel_pos[2],otro.mass);
                    }
                }
	       }
	      
	      Nacc[j][0] = Fx1; 
	      Nacc[j][1] = Fy1;
	      Nacc[j][2] = Fz1;

	      //Calculate new velocity 
	      for(int l=0; l<3; l++)
	      {
		este.vel[l] += 0.5*h*(acc[j][l] + Nacc[j][l]);
	      }

	   } 
	  //Write the Updated the values 
           Write_Pos(ofile,time);
           Write_Energy(E_output,time,type);

        }

	//stop clock and display time
	finish_VV=clock(); 
	cout<<"Total Time VV: "<<((double)(finish_VV - start_VV)/CLOCKS_PER_SEC)<<endl; 


	ofile.close();
	E_output.close();

 }

 //Calculates the position using the RK4 method. 
 void solver::RK4(planet &N,int type)
 {	
        double h = (N.t_f - N.t_i)/N.n;	
	double Fx,Fy,Fz; 
	double rel_pos[3];
 	double time;	

	//Setting up ks, the first [] is the total planets to be solved for and 
	//the second is the number of dimenstions 
	double k1_v[tot_p][3], k2_v[tot_p][3],k3_v[tot_p][3],k4_v[tot_p][3]; 
        double k1_x[tot_p][3], k2_x[tot_p][3],k3_x[tot_p][3],k4_x[tot_p][3];
	
		
	//Initializes Writting 
	char RK4[30];
        char Energy[30]; //contains the Kinetic energy, potential energy, and angular mom. 
        cout<<"Enter the name of the RK4 file: ";
        cin>>RK4;
        cout <<"Enter the name of the Energy file: ";
        cin>>Energy;

        std::ofstream ofile(RK4);
        std::ofstream E_output(Energy);
        Header_Pos(type,ofile);
        Header_Energy(type,E_output);

	//Write initial values to the file 
        time= 0.0;
        Write_Pos(ofile,time);
        Write_Energy(E_output,time,type);

	//Start Clock
	clock_t start_RK, finish_RK; 
	start_RK = clock();

	//Setting up the k values 
	for(int i=0; i<N.n;i++)
	{
	   time=(i+1)*h; 
  	
	   //Seting up K1	
	   for(int j=0; j<tot_p; j++)
	   {
		planet &este=all_planets[j];
		Fx=Fy=0.0;
				
                Force(Fx,Fy,Fz,este.pos[0],este.pos[1],este.pos[2],1);	
		
		if(type >0)
		{
		   for(int l=0; l<tot_p; l++)
		   {
		      if( j == l)
		      {
			Fx+=0.0;
			Fy+=0.0; 
			Fz+=0.0;   
		      }
		      else
		      {
		      	planet &otro=all_planets[l];
		      	for(int a=0; a<3;a++){rel_pos[a]=-(otro.pos[a]-este.pos[a]);}
	              	Force(Fx,Fy,Fz,rel_pos[0],rel_pos[1],rel_pos[2],otro.mass);
		      }
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
                Fx=Fy=0.0;
		
		for(int a=0; a<3;a++)
                {
                   rel_pos[a]=este.pos[a]+k1_x[j][a]/2.0;
                }
                Force(Fx,Fy,Fz,rel_pos[0],rel_pos[1],rel_pos[2],1);

		if(type >0)
		{
                   for(int l=0; l<tot_p; l++)
                   {   
		      if(j == l)
		      {
			Fx+=0.0;
			Fy+=0.0; 
			Fz+=0.0; 
		      }
		      else 
		      {
                      	planet &otro=all_planets[l];
                      	for(int a=0; a<3;a++)
                      	{
                        	rel_pos[a]=-((otro.pos[a]+k1_x[l][a]/2.0)-(este.pos[a]+k1_x[j][a]/2.0));
                      	}
                      	Force(Fx,Fy,Fz,rel_pos[0],rel_pos[1],rel_pos[2],otro.mass);
		      }
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
                Fx=Fy=0.0;
		
                
                for(int a=0; a<3;a++)
                {
                	rel_pos[a]=este.pos[a]+k2_x[j][a]/2.0;
                }
                Force(Fx,Fy,Fz,rel_pos[0],rel_pos[1],rel_pos[2],1);
		if(type>0)
		{
                    for(int l=0; l<tot_p; l++)
                    {
		       if(j == l)
		       {
			  Fx+=0.0; 
			  Fy+=0.0; 
			  Fz+=0.0;
		       }
		       else
		       {
                       	  planet &otro=all_planets[l];
                       	  for(int a=0; a<3;a++)
                       	  {
                         	 rel_pos[a]=-((otro.pos[a]+k2_x[l][a]/2.0)-(este.pos[a]+k2_x[j][a]/2.0));
                      	  }
                       	  Force(Fx,Fy,Fz,rel_pos[0],rel_pos[1],rel_pos[2],otro.mass);
		       }
                    }
		}

                k3_v[j][0] = h*Fx;
                k3_v[j][1] = h*Fy;
		k3_v[j][2] = h*Fz;
                for(int l= 0; l<2;l++)
                {
                   k3_x[j][l] = h*(este.vel[l]+ k2_v[j][l]/2.0);
                }
            }//End of loop

	    //Setting up K4
	    for(int j=0; j<tot_p; j++)
            {
                planet &este=all_planets[j];
                Fx=Fy=Fz=0.0;
               
                for(int a=0; a<3;a++)
                {
                   rel_pos[a]=este.pos[a]+k3_x[j][a];
                }
                Force(Fx,Fy,Fz,rel_pos[0],rel_pos[1],rel_pos[2],1);
		if(type>0)
		{
                   for(int l=0; l<tot_p; l++)
                   {
		      if(j == l)
		      {
			Fx+=0.0; 
			Fy+=0.0; 
			Fz+=0.0;
	              }
		      else
	 	      {	
                      	planet &otro=all_planets[l];
                      	for(int a=0; a<3;a++)
                      	{
                        	 rel_pos[a]=-((otro.pos[a]+k3_x[l][a])-(este.pos[a]+k3_x[j][a]));
                      	}
                      	Force(Fx,Fy,Fz,rel_pos[0],rel_pos[1],rel_pos[2],otro.mass);
		      }
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

	     //This updates the functions
             for(int j=0;j<tot_p;j++)
             {
                planet &este =all_planets[j];
		
		for(int l=0; l<3; l++)
		{
		   este.pos[l] += (k1_x[j][l] +2*(k2_x[j][l]+k3_x[j][l]) +k4_x[j][l])/6.0;
		   este.vel[l] += (k1_v[j][l] +2*(k2_v[j][l]+k3_v[j][l]) +k4_v[j][l])/6.0;
		}
             }
		
       	   //Write updated values to the file 
           Write_Pos(ofile,time);
           Write_Energy(E_output,time,type);
	
   	}

	finish_RK = clock(); 
	cout<<"Total Time RK: "<<((double)(finish_RK -start_RK)/CLOCKS_PER_SEC)<<endl; 

	ofile.close();
	E_output.close(); 
	
 }
 
 void solver::Force(double &Fx, double &Fy,double &Fz, double x, double y,double z, double m)
 {
	double r= sqrt(x*x+y*y+z*z); 

	Fx -= (G/(r*r*r))*x*m; 
	Fy -= (G/(r*r*r))*y*m; 	
	Fz -= (G/(r*r*r))*z*m;
 }
	
 void solver::Header_Pos(int type,std::ofstream& ofile)
 {

	ofile<<setiosflags(ios::showpoint | ios::uppercase);
	ofile<<"#Time(Yr)";
	
	if(type == 0)
	{	
	   ofile<<setw(20)<<"X_Ea"<<setw(20)<<"Y_Ea"<<setw(20)<<"Z_Ea";
	   ofile<<setw(20)<<"R_earth"<<endl;
	}
	else if(type == 1)
	{
	   ofile<<setw(20)<<"X_Ea"<<setw(20)<<"Y_Ea"<<setw(20)<<"Z_Ea";
	   ofile<<setw(20)<<"R_Earth";
	   ofile<<setw(20)<<"X_Ju"<<setw(20)<<"Y_Ju"<<setw(20)<<"Z_Ju";
	   ofile<<setw(20)<<"R_Jupiter"<<endl;	
	}
	else if(type == 2)
	{
	   ofile<<setw(20)<<"X_Su"<<setw(20)<<"Y_Su"<<setw(20)<<"Z_Su";
	   ofile<<setw(20)<<"R_Sun";
	   ofile<<setw(20)<<"X_Ea"<<setw(20)<<"Y_Ea"<<setw(20)<<"Z_Ea";
	   ofile<<setw(20)<<"R_Earth";
           ofile<<setw(20)<<"X_Ju"<<setw(20)<<"Y_Ju"<<setw(20)<<"Z_Ju";
	   ofile<<setw(20)<<"R_Jupiter"<<endl;
	}
	else if(type == 3)
	{
	   ofile<<setw(20)<<"X_Me"<<setw(20)<<"Y_Me"<<setw(20)<<"Z_Me";
	   ofile<<setw(20)<<"R_Mercury";
           ofile<<setw(20)<<"X_Ve"<<setw(20)<<"Y_Ve"<<setw(20)<<"Z_Ve";
	   ofile<<setw(20)<<"R_Venus";
           ofile<<setw(20)<<"X_Ea"<<setw(20)<<"Y_Ea"<<setw(20)<<"Z_Ea";
	   ofile<<setw(20)<<"R_Earth";
	   ofile<<setw(20)<<"X_Ma"<<setw(20)<<"Y_Ma"<<setw(20)<<"Z_Ma";
 	   ofile<<setw(20)<<"R_Mars";
           ofile<<setw(20)<<"X_Ju"<<setw(20)<<"Y_Ju"<<setw(20)<<"Z_Ju";
	   ofile<<setw(20)<<"R_Jupiter";
	   ofile<<setw(20)<<"X_Sa"<<setw(20)<<"Y_Sa"<<setw(20)<<"Z_Sa";
	   ofile<<setw(20)<<"R_Saturn";
           ofile<<setw(20)<<"X_Ur"<<setw(20)<<"Y_Ur"<<setw(20)<<"Z_Ur";
	   ofile<<setw(20)<<"R_Uranus";
           ofile<<setw(20)<<"X_Ne"<<setw(20)<<"Y_Ne"<<setw(20)<<"Z_Ne";
	   ofile<<setw(20)<<"R_Neptune";
           ofile<<setw(20)<<"X_Pu"<<setw(20)<<"Y_Pu"<<setw(20)<<"Z_Pu";
	   ofile<<setw(20)<<"R_Pluto"<<endl;
	}	
 }

void solver::Header_Energy(int type,std::ofstream& ofile)
 {
        ofile<<setiosflags(ios::showpoint | ios::uppercase);
        ofile<<"#Time(Yr)";

        if(type == 0)
        {
           ofile<<setw(20)<<"K_Ea"<<setw(20)<<"U_Ea"<<setw(20)<<"Etot_Ea";
	   ofile<<setw(20)<<"L_Ea"<<endl;
        }
        else if(type == 1)
        {
           ofile<<setw(20)<<"K_Ea"<<setw(20)<<"U_Ea"<<setw(20)<<"Etot_Ea";
	   ofile<<setw(20)<<"L_Ea";
           ofile<<setw(20)<<"K_Ju"<<setw(20)<<"U_Ju"<<setw(20)<<"Etot_Ju";
	   ofile<<setw(20)<<"L_Ju"<<endl;
        }
        else if(type == 2)
        {
           ofile<<setw(20)<<"K_Su"<<setw(20)<<"U_Su"<<setw(20)<<"Etot_Su";
	   ofile<<setw(20)<<"L_Su";
           ofile<<setw(20)<<"K_Ea"<<setw(20)<<"U_Ea"<<setw(20)<<"Etot_Ea";
	   ofile<<setw(20)<<"L_Ea";
           ofile<<setw(20)<<"K_Ju"<<setw(20)<<"U_Ju"<<setw(20)<<"Etot_Ju";
	   ofile<<setw(20)<<"L_Ju"<<endl;
        }
        else if(type == 3)
        {
           ofile<<setw(20)<<"K_Me"<<setw(20)<<"U_Me"<<setw(20)<<"Etot_Me"; 
	   ofile<<setw(20)<<"L_Me";
           ofile<<setw(20)<<"K_Ve"<<setw(20)<<"U_Ve"<<setw(20)<<"Etot_Ve";
	   ofile<<setw(20)<<"L_Ve";
           ofile<<setw(20)<<"K_Ea"<<setw(20)<<"U_Ea"<<setw(20)<<"Etot_Ea";
	   ofile<<setw(20)<<"L_Ea";
           ofile<<setw(20)<<"K_Ma"<<setw(20)<<"U_Ma"<<setw(20)<<"Etot_Ma";
	   ofile<<setw(20)<<"L_Ma";
           ofile<<setw(20)<<"K_Ju"<<setw(20)<<"U_Ju"<<setw(20)<<"Etot_Ju";
	   ofile<<setw(20)<<"L_Ju";
           ofile<<setw(20)<<"K_Sa"<<setw(20)<<"U_Sa"<<setw(20)<<"Etot_Sa";
	   ofile<<setw(20)<<"L_Sa";
           ofile<<setw(20)<<"K_Ur"<<setw(20)<<"U_Ur"<<setw(20)<<"Etot_Ur";
	   ofile<<setw(20)<<"L_Ur";
           ofile<<setw(20)<<"K_Ne"<<setw(20)<<"U_Ne"<<setw(20)<<"Etot_Ne";
	   ofile<<setw(20)<<"L_Ne";
           ofile<<setw(20)<<"K_Pl"<<setw(20)<<"U_Pl"<<setw(20)<<"Etot_pl";
	   ofile<<setw(20)<<"L_Pl"<<endl;
        }
 }

 void solver::Write_Pos(std::ofstream& ofile, double time)  
 {
     double R[tot_p]; 
   
     ofile<<time; 
     for(int i=0; i<tot_p; i++)
     {
	 planet &este = all_planets[i];
	 for(int j=0; j<3; j++)
	 {
	    ofile<<setw(20)<<este.pos[j];  
	 }

	 R[i] = sqrt(este.pos[0]*este.pos[0] + este.pos[1]*este.pos[1]+este.pos[2]*este.pos[2]);
         ofile<<setw(20)<<R[i];	
     }
     ofile<<endl; 
 }

 void solver::Write_Energy(std::ofstream& ofile, double time,int type)
 {
	double pot;
	double Etot; 
	double rel_p[3]; 
	ofile<<time; 

	for(int i=0; i<tot_p; i++)
	{
	   pot= 0.0;
	   Etot= 0.0;  
	   planet &este = all_planets[i];

	   //Gets the Kinetic Energy 
	   ofile<<setw(20)<<este.Kin_E();  

	   //Calcualted the Potential energy 
	   Potential(pot,este.pos[0],este.pos[1],este.pos[2],este.mass, 1); 
	   
	   if(type == 0)
	   {
		ofile<<setw(20)<<pot; 
	   }

	   if(type > 0)
	   {
		for(int j=0; j<tot_p; j++)
		{
		    planet &otro = all_planets[j]; 
	 	  
		    if(i == j)
		    {
			pot+=0.0; 
                    }  
		    else 
		    {
			for(int l=0; l<3; l++)
			{
			    rel_p[l] = este.pos[l] - otro.pos[l];
			}
			
			Potential(pot,rel_p[0],rel_p[1],rel_p[2],este.mass, otro.mass); 
	  	    }
		}
		ofile<<setw(20)<<pot;
	   }	
	   
           //Writes total energy 
	   Etot = este.Kin_E() + pot; 
	   ofile<<setw(20)<<Etot; 	
		
	   //Calculated the angular momentum 
           ofile<<setw(20)<<setprecision(7)<<este.Ang_M();
	}
	ofile<<endl; 
 }

 void solver::Potential(double &pot, double x, double y,double z, double m1, double m2)
 {
	double r = sqrt(x*x+y*y+z*z);
	pot -= (G*m1*m2)/r;
 }

