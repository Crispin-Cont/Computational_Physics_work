#include "/home/quetzalcoatl/Computational_Physics_work/Project3/Code/planet.h"
#include <iostream>

using namespace std; 

 planet::planet(){}

 planet::planet(double M, int step)
 {
	n=step; 
	mass=M; 

	if(M == 1.0)
	{
	   cout<<"Enter information for the SUN"<<endl;
	}
	else if(M ==1.6605e-7)
	{
	   cout<<"Enter infomation for MERCURY"<<endl; 	
	} 
	else if(M == 2.4483e-6)
        {
           cout<<"Enter information for the VENUS"<<endl;
        }
        else if(M ==3.0034e-6)
        {
           cout<<"Enter infomation for EARTH"<<endl;
        }
	 else if(M ==3.2278e-7)
        {
           cout<<"Enter infomation for MARS"<<endl;
        }
	else if(M == 9.5449e-4)
        {
           cout<<"Enter information for the JUPITER"<<endl;
        }
	else if(M == 9.5449e-1)
        {
           cout<<"Enter information for the JUPITER"<<endl;
        }
	else if(M == 9.5449e-3)
        {
           cout<<"Enter information for the JUPITER"<<endl;
        }
        else if(M == 2.8580e-4)
        {
           cout<<"Enter information for the SATURN"<<endl;
        }
        else if(M == 4.3656e-5)
        {
           cout<<"Enter infomation for URANUS"<<endl;
        }
	 else if(M ==5.1506e-5 )
        {
           cout<<"Enter infomation for NEPTUNE"<<endl;
        }
        else if(M == 6.5728e-8)
        {
           cout<<"Enter information for the PLUTO"<<endl;
        }
	
	int sepa; 
	cout<<"Do you want to input intial conditions,no(0) or yes(1): ";
	cin>>sepa;
	if(sepa == 1)
	{	
		cout<<"Initial position x_i: "; 
		cin>>pos[0]; 
		cout<<"Intial position y_i: "; 
		cin>>pos[1];
		cout<<"Initial position z_i: ";
 		cin>>pos[2]; 
		cout<<"Initial velocity Vx_i: ";
		cin>>vel[0]; 
		cout<<"Initial velocity Vy_i: "; 
		cin>>vel[1]; 
		cout<<"Initial velocity Vz_i: ";
		cin>>vel[2];
	}
	else
 	{
		if(M == 1.0)
        	{
			//Sun
                	pos[0]=-6.6275e-3;
                	pos[1]=-3.42121e-3;
                	pos[2]=1.96822e-4;
                	vel[0]=(6.12836e-6)*365;
                	vel[1]=(-6.7365e-6)*365;
                	vel[2]=(-1.17571e-7)*365;
	
        	}
	
        	if(M ==1.6605e-7)
       		{
			//Mercury
			pos[0]=-1.25329e-1;
                        pos[1]=-4.5394e-1;
                        pos[2]=-2.571179e-2;
                        vel[0]=(2.15648e-2)*365;
                        vel[1]=(-5.76037e-3)*365;               		       
			vel[2]=(-2.44891e-3)*365;
        	}
       		if(M == 2.4483e-6)
        	{
			//Venus
                        pos[0]=5.7835e-1;
                        pos[1]=-4.3512e-1;
                        pos[2]=-3.94689e-2;
                        vel[0]=(1.18909e-2)*365;
                        vel[1]=(1.61875e-2)*365;
                        vel[2]=(-4.64774e-4)*365;
        	}

       		if(M ==3.0034e-6)
        	{
			//Earth
			pos[0]=-9.9140e-1;
                        pos[1]=-1.6994e-1;
                        pos[2]=1.98187e-4;
                        vel[0]=(2.5895e-3)*365;
                        vel[1]=(-1.70385e-2)*365;
                        vel[2]=(5.1364e-7)*365;
       		}
        	else if(M ==3.2278e-7)
       		{
			//Mars
                        pos[0]=9.0233e-1;
                        pos[1]=1.1606;
                        pos[2]=2.2383e-3;
                        vel[0]=(-1.0490e-2)*365;
                        vel[1]=(9.797466e-3)*365;
                        vel[2]=(4.6327e-4)*365;
	        }
		if(M == 9.5449e-4)
        	{
           		//Jupiter
			pos[0]=3.55378;
                        pos[1]=3.47587;
                        pos[2]=-9.39615e-2;
                        vel[0]=(-5.3695e-3)*365;
                        vel[1]=(5.7509e-3)*365;
                        vel[2]=(9.6353e-5)*365;

        	}
        	if(M == 9.5449e-1)
        	{
			 //Jupiter
                        pos[0]=3.55378;
                        pos[1]=3.47587;
                        pos[2]=-9.39615e-2;
                        vel[0]=(-5.3695e-3)*365;
                        vel[1]=(5.7509e-3)*365;
                        vel[2]=(9.6353e-5)*365;
        	}
        	if(M == 9.5449e-3)
        	{
			 //Jupiter
                        pos[0]=3.55378;
                        pos[1]=3.47587;
                        pos[2]=-9.39615e-2;
                        vel[0]=(-5.3695e-3)*365;
                        vel[1]=(5.7509e-3)*365;
                        vel[2]=(9.6353e-5)*365;
        	}
        	if(M == 2.8580e-4)
        	{
			//Saturn
                        pos[0]=6.01042;
                        pos[1]=6.9007;
                        pos[2]=-3.59215e-1;
                        vel[0]=(-4.49875e-3)*365;
                        vel[1]=(3.65260e-3)*365;
                        vel[2]=(1.15655e-4)*365;
        	}
        	if(M == 4.3656e-5)
        	{
			//Uranus
                        pos[0]=1.4660e1;
                        pos[1]=-1.3499e1;
                        pos[2]=-2.40100e-1;
                        vel[0]=(2.635377e-3)*365;
                        vel[1]=(2.71037e-3)*365;
                        vel[2]=(-2.40986e-5)*365;	
        	}
         	if(M ==5.1506e-5 )
        	{
			//Neptune
                        pos[0]=1.70329e1;
                        pos[1]=-2.4836e1;
                        pos[2]=1.18927e-1;
                        vel[0]=(2.567966e-3)*365;
                        vel[1]=(1.79313e-3)*365;
                        vel[2]=(-9.6080e-5)*365;
        	}		
        	if(M == 6.5728e-8)
        	{
			//Pluto
                        pos[0]=-9.6135;
                        pos[1]=-2.80971e1;
                        pos[2]=5.78738;
                        vel[0]=(3.0528e-3)*365;
                        vel[1]=(-1.50824e-3)*365;
                        vel[2]=(-7.17669e-4)*365;
        	}


	}

	cout<<"Initial time t_i: "; 
	cin>>t_i; 
	cout<<"Final time t_f: "; 
	cin>>t_f;

 }

 double planet::Kin_E()
 {
	double Velocity=this->vel[0]*this->vel[0]+this->vel[1]*this->vel[1]+this->vel[2]*this->vel[2];
	return 0.5*this->mass*Velocity; 
 } 	
 
 double planet::Ang_M()
 {
	L[0]=this->pos[1]*this->vel[2]-this->pos[2]*this->vel[1];
	L[1]=this->pos[2]*this->vel[0]-this->pos[0]*this->vel[2];
	L[2]=this->pos[0]*this->vel[1]- this->pos[1]*this->vel[0];
	return this->mass*this->mass*(L[0]*L[0]+L[1]*L[1]+L[2]*L[2]); 
 }	
