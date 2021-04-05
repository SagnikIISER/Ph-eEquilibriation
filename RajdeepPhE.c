
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CODE for SOLVING DYSON EQUATION

First note down all the relevant formula plus diagram
Only true unknown in the whole buisness is the D and BarD
There are two modules written for the Retarded and Keldysh Component:
			a. Euler
			b. Self Consistent Field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include <stdio.h> 		      /*Import Standard Module*/
#include <stdlib.h> 		    /*Import Standard Library*/
#include <math.h>		        /*Import Math Module*/
#include <complex.h>        /*Import Modules for Complex Numbers*/
#include "Dawson.h"		      /*Import Module for Dawson Function*/
#include "green's.h"	     	/*Import Green's Functions and Baths*/


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Defining Iteration variables and Outputs
Greens Functions, Self Energies and Increments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


/*%%%%%%%%%%%%%%%%%%%%%%%%
Electron Sector
**************************/

double complex GR[2510][2510];
double complex IER[2510][2510];

double complex GA[2510][2510];

double complex GK[2510][2510];

double complex IEK1[2510][2510];
double complex IEK2[2510][2510];

double complex SigElR[2510][2510];
double complex SigElK[2510][2510];


/*%%%%%%%%%%%%%%%%%%%%%%%%%
Phonon Sector
**************************/

double complex DR[2510][2510];
double complex BarDR[2510][2510];
double complex IPR[2510][2510];

double complex DA[2510][2510];

double complex DK[2510][2510];
double complex BarDK[2510][2510];

double complex IPK1[2510][2510];
double complex IPK2[2510][2510];

double complex SigPhR[2510][2510];
double complex SigPhK[2510][2510];


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Main Module
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

double main() {



/*Definitions*/

  double t, tprime;				       /*time parameters*/
  double a,b; 					         /*to be used for lower and upper bound on time*/
  double tone, ttwo;				     /*extra variables*/
  double h;					             /*time increment*/

  double A;					             /*Lattice Constant*/
  int    n, klevel, ktot;				 /*Array Dimension*/
  double K;					             /*Momentum Formula*/
  double omega, epsilon;				 /*Initial Phononic, Electronic Energy Levels respectively*/
  double lambda;					       /*Perturbation Strength*/

  double Tphonon, Telectron;     /*Phonon and Electron Temp*/
  double nu;                     /*Chemical Potential*/
  double complex P;

    int i;
    int j;
    int l;

/*Working around the input variable calls*/

    a= 0.000000;					       /*Read value of a*/
    b= 250.000000;				       /*Read value of b*/
    h= 0.100000;					       /*Read value of h*/

    lambda=1.00000;				       /*Read value of lambda*/
    n=(int)((b-a)/h)+1;			     /*array Dimension*/

/*Phononic Lattice Parameters*/

    ktot= 19;					           /*Phononic Lattice Dimension*/
    A=1;						             /*Read value of Phononic Lattice Constant*/
    klevel=1; 			    		     /*Dummy Array for Phononic momentum*/



/*Reading the Initial Temp*/
    Tphonon=  1.000000;
    Telectron=1.000000;


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Energy Levels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*Momentum values at Dummy Indicies*/
   K=((-3.1418/A)+klevel*(6.2836/A)/(ktot-1.0));

/*Phononic Dispersion Relation*/
   omega=3.0*sqrt(sin((K*A/2.0))*sin((K*A/2.0)));

/*Electronic Dispersion Relation*/
   epsilon= K*K/2.00000;





/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Initial Condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


for (i=0; i<n; i++){


  /*Electronic Sector*/

        IER[i][i]= 0;
        GR[i][i]= -I;
        GK[i][i]= I*tanh((epsilon-nu)/(2.0*Telectron));

        DR[i+1][i]= GzeroR(epsilon,t+h,t);
        IER[i+1][i]= 0;

        GK[i+1][i]= GzeroK(epsilon,Telectron,nu,t+h,t);
        GK[i][i+1]=-conjf(DK[i+1][i]);


        IK1[i][i]=0;
        IK2[i][i]=0;
        //IK1[i+1][i]=(h/2.0)*lambda*lambda*SigmaR((i*h)+h,(i*h))*DK[i][i];
        //IK2[i+1][i]=(h/2.0)*lambda*lambda*SigK[i+1][i+1]*DR[i+1][i];

/*Phononic Sector*/

        IPR[i][i]=0;
        DR[i][i]=0;
        BarDR[i][i]= -1/2.0;

        DK[i][i]=-(I/(2.0*omega))*(1.0/(tanh((omega)/(2*Tsyst))));
        BarDK[i][i]= 0;

        DR[i+1][i]= -2.0*DzeroR(omega,t+h,t)*BarDR[i][i];
        IPR[i+1][i]= 0;

        DK[i+1][i]= -2.0*BarDzeroR(omega,(h*i)+h,(h*i))*DK[i][i];
        DK[i][i+1]= -conjf(DK[i+1][i]);


        IPK1[i][i]=0;
        IPK2[i][i]=0;
        //IPK1[i+1][i]=(h/2.0)*SigmaRP((i*h)+h,(i*h))*DK[i][i];
        //IPK2[i+1][i]=(h/2.0)*SigK[i+1][i+1]*DR[i+1][i];


   }

}
