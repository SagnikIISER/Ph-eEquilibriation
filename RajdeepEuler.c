
/*CODE for SOLVING DYSON EQUATION
First note down all the relevant formula plus diagram
Only true unknown in the whole buisness is the D and BarD
work out the equation for Keldysh Component
There are two modules written for the Retarded Component
			a. Euler
			b. Self Consistent Field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include <stdio.h> 						/*Import Standard Module*/
#include <math.h>							/*Import Math Module*/
#include <complex.h>          /*Import Modules for Complex Numbers*/
#include "Bessel.h"
#include "green's.h"				  /*Import Green's Functions and Baths*/



float main() {

/*Definitions*/

		float t, tprime;					/*time parameters*/
  	float a,b; 								/*to be used for lower and upper bound on time*/
  	float tone, ttwo, tthree;	/*extra variables*/
  	float h;									/*gap in time, "the epsilon"*/

		float A;									/*Lattice Constant*/
		int 	n, k, ktot;					/*array dimension*/
		float K;									/*Momentum formula*/
  	float omega;							/*omega, k; materials for dispersion relation*/
		float lambda;							/*perturbation strenght parameter*/


/*Working around the input function call*/

	a= 0.0;										  /*Read value of a*/
	b= 35.0;										/*Read value of b*/
	h= 0.1;									   	/*Read value of h*/

	lambda=1;									  /*Read value of lambda*/

  n=(b-a)/h+1;							    /*array Dimension*/
	ktot=5;									  /*array Dimension*/

	A=1;										    /*Read value of Lattice Constant*/
  k=1;											  /*Dummy Array for momentum*/

  tprime = a;



/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Defining Iteration variables and outputs
DR, BarDR and the integration I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

float complex DR[n][n];
float complex BarDR[n][n];
float complex IR[n][n];

/*
float complex DA[ktot][n][n];
float complex BarDA[ktot][n][n];


float complex DK[ktot][n][n];
float complex BarDK[ktot][n][n];

float	complex DUMMY[n][n];
float complex IK[n][n];

float N[ktot][n]; 							/*Occupation Number at time t and Momentum K*/


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The Module for Dyson Iteration using Euler Method:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


for (k=0; k<ktot; k=k+1){

/*Momentum values at Dummy Indicies*/
K=((-3.1418/A)+k*(6.2836/A)/(ktot-1));

/*Dispersion Relation*/
omega=2*sqrt(sin((K*A/2))*sin((K*A/2)));

/*To Plot the Dispersion Relation*/
//	printf("%f\t%f\n", K, omega );




/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*Initial Condition*/


for (t=a; t<=b; t=t+h){


IR[(int)(t/h)][(int)(t/h)]=0;
DR[(int)(t/h)][(int)(t/h)]=0;
BarDR[(int)(t/h)][(int)(t/h)]= 1;

/*
IK[(int)(t/h)][(int)(t/h)]=0;
DUMMY[(int)(t/h)][(int)(t/h)]=0;
*/

DR[(int)((t+h)/h)][(int)(t/h)]= DzeroR(omega,t+h,t)*I*BarDR[(int)(t/h)][(int)(t/h)];
IR[(int)((t+h)/h)][(int)(t/h)]=0;

}

/*The Retarded Part*/


float complex P; P=0.0;
int i;

for (t=tprime+h; t<=b; t=t+h)
	for (tprime=a; tprime<=b; tprime=tprime+h){
	 {BarDR[(int)(t/h)][(int)(tprime/h)]=  BarDzeroR(omega,t,t-h)*BarDR[(int)((t-h)/h)][(int)(tprime/h)]-omega*omega*DzeroR(omega,t,t-h)*DR[(int)((t-h)/h)][(int)(tprime/h)]+(h/2)*BarDzeroR(omega,t,t-h)*IR[(int)((t-h)/h)][(int)(tprime/h)];
		DR[(int)((t+h)/h)][(int)(tprime/h)]= DzeroR(omega,t+h,t)*BarDR[(int)(t/h)][(int)(tprime/h)]+BarDzeroR(omega,t+h,t)*DR[(int)(t/h)][(int)(tprime/h)]+(h/2)*DzeroR(omega,t+h,t)*IR[(int)(t/h)][(int)(tprime/h)];


				for (i = (int)((tprime+h)/h); i <= (int)((t)/h); i++)
				{
					P=P+h*SigmaR(omega,t+h,(h*i))*DR[i][(int)(tprime/h)];
				}
		IR[(int)((t+h)/h)][(int)(tprime/h)]=P;
		P=0.0;
	}
	}

/*The Advanced Part*/
/*
for (tprime=a; tprime<=b; tprime=tprime+h){
for (t=tprime+h; t<=b; t=t+h){
		DA[k][(int)((tprime)/h)][(int)(t/h)]=DR[(int)((t)/h)][(int)(tprime/h)];
		BarDA[k][(int)((tprime)/h)][(int)(t/h)]=BarDR[(int)((t)/h)][(int)(tprime/h)];
}
}

/*The Keldysh Part*/

/*

for ( t = a; t <=b ; t=t+h) {
for ( tprime = a; tprime <=b ; tprime=tprime+h) {

	for (i = (int)((a+h)/h); i <= (int)((t-h)/h); i++)
		{
			P=P+h*SigmaK(omega,t,(h*i))*DA[k][i][(int)(tprime/h)];
		}
	DUMMY[(int)((t)/h)][(int)(tprime/h)]=P;
	P=0.0;

	for (i = (int)((a+h)/h); i <= (int)((t-h)/h); i++)
		{
			P=P+h*DR[(int)(t/h)][i]*DUMMY[i][(int)(tprime/h)];
		}
	IK[(int)((t)/h)][(int)(tprime/h)]=P;
	P=0.0;


	DK[k][(int)((t)/h)][(int)((tprime)/h)]= DzeroK(omega,t,tprime)+IK[(int)((t)/h)][(int)(tprime/h)];

}
}


/*The Occupation Number*/
/*
		for (t=a; t<=b; t=t+h){
		N[k][(int)((t)/h)]=0.5*(DK[k][(int)((t)/h)][(int)((t)/h)]-1);
		}

}

/*The Statistics*/
/*
float VarE[n], EE[n], E[n], Ntot[n];
float EEPr, EPr, NPr;

EPr=0.0;
EEPr=0.0;
NPr=0.0;

, cimagf(DR[0][(int)(t/h)][0])


		for ( t = a; t <=b ; t=t+h) {
		for ( k = 0; k <ktot ; k=k+1) {
					NPr=NPr+N[k][(int)((t)/h)];
				}
				Ntot[(int)((t)/h)]=NPr;
				NPr=0.0;
			}


			for ( t = a; t <=b ; t=t+h) {
			for ( k = 0; k <ktot ; k=k+1) {
				K=((-3.1418/A)+k*(6.2836/A)/(ktot-1));
				omega=-2*2*cos(K);
				EPr=EPr+(omega)*N[k][(int)((t)/h)];
					}
					E[(int)((t)/h)]=EPr;
					EPr=0.0;
				}

			for ( t = a; t DzeroR(<=b ; t=t+h) {
			for ( k = 0; k <ktot ; k=k+1) {
						EEPr=EEPr+(omega)*(omega)*N[k][(int)((t)/h)];
						}
						EE[(int)((t)/h)]=EEPr;
						EEPr=0.0;
				}

			for ( t = a+h; t <=b ; t=t+h) {
				printf("%f\t%f\n", t, (E[(int)((t)/h)]/Ntot[(int)((t)/h)]) );
			}


//}

/*	for (t=a; t<=b; t=t+h){
			 printf("%f\t%f\n", t, DK[0][(int)((t)/h)][(int)((t)/h)] );
		 }
*/
}
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K loop ends here.
Put K independent printf statements after this.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

for (t=a; t<=b; t=t+h){
printf("%f\t%f\n", t ,cimagf(DR[(int)(t/h)][0]) );
}
}
