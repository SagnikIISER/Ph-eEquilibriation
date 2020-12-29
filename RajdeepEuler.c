
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
#include "Bessel.h"					/*Import Green's Functions and Baths*/
#include "green's.h"					/*Import Green's Functions and Baths*/

float main() {

/*Definitions*/

		float t, tprime;						/*time parameters*/
  	float a,b; 									/*to be used for lower and upper bound on time*/
  	float tone, ttwo, tthree;		/*extra variables*/
  	float h;										/*gap in time, "the epsilon"*/

		float A;									/*Lattice Constant*/
		int n, k, ktot;							/*array dimension*/
		float K;										/*Momentum formula*/
  	float omega;							  /*omega, k; materials for dispersion relation*/
		float lambda;								/*perturbation strenght parameter*/


/*Working around the input function call*/

	a= 0.0;										/*Read value of a*/
	b= 10.0;										/*Read value of b*/
	h= 0.1;										/*Read value of h*/
														/*Read value of omega*/
	lambda=1;									/*Read value of lambda*/

  n=(b-a)/h;							/*array Dimension*/
	ktot=20;									/*array Dimension*/

	A=1;										/*Read value of Lattice Constant*/
/*	k=1;											/*Dummy Array for momentum*/

/*	printf("the matrix dimension is:%d\n",n );   */

  tprime = a;

/*	printf("Value of time limits:%f\n",tprime ); */
/*	printf("Check for known values of Green's Functions:%f\n", DR(1.00, 2.00,3.0123) ); */


/*To Plot the input bare Green's Function and Bath*/


/*
		for ( t = a; t <= b ; t=t+h) {
				printf("%f\t%f\n",t, DzeroR(omega,t,10*h) );
				}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Defining Iteration variables and outputs
DR, BarDR and the integration I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

float DR[ktot][n][n];
float BarDR[ktot][n][n];
float IR[n][n];


float DA[ktot][n][n];
float BarDA[ktot][n][n];


float DK[ktot][n][n];
float BarDK[ktot][n][n];

float	DUMMY[n][n];
float IK[n][n];

float N[ktot][n]; 							/*Occupation Number at time t and Momentum K*/


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The Module for Dyson Iteration using Euler Method:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

for (k=0; k<ktot; k=k+1){

/*Momentum values at Dummy Indicies*/
//k=0;

K=((-3.1418/A)+k*(6.2836/A)/(ktot-1));
omega=-2*2*cos(K);

/*To Plot the Dispersion Relation*/

/*	printf("%d\t%f\t%f\n", k, K, omega );

/*Initial Condition*/


for (t=a; t<=b; t=t+h){


IR[(int)(t/h)][(int)(t/h)]=0;
DR[k][(int)(t/h)][(int)(t/h)]=0;
BarDR[k][(int)(t/h)][(int)(t/h)]= 1;

IK[(int)(t/h)][(int)(t/h)]=0;
DUMMY[(int)(t/h)][(int)(t/h)]=0;


DR[k][(int)((t+h)/h)][(int)(t/h)]= DzeroR(omega,t+h,t)*BarDR[k][(int)(t/h)][(int)(t/h)];
IR[(int)((t+h)/h)][(int)(t/h)]=0;

}

/*The Retarded Part*/


float P; P=0.0;
int i;

	for (tprime=a; tprime<=b; tprime=tprime+h){
	for (t=tprime+h; t<=b; t=t+h)
	 {BarDR[k][(int)(t/h)][(int)(tprime/h)]=  BarDzeroR(omega,t,t-h)*BarDR[k][(int)((t-h)/h)][(int)(tprime/h)]-omega*omega*DzeroR(omega,t,t-h)*DR[k][(int)((t-h)/h)][(int)(tprime/h)]+(h/2)*BarDzeroR(omega,t,t-h)*IR[(int)((t-h)/h)][(int)(tprime/h)];
		DR[k][(int)((t+h)/h)][(int)(tprime/h)]= DzeroR(omega,t+h,t)*BarDR[k][(int)(t/h)][(int)(tprime/h)]+BarDzeroR(omega,t+h,t)*DR[k][(int)(t/h)][(int)(tprime/h)]+(h/2)*DzeroR(omega,t+h,t)*IR[(int)(t/h)][(int)(tprime/h)];


				for (i = (int)((tprime+h)/h); i <= (int)((t)/h); i++)
				{
					P=P+h*SigmaR(omega,t+h,(h*i))*DR[k][i][(int)(tprime/h)];
				}
		IR[(int)((t+h)/h)][(int)(tprime/h)]=P;
		P=0.0;
	}
	}

/*The Advanced Part*/

for (tprime=a; tprime<=b; tprime=tprime+h){
for (t=tprime+h; t<=b; t=t+h){
		DA[k][(int)((tprime)/h)][(int)(t/h)]=DR[k][(int)((t)/h)][(int)(tprime/h)];
		BarDA[k][(int)((tprime)/h)][(int)(t/h)]=BarDR[k][(int)((t)/h)][(int)(tprime/h)];
}
}

/*The Keldysh Part*/



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
			P=P+h*DR[k][(int)(t/h)][i]*DUMMY[i][(int)(tprime/h)];
		}
	IK[(int)((t)/h)][(int)(tprime/h)]=P;
	P=0.0;


	DK[k][(int)((t)/h)][(int)((tprime)/h)]= DzeroK(omega,t,tprime)+IK[(int)((t)/h)][(int)(tprime/h)];

}
}


/*The Occupation Number*/

		for (t=a; t<=b; t=t+h){
		N[k][(int)((t)/h)]=0.5*(DK[k][(int)((t)/h)][(int)((t)/h)]-1);
		}

}

/*The Statistics*/

float VarE[n], EE[n], E[n], Ntot[n];
float EEPr, EPr, NPr;

EPr=0.0;
EEPr=0.0;
NPr=0.0;




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

			for ( t = a; t <=b ; t=t+h) {
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
