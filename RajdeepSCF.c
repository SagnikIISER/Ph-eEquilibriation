
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
#include "green's.h"					/*Import Green's Functions and Baths*/

float main() {

/*Definitions*/

		float t, tprime;						/*time parameters*/
  	float a,b; 									/*to be used for lower and upper bound on time*/
  	float tone, ttwo, tthree;		/*extra variables*/
  	float h;										/*gap in time, "the epsilon"*/

		int n;											/*array dimension*/
		float omega, k;							/*omega, k; materials for dispersion relation*/
		float lambda;								/*perturbation strenght parameter*/


/*Working around the input function call*/

	a= 0.0;										/*Read value of a*/
	b= 10;									/*Read value of b*/
	h= 0.1;										/*Read value of h*/
	omega= 1;									/*Read value of omega*/
	lambda=1;									/*Read value of lambda*/

  n=(b-a)/h;								/*array Dimension*/

/*	printf("the matrix dimension is:%d\n",n );        */

  tprime = a;

/*	printf("Value of time limits:%f\n",tprime ); */
/*	printf("Check for known values of Green's Functions:%f\n", DR(1.00, 2.00,3.0123) ); */


/*To Plot the input bare Green's Function and Bath*/


	/*
	for ( t = a; t <= b ; t=t+h) {
		printf("%f\t%f\n",t, DzeroR(omega,t,tprime) );
	}

	*/

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Defining Iteration variables and outputs
DR, BarDR and the integration I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

float DR[n][n];
float BarDR[n][n];
float I[n][n];

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The Module for Dyson Iteration using Euler Method:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*Initial Condition*/

I[0][0]=0;
DR[0][0]=0;
BarDR[0][0]= 1;

DR[(int)(h/h)][0]= DzeroR(omega,h,0)*BarDR[0][0];
I[(int)(h/h)][0]=0;


float P; P=0.0;
int i;

	for (t=a+h; t<=b; t=t+h)
	 {BarDR[(int)(t/h)][(int)(a/h)]=(DR[(int)(t/h)][(int)(a/h)]-DR[(int)((t-h)/h)][(int)(a/h)])/h;
		DR[(int)((t+h)/h)][(int)(a/h)]= DzeroR(omega,t+h,t)*BarDR[(int)(t/h)][(int)(a/h)]+BarDzeroR(omega,t+h,t)*DR[(int)(t/h)][(int)(a/h)]+(h/2)*DzeroR(omega,t+h,t)*I[(int)(t/h)][(int)(a/h)];


				for (i = (int)((a+h)/h); i <= (int)((t)/h); i++)
				{
					P=P+h*SigmaR(omega,t+h,(h*i))*DR[i][(int)(a/h)];
				}
		I[(int)((t+h)/h)][(int)(a/h)]=P;
		P=0.0;
	}

	/*To Plot the input bare Green's Function and Bath*/


		for ( t = a; t <= b ; t=t+h) {
			printf("%f\t%f\n",t, DR[(int)(t/h)][(int)(a/h)] );
		}


}
