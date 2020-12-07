
/*CODE for SOLVING DYSON EQUATION

First note down all the relevant formula plus diagram
Only true unknown in the whole buisness is the D and BarD
work out the equation for Keldysh Component
There are two modules written for the Retarded Component
			a. Euler
			b. Self Consistent Field

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include <stdio.h> 						/*Import Standard Module*/
#include "math.h"							/*Import Math Module*/

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

	b=0;											/*Read value of a*/
	a=10;											/*Read value of b*/
	h=0.1;										/*Read value of h*/

	omega=1;									/*Read value of omega*/
	lambda=1;									/*Read value of lambda*/

  n=(a-b)/h;								/*array Dimension*/

	printf("%d\n",n );

return 0.0;
}
