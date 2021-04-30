
#include<stdio.h>
#include<math.h>
#include"complex.h"
#include"green's.h"


int main()
{
    int itr, maxmitr;
    double l;
    double beta0, beta1;
    double nu0, nu1;
    double N,E;
    double a11, a12, a21, a22;
    double f1, f2;
    double det;


    double A;					            /*Lattice Constant*/
    int    klevel, ktot;  				/*Array Dimension*/
    double K;				             	/*Momentum Formula*/


    /*Lattice Parameters*/

        ktot= 200;					           /*Phononic Lattice Dimension*/
        A=1;						             /*Read value of Phononic Lattice Constant*/

        double epsilon[ktot];					        /*Initial Energy Level*/
        double n[ktot];					        /*Initial Energy Level*/

beta0 = 1.2;
nu0 = -1;

for (klevel=0; klevel<ktot; klevel=klevel+1){
    /*Momentum values at Dummy Indicies*/
           K=((-3.1418/A)+klevel*(6.2836/A)/(ktot-1.0));

    /*Electronic Dispersion Relation*/
           epsilon[klevel]= K*K/2.00000;

    /*Fermi Distribution*/
           n[klevel]= Fermi (epsilon[klevel], beta0, nu0 );

    /*Total Number Total Energy*/
           N=N+n[klevel];
           E=E+epsilon[klevel]*n[klevel];

//           printf("%d\t%f\n", klevel, n[klevel]);

}

/**********************
Beta Guess, Nu Guess
**********************/

beta0 = 1.3;
nu0 = -0.7;
l=0.000001;
maxmitr=100000;

/**********************
Newton Rhapson
**********************/
for (itr=1; itr<=maxmitr; itr++)
{

for (klevel=0; klevel<ktot; klevel=klevel+1){

a11 = a11 + FermiD1(epsilon[klevel], beta0, nu0 );
a12 = a12 + FermiD2(epsilon[klevel], beta0, nu0 );
a21 = a21 + epsilon[klevel]*FermiD1(epsilon[klevel], beta0, nu0 );
a22 = a22 + epsilon[klevel]*FermiD2(epsilon[klevel], beta0, nu0 );

f1=f1+Fermi (epsilon[klevel], beta0, nu0 );
f2=f2+epsilon[klevel]*Fermi (epsilon[klevel], beta0, nu0 );
}

det=a11*a22-a12*a21;
f1=f1-N;
f2=f2-E;  //Here is where the update comes

beta1 = beta0 - a22*f1/det+a12*f2/det;
nu1 = nu0 +a21*f1/det - a11*f2/det;

//printf("\n%f\t%f\n", beta1, nu1);

printf(" At Iteration no. %3d, beta = %9.6f, nu = %12.6f\n", itr, beta1, nu1);
if (fabs(beta1-beta0) < l)
{
if (fabs(nu1-nu0) < l)
{

    printf("After %3d iterations, beta = %8.6f, nu = %11.6f\n", itr, beta1, nu1);
    return 0;
}
}
nu0 = nu1;
beta0 = beta1;
}


    printf(" The required solution does not converge or iterations are insufficient\n");

}
