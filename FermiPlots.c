/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CODE for SOLVING DYSON EQUATION

FI1st note down all the relevant formula plus diagram
Only true unknown in the whole buisness is the D and BarD
There are two modules written for the Retarded and Keldysh Component:
			a. Euler
			b. Self Consistent Field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include <stdio.h>        /*Import Standard Module*/
#include <stdlib.h>       /*Import Standard Library*/
#include <omp.h>
#include <math.h>		      /*Import Math Module*/
#include <complex.h>      /*Import Modules for Complex Numbers*/
#include "Dawson.h"		    /*Import Module for Dawson Function*/
#include "green's.h"		  /*Import Green's Functions and Baths*/


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Defining Iteration variables and outputs
DR, BarDR and the integration I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

double complex GR[1010][1010];
double complex BarGR[1010][1010];

double complex SigK[1010][1010];
double complex GK[1010][1010];
double complex BarGK[1010][1010];
double complex DKthermal;

double complex I1, I2, I4, I5;
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Main Module
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

double main() {


/*Definitions*/

double t, tprime;				/*time parameters*/
double a,b; 					/*to be used for lower and upper bound on time*/
double tone, ttwo;				/*extra variables*/
double h;					/*gap in time, "the epsilon"*/

double A;					/*Lattice Constant*/
int    n, klevel, ktot;				/*Array Dimension*/
double K;					/*Momentum Formula*/
double epsilon;					/*Initial Energy Level*/
double lambda;					/*Perturbation Strength*/
double nomega;          /*Density of states*/

double Tbath, Telectron, nu;


double akka;
double complex P;


/*Working around the input variable calls*/

  a= 0.000000;					/*Read value of a*/
  b= 40.000000;				 /*Read value of b*/
  h= 0.050000;					/*Read value of h*/

  lambda=1.00000;				/*Read value of lambda*/
  n=(int)((b-a)/h)+1;			    	/*array Dimension*/
  ktot= 19;					/*array Dimension*/

  A= 1.0;						/*Read value of Lattice Constant*/
  klevel= 0; 			    		/*Dummy Array for momentum*/

  Telectron=0.000000;
  nu=0.5;

int i;
int j;
int k;
int l;


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The Module for Dyson Iteration using Euler Method:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*
for (klevel=0; klevel<ktot/2; klevel=klevel+1){

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Momentum values at Dummy Indicies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

 K=((-1/A)+klevel*(2/A)/(ktot-1.0));

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
System Dispersion Relation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/* epsilon = K*K/2;
 //nomega= 1/(exp(omega/Tsyst)-1);


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 To Plot the Dispersion Relation
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*
printf("%d\t%f\t%f\t%f\n", klevel, K, omega, nomega );
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Assigning Sigma K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

for (epsilon=0.01; epsilon<1.00; epsilon=epsilon+0.02) {


	for (i=0; i<n; i++){
    SigK[i][0]=SigmaK(i*h,0);
  }

  for (i=0; i<n; i++){
    SigK[0][i]=SigmaK(0,i*h);
  }

  for (j=1; j<n; j++){
    for (i=j; i<n; i++){
        SigK[i][j]=SigK[i-j][0];
    }
  }


  for (j=1; j<n; j++){
    for (i=j; i<n; i++){
        SigK[j][i]=SigK[0][i-j];
    }
  }

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Initial Condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/



for (i=0; i<n; i++){


            GR[i][i]= I;
            GK[i][i]= I*tanh((epsilon-nu)/(2.0*Telectron));

            GR[i+1][i]= GzeroR(epsilon,(i*h)+h,(i*h));

            GK[i+1][i]= GzeroK(epsilon,Telectron,nu,(i*h)+h,(i*h));
            GK[i][i+1]= -conjf(GK[i+1][i]);

}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The Retarded Part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


for (i=0; i<n; i++){
    	for (j=1; j<i; j++){
        for (l = j+1; l < i; l++)
        {
          P=P+h*lambda*lambda*SigmaR((i*h),(h*l))*GR[l][j];
        }

        I4=P;
        P=0.0;

        GR[i+1][j]= I*GzeroR(epsilon,(i*h)+h,(i*h))*GR[i][j]+(h/2.0)*GzeroR(epsilon,(i*h)+h,(i*h))*I4;//+(h/2.0)*GzeroR(epsilon,(i*h)+h,(i*h)+h)*I2;

    	}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*The Keldysh Part*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

    for (j=0; j<i; j++){


    for (l = 1; l <= i; l++)
      {
        P=P+h*lambda*lambda*SigmaR((i*h),(h*l))*GK[l][j];
      }

      I4=P+(h/2.0)*lambda*lambda*SigmaR((i*h),(j*h))*GK[j][j];;
      P=0.0;

          for (l = 1; l < j; l++)
          {
            P = P + h*lambda*lambda*SigK[i][l]*conjf(GR[j][l]);
          }

      I5= P + (h/2.0)*lambda*lambda*SigK[i][i]*conjf(GR[j][i]);
      P=0;

      GK[i+1][j]= I*GzeroR(epsilon,(i*h)+h,(i*h))*GK[i][j]+(h/2.0)*GzeroR(epsilon,(i*h)+h,(i*h))*(I4+I5);//+(h/2.0)*GzeroR(epsilon,(i*h)+h,(i*h)+h)*(I1+I2);
		  GK[j][i+1]=-conjf(GK[i+1][j]);
		}

  for (k=0; k<=i; k++)
  {
    for (l = 1; l <= k; l++)
    {
      P=P+h*lambda*lambda*SigmaR((k*h),(h*l))*GK[l][i+1];
    }

    I4=P+(h/2.0)*lambda*lambda*SigmaR((k*h),(i*h))*GK[i][i];
    P=0.0;

        for (l = 1; l < i; l++)
        {
          P = P + h*lambda*lambda*SigK[k][l]*conjf(GR[i+1][l]);
        }


    I5=P+(h/2.0)*lambda*lambda*SigK[k][k]*conjf(GR[i+1][k]);
    P=0.0;

    GK[k+1][i+1]= I*GzeroR(epsilon,(k*h)+h,(k*h))*GK[k][i+1]+(h/2.0)*GzeroR(epsilon,(k*h)+h,(k*h))*(I4+I5);//+(h/2.0)*GzeroR(epsilon,(k*h)+h,(k*h)+h)*(I4+I5);


  }
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K loop ends here.
Put K independent printf statements after this.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
//printf("\n%f\n%f\n%f", omega, crealf(DKthermal), -cimagf(DK[n][n])-crealf(DKthermal));

/*
for (i = 0; i < 151; i=i+25){
printf("%f\t%f\n", 1.00000 , 0.5*(cimagf(GK[i][i])+1) )  ;
}


for (i = 200; i < n; i=i+50){
printf("%f\t%f\n", 1.00000 , 0.5*(cimagf(GK[i][i])+1) )  ;
}
*/

printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", epsilon , 0.5*(cimagf(GK[0][0])+1), 0.5*(cimagf(GK[50][50])+1), 0.5*(cimagf(GK[100][100])+1), 0.5*(cimagf(GK[150][150])+1), 0.5*(cimagf(GK[200][200])+1), 0.5*(cimagf(GK[300][300])+1), 0.5*(cimagf(GK[400][400])+1), 0.5*(cimagf(GK[500][500])+1), 0.5*(cimagf(GK[600][600])+1), 0.5*(cimagf(GK[795][795])+1)  )  ;

}


}
