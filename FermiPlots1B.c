/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CODE for iterating DYSON EQUATION
FERMIONIC SECTOR

The module is meant to generatre nk/k plots at vearious instants for a
given range of energy levels:

Iterations are Self-Consistent and Exact 
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
double complex GK[1010][1010];

double complex SigK[1010][1010];


double complex I2, I3, I4, I5;


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Main Module
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

double main() {


/*Definitions*/

double a,b; 					   /*to be used for lower and upper bound on time*/
double h;					       /*gap in time, "the epsilon"*/

double epsilon;					 /*Initial Energy Level*/
double lambda;					 /*Perturbation Strength*/

double Telectron, nu;


double akka;
double complex P;


/*Working around the input variable calls*/

  a= 0.000000;					/*Read value of a*/
  b= 40.000000;				  /*Read value of b*/
  h= 0.050000;					/*Read value of h*/

  lambda=1.00000;				/*Read value of lambda*/
  n=(int)((b-a)/h)+1;	 	/*array Dimension*/



  Telectron=0.000000;
  nu=1.0;



int i;
int j;
int k;
int l;


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Assigning Sigma K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

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


for (epsilon=0.01; epsilon<2.00; epsilon=epsilon+0.02) {


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Initial Condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


for (i=0; i<n; i++){
            GR[i][i]= -I;
            GR[i+1][i]= GzeroR(epsilon,(i*h)+h,(i*h));
}


GK[0][0]= -I*tanh((epsilon-nu)/(2.0*Telectron));
GK[1][0]= GzeroK(epsilon,Telectron,nu, h , 0);
GK[0][1]= -conjf(GK[1][0]);



/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The Retarded Part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


for (i=0; i<n; i++){
    	for (j=1; j<i; j++){


        for (l = j+1; l < i; l++)
        {
          P=P+h*lambda*lambda*SigmaR((i*h),(h*l))*GR[l][j];
        }

        I4=P+(h/2.0)*SigmaR((i*h),(h*i))*GR[i][j]+(h/2.0)*SigmaR((i*h),(h*j))*GR[j][j];
        P=0.0;



        for (l = j+1; l < i+1; l++)
        {
          P=P+h*lambda*lambda*SigmaR((i*h)+h,(h*l))*GR[l][j];
        }

        I2=P+(h/2.0)*SigmaR((i*h)+h,(h*j))*GR[j][j];
        P=0.0;


        GR[i+1][j]= (I*GzeroR(epsilon,(i*h)+h,(i*h))*GR[i][j]+(h/2.0)*GzeroR(epsilon,(i*h)+h,(i*h))*I4-I*(h/2.0)*I2)/(1+I*(h/2.0)*SigmaR((i*h)+h,(h*i)+h));

    	}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*The Keldysh Part*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

    for (j=0; j<i; j++){


    for (l = 1; l <= i; l++)
      {
        P=P+h*lambda*lambda*SigmaR((i*h),(h*l))*GK[l][j];
      }

      I4=P+(h/2.0)*lambda*lambda*SigmaR((i*h),(i*h))*GK[i][j]+(h/2.0)*lambda*lambda*SigmaR((i*h),(j*h))*GK[j][j];
      P=0.0;

          for (l = 1; l < j; l++)
          {
            P = P + h*lambda*lambda*SigK[i][l]*conjf(GR[j][l]);
          }

      I5= P + (h/2.0)*lambda*lambda*SigK[i][j]*conjf(GR[j][j])+ (h/2.0)*lambda*lambda*SigK[i][i]*conjf(GR[j][i]);
      P=0;




      for (l = 1; l <= i; l++)
        {
          P=P+h*lambda*lambda*SigmaR((i*h)+h,(h*l))*GK[l][j];
        }

        I2=P+(h/2.0)*lambda*lambda*SigmaR((i*h)+h,(j*h))*GK[j][j];
        P=0.0;

            for (l = 1; l < j; l++)
            {
              P = P + h*lambda*lambda*SigK[i+1][l]*conjf(GR[j][l]);
            }

        I3= P + (h/2.0)*lambda*lambda*SigK[i+1][j]*conjf(GR[j][j])+ (h/2.0)*lambda*lambda*SigK[i+1][i+1]*conjf(GR[j][i+1]);
        P=0;



      GK[i+1][j]= (I*GzeroR(epsilon,(i*h)+h,(i*h))*GK[i][j]+(h/2.0)*GzeroR(epsilon,(i*h)+h,(i*h))*(I4+I5)-(h/2.0)*I*(I2+I3))/(1+I*(h/2.0)*lambda*lambda*SigmaR((i*h)+h,(i*h)+h));
		  GK[j][i+1]=-conjf(GK[i+1][j]);
		}

  for (k=0; k<=i; k++)
  {


          for (l = 1; l <= k; l++)
          {
            P=P+h*lambda*lambda*SigmaR((k*h),(h*l))*GK[l][i+1];
          }

          I4=P+(h/2.0)*lambda*lambda*SigmaR((k*h),(i*h))*GK[i][i+1]+(h/2.0)*lambda*lambda*SigmaR((k*h),(k*h))*GK[k][i+1];
          P=0.0;

              for (l = 1; l < i; l++)
              {
                P = P + h*lambda*lambda*SigK[k][l]*conjf(GR[i+1][l]);
              }


          I5=P+(h/2.0)*lambda*lambda*SigK[k][k]*conjf(GR[i+1][k])+(h/2.0)*lambda*lambda*SigK[k][i+1]*conjf(GR[i+1][i+1]);
          P=0.0;




          for (l = 1; l <= k; l++)
          {
            P=P+h*lambda*lambda*SigmaR((k*h)+h,(h*l))*GK[l][i+1];
          }

          I2=P+(h/2.0)*lambda*lambda*SigmaR((k*h)+h,(i*h))*GK[i][i+1];
          P=0.0;

              for (l = 1; l < i; l++)
              {
                P = P + h*lambda*lambda*SigK[k+1][l]*conjf(GR[i+1][l]);
              }


          I3=P+(h/2.0)*lambda*lambda*SigK[k+1][i+1]*conjf(GR[i+1][i+1])+(h/2.0)*lambda*lambda*SigK[k+1][k+1]*conjf(GR[i+1][k+1]);
          P=0.0;


      GK[k+1][i+1]= (I*GzeroR(epsilon,(k*h)+h,(k*h))*GK[k][i+1]+(h/2.0)*GzeroR(epsilon,(k*h)+h,(k*h))*(I4+I5)+(h/2.0)*GzeroR(epsilon,(k*h)+h,(k*h)+h)*(I4+I5))/(1+I*(h/2.0)*lambda*lambda*SigmaR((k*h)+h,(k*h)+h));


  }
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Print Statements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", epsilon , 0.5*(cimagf(GK[0][0])+1), 0.5*(cimagf(GK[50][50])+1), 0.5*(cimagf(GK[100][100])+1), 0.5*(cimagf(GK[150][150])+1), 0.5*(cimagf(GK[200][200])+1), 0.5*(cimagf(GK[300][300])+1), 0.5*(cimagf(GK[400][400])+1), 0.5*(cimagf(GK[500][500])+1), 0.5*(cimagf(GK[600][600])+1), 0.5*(cimagf(GK[795][795])+1)  )  ;

}


}
