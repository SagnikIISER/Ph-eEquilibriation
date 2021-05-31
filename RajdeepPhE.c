
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CODE for SOLVING DYSON EQUATION

FI1st note down all the relevant formula plus diagram
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

double complex GR[6][1010][1010];
double complex GK[6][1010][1010];
double complex SigElR[6][1010][1010];
double complex SigElK[6][1010][1010];


/*%%%%%%%%%%%%%%%%%%%%%%%%%
Phonon Sector
**************************/

double complex DR[6][1010][1010];
double complex BarDR[6][1010][1010];

double complex DK[6][1010][1010];
double complex BarDK[6][1010][1010];

double complex SigPhR[6][1010][1010];
double complex SigPhK[6][1010][1010];


double complex GKtt1[30];
double complex GKtt[30];

double nel[30][1010];
double nph[30][1010];

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Energy Spectra and Other Variables
**********************************/

double epsilon[30];
double omega[30];

double Beta[1010];
double Nu[1010];

double complex IPh1A[1010], IPh1B[1010], IPh2A[1010], IPh2B[1010], IPh3A[1010], IPh3B[1010], I2[1010], I3[1010];
double complex IEl1A[1010], IEl1B[1010], IEl2A[1010], IEl2B[1010], IEl3A[1010], IEl3B[1010];


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Main Module
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

double main() {



  /*Definitions*/

    double a,b; 					         /*to be used for lower and upper bound on time*/
    double h;					             /*time increment*/
    double A;					             /*Lattice Constant*/
    int    n, klevel, ktot;				 /*Array Dimension*/
    int    plevel;		        		 /*Array Dimension*/
    double K;					             /*Momentum Formula*/
    double lambda;					       /*Perturbation Strength*/

    double Tphonon, Telectron;     /*Phonon and Electron Temp*/
    double nu;                     /*Chemical Potential*/
    double complex P;


      int i;
      int j;
      int k;
      int l;

  /*Working around the input variable calls*/

      a= 0.000000;					       /*Read value of a*/
      b= 50.000000;				       /*Read value of b*/
      h= 0.100000;					       /*Read value of h*/

      lambda=1.00000;				       /*Read value of lambda*/
      n=(int)((b-a)/h)+1;			     /*array Dimension*/

  /*Phononic Lattice Parameters*/

      ktot= 1;					           /*Phononic Lattice Dimension*/
      A=1;						             /*Read value of Phononic Lattice Constant*/
  //  klevel=1; 			    		     /*Dummy Array for Phononic momentum*/



  /*Reading the Initial Temp*/
      Tphonon=  1.000000;
      Telectron=1.000000;
      nu = -1.0000;

  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Energy Levels
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

  for (klevel=0; klevel<ktot+1; klevel++) {

  /*Momentum values at Dummy Indicies*/
     K=((-3.1418/A)+klevel*(6.2836/A)/(ktot));

  /*Phononic Dispersion Relation*/
     omega[klevel]=3.5*sqrt(sin((K*A/2.0))*sin((K*A/2.0)));

  /*Electronic Dispersion Relation*/
     epsilon[klevel]= K*K/2.00000;

  }

  /*
  for (klevel=0; klevel< ktot+1; klevel++) {
  printf("%d\t%f\t%f\n", klevel, omega[klevel], epsilon[klevel] );
  }

  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Initial Condition
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


  for (i=0; i<n; i++){
  for (klevel=0; klevel< ktot+1; klevel++) {

    /*Electronic Sector*/

          GR[klevel][i][i]= -I;
          GK[klevel][i][i]= I*tanh((epsilon[klevel]-nu)/(2.0*Telectron));

          GR[klevel][i+1][i]= GzeroR(epsilon[klevel],(i*h)+h,(i*h));

          GK[klevel][i+1][i]= GzeroK(epsilon[klevel],Telectron,nu,(i*h)+h,(i*h));
          GK[klevel][i][i+1]=-conjf(GK[klevel][i+1][i]);


  /*Phononic Sector*/

  DR[klevel][i][i]=0;
  BarDR[klevel][i][i]= -1/2.0;

  DK[klevel][i][i]=-(I/(2.0*omega[klevel]))*(1.0/(tanh((omega[klevel])/(2*Tphonon))));
  BarDK[klevel][i][i]= 0;

  DR[klevel][i+1][i]= -2.0*DzeroR(omega[klevel],t+h,t)*BarDR[klevel][i][i];

  DK[klevel][i+1][i]= -2.0*BarDzeroR(omega[klevel],(h*i)+h,(h*i))*DK[klevel][i][i];
  DK[klevel][i][i+1]=-conjf([klevel]DK[i+1][i]);

    }
  }


  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Dyson Iteration
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

  for (i=1; i<n; i++){


    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Phononic Self Energies
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


        for (j=1; j<i; j++){
          for (klevel=0; klevel< ktot+1; klevel++) {
          for (plevel=0; plevel< ktot+1; plevel++) {

              if (klevel+plevel<=ktot)
                 {
                   P = P-lambda*lambda*I*(GR[klevel+plevel][i][j]*GK[plevel][j][i]+GK[klevel+plevel][i][j]*conjf(GR[plevel][i][j]));
                   }
              if (klevel+plevel>ktot)
                 {
                    P = P-lambda*lambda*I*(GR[klevel+plevel-ktot][i][j]*GK[plevel][j][i]+GK[klevel+plevel-ktot][i][j]*conjf(GR[plevel][i][j]));
                   }
                   }
              SigPhR[klevel][i][j]=-P;
              P=0;
                 }
                }


        for (j=1; j<i; j++){
           for (klevel=0; klevel< ktot+1; klevel++) {
           for (plevel=0; plevel< ktot+1; plevel++) {

               if (klevel+plevel<=ktot){
                  if (i>j){
                       P = P-lambda*lambda*I*(GK[klevel+plevel][i][j]*GK[plevel][j][i]+GR[klevel+plevel][i][j]*conjf(GR[plevel][i][j]));
                           }
                  else{
                       P = P-lambda*lambda*I*(GK[klevel+plevel][i][j]*GK[plevel][j][i]+GR[plevel][i][j]*conjf(GR[klevel+plevel][i][j]));
                           }
                        }
               if (klevel+plevel>ktot){
                  if (i>j){
                       P = P-lambda*lambda*I*(GK[klevel+plevel-ktot][i][j]*GK[plevel][j][i]+GR[klevel+plevel-ktot][i][j]*conjf(GR[plevel][i][j]));
                          }
                  else{
                       P = P-lambda*lambda*I*(GK[klevel+plevel-ktot][i][j]*GK[plevel][j][i]+GR[plevel][i][j]*conjf(GR[klevel+plevel-ktot][i][j])); 
                          }
                        }
                        }
               SigPhK[klevel][i][j]=P;
               P=0;
                 }
              }


         for (k=1; k<=i; k++){
            for (klevel=0; klevel< ktot+1; klevel++) {
            for (plevel=0; plevel< ktot+1; plevel++) {

              if (klevel+plevel<=ktot{
                  if (k>i){
                      P = P-lambda*lambda*I*(GK[klevel+plevel][k][i]*GK[plevel][i][k]+GR[klevel+plevel][k][i]*conjf(GR[plevel][k][i]));
                          }
                  else{
                      P = P-lambda*lambda*I*(GK[klevel+plevel][k][i]*GK[plevel][i][k]+GR[plevel][k][i]*conjf(GR[klevel+plevel][k][i]));
                          }
                      }
              if (klevel+plevel>ktot) {
                  if (k>i){
                      P = P-lambda*lambda*I*(GK[klevel+plevel-ktot][k][i]*GK[plevel][i][k]+GR[klevel+plevel-ktot][k][i]*conjf(GR[plevel][k][i]));
                          }
                  else{
                      P = P-lambda*lambda*I*(GK[klevel+plevel-ktot][k][i]*GK[plevel][i][k]+GR[plevel][k][i]*conjf(GR[klevel+plevel-ktot][k][i]));
                          }
                        }
                        }

              SigPhK[klevel][k][i]=P;
              P=0;
                      }
                    }


    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Electronic Self Energies
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


              for (j=1; j<i; j++){
                  for (klevel=0; klevel< ktot+1; klevel++) {
                  for (plevel=0; plevel< ktot+1; plevel++) {

                      if (klevel+plevel<=ktot)
                         {
                         P = P-lambda*lambda*I*(GR[klevel+plevel][i][j]*DK[plevel][j][i]+GK[klevel+plevel][i][j]*conjf(DR[plevel][i][j]));
                         }
                     if (klevel+plevel>ktot)
                         {
                         P = P-lambda*lambda*I*(GR[klevel+plevel-ktot][i][j]*DK[plevel][j][i]+GK[klevel+plevel-ktot][i][j]*conjf(DR[plevel][i][j]));
                         }
                       }
                     SigElR[klevel][i][j]=-P;
                     P=0;
                         }
                      }


               for (j=1; j<i; j++){
                     for (klevel=0; klevel< ktot+1; klevel++) {
                     for (plevel=0; plevel< ktot+1; plevel++) {

                       if (klevel+plevel<=ktot){
                          if (i>j){
                               P = P-lambda*lambda*I*(GK[klevel+plevel][i][j]*DK[plevel][j][i]+GR[klevel+plevel][i][j]*conjf(DR[plevel][i][j]));
                                  }
                         else{
                               P = P-lambda*lambda*I*(GK[klevel+plevel][i][j]*DK[plevel][j][i]+GR[plevel][i][j]*conjf(DR[klevel+plevel][i][j]));
                                  }
                              }
                       if (klevel+plevel>ktot){
                           if (i>j){
                               P = P-lambda*lambda*I*(GK[klevel+plevel-ktot][i][j]*DK[plevel][j][i]+GR[klevel+plevel-ktot][i][j]*conjf(DR[plevel][i][j]));
                                  }
                           else{
                               P = P-lambda*lambda*I*(GK[klevel+plevel-ktot][i][j]*DK[plevel][j][i]+GR[plevel][i][j]*conjf(DR[klevel+plevel-ktot][i][j]));
                                  }
                                }
                                }
                        SigElK[klevel][i][j]=P;
                        P=0;
                           }
                        }


            for (k=1; k<=i; k++){
                    for (klevel=0; klevel< ktot+1; klevel++) {
                    for (plevel=0; plevel< ktot+1; plevel++) {

                        if (klevel+plevel<=ktot{
                            if (k>i){
                              P = P-lambda*lambda*I*(GK[klevel+plevel][k][i]*GK[plevel][i][k]+GR[klevel+plevel][k][i]*conjf(GR[plevel][k][i]));
                                    }
                            else{
                              P = P-lambda*lambda*I*(GK[klevel+plevel][k][i]*GK[plevel][i][k]+GR[plevel][k][i]*conjf(GR[klevel+plevel][k][i]));
                                    }
                              }
                        if (klevel+plevel>ktot) {
                            if (k>i){
                              P = P-lambda*lambda*I*(GK[klevel+plevel-ktot][k][i]*GK[plevel][i][k]+GR[klevel+plevel-ktot][k][i]*conjf(GR[plevel][k][i]));
                                    }
                            else{
                              P = P-lambda*lambda*I*(GK[klevel+plevel-ktot][k][i]*GK[plevel][i][k]+GR[plevel][k][i]*conjf(GR[klevel+plevel-ktot][k][i]));
                                    }
                                }
                                }

                        SigElK[klevel][k][i]=P;
                        P=0;
                                }
                            }












        }


             for (i = 0; i < n; i++){
               printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\n", i*h , crealf(GR[0][i][2]), cimagf(GR[0][i][2]), crealf(DK[0][i][i]), cimagf(DK[1][i][i]), crealf(DR[0][i][2]), cimagf(DR[0][i][2]) )  ;
             }


}
