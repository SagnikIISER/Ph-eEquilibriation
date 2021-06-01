
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
double complex I4,I5;


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
      b= 50.000000;				         /*Read value of b*/
      h= 0.100000;					       /*Read value of h*/

      lambda=1.00000;				       /*Read value of lambda*/
      n=(int)((b-a)/h)+1;			     /*array Dimension*/

  /*Phononic Lattice Parameters*/

      ktot= 1;					           /*Phononic Lattice Dimension*/
      A=1;						             /*Read value of Phononic Lattice Constant*/
  //  klevel=1; 			    		     /*Dummy Array for Phononic momentum*/



  /*Reading the Initial Temp*/
      Tphonon=  0.500000;
      Telectron=0.500000;
      nu = 0.0;

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

  GR[klevel][i][i]= I;
  GK[klevel][i][i]= -I*tanh((epsilon[klevel]-nu)/(2.0*Telectron));

  GR[klevel][i+1][i]= GzeroR(epsilon[klevel],(i*h)+h,(i*h));

  GK[klevel][i+1][i]= GzeroK(epsilon[klevel],Telectron,nu,(i*h)+h,(i*h));
  GK[klevel][i][i+1]= -conjf(GK[klevel][i+1][i]);

  /*Phononic Sector*/

  DR[klevel][i][i]=0;
  BarDR[klevel][i][i]= -1/2.0;

  DK[klevel][i][i]=-(I/(2.0*omega[klevel]))*(1.0/(tanh((omega[klevel])/(2*Tphonon))));
  BarDK[klevel][i][i]= 0;

  DR[klevel][i+1][i]= -2.0*DzeroR(omega[klevel],(i*h)+h,(i*h))*BarDR[klevel][i][i];

  DK[klevel][i+1][i]= -2.0*BarDzeroR(omega[klevel],(h*i)+h,(h*i))*DK[klevel][i][i];
  DK[klevel][i][i+1]=-conjf(DK[klevel][i+1][i]);

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
                   P = P+lambda*lambda*(GR[klevel+plevel][i][j]*GK[plevel][j][i]+GK[klevel+plevel][i][j]*conjf(GR[plevel][i][j]));
                   }
              if (klevel+plevel>ktot)
                 {
                    P = P+lambda*lambda*(GR[klevel+plevel-ktot][i][j]*GK[plevel][j][i]+GK[klevel+plevel-ktot][i][j]*conjf(GR[plevel][i][j]));
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
                       P = P+lambda*lambda*(GK[klevel+plevel][i][j]*GK[plevel][j][i]+GR[klevel+plevel][i][j]*conjf(GR[plevel][i][j]));
                           }
                  else{
                       P = P+lambda*lambda*(GK[klevel+plevel][i][j]*GK[plevel][j][i]+GR[plevel][i][j]*conjf(GR[klevel+plevel][i][j]));
                           }
                        }
               if (klevel+plevel>ktot){
                  if (i>j){
                       P = P+lambda*lambda*(GK[klevel+plevel-ktot][i][j]*GK[plevel][j][i]+GR[klevel+plevel-ktot][i][j]*conjf(GR[plevel][i][j]));
                          }
                  else{
                       P = P+lambda*lambda*(GK[klevel+plevel-ktot][i][j]*GK[plevel][j][i]+GR[plevel][i][j]*conjf(GR[klevel+plevel-ktot][i][j]));
                          }
                        }
                        }
               SigPhK[klevel][i][j]=P;
               P=0;
                 }
              }


         for (k=1; k<i; k++){
            for (klevel=0; klevel< ktot+1; klevel++) {
            for (plevel=0; plevel< ktot+1; plevel++) {

              if (klevel+plevel<=ktot){
                  if (k>i){
                      P = P+lambda*lambda*(GK[klevel+plevel][k][i]*GK[plevel][i][k]+GR[klevel+plevel][k][i]*conjf(GR[plevel][k][i]));
                          }
                  else{
                      P = P+lambda*lambda*(GK[klevel+plevel][k][i]*GK[plevel][i][k]+GR[plevel][k][i]*conjf(GR[klevel+plevel][k][i]));
                          }
                      }
              if (klevel+plevel>ktot) {
                  if (k>i){
                      P = P+lambda*lambda*(GK[klevel+plevel-ktot][k][i]*GK[plevel][i][k]+GR[klevel+plevel-ktot][k][i]*conjf(GR[plevel][k][i]));
                          }
                  else{
                      P = P+lambda*lambda*(GK[klevel+plevel-ktot][k][i]*GK[plevel][i][k]+GR[plevel][k][i]*conjf(GR[klevel+plevel-ktot][k][i]));
                          }
                        }
                        }

              SigPhK[klevel][k][i]=P;
              P=0;
                      }
                    }

  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  The Phonons
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


        /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        The Retarded Part
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

        for (klevel=0; klevel< ktot+1; klevel++) {
      	for (j=1; j<i; j++){
            for (l = j+1; l < i; l++)
                {
                  P=P+h*SigPhR[klevel][i][l]*DR[klevel][l][j];
                }

                IPh1A[j]=P;
                P=0.0;

            BarDR[klevel][i][j] = -2.0*BarDzeroR(omega[klevel],(i*h),(i*h)-h)*BarDR[klevel][i-1][j]+2.0*omega[klevel]*omega[klevel]*DzeroR(omega[klevel],(i*h),(i*h)-h)*DR[klevel][i-1][j]+(h/2.0)*BarDzeroR(omega[klevel],(i*h),(i*h)-h)*IPh1B[j];
            DR[klevel][i+1][j]= -2.0*DzeroR(omega[klevel],(i*h)+h,(i*h))*BarDR[klevel][i][j]-2.0*BarDzeroR(omega[klevel],(i*h)+h,(i*h))*DR[klevel][i][j]+(h/2.0)*DzeroR(omega[klevel],(i*h)+h,(i*h))*IPh1A[j];

            IPh1B[j]=IPh1A[j];

                }
              }

        /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
        /*The Keldysh Part*/
        /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

        for (klevel=0; klevel< ktot+1; klevel++) {
        for (j=0; j<i; j++){
            for (l = 1; l <= i; l++)
                {
                  P=P+h*SigPhR[klevel][i][l]*DK[klevel][l][j];
                }

             IPh3A[j]=P+(h/2.0)*SigPhR[klevel][i][j]*DK[klevel][j][j];
             P=0.0;

             for (l = 1; l < j; l++)
                 {
                   P = P + h*SigPhK[klevel][i][l]*conjf(DR[klevel][j][l]);
                 }

             IPh2A[j]= P + (h/2.0)*SigPhK[klevel][i][i]*conjf(DR[klevel][j][i]);
             P=0;

             BarDK[klevel][i][j] = -2.0*BarDzeroR(omega[klevel],(i*h),(i*h)-h)*BarDK[klevel][i-1][j]+2.0*omega[klevel]*omega[klevel]*DzeroR(omega[klevel],(i*h),(i*h)-h)*DK[klevel][i-1][j]+(h/2.0)*BarDzeroR(omega[klevel],(i*h),(i*h)-h)*(IPh3B[j]+IPh2B[j]);
             DK[klevel][i+1][j]= -2.0*DzeroR(omega[klevel],(i*h)+h,(i*h))*BarDK[klevel][i][j]-2.0*BarDzeroR(omega[klevel],(i*h)+h,(i*h))*DK[klevel][i][j]+(h/2.0)*DzeroR(omega[klevel],(i*h)+h,(i*h))*(IPh3A[j]+IPh2A[j]);

             IPh3B[j]=IPh3A[j];
             IPh2B[j]=IPh2A[j];


             DK[klevel][j][i+1]=-conjf(DK[klevel][i+1][j]);
             BarDK[klevel][j][i]=conjf(BarDK[klevel][i][j]);

             }
            }

        for (klevel=0; klevel< ktot+1; klevel++) {
        for (k=1; k<=i; k++){
             for (l = 1; l <= k; l++)
                {
                  P=P+h*SigPhR[klevel][k][l]*DK[klevel][l][i];
                  }

             I3[k+1]=P+(h/2.0)*SigPhR[klevel][i][i]*DK[klevel][i][i];
             P=0.0;

             for (l = 1; l < i; l++)
                {
                P = P + h*SigPhK[klevel][k][l]*conjf(DR[klevel][i][l]);
                }


            I2[k+1]=P+(h/2.0)*SigPhK[klevel][k][k]*conjf(DR[klevel][i][k]);
            P=0.0;

            BarDK[klevel][k][i] = -2.0*BarDzeroR(omega[klevel],(k*h),(k*h)-h)*BarDK[klevel][k-1][i]+2.0*omega[klevel]*omega[klevel]*DzeroR(omega[klevel],(k*h),(k*h)-h)*DK[klevel][k-1][i]+(h/2.0)*BarDzeroR(omega[klevel],(k*h),(k*h)-h)*(I3[k-1]+I2[k-1]);
            DK[klevel][k+1][i]= -2.0*DzeroR(omega[klevel],(k*h)+h,(k*h))*BarDK[klevel][k][i]-2.0*BarDzeroR(omega[klevel],(k*h)+h,(k*h))*DK[klevel][k][i]+(h/2.0)*DzeroR(omega[klevel],(k*h)+h,(k*h))*(I3[k]+I2[k]);
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
                         P = P+lambda*lambda*(GR[klevel+plevel][i][j]*DK[plevel][j][i]+GK[klevel+plevel][i][j]*conjf(DR[plevel][i][j]));
                         }
                     if (klevel+plevel>ktot)
                         {
                         P = P+lambda*lambda*(GR[klevel+plevel-ktot][i][j]*DK[plevel][j][i]+GK[klevel+plevel-ktot][i][j]*conjf(DR[plevel][i][j]));
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
                               P = P+lambda*lambda*(GK[klevel+plevel][i][j]*DK[plevel][j][i]+GR[klevel+plevel][i][j]*conjf(DR[plevel][i][j]));
                                  }
                         else{
                               P = P+lambda*lambda*(GK[klevel+plevel][i][j]*DK[plevel][j][i]+GR[plevel][i][j]*conjf(DR[klevel+plevel][i][j]));
                                  }
                              }
                       if (klevel+plevel>ktot){
                           if (i>j){
                               P = P+lambda*lambda*(GK[klevel+plevel-ktot][i][j]*DK[plevel][j][i]+GR[klevel+plevel-ktot][i][j]*conjf(DR[plevel][i][j]));
                                  }
                           else{
                               P = P+lambda*lambda*(GK[klevel+plevel-ktot][i][j]*DK[plevel][j][i]+GR[plevel][i][j]*conjf(DR[klevel+plevel-ktot][i][j]));
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

                        if (klevel+plevel<=ktot){
                            if (k>i){
                              P = P+lambda*lambda*(GK[klevel+plevel][k][i]*GK[plevel][i][k]+GR[klevel+plevel][k][i]*conjf(GR[plevel][k][i]));
                                    }
                            else{
                              P = P+lambda*lambda*(GK[klevel+plevel][k][i]*GK[plevel][i][k]+GR[plevel][k][i]*conjf(GR[klevel+plevel][k][i]));
                                    }
                              }
                        if (klevel+plevel>ktot) {
                            if (k>i){
                              P = P+lambda*lambda*(GK[klevel+plevel-ktot][k][i]*GK[plevel][i][k]+GR[klevel+plevel-ktot][k][i]*conjf(GR[plevel][k][i]));
                                    }
                            else{
                              P = P+lambda*lambda*(GK[klevel+plevel-ktot][k][i]*GK[plevel][i][k]+GR[plevel][k][i]*conjf(GR[klevel+plevel-ktot][k][i]));
                                    }
                                }
                                }

                        SigElK[klevel][k][i]=P;
                        P=0;
                                }
                            }


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  The Electrons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


        /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         The Retarded Part
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

          for (klevel=0; klevel< ktot+1; klevel++) {
          for (j=1; j<i; j++){
             for (l = j+1; l < i; l++)
                {
                 P=P+h*SigElR[klevel][i][l]*GR[klevel][l][j];
                }

             I4=P;
             P=0.0;

           GR[klevel][i+1][j]= I*GzeroR(epsilon[klevel],(i*h)+h,(i*h))*GR[klevel][i][j]+(h/2.0)*GzeroR(epsilon[klevel],(i*h)+h,(i*h))*I4;

               }
             }


        /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
        /*The Keldysh Part*/
        /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

          for (klevel=0; klevel< ktot+1; klevel++) {
          for (j=0; j<i; j++){
              for (l = 1; l <= i; l++)
                {
                P=P+h*SigElR[klevel][i][l]*GK[klevel][l][j];
                }

          I4=P+(h/2.0)*SigElR[klevel][i][j]*GK[klevel][j][j];;
          P=0.0;

                for (l = 1; l < j; l++)
                {
                P = P + h*SigElK[klevel][i][l]*conjf(GR[klevel][j][l]);
                }

          I5= P + (h/2.0)*SigElK[klevel][i][i]*conjf(GR[klevel][j][i]);
          P=0;

          GK[klevel][i+1][j]= I*GzeroR(epsilon[klevel],(i*h)+h,(i*h))*GK[klevel][i][j]+(h/2.0)*GzeroR(epsilon[klevel],(i*h)+h,(i*h))*(I4+I5);
          GK[klevel][j][i+1]=-conjf(GK[klevel][i+1][j]);
                }
              }



          for (klevel=0; klevel< ktot+1; klevel++) {
          for (k=0; k<=i; k++){
               for (l = 1; l <= k; l++)
                {
                P=P+h*SigElR[klevel][k][l]*GK[klevel][l][i+1];
                }

          I4=P+(h/2.0)*SigElR[klevel][k][i]*GK[klevel][i][i];
          P=0.0;

               for (l = 1; l < i; l++)
                {
                P = P + h*SigElK[klevel][k][l]*conjf(GR[klevel][i+1][l]);
                }


          I5=P+(h/2.0)*SigElK[klevel][k][k]*conjf(GR[klevel][i+1][k]);
          P=0.0;

          GK[klevel][k+1][i+1]= I*GzeroR(epsilon[klevel],(k*h)+h,(k*h))*GK[klevel][k][i+1]+(h/2.0)*GzeroR(epsilon[klevel],(k*h)+h,(k*h))*(I4+I5);


          }
        }


 }

 /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 K loop ends here.
 Put K independent printf statements after this.
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

             for (i = 0; i < n; i++){
               printf("%f\t%f\t%f\n", i*h , crealf(GR[0][i][3]),  crealf(DR[0][i][3]) )  ;
             }


}
