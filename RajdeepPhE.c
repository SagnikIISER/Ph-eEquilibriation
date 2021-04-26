
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

double complex GR[6][1510][1510];
double complex GK[6][1510][1510];
double complex SigElR[6][1510][1510];
double complex SigElK[6][1510][1510];


/*%%%%%%%%%%%%%%%%%%%%%%%%%
Phonon Sector
**************************/

double complex DR[6][1510][1510];
double complex BarDR[6][1510][1510];

double complex DK[6][1510][1510];
double complex BarDK[6][1510][1510];

double complex SigPhR[6][1510][1510];
double complex SigPhK[6][1510][1510];

/*%%%%%%%%%%%%%%%%%%%%%%%%
Electron Sector
**************************/

double complex I1[1510];
double complex I2[1510];






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


  double complex partial_Sum, total_Sum;

    int i;
    int j;
    int l;

/*Working around the input variable calls*/

    a= 0.000000;					       /*Read value of a*/
    b= 15.000000;				       /*Read value of b*/
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
    nu = -1.0000;

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Energy Levels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*Momentum values at Dummy Indicies*/
   K=((-3.1418/A)+klevel*(6.2836/A)/(ktot-1.0));

/*Phononic Dispersion Relation*/
   omega=3.0*sqrt(sin((K*A/2.0))*sin((K*A/2.0)));

/*Electronic Dispersion Relation*/
   epsilon= K*K/2.00000;



//printf("%f\t%f\n", omega,epsilon );

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Initial Condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


for (i=0; i<n; i++){


  /*Electronic Sector*/

        GR[0][i][i]= -I;
        GK[0][i][i]= I*tanh((epsilon-nu)/(2.0*Telectron));

        GR[0][i+1][i]= GzeroR(epsilon,t+h,t);

        GK[0][i+1][i]= GzeroK(epsilon,Telectron,nu,t+h,t);
        GK[0][i][i+1]=-conjf(DK[0][i+1][i]);


/*Phononic Sector*/

        DR[0][i][i]=0;
        BarDR[0][i][i]= -1/2.0;

        DK[0][i][i]=-(I/(2.0*omega))*(1.0/(tanh((omega)/(2*Tphonon))));
        BarDK[0][i][i]= 0;

        DR[0][i+1][i]= -2.0*DzeroR(omega,t+h,t)*BarDR[0][i][i];


        DK[0][i+1][i]= -2.0*BarDzeroR(omega,(h*i)+h,(h*i))*DK[0][i][i];
        DK[0][i][i+1]= -conjf(DK[0][i+1][i]);



   }



/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dyson Iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


for (j=0; j<n; j++){


  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Electron Sector
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

     /*The Retarded Part*/

       	for (i=j+1; i<n; i++){

       	 	GR[0][i+1][j]= -I*GzeroR(epsilon,(i*h)+h,(i*h))*GR[0][i][j]+(h/2.0)*GzeroR(omega,(i*h)+h,(i*h))*I1[i];

       				for (l = j+1; l < i; l++)
       				{
       					P=P+h*SigElR[0][i+1][l]*GR[0][l][j];
       				}


       		I1[i+1]=P;
     		  P=0.0;

          }
     /*The Keldysh Part*/

           for (i=1; i<n ; i++){

       		 GK[0][i+1][j]= I*GzeroR(omega,(i*h)+h,(i*h))*GK[0][i][j]+(h/2.0)*DzeroR(omega,(i*h)+h,(i*h))*(I1[i]+I2[i]);

       				for (l = 1; l <= i; l++)
       				{
       					P=P+h*SigElR[0][i+1][l]*DK[0][l][j];
       				}

       		I1[i+1]=P+(h/2.0)*SigElR[0][i+1][j]*DK[0][j][j];
       		P=0.0;


               //#pragma omp parallel private(partial_Sum) shared(total_Sum)
               //{

                     partial_Sum = 0;
                     total_Sum = 0;

               //#pragma omp for

       						for (l = 1; l < j; l++)
       						{
       							partial_Sum = partial_Sum+h*SigElK[0][i+1][l]*conjf(GR[0][j][l]);
       						}


                   //Create thread safe region.
                   //#pragma omp critical
                   //{
                           //add each threads partial sum to the total sum
                           total_Sum = partial_Sum;
                   //}
                   I2[i+1]=total_Sum+(h/2.0)*SigElK[0][i+1][i+1]*conjf(GR[0][j][i+1]);

                   //}
                   partial_Sum = 0;
                   total_Sum = 0;

       		GK[0][j][i+1]=-conjf(GK[0][i+1][j]);

       	}

        for (i=0; i<n; i++){
          I1[i]=0;
          I2[i]=0;
        }

    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Phononic Self Energies
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

        for (i=j+1; i<n; i++){
            SigPhR[0][i][j]= -crealf(I*GK[0][i][j]*conjf(GR[0][i][j]));
        }

        for (i=1; i<n; i++){
          if (i>j) {
            SigPhK[0][i][j]= GR[0][i][j]*conjf(GR[0][i][j])+GK[0][i][j]*conjf(GK[0][i][j]) ;
                    }
                    else
                  {SigPhK[0][i][j]= GR[0][j][i]*conjf(GR[0][j][i])+GK[0][i][j]*conjf(GK[0][i][j]) ;
                  }

        }


     /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Phonon Sector
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/



     /*The Retarded Part*/

         	for (i=j+1; i<n; i++){
         	 	BarDR[0][i][j] = -2.0*BarDzeroR(omega,(i*h),(i*h)-h)*BarDR[0][i-1][j]+2.0*omega*omega*DzeroR(omega,(i*h),(i*h)-h)*DR[0][i-1][j]+(h/2.0)*BarDzeroR(omega,(i*h),(i*h)-h)*I1[i-1];
         	 	DR[0][i+1][j]= -2.0*DzeroR(omega,(i*h)+h,(i*h))*BarDR[0][i][j]-2.0*BarDzeroR(omega,(i*h)+h,(i*h))*DR[0][i][j]+(h/2.0)*DzeroR(omega,(i*h)+h,(i*h))*I1[i];

         				for (l = j+1; l < i; l++)
         				{
         					P=P+h*SigPhR[0][i+1][l]*DR[0][l][j];
         				}

         		I1[i+1]=P;
         		P=0.0;
         	}



     /*The Keldysh Part*/


     	for (i=1; i<n ; i++){
     	  BarDK[0][i][j] = -2.0*BarDzeroR(omega,(i*h),(i*h)-h)*BarDK[0][i-1][j]+2.0*omega*omega*DzeroR(omega,(i*h),(i*h)-h)*DK[0][i-1][j]+(h/2.0)*BarDzeroR(omega,(i*h),(i*h)-h)*(I1[i-1]+I2[i-1]);
     		DK[0][i+1][j]= -2.0*DzeroR(omega,(i*h)+h,(i*h))*BarDK[0][i][j]-2.0*BarDzeroR(omega,(i*h)+h,(i*h))*DK[0][i][j]+(h/2.0)*DzeroR(omega,(i*h)+h,(i*h))*(I1[i]+I2[i]);

     				for (l = 1; l <= i; l++)
     				{
     					P=P+h*SigPhR[0][i+1][l]*DK[0][l][j];
     				}

     		I1[i+1]=P+(h/2.0)*SigPhR[0][i+1][j]*DK[0][j][j];


             //#pragma omp parallel private(partial_Sum) shared(total_Sum)
             //{

                   partial_Sum = 0;
                   total_Sum = 0;

             //#pragma omp for

     						for (l = 1; l < j; l++)
     						{
     							partial_Sum = partial_Sum+h*SigPhK[0][i+1][l]*conjf(DR[0][j][l]);
     						}


                 //Create thread safe region.
                 //#pragma omp critical
                 //{
                         //add each threads partial sum to the total sum
                         total_Sum = partial_Sum;
                 //}
                 I2[i+1]=total_Sum+(h/2.0)*SigPhK[0][i+1][i+1]*conjf(DR[0][j][i]);

                 //}
                 partial_Sum = 0;
                 total_Sum = 0;

     		DK[0][j][i+1]=-conjf(DK[0][i+1][j]);
     		BarDK[0][j][i]=conjf(BarDK[0][i][j]);

     	}





    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Electronic Self Energies
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

        for (i=j+1; i<n; i++){
            SigElR[0][i][j]= I*GR[0][i][j]*DK[0][i][j]+I*GK[0][i][j]*DR[0][i][j];
        }

        for (i=1; i<n; i++){
          if (i>j) {
            SigElK[0][i][j]= -GR[0][i][j]*DR[0][i][j]-GK[0][i][j]*DK[0][i][j];
          }
          else{  SigElK[0][i][j]= -conjf(GR[0][j][i])*DR[0][j][i]-conjf(GK[0][j][i])*DK[0][j][i];}

        }


   }


   for (i = 0; i < n; i++){
   printf("%f\t%f\t%f\n", i*h , -cimagf(DK[0][i][i]), cimagf(GK[0][i][i]) )  ;
   }



}
