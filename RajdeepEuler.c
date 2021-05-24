
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CODE for SOLVING DYSON EQUATION

FI1st note down all the relevant formula plus diagram
Only true unknown in the whole buisness is the D and BarD
There are two modules written for the Retarded and Keldysh Component:
			a. Euler
			b. Self Consistent Field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include <stdio.h> 		/*Import Standard Module*/
#include <stdlib.h> 		/*Import Standard Library*/
#include <omp.h>
#include <math.h>		/*Import Math Module*/
#include <complex.h>          	/*Import Modules for Complex Numbers*/
#include "Dawson.h"		/*Import Module for Dawson Function*/
#include "green's.h"		/*Import Green's Functions and Baths*/


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Defining Iteration variables and outputs
DR, BarDR and the integration I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

double complex DR[2510][2510];
double complex BarDR[2510][2510];

double complex SigK[2510][2510];
double complex DK[2510][2510];
double complex BarDK[2510][2510];
double complex DKthermal;

double complex I1[2510][2510];
double complex I2[2510][2510];
double complex I3[2510][2510];

double complex I1A, I1B, I2A, I2B, I3A, I3B;
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
double omega;					/*Initial Energy Level*/
double lambda;					/*Perturbation Strength*/
double nomega;          /*Density of states*/

double sigma, Tbath, Tsyst;


double complex partial_Sum, total_Sum;
double akka;
double complex P;


/*Working around the input variable calls*/

  a= 0.000000;					/*Read value of a*/
  b= 25.000000;				 /*Read value of b*/
  h= 0.050000;					/*Read value of h*/

  lambda=1.00000;				/*Read value of lambda*/
  n=(int)((b-a)/h)+1;			    	/*array Dimension*/
  ktot= 19;					/*array Dimension*/

  A= 1.0;						/*Read value of Lattice Constant*/
  klevel= 0; 			    		/*Dummy Array for momentum*/

  sigma=10.000000;
  Tsyst=1.000000;
  Tbath=1.200000;

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

 K=((-3.1418/A)+klevel*(6.2836/A)/(ktot-1.0));

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
System Dispersion Relation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

 omega=2.0*sqrt(sin((K*A/2.0))*sin((K*A/2.0)));
 nomega= 1/(exp(omega/Tsyst)-1);

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 To Plot the Dispersion Relation
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*
printf("%d\t%f\t%f\t%f\n", klevel, K, omega, nomega );

}
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Computing the DK thermal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*



P=0.0;

for(akka=-100; akka <= 100; akka=akka+h){
P = P + h*(lambda*lambda*akka*exp(-(akka*akka/(sigma*sigma))))/((2.0*akka*akka-2.0*omega*omega+lambda*lambda*(sigma/sqrt(3.14159265))-lambda*lambda*(2.0/sqrt(3.14159265))*akka*Dawson(akka/sigma))*(2.0*akka*akka-2.0*omega*omega+lambda*lambda*(sigma/sqrt(3.14159265))-lambda*lambda*(2.0/sqrt(3.14159265))*akka*Dawson(akka/sigma))+(lambda*lambda*lambda*lambda*akka*akka*exp(-2.0*(akka*akka/(sigma*sigma)))))*(1.0/tanh((akka/(2.0*Tbath))));
}

DKthermal = P/(3.14159265);


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

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Initial Condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/



for (i=0; i<n; i++){

DR[i][i]=0;

BarDR[i][i]= -1/2.0;

DK[i][i]=-(I/(2.0*omega))*(1.0/(tanh((omega)/(2*Tsyst))));
BarDK[i][i]= 0;

DR[i+1][i]= -2.0*DzeroR(omega,t+h,t)*BarDR[i][i];

DK[i+1][i]= -2.0*BarDzeroR(omega,(h*i)+h,(h*i))*DK[i][i];
DK[i][i+1]=-conjf(DK[i+1][i]);


}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The Retarded Part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


for (i=1; i<n; i++){
    	for (j=1; j<i; j++){
    	 	BarDR[i][j] = -2.0*BarDzeroR(omega,(i*h),(i*h)-h)*BarDR[i-1][j]+2.0*omega*omega*DzeroR(omega,(i*h),(i*h)-h)*DR[i-1][j]+(h/2.0)*BarDzeroR(omega,(i*h),(i*h)-h)*I1[i-1][j];
    	 	DR[i+1][j]= -2.0*DzeroR(omega,(i*h)+h,(i*h))*BarDR[i][j]-2.0*BarDzeroR(omega,(i*h)+h,(i*h))*DR[i][j]+(h/2.0)*DzeroR(omega,(i*h)+h,(i*h))*I1[i][j];

    				for (l = j+1; l < i; l++)
    				{
    					P=P+h*lambda*lambda*SigmaR((i*h)+h,(h*l))*DR[l][j];
    				}

    		I1[i+1][j]=P+(1.0/(2.0*(sqrt(3.14159265))))*lambda*lambda*sigma*DR[i+1][j];
    		P=0.0;
    	}
    	}





/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*The Keldysh Part*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

    for (i=1; i<n; i++){
    for (j=0; j<i; j++){
    BarDK[i][j] = -2.0*BarDzeroR(omega,(i*h),(i*h)-h)*BarDK[i-1][j]+2.0*omega*omega*DzeroR(omega,(i*h),(i*h)-h)*DK[i-1][j]+(h/2.0)*BarDzeroR(omega,(i*h),(i*h)-h)*(I3[i-1][j]+I2[i-1][j]);
		DK[i+1][j]= -2.0*DzeroR(omega,(i*h)+h,(i*h))*BarDK[i][j]-2.0*BarDzeroR(omega,(i*h)+h,(i*h))*DK[i][j]+(h/2.0)*DzeroR(omega,(i*h)+h,(i*h))*(I3[i][j]+I2[i][j]);

				for (l = 1; l <= i; l++)
				{
					P=P+h*lambda*lambda*SigmaR((i*h)+h,(h*l))*DK[l][j];
				}

		I3[i+1][j]=P+(h/2.0)*lambda*lambda*SigmaR((i*h)+h,(j*h))*DK[j][j]+(1.0/(2.0*sqrt(3.14159265)))*lambda*lambda*sigma*DK[i+1][j];
		P=0.0;

              partial_Sum = 0;
              total_Sum = 0;

						for (l = 1; l < j; l++)
						{
						}

            I2[i+1][j]= P + (h/2.0)*lambda*lambda*SigK[i+1][i+1]*conjf(DR[j][i+1]);
            P = P + h*lambda*lambda*SigK[i+1][l]*conjf(DR[j][l]);
            P=0;

		DK[j][i+1]=-conjf(DK[i+1][j]);
		BarDK[j][i]=conjf(BarDK[i][j]);

	}

  for (k=1; k<=i; k++)

  {
  BarDK[k][i] = -2.0*BarDzeroR(omega,(k*h),(k*h)-h)*BarDK[k-1][i]+2.0*omega*omega*DzeroR(omega,(k*h),(k*h)-h)*DK[k-1][i]+(h/2.0)*BarDzeroR(omega,(k*h),(k*h)-h)*(I3[k-1][i]+I2[k-1][i]);
  DK[k+1][i]= -2.0*DzeroR(omega,(k*h)+h,(k*h))*BarDK[k][i]-2.0*BarDzeroR(omega,(k*h)+h,(k*h))*DK[k][i]+(h/2.0)*DzeroR(omega,(k*h)+h,(k*h))*(I3[k][i]+I2[k][i]);

      for (l = 1; l <= k; l++)
      {
        P=P+h*lambda*lambda*SigmaR((k*h)+h,(h*l))*DK[l][i];
      }

      I3[k+1][i]=P+(h/2.0)*lambda*lambda*SigmaR((k*h)+h,(i*h))*DK[i][i]+(1.0/(2.0*sqrt(3.14159265)))*lambda*lambda*sigma*DK[k+1][i];

      P=0.0;

          for (l = 1; l < i; l++)
          {
            P = P + h*lambda*lambda*SigK[k+1][l]*conjf(DR[i][l]);
          }


          I2[k+1][i]=P+(h/2.0)*lambda*lambda*SigK[k+1][k+1]*conjf(DR[i][k+1]);

          P=0.0;
  }
	}






/*The Occupation Number*/
/*
		for (i=0; i<n ; i++){
		N[i]=0.5*(-2*omega*cimagf(DK[i][i])-1);
		}



/*The Statistics*/
/*
double VarE[n], EE[n], E[n], Ntot[n];-2*omega*Dawson(omega/(sqrt(2)*sigma))/sqrt(3.14159265)
double EEPr, EPr, NPr;

EPr=0.0;
EEPr=0.0;
NPr=0.0;

, cimagf(DR[0][i][0])


		for ( t = a; t <=b ; t=t+h) {-2*omega*Dawson(omega/(sqrt(2)*sigma))/sqrt(3.14159265)
		for ( k = 0; k <ktot ; k=k+1) {
					NPr=NPr+N[k][i];
				}
				Ntot[i]=NPr;
				NPr=0.0;
			}


			for ( t = a; t <=b ; t=t+h) {
			for ( k = 0; k <ktot ; k=k+1) {
				K=((-3.1418/A)+k*(6.2836/A)/(ktot-1));
				omega=-2*2*cos(K);
				EPr=EPr+(omega)*N[k][i];
					}
					E[i]=EPr;
					EPr=0.0;
				}

			for ( t = a; t DzeroR(<=b ; t=t+h) {
			for ( k = 0; k <ktot ; k=k+1) {
						EEPr=EEPr+(omega)*(omega)*N[k][i];
						}
						EE[i]=EEPr;
						EEPr=0.0;
				}

			for ( t = a+h; t <=b ; t=t+h) {
				printf("%f\t%f\n", t, (E[i]/Ntot[i]) );
			}


//}

*/
/*	for (i=0; i<n ; i++){
			 printf("%f\t%f\n", t, DK[0][i][i] );
		 }

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K loop ends here.
Put K independent printf statements after this.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
//printf("\n%f\n%f\n%f", omega, crealf(DKthermal), -cimagf(DK[n][n])-crealf(DKthermal));


for (i = 0; i < n; i++){
printf("%f\t%f\t%f\n", i*h , crealf(DR[i][12]), -cimagf(DK[i][i]))  ;
}

}
