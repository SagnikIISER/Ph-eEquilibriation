
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CODE for SOLVING DYSON EQUATION

First note down all the relevant formula plus diagram
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
double complex IR[2510][2510];

double complex SigK[2510][2510];
double complex DA[2510][2510];

double complex DK[2510][2510];
double complex BarDK[2510][2510];
double complex DKthermal;

double complex IK1[2510][2510];
double complex IK2[2510][2510];


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

double sigma, Tbath, Tsyst;


double complex partial_Sum, total_Sum;
double akka;
double complex P;


/*Working around the input variable calls*/

  a= 0.000000;					/*Read value of a*/
  b= 250.000000;				 /*Read value of b*/
  h= 0.100000;					/*Read value of h*/

  lambda=0.50000;				/*Read value of lambda*/
  n=(int)((b-a)/h)+1;			    	/*array Dimension*/
  ktot= 19;					/*array Dimension*/

  A=1;						/*Read value of Lattice Constant*/
  klevel=1; 			    		/*Dummy Array for momentum*/

  sigma=5.000000;
  Tsyst=1.000000;
  Tbath=0.800000;

int i;
int j;
int k;
int l;


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The Module for Dyson Iteration using Euler Method:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


//for (klevel=0; klevel<ktot; klevel=klevel+1){

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Momentum values at Dummy Indicies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

 K=((-3.1418/A)+klevel*(6.2836/A)/(ktot-1.0));

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
System Dispersion Relation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

 omega=2.0*sqrt(sin((K*A/2.0))*sin((K*A/2.0)));


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 To Plot the Dispersion Relation
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

//printf("%f\t%f\n", K, omega );


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


IR[i][i]=0;
DR[i][i]=0;

BarDR[i][i]= -1/2.0;

DK[i][i]=-(I/(2.0*omega))*(1.0/(tanh((omega)/(2*Tsyst))));
BarDK[i][i]= 0;

DR[i+1][i]= -2.0*DzeroR(omega,t+h,t)*BarDR[i][i];
IR[i+1][i]= +1.0/(2.0*sqrt(3.14159265))*lambda*lambda*sigma*DR[i+1][i];


DK[i+1][i]= -2.0*BarDzeroR(omega,(h*i)+h,(h*i))*DK[i][i];
DK[i][i+1]=-conjf(DK[i+1][i]);


IK1[i][i]=0;
IK2[i][i]=0;
IK1[i+1][i]=(h/2.0)*lambda*lambda*SigmaR((i*h)+h,(i*h))*DK[i][i]+(1.0/(2.0*sqrt(3.14159265)))*lambda*lambda*sigma*DK[i+1][i];
IK2[i+1][i]=(h/2.0)*lambda*lambda*SigK[i+1][i+1]*DR[i+1][i];
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The Retarded Part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/



    	for (j=0; j<n; j++){
    	for (i=j+1; i<n; i++){
    	 	BarDR[i][j] = -2.0*BarDzeroR(omega,(i*h),(i*h)-h)*BarDR[i-1][j]+2.0*omega*omega*DzeroR(omega,(i*h),(i*h)-h)*DR[i-1][j]+(h/2.0)*BarDzeroR(omega,(i*h),(i*h)-h)*IR[i-1][j];
    	 	DR[i+1][j]= -2.0*DzeroR(omega,(i*h)+h,(i*h))*BarDR[i][j]-2.0*BarDzeroR(omega,(i*h)+h,(i*h))*DR[i][j]+(h/2.0)*DzeroR(omega,(i*h)+h,(i*h))*IR[i][j];

    				for (l = j+1; l < i; l++)
    				{
    					P=P+h*lambda*lambda*SigmaR((i*h)+h,(h*l))*DR[l][j];
    				}

    		IR[i+1][j]=P+(1.0/(2.0*(sqrt(3.14159265))))*lambda*lambda*sigma*DR[i+1][j];
    		P=0.0;
    	}
    	}





/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The Advanced Part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


for (j=1; j<n ; j++){
for (i=0; i<n ; i++){
		DA[j][i]=conjf(DR[i][j]);
}
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*The Keldysh Part*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


	for (j=1; j<n ; j++){
	for (i=1; i<n ; i++){
	  BarDK[i][j] = -2.0*BarDzeroR(omega,(i*h),(i*h)-h)*BarDK[i-1][j]+2.0*omega*omega*DzeroR(omega,(i*h),(i*h)-h)*DK[i-1][j]+(h/2.0)*BarDzeroR(omega,(i*h),(i*h)-h)*(IK1[i-1][j]+IK2[i-1][j]);
		DK[i+1][j]= -2.0*DzeroR(omega,(i*h)+h,(i*h))*BarDK[i][j]-2.0*BarDzeroR(omega,(i*h)+h,(i*h))*DK[i][j]+(h/2.0)*DzeroR(omega,(i*h)+h,(i*h))*(IK1[i][j]+IK2[i][j]);

				for (l = 1; l <= i; l++)
				{
					P=P+h*lambda*lambda*SigmaR((i*h)+h,(h*l))*DK[l][j];
				}

		IK1[i+1][j]=P+(h/2.0)*lambda*lambda*SigmaR((i*h)+h,(j*h))*DK[j][j]+(1.0/(2.0*sqrt(3.14159265)))*lambda*lambda*sigma*DK[i+1][j];
		P=0.0;


        //#pragma omp parallel private(partial_Sum) shared(total_Sum)
        //{

              partial_Sum = 0;
              total_Sum = 0;

        //#pragma omp for

						for (l = 1; l < j; l++)
						{
							partial_Sum = partial_Sum+h*lambda*lambda*SigK[i+1][l]*DA[l][j];
						}


            //Create thread safe region.
            //#pragma omp critical
            //{
                    //add each threads partial sum to the total sum
                    total_Sum = partial_Sum;
            //}
            IK2[i+1][j]=total_Sum+(h/2.0)*lambda*lambda*SigK[i+1][i+1]*DA[i+1][j];

            //}
            partial_Sum = 0;
            total_Sum = 0;

		DK[j][i+1]=-conjf(DK[i+1][j]);
		BarDK[j][i]=conjf(BarDK[i][j]);

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
printf("%f\t%f\n", i*h , -cimagf(DK[i][i]))  ;
}
}
