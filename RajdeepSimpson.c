
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
#include <math.h>		/*Import Math Module*/
#include <complex.h>          	/*Import Modules for Complex Numbers*/
#include "Dawson.h"		/*Import Module for Dawson Function*/
#include "green's.h"		/*Import Green's Functions and Baths*/


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Defining Iteration variables and outputs
DR, BarDR and the integration I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

float complex DR[2000][2000];
float complex BarDR[2000][2000];
float complex IR[2000][2000];

float complex DA[2000][2000];

float complex DK[2000][2000];
float complex BarDK[2000][2000];
float complex DKthermal;
float complex DKthermal2;

float complex IK1[2000][2000];
float complex IK2[2000][2000];

float N[2000];

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Main Module
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

float main() {


/*Definitions*/

float t, tprime;				/*time parameters*/
float a,b; 					/*to be used for lower and upper bound on time*/
float tone, ttwo;				/*extra variables*/
float h;					/*gap in time, "the epsilon"*/

float A;					/*Lattice Constant*/
int   n, k, ktot;				/*Array Dimension*/
float K;					/*Momentum Formula*/
float omega;					/*Initial Energy Level*/
float lambda;					/*Perturbation Strength*/

float sigma, Tbath, Tsyst;



/*Working around the input variable calls*/

  a= 0.0;					/*Read value of a*/
  b= 50.0;					/*Read value of b*/
  h= 0.05;					/*Read value of h*/

  lambda=2.0;					/*Read value of lambda*/

  n=(b-a)/h+1;			    		/*array Dimension*/
  ktot= 35;					/*array Dimension*/

  A=1;						/*Read value of Lattice Constant*/
  k=5; 			    			/*Dummy Array for momentum*/

  sigma=10.0;
  Tsyst=1.0;




/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The Module for Dyson Iteration using Euler Method:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


//for (k=0; k<ktot; k=k+1){

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Momentum values at Dummy Indicies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

 K=((-3.1418/A)+k*(6.2836/A)/(ktot-1.0));

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


float akka;
float complex P;

/*
for (Tbath=0.1 ; Tbath <= 10 ;  Tbath= Tbath+h) {

P=0.0;

for(akka=-100; akka <= 100; akka=akka+h){
P = P + h*(-lambda*lambda*akka*exp(-(akka*akka/(sigma*sigma))))/((2.0*akka*akka-2.0*omega*omega-lambda*lambda*(sigma/sqrt(3.14159265))+lambda*lambda*(2.0/sqrt(3.14159265))*akka*Dawson(akka/sigma))*(2.0*akka*akka-2.0*omega*omega-lambda*lambda*(sigma/sqrt(3.14159265))+lambda*lambda*(2.0/sqrt(3.14159265))*akka*Dawson(akka/sigma))+(lambda*lambda*lambda*lambda*akka*akka*exp(-2.0*(akka*akka/(sigma*sigma)))))*(1.0/tanh((akka/(2.0*Tbath))));
}
DKthermal = P/(3.14159265);




printf("%f\t%f\n", Tbath, -crealf(DKthermal));


}

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*Initial Condition*/

for (t=a; t<=b; t=t+2.0*h){


IR[(int)(t/h)][(int)(t/h)]=0;
DR[(int)(t/h)][(int)(t/h)]=0;


IR[(int)((t+h)/h)][(int)(t/h)]=0;
DR[(int)((t+h)/h)][(int)(t/h)]=-2.0*DzeroR(omega,t+h,t)*BarDR[(int)(t/h)][(int)(t/h)];

BarDR[(int)(t/h)][(int)(t/h)]= -1/2.0;
BarDR[(int)((t+h)/h)][(int)(t/h)]= BarDzeroR(omega,t+h,t);


DK[(int)(t/h)][(int)(t/h)]=-(I/(2.0*omega))*(1.0/(tanh((omega)/(2*Tsyst))));
BarDK[(int)(t/h)][(int)(t/h)]= 0;

DR[(int)((t+2.0*h)/h)][(int)(t/h)]= -2.0*DzeroR(omega,t+2.0*h,t)*BarDR[(int)(t/h)][(int)(t/h)];
DR[(int)((t+3.0*h)/h)][(int)(t/h)]= -2.0*DzeroR(omega,t+3.0*h,t+h)*BarDR[(int)((t+h)/h)][(int)(t/h)]-2.0*BarDzeroR(omega,t+3.0*h,t+h)*DR[(int)((t+h)/h)][(int)(t/h)];


IR[(int)((t+2.0*h)/h)][(int)(t/h)]=h*lambda*lambda*SigmaR(t+2.0*h,t+h)*DR[(int)((t+h)/h)][(int)(t/h)]-(1.0/(2.0*sqrt(3.14159265)))*lambda*lambda*sigma*DR[(int)((t+2.0*h)/h)][(int)(t/h)];
IR[(int)((t+3.0*h)/h)][(int)(t/h)]=h*lambda*lambda*SigmaR(t+3.0*h,t+2.0*h)*DR[(int)((t+2.0*h)/h)][(int)(t/h)]+h*lambda*lambda*SigmaR(t+3.0*h,t+h)*DR[(int)((t+h)/h)][(int)(t/h)]-(1.0/(2.0*sqrt(3.14159265)))*lambda*lambda*sigma*DR[(int)((t+3.0*h)/h)][(int)(t/h)];


DK[(int)((t+2.0*h)/h)][(int)(t/h)]= -2.0*BarDzeroR(omega,t+2.0*h,t)*DK[(int)(t/h)][(int)(t/h)];
DK[(int)((t)/h)][(int)((t+2.0*h)/h)]=-conjf(DK[(int)((t+2.0*h)/h)][(int)(t/h)]);


IK1[(int)(t/h)][(int)(t/h)]=0;
IK2[(int)(t/h)][(int)(t/h)]=0;
//IK1[(int)((t+2.0*h)/h)][(int)(t/h)]=(h/2.0)*lambda*lambda*SigmaR(t+2.0*h,t)*DK[(int)(t/h)][(int)(t/h)]-(1.0/(2.0*sqrt(3.14159265)))*lambda*lambda*sigma*DK[(int)((t+2.0*h)/h)][(int)(t/h)];
//IK2[(int)((t+2.0*h)/h)][(int)(t/h)]=(h/2.0)*lambda*lambda*SigmaK(t+2.0*h,t+2.0*h)*DR[(int)((t)/h)][(int)(t+2.0*h/h)];
}

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The Retarded Part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

P=0.0;
int i;


	for (tprime=a; tprime<b; tprime=tprime+h){
	for (t=tprime+2.0*h; t<b; t=t+2.0*h){


	 	BarDR[(int)(t/h)][(int)(tprime/h)] = -2.0*BarDzeroR(omega,t,t-2.0*h)*BarDR[(int)((t-2.0*h)/h)][(int)(tprime/h)]+2.0*omega*omega*DzeroR(omega,t,t-2.0*h)*DR[(int)((t-2.0*h)/h)][(int)(tprime/h)]+(h/6.0)*(2.0*BarDzeroR(omega,t,t-h)*IR[(int)((t-h)/h)][(int)(tprime/h)]);
	 	DR[(int)((t+2.0*h)/h)][(int)(tprime/h)]= -2.0*DzeroR(omega,t+2.0*h,t)*BarDR[(int)(t/h)][(int)(tprime/h)]-2.0*BarDzeroR(omega,t+2.0*h,t)*DR[(int)(t/h)][(int)(tprime/h)]+(h/6.0)*(2.0*DzeroR(omega,t+2.0*h,t+h)*IR[(int)((t+h)/h)][(int)(tprime/h)]);

    BarDR[(int)((t+h)/h)][(int)(tprime/h)] = -2.0*BarDzeroR(omega,(t+h),(t-h))*BarDR[(int)((t-h)/h)][(int)(tprime/h)]+2.0*omega*omega*DzeroR(omega,(t+h),(t-h))*DR[(int)((t-h)/h)][(int)(tprime/h)]+(h/6.0)*(4.0*BarDzeroR(omega,(t+h),(t))*IR[(int)((t)/h)][(int)(tprime/h)]+BarDzeroR(omega,(t+h),(t-h))*IR[(int)((t-h)/h)][(int)(tprime/h)]);
    DR[(int)((t+3.0*h)/h)][(int)(tprime/h)]= -2.0*DzeroR(omega,t+3.0*h,(t+h))*BarDR[(int)((t+h)/h)][(int)(tprime/h)]-2.0*BarDzeroR(omega,t+3.0*h,(t+h))*DR[(int)((t+h)/h)][(int)(tprime/h)]+(h/6.0)*(4.0*DzeroR(omega,t+3.0*h,(t+2.0*h))*IR[(int)((t+2.0*h)/h)][(int)(tprime/h)]+DzeroR(omega,t+3.0*h,(t+h))*IR[(int)((t+h)/h)][(int)(tprime/h)]);


				for (i = (int)((tprime+2.0*h)/h); i <= (int)((t)/h); i++)
				{
					P=P+h*lambda*lambda*SigmaR(t+2.0*h,(h*i))*DR[i][(int)(tprime/h)];
				}

		IR[(int)((t+2.0*h)/h)][(int)(tprime/h)]=P+(1.0/(sqrt(3.14159265)))*lambda*lambda*sigma*DR[(int)((t+2.0*h)/h)][(int)(tprime/h)];
		P=0.0;

  				for (i = (int)((tprime+3.0*h)/h); i <= (int)((t+h)/h); i++)
  				{
  					P=P+h*lambda*lambda*SigmaR(t+3.0*h,(h*i))*DR[i][(int)(tprime/h)];
  				}

  		IR[(int)((t+3.0*h)/h)][(int)(tprime/h)]=P+(1.0/(sqrt(3.14159265)))*lambda*lambda*sigma*DR[(int)((t+3.0*h)/h)][(int)(tprime/h)];
  		P=0.0;
  	}
  	}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The Advanced Part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*

for (tprime=a; tprime<=b; tprime=tprime+2.0*h){
for (t=a; t<=b; t=t+2.0*h){
		DA[(int)((tprime)/h)][(int)(t/h)]=conjf(DR[(int)(t/h)][(int)(tprime/h)]);
}
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*The Keldysh Part*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*

	for (tprime=a; tprime<b; tprime=tprime+2.0*h){
	for (t=a+h; t<b; t=t+2.0*h){
	  BarDK[(int)(t/h)][(int)(tprime/h)] = -2.0*BarDzeroR(omega,t,t-2.0*h)*BarDK[(int)((t-2.0*h)/h)][(int)(tprime/h)]+2.0*omega*omega*DzeroR(omega,t,t-2.0*h)*DK[(int)((t-2.0*h)/h)][(int)(tprime/h)]+(h/2.0)*BarDzeroR(omega,t,t-2.0*h)*(IK1[(int)((t-2.0*h)/h)][(int)(tprime/h)]+IK2[(int)((t-2.0*h)/h)][(int)(tprime/h)]);
		DK[(int)((t+2.0*h)/h)][(int)(tprime/h)]= -2.0*DzeroR(omega,t+2.0*h,t)*BarDK[(int)(t/h)][(int)(tprime/h)]-2.0*BarDzeroR(omega,t+2.0*h,t)*DK[(int)(t/h)][(int)(tprime/h)]+(h/2.0)*DzeroR(omega,t+2.0*h,t)*(IK1[(int)(t/h)][(int)(tprime/h)]+IK2[(int)(t/h)][(int)(tprime/h)]);


				for (i = (int)((a+h)/h); i <= (int)((t)/h); i++)
				{
					P=P+h*lambda*lambda*SigmaR(t+2.0*h,(h*i))*DK[i][(int)(tprime/h)];
				}

		IK1[(int)((t+2.0*h)/h)][(int)(tprime/h)]=P+(h/2.0)*lambda*lambda*SigmaR(t+2.0*h,tprime)*DK[(int)(tprime/h)][(int)(tprime/h)]+(1.0/(sqrt(3.14159265)))*lambda*lambda*sigma*DK[(int)((t+2.0*h)/h)][(int)(tprime/h)];
		P=0.0;


						for (i = (int)((a+h)/h); i < (int)((tprime)/h); i++)
						{
							P=P+h*lambda*lambda*SigmaK(t+2.0*h,(h*i))*DA[i][(int)(tprime/h)];
						}

		IK2[(int)((t+2.0*h)/h)][(int)(tprime/h)]=P+(h/2.0)*lambda*lambda*SigmaK(t+2.0*h,t+2.0*h)*DA[(int)((t+2.0*h)/h)][(int)(tprime/h)];
		P=0.0;

		DK[(int)((tprime)/h)][(int)((t+2.0*h)/h)]=-conjf(DK[(int)((t+2.0*h)/h)][(int)(tprime/h)]);
		BarDK[(int)((tprime)/h)][(int)(t/h)]=conjf(BarDK[(int)(t/h)][(int)(tprime/h)]);

	}
	}






/*The Occupation Number*/
/*
		for (t=a; t<=b; t=t+2.0*h){
		N[(int)((t)/h)]=0.5*(-2*omega*cimagf(DK[(int)(t/h)][(int)(t/h)])-1);
		}



/*The Statistics*/
/*
float VarE[n], EE[n], E[n], Ntot[n];-2*omega*Dawson(omega/(sqrt(2)*sigma))/sqrt(3.14159265)
float EEPr, EPr, NPr;

EPr=0.0;
EEPr=0.0;
NPr=0.0;

, cimagf(DR[0][(int)(t/h)][0])


		for ( t = a; t <=b ; t=t+2.0*h) {-2*omega*Dawson(omega/(sqrt(2)*sigma))/sqrt(3.14159265)
		for ( k = 0; k <ktot ; k=k+1) {
					NPr=NPr+N[k][(int)((t)/h)];
				}
				Ntot[(int)((t)/h)]=NPr;
				NPr=0.0;
			}


			for ( t = a; t <=b ; t=t+2.0*h) {
			for ( k = 0; k <ktot ; k=k+1) {
				K=((-3.1418/A)+k*(6.2836/A)/(ktot-1));
				omega=-2*2*cos(K);
				EPr=EPr+(omega)*N[k][(int)((t)/h)];
					}
					E[(int)((t)/h)]=EPr;
					EPr=0.0;
				}

			for ( t = a; t DzeroR(<=b ; t=t+2.0*h) {
			for ( k = 0; k <ktot ; k=k+1) {
						EEPr=EEPr+(omega)*(omega)*N[k][(int)((t)/h)];
						}
						EE[(int)((t)/h)]=EEPr;
						EEPr=0.0;
				}

			for ( t = a+h; t <=b ; t=t+2.0*h) {
				printf("%f\t%f\n", t, (E[(int)((t)/h)]/Ntot[(int)((t)/h)]) );
			}


//}

*/
/*	for (t=a; t<=b; t=t+2.0*h){
			 printf("%f\t%f\n", t, DK[0][(int)((t)/h)][(int)((t)/h)] );
		 }

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K loop ends here.
Put K independent printf statements after this.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
//printf("\n%f\t%f", omega, crealf(DKthermal))  ;


for (t = a+h; t < (b-h); t=t+2*h){
printf("%f\t%f\n", t , crealf(DR[(int)(t/h)][(int)(a/h)]))  ;
}
}
