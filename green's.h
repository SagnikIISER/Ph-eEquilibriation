/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Definition for Function Calls
   for all Analytical Green's functions
   in Phononic and Electronic Sectors,
   Baths and Self Energies for Test Module
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/



/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Phononic Bare Green's Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/****************************
   Function call for DzeroR
****************************/

    double complex DzeroR (double omega, double t, double tprime)
        {
        double complex D;
        if (t<tprime)
        {
            D=0;
        }

        else{
          D = -(1/(2.0*omega))*sin(omega*(t-tprime));
        }

   		  return D;
   	}


/*******************************
   Function call for BarDzeroR
*******************************/

     double complex BarDzeroR(double omega, double t, double tprime)
        {
        double complex D;
        if (t<tprime)
        {
           D=0;
        }
        else{
           D = -cos(omega*(t-tprime))/2.0;
         }

       return D;
       }


/*******************************
    Function call for DzeroK
 ******************************/

     double complex DzeroK(double omega, double Tsyst, double t, double tprime) {

     		double complex D;
        	D = -(I/(2.0*omega))*((cos((omega*(t-tprime)))/(tanh(omega/(2.0*Tsyst)))));
        	return D;
     	}




/*******************************
   Function call for BarDzeroK
********************************/


      double complex BarDzeroK(double omega, double Tsyst, double t, double tprime) {

          double complex D;
          D = (I/2.0)*((sin((omega*(t-tprime)))/(tanh((omega)/(2*Tsyst)))));
          return D;
      }



/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Electronic Bare Green's Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/****************************
   Function call for GzeroR
****************************/

      double complex GzeroR (double omega, double t, double tprime)
          {
          double complex G;
          if (t<tprime)
          {
             G=0;
          }

          else{
             G = -I*cos(omega*(t-tprime))-sin(omega*(t-tprime));
              }

         		return G;
         	}


/*******************************
  Function call for GzeroK
 ******************************/

       double complex GzeroK(double omega, double Tsyst, double nu, double t, double tprime) {

        		double complex G;
              	G =(I*cos(omega*(t-tprime))+sin(omega*(t-tprime)))*(tanh((omega-nu)/(2.0*Tsyst)));
              	return G;
           	}




/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Function Call For Self Energies
for the Test Module
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/



/*******************************
   Function call for SigmaR
********************************/

     double complex SigmaR( double t, double tprime)
     {

              double complex D;
              double sigma;

              sigma = 2.0;

              if (t<=tprime) {
                  D=0;
                }
         		  else{
                //D= -(1.0/(2.0*sqrt(3.14159265)))*sigma*sigma*sigma*(t-tprime)*exp(-sigma*sigma*(t-tprime)*(t-tprime)/4.0);
                D = (2*I/(3.14159265*sigma))*j1(-2.0*sigma*(t-tprime))/(t-tprime);
              }

              return D;
         	}



/*********************************
   Function call for SigmaK
*********************************/

    double complex SigmaK( double t, double tprime) {


                                 double complex D,P;
                                 double sigma,Tbath, nu;
                                 double akka;
                                 double h;

                                 h=0.02;
                                 P=0.0;
                                 sigma = 2.0;
                                 Tbath = 0.0;
                                 nu = 0.5;

                                 for(akka=-2.0*sigma; akka <= 2.0*sigma; akka=akka+h){
                                 P = P+ I*h*((2.0/sigma)*sqrt(1-(akka*akka)/(4*sigma*sigma)))*(tanh(((akka-nu)/(2.0*Tbath))))*((cos(akka*(t-tprime)))-I*sin(akka*(t-tprime)));
                                 }
                                 D = P/(3.14159265*3.14159265);
                                 P = 0.0;

             return D;
             }



/*%%%%%%%%%%%%%%%%%%%%%%
    Fermi functions
***********************/

double Fermi( double epsilon, double beta, double nu) {
float n;
n = 1/(exp(beta*(epsilon-nu))+1);
return n;
}



double FermiD1( double epsilon, double beta, double nu) {
float n;
n = -(epsilon-nu)*exp(beta*(epsilon-nu))/((exp(beta*(epsilon-nu))+1)*(exp(beta*(epsilon-nu))+1));
return n;
}



double FermiD2( double epsilon, double beta, double nu) {
float n;
n = beta*exp(beta*(epsilon-nu))/((exp(beta*(epsilon-nu))+1)*(exp(beta*(epsilon-nu))+1));
return n;
}
