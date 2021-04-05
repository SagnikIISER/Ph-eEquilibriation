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

              sigma = 10.0;

              if (t<tprime) {
                  D=0;
         		  else{
                D= -(1.0/(2.0*sqrt(3.14159265)))*sigma*sigma*sigma*(t-tprime)*exp(-sigma*sigma*(t-tprime)*(t-tprime)/4.0);
              }

              return D;
         	}

        }

/*********************************
   Function call for SigmaRP
*********************************/

    double complex SigmaRP( double t, double tprime) {


             double complex D;
             double sigma;

             sigma = 10.0;

             if (t<tprime) {
             D=0;
             }
        		  else{
                D = (1.0/(2.0*sqrt(3.14159265)))*sigma*sigma*sigma*(sigma*sigma*(t-tprime)*(t-tprime)/2.0-1.0)*exp(-sigma*sigma*(t-tprime)*(t-tprime)/4.0);
            }

            return D;
        }


/*********************************
   Function call for SigmaK
*********************************/

    double complex SigmaK( double t, double tprime) {


                                 double complex D,P;
                                 double sigma,Tbath;
                                 double akka;
                                 double h;

                                 h=0.02;
                                 P=0.0;
                                 sigma =10.0;
                                 Tbath = 1.2;

                                 for(akka=-40.0; akka <= 40.0; akka=akka+h){
                                 P = P+ h*(-2.0*I*akka*exp(-(akka*akka)/(sigma*sigma)))*(1.0/tanh((akka/(2.0*Tbath))))*((cos(akka*(t-tprime)))-I*sin(akka*(t-tprime)));
                                 }
                                 D = P/(2.0*3.14159265);
                                 P = 0.0;

             return D;
             }
