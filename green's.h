/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Definition for function calls
   for all Analytical Green's functions
   in Phononic and Electronic Sectors,
   Baths and Self Energies

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/



   /*Function call for DzeroR*/

   float DzeroR(float omega, float t, float tprime) {

   		float D;
   		D = (1/(omega))*sin(omega*(t-tprime));
   		return D;
   	}



   /*Function call for BarDzeroR*/

   float BarDzeroR(float omega, float t, float tprime) {

   		float D;
   		D = cos(omega*(t-tprime));
   		return D;
   	}




    /*Function call for DzeroK*/

     float DzeroK(float omega, float t, float tprime) {

     		float D;

        if (t==tprime) {
        	D = (1/tanh((omega+25)/(2*100)));
        }
        if (t>tprime) {
        	D = (1/omega)*(1/tanh(omega+25/(2*100)))*sin(omega*(t-tprime));
        }
        if (t<tprime) {
          D = -(1/omega)*(1/tanh(omega+25/(2*100)))*sin(omega*(t-tprime));
        }

     		return D;
     	}



    /*Function call for SigmaR*/

    float SigmaR(float omega, float t, float tprime) {

          float D;
          D = (1/(2*10))*Bessel(2*10*(t-tprime))/(t-tprime);
          return D;
          }

    /*Function call for SigmaK*/

    float SigmaK(float omega, float t, float tprime) {


      float  P, Q, tbath, w;
      P=0.0;
      tbath=10;

        for (w = (-2*tbath); w < (2*tbath); w=w+0.1) {
        Q=(1-(w*w/(4*tbath*tbath)));
        P=P+(0.1)*(2/tbath)*sqrt(Q)*(sin((w)*(t-tprime)))*tanh((w+25)/(2*25));
      }

      return P;

      }
