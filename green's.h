/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Definition for function calls
   for all Analytical Green's functions
   in Phononic and Electronic Sectors,
   Baths and Self Energies

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/



   /*Function call for DzeroR*/

    float complex DzeroR(float omega, float t, float tprime) {
        float complex D;
   		  D = (1/omega)*sin(omega*(t-tprime));
   		  return D;
   	}



   /*Function call for BarDzeroR*/

     float complex BarDzeroR(float omega, float t, float tprime) {

   		  float complex D;
   		  D = cos(omega*(t-tprime));
   		  return D;
   	}




    /*Function call for DzeroK*/

     float complex DzeroK(float omega, float t, float tprime) {

     		float complex D;
        float Tsyst;
        Tsyst=100;
        	D = (1/tanh((omega)/(2*Tsyst)));
        	return D;
     	}



    /*Function call for SigmaR*/

     float complex SigmaR(float omega, float t, float tprime) {

              float complex D;
              float sigma;
              sigma = 8.0;
              D = 2*omega*Dawson(omega/(sqrt(2)*sigma))/sqrt(3.14159265)-I*omega*exp(-(omega*omega)/(sigma*sigma));
              return D;
              }

    /*Function call for SigmaK*/

    float complex SigmaK(float omega, float t, float tprime) {

             float complex D;
             float sigma, Tbath;
             sigma = 8.0;
             Tbath = 10;
             D = -2*I*omega*exp(-(omega*omega)/(sigma*sigma))/tanh((omega/(2*Tbath)));
             return D;
             }
