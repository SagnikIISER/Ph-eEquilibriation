/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Definition for function calls
   for all Analytical Green's functions
   in Phononic and Electronic Sectors,
   Baths and Self Energies

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/



   /*Function call for DzeroR*/

    float complex DzeroR(float omega, float t, float tprime) {
        float complex D;
        if (t<tprime) {
            D=0;
        }
   		  else{
          D = (1/omega)*sin(omega*(t-tprime));
        }

   		  return D;
   	}


   /*Function call for BarDzeroR*/

     float complex BarDzeroR(float omega, float t, float tprime) {

       float complex D;
       if (t<tprime) {
           D=0;
       }
       else{
           D = cos(omega*(t-tprime));
       }

       return D;
       }


   /*Function call for DzeroR*/

        float complex DzeroA(float omega, float t, float tprime) {
            float complex D;
            if (t>tprime) {
                D=0;
            }
       		  else{
              D = (1/omega)*sin(omega*(tprime-t));
            }

       		  return D;
       	}


    /*Function call for DzeroK*/

     float complex DzeroK(float omega, float t, float tprime) {

     		float complex D;
        float Tsyst;
        Tsyst= 0.5;
        	D = -(I/2*omega)*((cos((-omega*(t-tprime)))/(tanh((omega)/(2*Tsyst)))));//-I*(1/(2*omega))*(1+exp(-(omega/Tsyst))*(cos(omega*(t-tprime))+I*sin(omega*(t-tprime))))
        	return D;
     	}



    /*Function call for SigmaR*/

     float complex SigmaR( float t, float tprime) {


              float complex D,P;
              float sigma;
              float akka;
              float h;

              h=0.1;
              P=0.0;
              sigma = 1.0;

              if (t<tprime) {
                  D=0;
              }
         		  else{
                /*for(akka=-5; akka <= 5; akka=akka+h){
                P = P+ h*(1/(2*3.14159265))*(2*akka*Dawson(akka/(sqrt(2)*sigma))/sqrt(3.14159265)+I*akka*exp(-(akka*akka)/(sigma*sigma)))*(cos(akka*(t-tprime))+I*sin(akka*(t-tprime)));
                }
                D = P;
                P = 0.0;*/
                D= 2*sigma*sigma*sigma*(t-tprime)*exp(-sigma*sigma*(t-tprime)*(t-tprime)/2);
              }

         		  return D;
         	}


          /*Function call for SigmaA*/

           float complex SigmaA( float t, float tprime) {

                           float complex D,P;
                           float sigma;
                           float akka;
                           float h;

                           h=0.1;
                           P=0.0;
                           sigma = 1.0;

                    if (t>tprime) {
                        D=0;
                    }
               		  else{
/*                      for(akka=-5; akka <= 5; akka=akka+h){
                      P = P+ h*(1/(2*3.14159265))*(2*akka*Dawson(akka/(sqrt(2)*sigma))/sqrt(3.14159265)+I*akka*exp(-(akka*akka)/(sigma*sigma)))*(cos(akka*(t-tprime))-I*sin(akka*(t-tprime)));
                      }
                      D = P;
                      P = 0.0;
*/
                      D= 2*sigma*sigma*sigma*(t-tprime)*exp(-sigma*sigma*(t-tprime)*(t-tprime)/2);;
                    }

               		  return D;
               	}

    /*Function call for SigmaK*/

    float complex SigmaK( float t, float tprime) {


                                 float complex D,P;
                                 float sigma,Tbath;
                                 float akka;
                                 float h;

                                 h=0.1;
                                 P=0.0;
                                 sigma = 1.0;
                                 Tbath = 2.0;

                                 for(akka=-5; akka <= 5; akka=akka+h){
                                 P = P+ h*(-2*I*akka*exp(-(akka*akka)/(sigma*sigma)))*(1/tanh((akka/(2*Tbath))))*(cos(akka*(t-tprime))-I*sin(akka*(t-tprime)));
                                 }
                                 D =(P/(2*3.14159265));
                                 P = 0.0;

             return D;
             }
