/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Definition for function calls
   for all Analytical Green's functions
   in Phononic and Electronic Sectors,
   Baths and Self Energies

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

   /*Function call for DzeroR*/

   float DzeroR(float omega, float t, float tprime) {

   		float D;

   		D = (1/(2*omega))*sin(omega*(t-tprime));
   		return D;
   	}

   /*Function call for BarDzeroR*/

   float BarDzeroR(float omega, float t, float tprime) {

   		float D;

   		D = (1/(2))*cos(omega*(t-tprime));
   		return D;
   	}

   /*Function call for SigmaR*/

   float SigmaR(float omega, float t, float tprime) {

       float D;
       D = 0;
       return 0;
     }
