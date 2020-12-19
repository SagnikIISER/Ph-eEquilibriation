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

   /*Function call for SigmaR*/

   float SigmaR(float omega, float t, float tprime) {

       float D;
       if (t==tprime) {
         D = 1;
       }
       else {D=0;}
       return 0;
     }
