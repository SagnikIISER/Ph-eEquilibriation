float Dawson( float x )
  /*------------------------------------------------------------*/
  /* PURPOSE: Evaluate Bessel function of first kind and order  */
  /*          1 at input x                                      */
  /*------------------------------------------------------------*/
  {
     float h , P;
     int i;
     h=0.1;
     P=0.00;

     for (i=0; i < (x/h); i++) {
        P=P+h*exp(i*i*h*h);
     }
     P=P+h+(h*2)*exp(x*x);
     P=P*exp(-x*x);
     return P;
  }
