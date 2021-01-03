


float Dawson( float x )

  {
     float h , P;
     float a;
     h=0.1;
     P=0.0;
     for (a=h; a < fabs(x); a=a+h) {
        P=P+h*exp((a*a));
     }

     if (x>0) {
        P=P+(h/2)+(h/2)*(exp(x*x));
        P=P*exp(-x*x);
     }

      if (x<0) {
        P=P+(h/2)-(h/2)*(exp(x*x));
        P=P*exp(-x*x);
      }


      return P;


  }

/*
  double dawsonf(float x){
          int j,N;
          double r,y,delx=pow(10.0,-4);
          double complex *term;
          double sum=0.0;

          N=(int)(fabs(x/delx))+1;
          term=(double complex *)malloc(N*sizeof(double complex));

          for(j=0;j<N;j++){
                  y=j*delx;
                  r=y*y-x*x;
                  term[j]=exp(r)+I*0.0;
          }// j loop ends
          sum=creal(simpson_comb(term,N,delx));
          if(x<0.0) sum*=-1.0;

          free(term);
          return sum;
        }
*/
