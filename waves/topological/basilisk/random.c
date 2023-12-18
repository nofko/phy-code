#include "grid/multigrid1D.h"
#include "green-naghdi.h"


double h0 = -0.055;
double h1 = -0.015;

double LL = 0.25;
double d  = 0.1;
double pd = 0.0205;
  
double A = 0.0;
double f = 0.0;

double x0 = 0.521;
double x1 = 0.782;
double x2 = 1.011;
double x3 = 1.271;
double x4 = 1.530;
double x5 = 1.831;
double x6 = 1.990;
double x7 = 2.206;
double x8 = 2.575;


event init (i = 0)
{

  double slope = (h1-h0)/pd;

  u.n[left]  = - radiation (A*sin(2.*pi*t*f));
  u.n[right] = + radiation (0);

  foreach() {
    
    zb[] = (x < x0        ? h0 :
	    x < x0+pd     ? h0 + (x-x0)*slope :
	    x < x0+d-pd   ? h1 :
	    x < x0+d      ? h0 - (x-x0-d)*slope :
	    
	    x < x1        ? h0 :
	    x < x1+pd     ? h0 + (x-x1)*slope :
	    x < x1+d-pd   ? h1 :
	    x < x1+d      ? h0 - (x-d-x1)*slope :

	    x < x2        ? h0 :
	    x < x2+pd     ? h0 + (x-x2)*slope :
	    x < x2+d-pd   ? h1 :
	    x < x2+d      ? h0 - (x-d-x2)*slope :

	    x < x3        ? h0 :
	    x < x3+pd     ? h0 + (x-x3)*slope :
	    x < x3+d-pd   ? h1 :
	    x < x3+d      ? h0 - (x-d-x3)*slope :

	    x < x4        ? h0 :
	    x < x4+pd     ? h0 + (x-x4)*slope :
	    x < x4+d-pd   ? h1 :
	    x < x4+d      ? h0 - (x-d-x4)*slope :

	    x < x5        ? h0 :
	    x < x5+pd     ? h0 + (x-x5)*slope :
	    x < x5+d-pd   ? h1 :
	    x < x5+d      ? h0 - (x-d-x5)*slope :

	    x < x6        ? h0 :
	    x < x6+pd     ? h0 + (x-x6)*slope :
	    x < x6+d-pd   ? h1 :
	    x < x6+d      ? h0 - (x-d-x6)*slope :

	    x < x7        ? h0 :
	    x < x7+pd     ? h0 + (x-x7)*slope :
	    x < x7+d-pd   ? h1 :
	    x < x7+d      ? h0 - (x-d-x7)*slope :

	    x < x8        ? h0 :
	    x < x8+pd     ? h0 + (x-x8)*slope :
	    x < x8+d-pd   ? h1 :
	    x < x8+d      ? h0 - (x-d-x8)*slope :

	    
	    -0.055);
    
    h[] = - zb[];
  }

}


event friction (i++) {
  foreach() {
    double a = h[] < dry ? HUGE : 1. + 9e-3*dt*norm(u)/h[];
    foreach_dimension()
      u.x[] /= a;
  }
}

int main(int argc, char *argv[]) {

  sscanf(argv[1],"%lf",&A);
  sscanf(argv[2],"%lf",&f);
  
  N = 1024;
  L0 = 4;
  G = 9.81;
  alpha_d = 1.;
  
  run();

}


event output (i += 1; t <= 15) {
  
  foreach()
    fprintf (stdout, "%g ", eta[]);
    fprintf (stdout, "\n");
  
}
