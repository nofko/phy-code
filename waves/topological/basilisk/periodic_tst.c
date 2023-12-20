
#include "grid/multigrid1D.h"
#include "green-naghdi.h"


double h0 = -0.055;
double h1 = -0.015;

double LL = 0.25;
double d  = 0.1;
double x0 = 0.5;
double pd = 0.0205;
  
double A = 0.0;
double f = 0.0;

  
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
	    
	    x < x0+LL        ? h0 :
	    x < x0+pd+LL     ? h0 + (x-x0-LL)*slope :
	    x < x0+d-pd+LL   ? h1 :
	    x < x0+d+LL      ? h0 - (x-x0-d-LL)*slope :

	    x < x0+LL*2        ? h0 :
	    x < x0+pd+LL*2     ? h0 + (x-x0-LL*2)*slope :
	    x < x0+d-pd+LL*2   ? h1 :
	    x < x0+d+LL*2      ? h0 - (x-x0-d-LL*2)*slope :

	    x < x0+LL*3        ? h0 :
	    x < x0+pd+LL*3     ? h0 + (x-x0-LL*3)*slope :
	    x < x0+d-pd+LL*3   ? h1 :
	    x < x0+d+LL*3      ? h0 - (x-x0-d-LL*3)*slope :

	    x < x0+LL*4        ? h0 :
	    x < x0+pd+LL*4     ? h0 + (x-x0-LL*4)*slope :
	    x < x0+d-pd+LL*4   ? h1 :
	    x < x0+d+LL*4      ? h0 - (x-x0-d-LL*4)*slope :

	    x < x0+LL*5        ? h0 :
	    x < x0+pd+LL*5     ? h0 + (x-x0-LL*5)*slope :
	    x < x0+d-pd+LL*5   ? h1 :
	    x < x0+d+LL*5      ? h0 - (x-x0-d-LL*5)*slope :

	    x < x0+LL*6        ? h0 :
	    x < x0+pd+LL*6     ? h0 + (x-x0-LL*6)*slope :
	    x < x0+d-pd+LL*6   ? h1 :
	    x < x0+d+LL*6      ? h0 - (x-x0-d-LL*6)*slope :

	    x < x0+LL*7        ? h0 :
	    x < x0+pd+LL*7     ? h0 + (x-x0-LL*7)*slope :
	    x < x0+d-pd+LL*7   ? h1 :
	    x < x0+d+LL*7      ? h0 - (x-x0-d-LL*7)*slope :

	    x < x0+LL*8        ? h0 :
	    x < x0+pd+LL*8     ? h0 + (x-x0-LL*8)*slope :
	    x < x0+d-pd+LL*8   ? h1 :
	    x < x0+d+LL*8      ? h0 - (x-x0-d-LL*8)*slope :

	    
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


event output (i += 1; t <= 12) {
  
  foreach()
    fprintf (stdout, "%g ", eta[]);
    fprintf (stdout, "\n");
  
}

