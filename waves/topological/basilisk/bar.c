#include "grid/multigrid1D.h"
#include "green-naghdi.h"


double h0 = -0.055;
double h1 = -0.015;

double LL = 0.25;
double d  = 0.08;
double x0 = 0.5;
double pd = 0.01;
  
double A = 0.014;

  
event init (i = 0)
{

  double slope = (h1-h0)/pd;

  u.n[left]  = - radiation (A*sin(2.*pi*t*0.5));
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


int main() {

  
  N = 2048;
  L0 = 6;
  G = 9.81;
  
  run();

}


event output (i += 2; t <= 20) {
 
  foreach()
      fprintf (stdout, "%g ", eta[]);
    fprintf (stdout, "\n");
  
}
