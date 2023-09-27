#include "grid/multigrid1D.h"
#include "green-naghdi.h"


double h0 = -0.055;
double h1 = -0.015;

double LL = 0.25;
double d  = 0.08;
double x0 = 0.5;
double pd = 0.01;
  
double A = 0.0;
double f = 0.0;

  
event init (i = 0)
{

  
  u.n[left]  = - radiation (A*sin(2.*pi*t*f));
  u.n[right] = + radiation (0);

  foreach() {
    
    zb[] = (h0);
    
    h[] = - zb[];
  }

}


int main(int argc, char *argv[]) {

  sscanf(argv[1],"%lf",&A);
  sscanf(argv[2],"%lf",&f);
  
  /*A = atof(argv[0]);
  printf("%4.8f\n", A);*/
    
  N = 2048;
  L0 = 12;
  G = 9.81;
  
  run();

}

void plot_profile (double t, FILE * fp)
{
  fprintf (fp,
	   "set title 't = %.2f'\n"
	   "p [0:25][-0.12:0.04]'-' u 1:3:2 w filledcu lc 3 t ''\n", t);
  foreach()
    fprintf (fp, "%g %g %g\n", x, eta[], zb[]);
  fprintf (fp, "e\n\n");
  fflush (fp);
}

event output (i += 1; t <= 20) {
 
  foreach()
      fprintf (stdout, "%g ", eta[]);
    fprintf (stdout, "\n");
  
}
