@if _XOPEN_SOURCE < 700
  @undef _XOPEN_SOURCE
  @define _XOPEN_SOURCE 700
@endif
@if _GNU_SOURCE
@include <stdint.h>
@include <string.h>
@include <fenv.h>
@endif
#define _CATCH
#define dimension 1
#include "common.h"
#ifndef BASILISK_HEADER_0
#define BASILISK_HEADER_0
#line 1 "bar.c"
#include "grid/multigrid1D.h"
#include "green-naghdi.h"


double h0 = 0.055;
double h1 = 0.015

double LL = 0.25;
double d  = 0.08;
double x0 = 0.5;
double pd = 0.01;
double slope = (h1-h0)/pd
  
double A = 0.003;

event init (i = 0)
{

  u.n[left]  = - radiation (A*sin(2.*pi*t*1.5));
  u.n[right] = + radiation (0);

  foreach() {
    
    zb[] = (x < x0        ? -0.055 :
	    x < x0+pd     ? -0.055 + (x-x0)*slope :
	    x < x0+d-2*pd ? -0.015 :
	    x < x0+d      ? -0.055 - (x-x0+d)*slop :
	    
	    x < x0+L        ? -0.055 :
	    x < x0+pd+L     ? -0.055 + (x-x0+L)*slope :
	    x < x0+d-2*pd+L ? -0.015 :
	    x < x0+d+L      ? -0.055 - (x-x0+d+L)*slope :

	    
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

#endif
