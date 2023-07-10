#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/check_eta.h"
#include "layered/perfs.h"



int main() {
  
  N = 2048;
  L0 = 50;
  G = 9.81;

  nl = 2;  
  breaking = 0.1;
  CFL_H = 0.5;

  run();
}



event init (i = 0)
{

  u.n[left]  = - radiation (0.03*sin(2.*pi*t/2.02));
  u.n[right] = + radiation (0);

    foreach() {
    zb[] = (x < 6 ? -0.4 :
	    x < 12 ? -0.4 + (x - 6.)/6.*0.3 :
	    x < 14 ? -0.1 :
	    x < 17 ? -0.1 - (x - 14.)/3.*0.3 :
	    -0.4);
    
    foreach_layer()
      h[] = max(- zb[], 0.)/nl;

    }
}


Gauge gauges[] = {
  {"WG4",  10.5},
  {"WG5",  12.5},
  {"WG6",  13.5},
  {"WG7",  14.5},
  {"WG8",  15.7},
  {"WG9",  17.3},
  {"WG10", 19},
  {"WG11", 21},
  {NULL}
};

event output (i += 2; t <= 40)
  output_gauges (gauges, {eta});
