#include "pluto.h"
#include "stdlib.h"

#define real double
#define REAL double


void Cool_function(REAL ei, REAL dens, REAL &cool_sum); 

void StepMeSingle(REAL &y, REAL &dydx, REAL &x,REAL htry,REAL &hdid,REAL &hnext,
                  REAL &yscal, void (derivs)(REAL,REAL,REAL&),
                  int &jcounts,REAL dens, REAL escale);

  void CoolDriver(double &ei_in, double dens, double dt_in);

  double Cool_gas(REAL ei, REAL dens, REAL dt);
