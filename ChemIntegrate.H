#include <stdlib.h>
#include <iomanip>
#include <math.h>



// Solves matrix equation:
//        a * x = b
// and puts solution in b
void leqs(REAL (*a)[NINT],REAL* b);

void ChemIntegrate(REAL start,REAL firstDT,REAL tmin,REAL dt,REAL* eosIn, REAL* ymass, void (*derivs)(REAL,REAL*,REAL*,REAL*, double* ), double* rpar, REAL* RatRaw, int &jcounts);

void StepMe(REAL* y,REAL* dydx,REAL &x,REAL htry,REAL &hdid,REAL &hnext,REAL* yscal, void (*derivs)(REAL,REAL*,REAL*,REAL*,double*), int &jcounts, REAL* ratraw, double* rpar);



// Rosenbrock ODE semi-implicit method for a single variable 
// (specific internal energy in our case)
void StepMeSingle(REAL &y, REAL &dydx, REAL &x,REAL htry,REAL &hdid,REAL &hnext,
                  REAL &yscal, void (derivs)(REAL,REAL*,REAL*,REAL*,REAL&),
                  int &jcounts,REAL* ytot, REAL* ratraw, double* rpar, REAL escale);
