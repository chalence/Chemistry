
#include <math.h>
//using namespace ChemNS;

//ARS adding include statements for high-density molecular cooling rates
#if DO_CHEMISTRY == PRIMORDIAL
//  #include "mol_data.h"
//  #include "atomic_data.h"
//  #include "h2_opac.cpp"
  void compute_h2_opacity(double temp, double N_H2_eff, double &H2_opacity_correction);
#endif
////////////////////////////////

// external fortran subroutines:
extern "C" {
  void co_cooling_(const double*, const double*, const double*, const double*, const double*, const double*, const double*, const double*, double*);

  void co_shielding_(const double*, const double*, double*, double*);

  void fine_structure_cooling_(const double* ,const double* ,const double* ,const double* ,const double* ,const double* ,const double* ,const double* ,const double* ,const double*, double* );
}

void Chem_EOS(int mode, double* eosIn, double* y, int nSpec, int loc);

void AbortChem();


// determine electron abunance from charge conservation:
// network dependent
void fill_species(double* ymass, double* ytot);


void Cool_function(double ei, double* ytot, double* ratraw, double* rpar, double &cool_sum);  

void Cool_gas(double* eosIn, double* ymass, double* ratraw, double* rpar, double dt);


void Chem_Burner(double dt, double* eosIn, double* yIn, double* yOut, double* rpar);
