

void AbortChem();

void ChemDriver(const double &dt_in, const double &den_in, double &ei_in, double &gamma, double* y, double* rpar);

void Chem_EOS(int mode, double* eosIn, double* y, int nSpec, int loc);
