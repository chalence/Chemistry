
#include <iostream>
#include <iomanip>
using namespace std;

#include "ChemGlobals.H"

#include "ChemDriver.H"
#include "ChemBurner.H"
using namespace ChemNS;



double yin_original[NINT];
double rpar_original[nrpar];
double ei_original;
double dt_original;
double dens_original;


extern double yin_original[NINT];
extern double rpar_original[nrpar];
extern double ei_original;
extern double dt_original;
extern double dens_original;

void ChemDriver(const double &dt_in, const double &den_in, double &ei_in, double &gamma_here, double* y_fromOrion, double* rpar)
{
    int i;
    double smallx = 1.0e-20;
    
    // Initialize arrays (xIN will be read in, actually)

double xIn[NINT] = {0.0};
double yIn[NINT] = {0.0};
double yOut[NINT] = {0.0};
double ytot[NSPECIES];

    double eosIn[6];  // Density, Temperature, EngDen, Pressure, Metal Fraction, gamma

    // save original values for use in AbortChem
    for(i=0;i<NINT;i++) yin_original[i] = y_fromOrion[i];
    for(i=0;i<nrpar;i++) rpar_original[i] = rpar[i];
    ei_original = ei_in;
    dt_original = dt_in;
    dens_original = den_in;
    //

    
    // Initialization parameters
    // and other parameters that might become read in
    // CTSS: deprecated, J21 is now in rpar array
    // and shielding factors should be determined consistently
    //bool chem_cc_case = 1;
    //double fssh2 = 1.0;
    //double fsshd = 1.0;
    //double J21 = 0.0;  // J21 radiation
    //double RadParams[3] = {J21,fssh2,fsshd};
    
    
    // Load from Orion2
    for(i=0;i<NINT;i++) xIn[i] = y_fromOrion[i];
    double dt = dt_in;
    eosIn[0] = den_in;
    eosIn[1] = 0.0; //Temp
    eosIn[2] = ei_in;
    eosIn[3] = 0.0; //P
    eosIn[4] = 0.0; // mass frac
    eosIn[5] = gamma_here;
    rpar[gamma_par] = gamma_here;
    rpar[rho_par] = den_in;
    rpar[e0_par] = ei_in;
    rpar[temp_par] = -1.0;

    // check initial abundances
    for(i=0;i<NINT;i++)
    {
      if(xIn[i] < -1.0){
	cout << "big negative abundance upon entering chemistry! " << i << " " << xIn[i] << endl;
	AbortChem();
      }
    }
    if(xIn[iH] > 1.0){
      if(xIn[iH] > 2.0){
	cout << "H abundance too high upon entering chemistry" << endl;
	AbortChem();
      }
      xIn[iH] = 1.0 - xIn[iHP] - 2.0 * xIn[iH2];
    }
    if(xIn[iH2] > 0.5){
      if(xIn[iH2] > 1.0){
	cout << "too H2 abundance too high upon entering chemistry" << endl;
	AbortChem();
      }
      xIn[iH2] = 0.5 * (1.0 - xIn[iH] - xIn[iHP]);
    }
    
      
    

 
    // yIn = abundance array of tracked species (size NINT)
    // ytot = abundance array of ALL species (including electrons, size NINT+1)
    //double ytot[NSPECIES];
    for(i=0;i<NINT;i++)
    {
      yIn[i] = xIn[i] ; 
	ytot[i] = yIn[i];
        if(yIn[i]< smallx ) yIn[i] = smallx;
    }
 

    fill_species(yIn,ytot);
    Chem_EOS(2,eosIn,ytot,NSPECIES,1);
    rpar[temp_par] = eosIn[1];
    rpar_original[temp_par] = rpar[temp_par];
    
 
    // Call BURNER
    Chem_Burner(dt,eosIn,yIn,yOut,rpar);
    
 
    fill_species(yIn,ytot);
    Chem_EOS(1,eosIn,ytot,NSPECIES,2);
        
    for(i=0;i<NINT;i++) y_fromOrion[i] = yOut[i]; 
    //den_in = eosIn[0];
    ei_in = eosIn[2];    
   
    //ARS update gamma_here to new value! 
    gamma_here = eosIn[5];
    

// check final abundances
    for(i=0;i<NINT;i++)
    {
      if(yOut[i] < 0.0){
	cout << "negative abundance upon leaving chemistry! " << i << " " << xIn[i] << endl;
	AbortChem();
      }
    }


    if(yOut[iH] > 1.5){
      cout << "H abundance too high upon leaving chemistry" << endl;
      cout << "current abundances:" << endl;
      for(i=0;i<NINT;i++){
	cout << i << " " << yOut[i] << endl;
      }
      AbortChem();
    }
    if(yOut[iH2] > 1.0){
      cout << "H2 abundance too high upon leaving chemistry" << endl;
      AbortChem();
    }
    
    
    if(eosIn[1] > 1.0e9){
      cout << "leaving chemistry with temp too high, > 1.0e9. Temp = " << eosIn[1] << endl;
      AbortChem();
    }
    

}



void Chem_EOS(int mode, double* eosIn, double* y, int nSpec, int loc)
{
    //double mh = 1.67e-24;
    //double kb = 1.38e-16;
    double yhetot, fac;
    
    double den,temp,ei,pres,gc;
    den=eosIn[0];
    temp=eosIn[1];
    ei=eosIn[2];
    pres=eosIn[3];
    gc=eosIn[5];
    
    // CTSS: P = ntot * kb * T = rho * kb * T / (mu * mh)
    //       P = E * (gamma-1), where E = volumetric energy density
    //       P = ei / rho * (gamma-1), where ei = specific energy density
    //       To a very good approximation, ntot = (1 + yhe + ye - yh2) * nh
    //       where nh is the number density of hydrogen nuclei given by
    //       nh = rho / (mh * (1.0 + 4*yhe)) 
    //       ...which basically describes what's being done below:

    yhetot = y[iHE] + y[iHEP];
    fac = (1.0 + yhetot + y[iELEC] - y[iH2]) * kboltz / ((1.0 + 4.0*yhetot) * mh * (gc-1.0));

    if(mode==1) // determines ei and pressure from rho and T
    {
      ei = fac * temp;
      pres = (gc-1.0) * den * ei;
      
      eosIn[0]=den;
      eosIn[1]=temp;
      eosIn[2]=ei;
      eosIn[3]=pres;
    }
    else if(mode==2) //determines T and P from rho and ei
    {
      temp = ei / fac;
      pres = (gc-1.0) * den * ei;
      
      eosIn[0]=den;
      eosIn[1]=temp;
      eosIn[2]=ei;
      eosIn[3]=pres;
    }
    else if(mode==3) // determines T and ei from rho and P
    {
      // seems like this mode is never used...
      cout << "EOS mode 3 not done, and never really needed..." << endl;
      AbortChem();
    }
    else // doesn't change anything... not good.
    {
      cout << "Bad option in Chem_EOS " << mode << endl;
      AbortChem();
        
    }
    
    if(eosIn[1] != eosIn[1] || eosIn[2] != eosIn[2]) {
      cout << "nan after EOS! mode = " << mode << endl;
      cout << "loc = " << loc << endl;
      cout << "eos[0] = " << eosIn[0] << endl;
      cout << "eos[1] = " << eosIn[1] << endl;
      cout << "eos[2] = " << eosIn[2] << endl;
      cout << "gc = " << gc << endl;
      AbortChem();
    }
    
    if(temp > 1.0e9){
      cout << "temp too high after EOS. temp = " << temp << endl;
      cout << "mode = " << mode << endl;
      cout << "loc = " << loc << endl;
      cout << "eos[0] = " << eosIn[0] << endl;
      cout << "eos[1] = " << eosIn[1] << endl;
      cout << "eos[2] = " << eosIn[2] << endl;
      cout << "yh2 = " << y[iH2] << endl;
      cout << "yie = " << y[iELEC] << endl;
      AbortChem();
    
    }

}



void AbortChem()
{ 
  int i;
  cout << "Aborting chemistry!" << endl;
  cout << "original rpar: " << endl;
  for(i=0;i<nrpar;i++) cout << i << " " << rpar_original[i] << endl;
  cout << "original abundances:" << endl;
  for(i=0;i<NINT;i++) cout << i << " " << yin_original[i] << endl;
  cout << "original ei = " << ei_original << endl;
  cout << "original dt = " << dt_original << endl;
  cout << "original density = " << dens_original << endl;
  cout << "goodbye..."<<endl; 
  exit(0);
}
