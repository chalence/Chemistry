
#include <math.h> 
#include <iostream> 
#include <fstream>
#include <string>
#include "ChemDriver.H"
#include "ChemGlobals.H"


using namespace std;

void fill_rpar(double* rpar);
void fill_abundances(double* yIn);
		   
int main(){

  cout << "Testing chemistry..." << endl;

  double dt = 1.0e10;
  double den_in = 1.0e-24;
  double ei = 1.0e9;
  
  double gamma = 5.0/3.0;
  double rpar[nrpar];
  // starting mass fractions;
  double yIn[NINT];

  fill_abundances(yIn);
  fill_rpar(rpar);



  ofstream dump;
  string fout = "ff.dat";
  dump.open((fout).c_str());

  // let's set up a free-fall collapse

  double tff = 0.000476;
  double rho = 1.0e-26;
  double rho_0 = rho;
  double rho_prev = rho_0;
  double rho_dot;
  double step;
  double t = 1.0;
  double rho_stop = 1.0e-18;
  double temperature = 200.0;
  double nh = rho / mh / (1.0 + 4.0*abhe);
  double ntot = nh * (abhe + 1.0);
  ei = temperature / (gamma - 1.0) * ntot * kboltz / rho;
  while(rho < rho_stop) // Interates until dt is advanced
    {
      rho = pow(pow(rho_0,-0.5) - tff * t / 2.0 ,-2.0);
      rho_dot = tff * pow(rho,1.5);
      step = 0.01 * rho / rho_dot;
      
      // adiabatically heat gas here
      ei = ei * pow(rho/rho_prev,gamma-1.0);
      cout << "input ei = " << ei << endl;

      rho_prev = rho;
      
      ChemDriver(step,rho,ei,gamma,yIn,rpar);
      temperature = ei * (gamma - 1.0) /  kboltz * mh * 1.3 / (1.0 + abhe - 0.5);
      cout << "t = " << t << " rho = " << rho << ", temp = " << temperature << " h2 =" << yIn[iH2] << endl;
      dump << rho/1.67e-24 << " " << t << " " << temperature << " " << yIn[iH2] /2.0 * 1.3 << " " << yIn[iHP] << " " << yIn[iH] <<  endl;
      t = t + step;
    }

  dump.close(); 
  cout << "ending mass fractions = " << endl;
  for(int i=0;i<NINT;i++) cout << i << " " << yIn[i] << endl;
  cout << " " << endl;
  cout << "ending T = " << ei * (gamma - 1.0) /  kboltz * mh * 1.2 << endl;
  cout << "ending time = " << t << endl;


  cout << "done with chemdriver!" << endl;

}





void fill_rpar(double* rpar){

  for(int i=0;i<nrpar;i++) rpar[i] = 0.0;
  rpar[g0_par]= 0.000000036;
  rpar[NH_par]= 1e+10;
  rpar[NH2_par]= 1e+10;
  rpar[NCO_par]= 1e+10;
  rpar[tdust_par]= 10;
  rpar[zeta_par]= 1e-25;
  rpar[dx_par]= 4.70837e+13;
  rpar[divv_par]= 1.08538e-08;
  rpar[divv1_par]= 7.07356e-09;
  rpar[rho_par]= 0;
  rpar[vturb_par]= 285909;
  rpar[gradrho_par]= 6.36486e-25;
  rpar[gamma_par]= 0;
  rpar[temp_par]= 1510.46;

}

void fill_abundances(double* yIn){

#if DO_CHEMISTRY == CONTEMPORARY
  double yh, yhm, yh2, yhp, yh2p, yh3p;
  double yhe, yhep, yc, ycp, yo, ychx;
  double yohx, yhcop, yco, ym, ymp;
  yhm = 9.5e-13;
  yh2 = 1.0e-6;
  yhp = 7.5e-5;
  yh2p = 2.92e-13;
  yh3p = 5.92e-14;
  yhe = abhe;
  yhep = 1.46e-5;
  yc = 1.6e-4;
  ycp = 1.0e-7;
  yo = 3.2e-4;
  ychx = 1.3e-12;
  yohx = 4.3e-16;
  yhcop = 4.56e-19;
  yco = 6.0e-10;
  ym = 3.9e-10;
  ymp = 1.99e-7;

  yh = 1.0 - yhp - 2.0 * yh2;

  yIn[iH] = yh;
  yIn[iHM] = yhm ;
  yIn[iHP] = yhp;
  yIn[iH2] = yh2;
  yIn[iH2P] = yh2p ;
  yIn[iH3P] = yh3p;
  yIn[iHE] = yhe;
  yIn[iHEP] = yhep;
  yIn[iC] = yc;
  yIn[iCP] = ycp;
  yIn[iO] = yo;
  yIn[iCHX] = ychx;
  yIn[iOHX] = yohx;
  yIn[iCO] = yco;
  yIn[iHCOP] = yhcop;
  yIn[iM] = ym;
  yIn[iMP] = ymp;
#endif

#if DO_CHEMISTRY == PRIMORDIAL
  yIn[iHP] = 1.0e-4;
  yIn[iH] = 0.99;
  yIn[iHM] = 1.0e-10;
  yIn[iH2P] = 1.0e-10;
  yIn[iH2] = 1.0e-6;
  yIn[iDP] = DeutAbund * 1.0e-4;
  yIn[iD] = DeutAbund;
  yIn[iDM] = 1.0e-10;
  yIn[iHDP] = 1.0e-10;
  yIn[iHD] = 1.0e-10;
  yIn[iD2P] = 1.0e-10;
  yIn[iD2] = 1.0e-10;
  yIn[iHEP] = 1.0e-8;
  yIn[iHE] = abhe;
#endif

}
