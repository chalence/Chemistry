
#include <math.h> 
#include <iostream> 
#include <fstream>
#include <string>
#include "ChemGlobals.H"
//#include "ChemExternal.H"

using namespace std;
using namespace ChemNS;
		   
void ChemDriver(const double&, const double&, double&, double&, double*, double*);
	
int main()
{

  cout << "Starting chemistry" <<  endl;
  
  REAL xhe = 0.08;
  REAL gamma = 5.0/3.0;
  
  REAL dt = 1.0e15;
  REAL rho = 1.0e-25;
  REAL temperature = 200.0;
  REAL rho_0 = rho;

  REAL ei, nh, ntot;
  REAL rpar[nrpar];
  
  // starting mass fractions;
  REAL yIn[NINT], ysav[NINT];
  
  for(int i=0;i<nrpar;i++) rpar[i] = 0.0;

#if DO_CHEMISTRY == PRIMORDIAL

  double yhp, yh, yhm, yh2p, yh2, ydp, yd, ydm, yhdp, yhd, yd2p, yd2, yhep, yhe, yhepp;
  yhp = 1.0e-4;
  yhm = 1.0e-15;
  yh2p = 1.0e-15;
  yh2 = 1.0e-6;
  ydp = 1.0e-8;
  yd = 4.0e-5;
  ydm = 1.0e-15;
  yhdp = 1.0e-15;
  yhd = 1.0e-10;
  yh2p = 1.0e-15;
  yd2 = 1.0e-15;
  yhep = 1.0e-10;
  yhe = 0.08;
  yhepp = 1.0e-15;

  yh = 1.0 - yhp - 2.0 * yh2;
  
  yIn[iHP] = yhp ;
  yIn[iH] = yh;
  yIn[iHM] = yhm;
  yIn[iH2P] = yh2p;
  yIn[iH2] = yh2;
  yIn[iDP] = ydp;
  yIn[iD] = yd;
  yIn[iDM] = ydm;
  yIn[iHDP] = yhdp;
  yIn[iHD] = yhd;
  yIn[iD2P] = yd2p;
  yIn[iD2] = yd2;
  yIn[iHEP] = yhep;
  yIn[iHE] = yhe;
  yIn[iHEPP] = yhepp;


  rpar[0]= 0;
  rpar[1]= 0.36;
  rpar[2]= 1e+10;
  rpar[3]= 1e+10;
  rpar[4]= 1e+10;
  rpar[5]= 10;
  rpar[6]= 1e-17;
  rpar[7]= 4.70837e+13;
  rpar[8]= 1.08538e-08;
  rpar[9]= 7.07356e-09;
  rpar[10]= 0;
  rpar[11]= 285909;
  rpar[12]= 6.36486e-25;
  rpar[13]= 0;
  rpar[14]= 1510.46;


  yIn[0]= 3.1726e-18;
  yIn[1]= 4.91182e-11;
  yIn[2]= 4.01178e-17;
  yIn[3]= 1.93232e-15;
  yIn[4]= 0.457598;
  yIn[5]= 9.15196e-21;
  yIn[6]= 9.15196e-21;
  yIn[7]= 9.15196e-21;
  yIn[8]= 1.79413e-16;
  yIn[9]= 1.9925e-05;
  yIn[10]= 2.0126e-16;
  yIn[11]= 1.42614e-10;
  yIn[12]= 9.15196e-21;
  yIn[13]= 0.0553155;
  yIn[14]= 9.15196e-21;


  
#endif

#if DO_CHEMISTRY == CONTEMPORARY

 double yh, yhm, yh2, yhp, yh2p, yh3p, yhe, yhep, yc, ycp, yo, ychx, yohx, yhcop, yco, ym, ymp;
  yhm = 9.5e-13;
  yh2 = 1.0e-6;
  yhp = 7.5e-5;
  yh2p = 2.92e-13;
  yh3p = 5.92e-14;
  yhe = 0.08;
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

 
     //rpar[g0_par] = 0.36;
     //rpar[tdust_par] = 10.0;
     //rpar[NH_par] = 1.0e10;
     //rpar[NH2_par] = 1.0e10;
     //rpar[NCO_par] = 1.0e10;
     //rpar[zeta_par] = 1.0e-17;
     //rpar[divv_par] = 1.3506e-13;
  
  nh = rho / mh / (1.0 + 4.0*xhe);
  ntot = nh * (xhe + 1.0);
  ei = temperature / (gamma - 1.0) * ntot * kboltz / rho;
  cout << "ei = " << ei << endl;
  
  cout << "starting T = " << temperature << endl;

  REAL ei0;
  REAL step;
  REAL tff = 0.000476;
  REAL t=1.0;
  REAL duration=3.0e50; 
  REAL t_start;
  REAL rho_dot;
  REAL rho_prev;
  double gamma_here;
  gamma_here = 1.6666667;
  ofstream dump;
  string fout = "evol.dat";
  dump.open((fout).c_str());
  
  step = 6.03535e+06;
  // hack: commented out density is right one for popIII test
  //rho = 1.30563e-10;
  rho = 1.30563e-13;
  ei = 8.84501e+10;
  ei0 = ei;


  
  cout << "beginning abundances:" << endl;
  for(int i=0;i<NINT;i++){ 
    ysav[i] = yIn[i];
    cout << i << " " << yIn[i] << endl;
  }

  //step = 2.98783e9;
  //rho =  4.95657e-21;
  //ei = 2.39114e+13;
  gamma_here = 5.0/3.0;
  ChemDriver(step,rho,ei,gamma_here,yIn,rpar);  

  cout << "done with evol!" << endl;
  exit(0);

  // for equilibrium temperature as function of density
  /*
  ofstream dump1;
  string fout1 = "temp_eq_simp.dat";
  dump1.open((fout1).c_str());

  int ndens = 200;
  double rhol = 2.0e-27;
  double rhoh = 3.0e-18;
  duration = 1.0e17;
  

  for(int i=0;i<ndens;i++) {

    //for(int s=0;s<NINT;s++) yIn[s] = ysav[s];

    rho = rhol * pow(10.0, log10(rhoh/rhol) * double(i) / double(ndens-1)   );
    cout << " rho = " << rho  <<endl;
    
    ei = 3.0e9;
    step = 1.0;
    t = 1.0;
    while(t<duration){
      step = step * 1.1;
      ChemDriver(step,rho,ei,gamma_here,yIn,rpar);  
      t = t + step;
    }
    
    temperature = ei * (gamma_here - 1.0) /  kb * mh * 1.3 / (1.0 + 0.08 - yIn[iH2]);
    dump1 << rho << " " << temperature << " " << yIn[iCO] << " " << yIn[iC] << " " << yIn[iCP] << " " << yIn[iCHX] << " " << yIn[iHCOP] << " " << yIn[iOHX] << " " <<  yIn[iO] << " " << yIn[iH3P] << endl;
    
  }
  
  dump1.close(); 
  return 0;
  */

  //
  
  rho = 1.0e-26;
  rho_0 = rho;
  rho_prev = rho_0;
  step = 1.0;
  t = 1.0;
  while(t < duration) // Interates until dt is advanced
    {
      rho = pow(pow(rho_0,-0.5) - tff * t / 2.0 ,-2.0);
      rho_dot = tff * pow(rho,1.5);
      //rho = rho_0;
      //step = step * 1.1;
      step = 0.01 * rho / rho_dot;
      
      // adiabatically heat gas here
      ei = ei * pow(rho/rho_prev,1.66667-1);
      cout << "input ei = " << ei << endl;

      rho_prev = rho;
      
      ChemDriver(step,rho,ei,gamma_here,yIn,rpar);
      temperature = ei * (gamma - 1.0) /  kboltz * mh * 1.3 / (1.0 + 0.08 - 0.5);
      cout << "t = " << t << " rho = " << rho << ", temp = " << temperature << " h2 =" << yIn[iH2] << endl;
      dump << rho/1.67e-24 << " " << t << " " << temperature << " " << yIn[iH2] /2.0 * 1.3 << " " << yIn[iHP] << " " << yIn[iH] <<  endl;
      t = t + step;
      
      if(rho > 1.0e-10) break;

    }
  dump.close(); 
  cout << "ending mass fractions = " << endl;
  for(int i=0;i<NINT;i++) cout << i << " " << yIn[i] << endl;
  cout << " " << endl;
  cout << "ending T = " << ei * (gamma - 1.0) /  kboltz * mh * 1.2 << endl;
  cout << "ending time = " << t << endl;
  

  return 0;
		 
		     
}
