#define nrpar       30
#define g0_par      1
#define NH_par      2
#define NH2_par     3
#define NCO_par     4
#define tdust_par   5
#define zeta_par    6
#define dx_par      7
#define divv_par    8
#define divv1_par   9
#define rho_par     10
#define vturb_par   11
#define gradrho_par 12
#define gamma_par   13
#define temp_par    14
#define e0_par      15
#define temp0_par   16
#define h2var0_par  17  //ARS added rpar definition here which is passed from GetNetworkRates to Cool_function
#define J21_par  18

//ARS: dfac can be set to "den" if tracer fields are densities, or "1.0" if they are abundances.  This is used in determineGamma.cpp and AMRLevelOrion.cpp
//#define dfac        1.0  
#define dfac        den 


//ARS added a few global values as well as the species array defintions
#define MRHE2H 0.24
#define MRH    0.76
#define DeutAbund 2.6e-5
#define abhe ((1. - MRH)/4./MRH)
#define h2conv (2.0/1.3158)
#define heconv (4.0/1.3158)

#if DO_CHEMISTRY == PRIMORDIAL
  #define        iHP    0
  #define        iH     1
  #define        iHM    2
  #define        iH2P   3
  #define        iH2    4
  #define        iDP    5
  #define        iD     6
  #define        iDM    7
  #define        iHDP   8
  #define        iHD    9
  #define        iD2P   10
  #define        iD2    11
  #define        iHEP   12
  #define        iHE    13
  #define        iHEPP  14
  #define        iELEC  15
//#endif
//ARS slight modified the DO_CHEMISTRY if statement to be if/elif instead of if/then-if/then

#elif DO_CHEMISTRY == CONTEMPORARY
  #define        iH    0
  #define        iHM     1
  #define        iH2    2
  #define        iHP    3
  #define        iH2P   4
  #define        iH3P   5
  #define        iHE    6
  #define        iHEP   7
  #define        iC    8
  #define        iCP   9
  #define        iO    10
  #define       iCHX   11
  #define        iOHX    12
  #define        iHCOP   13
  #define       iCO   14
  #define        iM    15
  #define        iMP  16
  #define        iELEC  17
#endif  //DO_CHEMISTRY




