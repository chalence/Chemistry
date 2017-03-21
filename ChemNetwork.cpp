
#include <iostream>
#include <math.h>
#include "ChemGlobals.H"
#include "ChemNetwork.H"
#include "ChemDriver.H"
#include "ChemBurner.H"
using namespace ChemNS;
using namespace std;

// yt = tracked species of size NINT
// ys = total species of size NSPECIES
// GetNetwork needs to reconstruct ys from yt by calling fill_species
void GetNetwork(REAL tt,REAL* dydt, REAL* ratraw, REAL* yt, double* rpar)
{
    int i;
    
    for(i=0;i<NINT;i++) dydt[i]=0.0;
    
    REAL ys[NSPECIES];
    fill_species(yt,ys);


#if DO_CHEMISTRY == PRIMORDIAL
    
    dydt[iHP] = -ratraw[iR003]*ys[iH]*ys[iHP]+ratraw[iR004]*ys[iH]*ys[iH2P] -ratraw[iR005]*ys[iHM]*ys[iHP]-ratraw[iR007]*ys[iH2]*ys[iHP] +ratraw[iR012]*ys[iH]*ys[iELEC]-ratraw[iR013]*ys[iHP]*ys[iELEC] -ratraw[iR016]*ys[iHP]*ys[iHM]+ratraw[iR024]*ys[iH2]*ys[iHEP] +ratraw[iR026]*ys[iHEP]*ys[iH]-ratraw[iR027]*ys[iHE]*ys[iHP] -ratraw[iR034]*ys[iD]*ys[iHP]+ratraw[iR035]*ys[iH]*ys[iDP] +ratraw[iR038]*ys[iHDP]*ys[iH]+ratraw[iR039]*ys[iH2]*ys[iDP] -ratraw[iR041]*ys[iHD]*ys[iHP]-ratraw[iR042]*ys[iD]*ys[iHP] -ratraw[iR060]*ys[iHP]*ys[iDM]-ratraw[iR067]*ys[iHP]*ys[iDM] +ratraw[iR082]*ys[iH2P]*ys[iD]+ratraw[iR085]*ys[iHDP]*ys[iD] +ratraw[iR087]*ys[iH]*ys[iD2P]-ratraw[iR092]*ys[iHD]*ys[iHP] -ratraw[iR093]*ys[iHD]*ys[iHP]+ratraw[iR095]*ys[iHD]*ys[iDP] -ratraw[iR097]*ys[iD2]*ys[iHP]-ratraw[iR098]*ys[iD2]*ys[iHP] -ratraw[iR099]*ys[iD2]*ys[iHP]+ratraw[iR102]*ys[iHD]*ys[iHEP] +ratraw[iR118]*ys[iH2P]+ratraw[iR120]*ys[iHDP];
    
    
    dydt[iH] =  -ratraw[iR001]*ys[iH]*ys[iELEC]-ratraw[iR002]*ys[iHM]*ys[iH]-ratraw[iR003]*ys[iH]*ys[iHP]-ratraw[iR004]*ys[iH]*ys[iH2P]+ratraw[iR005]*ys[iHM]*ys[iHP]+ratraw[iR005]*ys[iHM]*ys[iHP]+ratraw[iR006]*ys[iH2P]*ys[iELEC]+ratraw[iR006]*ys[iH2P]*ys[iELEC]+ratraw[iR007]*ys[iH2]*ys[iHP]+ratraw[iR008]*ys[iH2]*ys[iELEC]+ratraw[iR008]*ys[iH2]*ys[iELEC]-ratraw[iR009]*ys[iH2]*ys[iH]+ratraw[iR009]*ys[iH2]*ys[iH]+ratraw[iR009]*ys[iH2]*ys[iH]+ratraw[iR009]*ys[iH2]*ys[iH]+ratraw[iR010]*ys[iH2]*ys[iH2] +ratraw[iR010]*ys[iH2]*ys[iH2]+ratraw[iR011]*ys[iH2]*ys[iHE]    +ratraw[iR011]*ys[iH2]*ys[iHE]-ratraw[iR012]*ys[iH]*ys[iELEC]    +ratraw[iR013]*ys[iHP]*ys[iELEC]+ratraw[iR014]*ys[iHM]*ys[iELEC]    -ratraw[iR015]*ys[iHM]*ys[iH]+ratraw[iR015]*ys[iHM]*ys[iH]    +ratraw[iR015]*ys[iHM]*ys[iH]+ratraw[iR021]*ys[iHM]*ys[iH2P]    +ratraw[iR022]*ys[iHM]*ys[iH2P]+ratraw[iR022]*ys[iHM]*ys[iH2P]    +ratraw[iR022]*ys[iHM]*ys[iH2P]+ratraw[iR023]*ys[iH2]*ys[iELEC]    +ratraw[iR024]*ys[iH2]*ys[iHEP]-ratraw[iR026]*ys[iHEP]*ys[iH]    +ratraw[iR027]*ys[iHE]*ys[iHP]+ratraw[iR028]*ys[iHEP]*ys[iHM]    +ratraw[iR029]*ys[iHE]*ys[iHM]-ratraw[iR030]*ys[iH]*ys[iH]*ys[iH]    -ratraw[iR030]*ys[iH]*ys[iH]*ys[iH]-ratraw[iR030]*ys[iH]*ys[iH]*ys[iH]    +ratraw[iR030]*ys[iH]*ys[iH]*ys[iH]-ratraw[iR031]*ys[iH]*ys[iH]*ys[iH2]    -ratraw[iR031]*ys[iH]*ys[iH]*ys[iH2]-ratraw[iR032]*ys[iH]*ys[iH]*ys[iHE]    -ratraw[iR032]*ys[iH]*ys[iH]*ys[iHE]+ratraw[iR034]*ys[iD]*ys[iHP]    -ratraw[iR035]*ys[iH]*ys[iDP]-ratraw[iR036]*ys[iH]*ys[iD]    +ratraw[iR037]*ys[iH2]*ys[iD]-ratraw[iR038]*ys[iHDP]*ys[iH]    -ratraw[iR040]*ys[iHD]*ys[iH]-ratraw[iR043]*ys[iH]*ys[iDP]    +ratraw[iR044]*ys[iHDP]*ys[iELEC]+ratraw[iR048]*ys[iH2P]*ys[iD]    -ratraw[iR050]*ys[iHDP]*ys[iH]-ratraw[iR052]*ys[iH]*ys[iDM]    +ratraw[iR053]*ys[iD]*ys[iHM]-ratraw[iR055]*ys[iH]*ys[iDM]    +ratraw[iR057]*ys[iHD]*ys[iELEC]-ratraw[iR064]*ys[iDM]*ys[iH]    +ratraw[iR064]*ys[iDM]*ys[iH]+ratraw[iR066]*ys[iDP]*ys[iHM]    +ratraw[iR067]*ys[iHP]*ys[iDM]+ratraw[iR070]*ys[iH2P]*ys[iDM]    +ratraw[iR070]*ys[iH2P]*ys[iDM]+ratraw[iR071]*ys[iHDP]*ys[iHM]    +ratraw[iR072]*ys[iHDP]*ys[iHM]+ratraw[iR072]*ys[iHDP]*ys[iHM]    +ratraw[iR074]*ys[iHDP]*ys[iDM]+ratraw[iR075]*ys[iD2P]*ys[iHM]    +ratraw[iR076]*ys[iD2P]*ys[iHM]-ratraw[iR083]*ys[iHDP]*ys[iH]    +ratraw[iR084]*ys[iHDP]*ys[iD]-ratraw[iR087]*ys[iH]*ys[iD2P]    -ratraw[iR088]*ys[iD2P]*ys[iH]-ratraw[iR089]*ys[iD2P]*ys[iH]    +ratraw[iR091]*ys[iH2]*ys[iDP]+ratraw[iR092]*ys[iHD]*ys[iHP]    +ratraw[iR096]*ys[iHD]*ys[iDP]+ratraw[iR099]*ys[iD2]*ys[iHP]    +ratraw[iR103]*ys[iHD]*ys[iHEP]+ratraw[iR106]*ys[iHD]*ys[iD]    -ratraw[iR107]*ys[iD2]*ys[iH]-ratraw[iR108]*ys[iHD]*ys[iH]    +ratraw[iR108]*ys[iHD]*ys[iH]+ratraw[iR108]*ys[iHD]*ys[iH]    +ratraw[iR109]*ys[iHD]*ys[iH2]+ratraw[iR110]*ys[iHD]*ys[iHE]    +ratraw[iR111]*ys[iHD]*ys[iELEC]-ratraw[iR112]*ys[iD2]*ys[iH]    +ratraw[iR112]*ys[iD2]*ys[iH]+ratraw[iR116]*ys[iHM]    +ratraw[iR118]*ys[iH2P]+ratraw[iR119]*ys[iHDP]    +ratraw[iR122]*ys[iH2]+ratraw[iR122]*ys[iH2]+ratraw[iR123]*ys[iHD];
    
    dydt[iHM] =  ratraw[iR001]*ys[iH]*ys[iELEC]-ratraw[iR002]*ys[iHM]*ys[iH]    -ratraw[iR005]*ys[iHM]*ys[iHP]-ratraw[iR014]*ys[iHM]*ys[iELEC]    -ratraw[iR015]*ys[iHM]*ys[iH]-ratraw[iR016]*ys[iHP]*ys[iHM]    -ratraw[iR021]*ys[iHM]*ys[iH2P]-ratraw[iR022]*ys[iHM]*ys[iH2P]    +ratraw[iR023]*ys[iH2]*ys[iELEC]-ratraw[iR028]*ys[iHEP]*ys[iHM]    -ratraw[iR029]*ys[iHE]*ys[iHM]+ratraw[iR052]*ys[iH]*ys[iDM]    -ratraw[iR053]*ys[iD]*ys[iHM]-ratraw[iR054]*ys[iD]*ys[iHM]    +ratraw[iR058]*ys[iHD]*ys[iELEC]-ratraw[iR061]*ys[iDP]*ys[iHM]    -ratraw[iR066]*ys[iDP]*ys[iHM]-ratraw[iR071]*ys[iHDP]*ys[iHM]    -ratraw[iR072]*ys[iHDP]*ys[iHM]-ratraw[iR075]*ys[iD2P]*ys[iHM]    -ratraw[iR076]*ys[iD2P]*ys[iHM]-ratraw[iR116]*ys[iHM];
    
    dydt[iH2P] =  ratraw[iR003]*ys[iH]*ys[iHP]-ratraw[iR004]*ys[iH]*ys[iH2P]    -ratraw[iR006]*ys[iH2P]*ys[iELEC]+ratraw[iR007]*ys[iH2]*ys[iHP]    +ratraw[iR016]*ys[iHP]*ys[iHM]-ratraw[iR021]*ys[iHM]*ys[iH2P]    -ratraw[iR022]*ys[iHM]*ys[iH2P]+ratraw[iR025]*ys[iH2]*ys[iHEP]    -ratraw[iR048]*ys[iH2P]*ys[iD]+ratraw[iR050]*ys[iHDP]*ys[iH]    -ratraw[iR069]*ys[iH2P]*ys[iDM]-ratraw[iR070]*ys[iH2P]*ys[iDM]    -ratraw[iR081]*ys[iD]*ys[iH2P]-ratraw[iR082]*ys[iH2P]*ys[iD]    +ratraw[iR090]*ys[iH2]*ys[iDP]+ratraw[iR093]*ys[iHD]*ys[iHP]    -ratraw[iR118]*ys[iH2P];
    
    dydt[iH2] =  ratraw[iR002]*ys[iHM]*ys[iH]+ratraw[iR004]*ys[iH]*ys[iH2P]    -ratraw[iR007]*ys[iH2]*ys[iHP]-ratraw[iR008]*ys[iH2]*ys[iELEC]    -ratraw[iR009]*ys[iH2]*ys[iH]-ratraw[iR010]*ys[iH2]*ys[iH2]    -ratraw[iR010]*ys[iH2]*ys[iH2]+ratraw[iR010]*ys[iH2]*ys[iH2]    -ratraw[iR011]*ys[iH2]*ys[iHE]+ratraw[iR021]*ys[iHM]*ys[iH2P]    -ratraw[iR023]*ys[iH2]*ys[iELEC]-ratraw[iR024]*ys[iH2]*ys[iHEP]    -ratraw[iR025]*ys[iH2]*ys[iHEP]+ratraw[iR030]*ys[iH]*ys[iH]*ys[iH]    -ratraw[iR031]*ys[iH]*ys[iH]*ys[iH2]+ratraw[iR031]*ys[iH]*ys[iH]*ys[iH2]    +ratraw[iR031]*ys[iH]*ys[iH]*ys[iH2]+ratraw[iR032]*ys[iH]*ys[iH]*ys[iHE]    -ratraw[iR037]*ys[iH2]*ys[iD]-ratraw[iR039]*ys[iH2]*ys[iDP]    +ratraw[iR040]*ys[iHD]*ys[iH]+ratraw[iR041]*ys[iHD]*ys[iHP]    +ratraw[iR069]*ys[iH2P]*ys[iDM]+ratraw[iR081]*ys[iD]*ys[iH2P]    +ratraw[iR083]*ys[iHDP]*ys[iH]-ratraw[iR090]*ys[iH2]*ys[iDP]    -ratraw[iR091]*ys[iH2]*ys[iDP]-ratraw[iR109]*ys[iHD]*ys[iH2]    +ratraw[iR109]*ys[iHD]*ys[iH2]-ratraw[iR113]*ys[iD2]*ys[iH2]    +ratraw[iR113]*ys[iD2]*ys[iH2]-ratraw[iR122]*ys[iH2];
    
    dydt[iDP] =  -ratraw[iR033]*ys[iDP]*ys[iELEC]+ratraw[iR034]*ys[iD]*ys[iHP]    -ratraw[iR035]*ys[iH]*ys[iDP]-ratraw[iR039]*ys[iH2]*ys[iDP]    +ratraw[iR041]*ys[iHD]*ys[iHP]-ratraw[iR043]*ys[iH]*ys[iDP]    +ratraw[iR045]*ys[iD]*ys[iELEC]+ratraw[iR046]*ys[iHEP]*ys[iD]    -ratraw[iR047]*ys[iHE]*ys[iDP]+ratraw[iR049]*ys[iHDP]*ys[iD]    -ratraw[iR061]*ys[iDP]*ys[iHM]-ratraw[iR062]*ys[iDP]*ys[iDM]    -ratraw[iR066]*ys[iDP]*ys[iHM]-ratraw[iR068]*ys[iDP]*ys[iDM]    -ratraw[iR080]*ys[iD]*ys[iDP]+ratraw[iR081]*ys[iD]*ys[iH2P]    +ratraw[iR083]*ys[iHDP]*ys[iH]+ratraw[iR086]*ys[iD]*ys[iD2P]    +ratraw[iR089]*ys[iD2P]*ys[iH]-ratraw[iR090]*ys[iH2]*ys[iDP]    -ratraw[iR091]*ys[iH2]*ys[iDP]-ratraw[iR094]*ys[iHD]*ys[iDP]    -ratraw[iR095]*ys[iHD]*ys[iDP]-ratraw[iR096]*ys[iHD]*ys[iDP]    +ratraw[iR097]*ys[iD2]*ys[iHP]-ratraw[iR100]*ys[iD2]*ys[iDP]    +ratraw[iR103]*ys[iHD]*ys[iHEP]+ratraw[iR105]*ys[iD2]*ys[iHEP]    +ratraw[iR119]*ys[iHDP]+ratraw[iR121]*ys[iD2P];
    
    dydt[iD] =  ratraw[iR033]*ys[iDP]*ys[iELEC]-ratraw[iR034]*ys[iD]*ys[iHP]    +ratraw[iR035]*ys[iH]*ys[iDP]-ratraw[iR036]*ys[iH]*ys[iD]    -ratraw[iR037]*ys[iH2]*ys[iD]+ratraw[iR040]*ys[iHD]*ys[iH]    -ratraw[iR042]*ys[iD]*ys[iHP]+ratraw[iR044]*ys[iHDP]*ys[iELEC]    -ratraw[iR045]*ys[iD]*ys[iELEC]-ratraw[iR046]*ys[iHEP]*ys[iD]    +ratraw[iR047]*ys[iHE]*ys[iDP]-ratraw[iR048]*ys[iH2P]*ys[iD]    -ratraw[iR049]*ys[iHDP]*ys[iD]+ratraw[iR050]*ys[iHDP]*ys[iH]    -ratraw[iR051]*ys[iD]*ys[iELEC]+ratraw[iR052]*ys[iH]*ys[iDM]    -ratraw[iR053]*ys[iD]*ys[iHM]-ratraw[iR054]*ys[iD]*ys[iHM]    -ratraw[iR056]*ys[iD]*ys[iDM]+ratraw[iR058]*ys[iHD]*ys[iELEC]    +ratraw[iR059]*ys[iD2]*ys[iELEC]+ratraw[iR063]*ys[iDM]*ys[iELEC]    +ratraw[iR064]*ys[iDM]*ys[iH]+ratraw[iR065]*ys[iDM]*ys[iHE]    +ratraw[iR066]*ys[iDP]*ys[iHM]+ratraw[iR067]*ys[iHP]*ys[iDM]+ratraw[iR068]*ys[iDP]*ys[iDM]+ratraw[iR068]*ys[iDP]*ys[iDM]+ratraw[iR069]*ys[iH2P]*ys[iDM]+ratraw[iR070]*ys[iH2P]*ys[iDM]    +ratraw[iR072]*ys[iHDP]*ys[iHM]+ratraw[iR073]*ys[iHDP]*ys[iDM]+ratraw[iR074]*ys[iHDP]*ys[iDM]+ratraw[iR074]*ys[iHDP]*ys[iDM]+ratraw[iR076]*ys[iD2P]*ys[iHM]+ratraw[iR076]*ys[iD2P]*ys[iHM]+ratraw[iR077]*ys[iD2P]*ys[iDM]+ratraw[iR078]*ys[iD2P]*ys[iDM]    +ratraw[iR078]*ys[iD2P]*ys[iDM]+ratraw[iR078]*ys[iD2P]*ys[iDM]    +ratraw[iR079]*ys[iHEP]*ys[iDM]-ratraw[iR080]*ys[iD]*ys[iDP]    -ratraw[iR081]*ys[iD]*ys[iH2P]-ratraw[iR082]*ys[iH2P]*ys[iD]    -ratraw[iR084]*ys[iHDP]*ys[iD]-ratraw[iR085]*ys[iHDP]*ys[iD]    -ratraw[iR086]*ys[iD]*ys[iD2P]+ratraw[iR088]*ys[iD2P]*ys[iH]    +ratraw[iR090]*ys[iH2]*ys[iDP]+ratraw[iR093]*ys[iHD]*ys[iHP]    +ratraw[iR094]*ys[iHD]*ys[iDP]+ratraw[iR098]*ys[iD2]*ys[iHP]    +ratraw[iR100]*ys[iD2]*ys[iDP]+ratraw[iR102]*ys[iHD]*ys[iHEP]    +ratraw[iR105]*ys[iD2]*ys[iHEP]-ratraw[iR106]*ys[iHD]*ys[iD]    +ratraw[iR107]*ys[iD2]*ys[iH]+ratraw[iR108]*ys[iHD]*ys[iH]    +ratraw[iR109]*ys[iHD]*ys[iH2]+ratraw[iR110]*ys[iHD]*ys[iHE]    +ratraw[iR111]*ys[iHD]*ys[iELEC]+ratraw[iR112]*ys[iD2]*ys[iH]    +ratraw[iR112]*ys[iD2]*ys[iH]+ratraw[iR113]*ys[iD2]*ys[iH2]    +ratraw[iR113]*ys[iD2]*ys[iH2]+ratraw[iR114]*ys[iD2]*ys[iHE]    +ratraw[iR114]*ys[iD2]*ys[iHE]+ratraw[iR115]*ys[iD2]*ys[iELEC]    +ratraw[iR115]*ys[iD2]*ys[iELEC]+ratraw[iR117]*ys[iDM]    +ratraw[iR120]*ys[iHDP]+ratraw[iR121]*ys[iD2P]+ratraw[iR123]*ys[iHD]    +ratraw[iR124]*ys[iD2]+ratraw[iR124]*ys[iD2];
    
    dydt[iDM] =
    ratraw[iR051]*ys[iD]*ys[iELEC]-ratraw[iR052]*ys[iH]*ys[iDM]
    +ratraw[iR053]*ys[iD]*ys[iHM]-ratraw[iR055]*ys[iH]*ys[iDM]
    -ratraw[iR056]*ys[iD]*ys[iDM]+ratraw[iR057]*ys[iHD]*ys[iELEC]
    +ratraw[iR059]*ys[iD2]*ys[iELEC]-ratraw[iR060]*ys[iHP]*ys[iDM]
    -ratraw[iR062]*ys[iDP]*ys[iDM]-ratraw[iR063]*ys[iDM]*ys[iELEC]
    -ratraw[iR064]*ys[iDM]*ys[iH]-ratraw[iR065]*ys[iDM]*ys[iHE]
    -ratraw[iR067]*ys[iHP]*ys[iDM]-ratraw[iR068]*ys[iDP]*ys[iDM]
    -ratraw[iR069]*ys[iH2P]*ys[iDM]-ratraw[iR070]*ys[iH2P]*ys[iDM]
    -ratraw[iR073]*ys[iHDP]*ys[iDM]-ratraw[iR074]*ys[iHDP]*ys[iDM]
    -ratraw[iR077]*ys[iD2P]*ys[iDM]-ratraw[iR078]*ys[iD2P]*ys[iDM]
    -ratraw[iR079]*ys[iHEP]*ys[iDM]-ratraw[iR117]*ys[iDM];
    
    dydt[iHDP] =  -ratraw[iR038]*ys[iHDP]*ys[iH]+ratraw[iR042]*ys[iD]*ys[iHP]    +ratraw[iR043]*ys[iH]*ys[iDP]-ratraw[iR044]*ys[iHDP]*ys[iELEC]    +ratraw[iR048]*ys[iH2P]*ys[iD]-ratraw[iR049]*ys[iHDP]*ys[iD]    -ratraw[iR050]*ys[iHDP]*ys[iH]+ratraw[iR060]*ys[iHP]*ys[iDM]    +ratraw[iR061]*ys[iDP]*ys[iHM]-ratraw[iR071]*ys[iHDP]*ys[iHM]    -ratraw[iR072]*ys[iHDP]*ys[iHM]-ratraw[iR073]*ys[iHDP]*ys[iDM]    -ratraw[iR074]*ys[iHDP]*ys[iDM]-ratraw[iR083]*ys[iHDP]*ys[iH]    -ratraw[iR084]*ys[iHDP]*ys[iD]-ratraw[iR085]*ys[iHDP]*ys[iD]    +ratraw[iR088]*ys[iD2P]*ys[iH]+ratraw[iR091]*ys[iH2]*ys[iDP]    +ratraw[iR092]*ys[iHD]*ys[iHP]+ratraw[iR094]*ys[iHD]*ys[iDP]    +ratraw[iR098]*ys[iD2]*ys[iHP]+ratraw[iR101]*ys[iHD]*ys[iHEP]    -ratraw[iR119]*ys[iHDP]-ratraw[iR120]*ys[iHDP];
    
    dydt[iHD] =  ratraw[iR036]*ys[iH]*ys[iD]+ratraw[iR037]*ys[iH2]*ys[iD]    +ratraw[iR038]*ys[iHDP]*ys[iH]+ratraw[iR039]*ys[iH2]*ys[iDP]    -ratraw[iR040]*ys[iHD]*ys[iH]-ratraw[iR041]*ys[iHD]*ys[iHP]    +ratraw[iR049]*ys[iHDP]*ys[iD]+ratraw[iR054]*ys[iD]*ys[iHM]    +ratraw[iR055]*ys[iH]*ys[iDM]-ratraw[iR057]*ys[iHD]*ys[iELEC]    -ratraw[iR058]*ys[iHD]*ys[iELEC]+ratraw[iR071]*ys[iHDP]*ys[iHM]    +ratraw[iR073]*ys[iHDP]*ys[iDM]+ratraw[iR082]*ys[iH2P]*ys[iD]    +ratraw[iR089]*ys[iD2P]*ys[iH]-ratraw[iR092]*ys[iHD]*ys[iHP]    -ratraw[iR093]*ys[iHD]*ys[iHP]-ratraw[iR094]*ys[iHD]*ys[iDP]    -ratraw[iR095]*ys[iHD]*ys[iDP]-ratraw[iR096]*ys[iHD]*ys[iDP]    +ratraw[iR097]*ys[iD2]*ys[iHP]-ratraw[iR101]*ys[iHD]*ys[iHEP]    -ratraw[iR102]*ys[iHD]*ys[iHEP]-ratraw[iR103]*ys[iHD]*ys[iHEP]    -ratraw[iR106]*ys[iHD]*ys[iD]+ratraw[iR107]*ys[iD2]*ys[iH]    -ratraw[iR108]*ys[iHD]*ys[iH]-ratraw[iR109]*ys[iHD]*ys[iH2]    -ratraw[iR110]*ys[iHD]*ys[iHE]-ratraw[iR111]*ys[iHD]*ys[iELEC]-ratraw[iR123]*ys[iHD];
    
    dydt[iD2P] =  ratraw[iR062]*ys[iDP]*ys[iDM]-ratraw[iR075]*ys[iD2P]*ys[iHM]    -ratraw[iR076]*ys[iD2P]*ys[iHM]-ratraw[iR077]*ys[iD2P]*ys[iDM]    -ratraw[iR078]*ys[iD2P]*ys[iDM]+ratraw[iR080]*ys[iD]*ys[iDP]    +ratraw[iR084]*ys[iHDP]*ys[iD]-ratraw[iR086]*ys[iD]*ys[iD2P]    -ratraw[iR087]*ys[iH]*ys[iD2P]-ratraw[iR088]*ys[iD2P]*ys[iH]    -ratraw[iR089]*ys[iD2P]*ys[iH]+ratraw[iR096]*ys[iHD]*ys[iDP]    +ratraw[iR099]*ys[iD2]*ys[iHP]+ratraw[iR100]*ys[iD2]*ys[iDP]    +ratraw[iR104]*ys[iD2]*ys[iHEP]-ratraw[iR121]*ys[iD2P];
    
    dydt[iD2] =  ratraw[iR056]*ys[iD]*ys[iDM]-ratraw[iR059]*ys[iD2]*ys[iELEC]+ratraw[iR075]*ys[iD2P]*ys[iHM]+ratraw[iR077]*ys[iD2P]*ys[iDM]    +ratraw[iR085]*ys[iHDP]*ys[iD]+ratraw[iR086]*ys[iD]*ys[iD2P]+ratraw[iR087]*ys[iH]*ys[iD2P]+ratraw[iR095]*ys[iHD]*ys[iDP]    -ratraw[iR097]*ys[iD2]*ys[iHP]-ratraw[iR098]*ys[iD2]*ys[iHP]-ratraw[iR099]*ys[iD2]*ys[iHP]-ratraw[iR100]*ys[iD2]*ys[iDP]    -ratraw[iR104]*ys[iD2]*ys[iHEP]-ratraw[iR105]*ys[iD2]*ys[iHEP]+ratraw[iR106]*ys[iHD]*ys[iD]-ratraw[iR107]*ys[iD2]*ys[iH]    -ratraw[iR112]*ys[iD2]*ys[iH]-ratraw[iR113]*ys[iD2]*ys[iH2]-ratraw[iR114]*ys[iD2]*ys[iHE]-ratraw[iR115]*ys[iD2]*ys[iELEC]    -ratraw[iR124]*ys[iD2];
    
    dydt[iHEP] =   ratraw[iR017]*ys[iHE]*ys[iELEC]-ratraw[iR018]*ys[iHEP]*ys[iELEC]-ratraw[iR019]*ys[iHEP]*ys[iELEC]+ratraw[iR020]*ys[iHEPP]*ys[iELEC]-ratraw[iR024]*ys[iH2]*ys[iHEP]-ratraw[iR025]*ys[iH2]*ys[iHEP]-ratraw[iR026]*ys[iHEP]*ys[iH]+ratraw[iR027]*ys[iHE]*ys[iHP]-ratraw[iR028]*ys[iHEP]*ys[iHM]-ratraw[iR046]*ys[iHEP]*ys[iD]+ratraw[iR047]*ys[iHE]*ys[iDP]-ratraw[iR079]*ys[iHEP]*ys[iDM]-ratraw[iR101]*ys[iHD]*ys[iHEP]-ratraw[iR102]*ys[iHD]*ys[iHEP]-ratraw[iR103]*ys[iHD]*ys[iHEP]-ratraw[iR104]*ys[iD2]*ys[iHEP]-ratraw[iR105]*ys[iD2]*ys[iHEP];
                                                    
    dydt[iHE] =  -ratraw[iR011]*ys[iH2]*ys[iHE]+ratraw[iR011]*ys[iH2]*ys[iHE]-ratraw[iR017]*ys[iHE]*ys[iELEC]+ratraw[iR019]*ys[iHEP]*ys[iELEC]+ratraw[iR024]*ys[iH2]*ys[iHEP]+ratraw[iR025]*ys[iH2]*ys[iHEP]+ratraw[iR026]*ys[iHEP]*ys[iH]-ratraw[iR027]*ys[iHE]*ys[iHP]+ratraw[iR028]*ys[iHEP]*ys[iHM]-ratraw[iR029]*ys[iHE]*ys[iHM]+ratraw[iR029]*ys[iHE]*ys[iHM]-ratraw[iR032]*ys[iH]*ys[iH]*ys[iHE]+ratraw[iR032]*ys[iH]*ys[iH]*ys[iHE]+ratraw[iR046]*ys[iHEP]*ys[iD]-ratraw[iR047]*ys[iHE]*ys[iDP]-ratraw[iR065]*ys[iDM]*ys[iHE]        +ratraw[iR065]*ys[iDM]*ys[iHE]+ratraw[iR079]*ys[iHEP]*ys[iDM]        +ratraw[iR101]*ys[iHD]*ys[iHEP]+ratraw[iR102]*ys[iHD]*ys[iHEP]        +ratraw[iR103]*ys[iHD]*ys[iHEP]+ratraw[iR104]*ys[iD2]*ys[iHEP]        +ratraw[iR105]*ys[iD2]*ys[iHEP]-ratraw[iR110]*ys[iHD]*ys[iHE]        +ratraw[iR110]*ys[iHD]*ys[iHE]-ratraw[iR114]*ys[iD2]*ys[iHE]        +ratraw[iR114]*ys[iD2]*ys[iHE];
        
    dydt[iHEPP] =  ratraw[iR018]*ys[iHEP]*ys[iELEC]-ratraw[iR020]*ys[iHEPP]*ys[iELEC];
                                                                      

#endif

#if DO_CHEMISTRY == CONTEMPORARY
    
    dydt[iH] = -ys[iH] * ys[iELEC] * ratraw[iR001]
      - ys[iHM] * ys[iH] * ratraw[iR002]
      - ys[iH] * ys[iHP] * ratraw[iR003]
      - ys[iH] * ys[iH2P] * ratraw[iR004]
      + 2.0 * ys[iHM] * ys[iHP] * ratraw[iR005]
      + 2.0 * ys[iH2P] * ys[iELEC] * ratraw[iR006]
      + ys[iH2] * ys[iHP] * ratraw[iR007]
      + 2.0 * ys[iH2] * ys[iELEC] * ratraw[iR008]
      + 2.0 * ys[iH2] * ys[iH] * ratraw[iR009]
      + 2.0 * ys[iH2] * ys[iH2] * ratraw[iR010]
      - ys[iH] * ys[iELEC] * ratraw[iR011]
      + ys[iHP] * ys[iELEC] * ratraw[iR012]
      + ys[iHM] * ys[iELEC] * ratraw[iR013]
      + ys[iHM] * ys[iH] * ratraw[iR014]
      - ys[iHEP] * ys[iH] * ratraw[iR018]
      + ys[iHE] * ys[iHP] * ratraw[iR019]
      + 2.0 * ys[iH2] * ys[iHE] * ratraw[iR020]
      + ys[iH2] * ys[iHEP] * ratraw[iR022]
      + ys[iHEP] * ys[iHM] * ratraw[iR023]
      + ys[iH3P] * ys[iELEC] * ratraw[iR024]
      + 3.0 * ys[iH3P] * ys[iELEC] * ratraw[iR025]
      - 2.0 * ys[iH] * ys[iH] * ys[iH] * ratraw[iR027]
      - 2.0 * ys[iH] * ys[iH] * ys[iH2] * ratraw[iR028] 
      - 2.0 * ys[iH] * ys[iH] * ys[iHE] * ratraw[iR029] 
      - 2.0 * ys[iH] * ratraw[iR030]
      + ys[iCP] * ys[iH2] * ratraw[iR035]
      + ys[iO] * ys[iCHX] * ratraw[iR037]
      + ys[iC] * ys[iOHX] * ratraw[iR038]
      + ys[iHCOP] * ys[iELEC] * ratraw[iR040]
      + ys[iH2P] * ys[iH2] * ratraw[iR043]
      - ys[iH3P] * ys[iH] * ratraw[iR044]
      + ys[iOHX] * ys[iO] * ratraw[iR045]
      + ys[iCHX] * ratraw[iP001]
      + ys[iOHX] * ratraw[iP002]
      + ys[iHCOP] * ratraw[iP004]
      + ys[iHM] * ratraw[iP005]
      + ys[iH2P] * ratraw[iP006]
      + 2.0 * ys[iH2] * ratraw[iP007]
      + ys[iH3P] * ratraw[iP009]
      - ys[iH] * ratraw[iP012]
      + 2.0 * ys[iH2] * ratraw[iP016] ;
    
    dydt[iHM] = ys[iH] * ys[iELEC] * ratraw[iR001]
      - ys[iHM] * ys[iH] * ratraw[iR002]
      - ys[iHM] * ys[iHP] * ratraw[iR005]
      - ys[iHM] * ys[iELEC] * ratraw[iR013] 
      - ys[iHM] * ys[iH] * ratraw[iR014]
      - ys[iHM] * ys[iHP] * ratraw[iR015]
      - ys[iHEP] * ys[iHM] * ratraw[iR023]
      - ys[iHM] * ratraw[iP005] ;

    dydt[iH2] = ys[iHM] * ys[iH] * ratraw[iR002]
      + ys[iH] * ys[iH2P] * ratraw[iR004] 
      - ys[iH2] * ys[iHP] * ratraw[iR007]
      - ys[iH2] * ys[iELEC] * ratraw[iR008] 
      - ys[iH2] * ys[iH] * ratraw[iR009]
      - ys[iH2] * ys[iH2] * ratraw[iR010]
      - ys[iH2] * ys[iHE] * ratraw[iR020]
      - ys[iH2] * ys[iHEP] * ratraw[iR021]
      - ys[iH2] * ys[iHEP] * ratraw[iR022]
      + ys[iH3P] * ys[iELEC] * ratraw[iR024]
      - ys[iH2] * ys[iHP] * ratraw[iR026]
      + ys[iH] * ys[iH] * ys[iH] * ratraw[iR027] 
      + ys[iH] * ys[iH] * ys[iH2] * ratraw[iR028]
      + ys[iH] * ys[iH] * ys[iHE] * ratraw[iR029]
      + ys[iH] * ratraw[iR030] // dust grain formation
      + ys[iH3P] * ys[iC] * ratraw[iR031] 
      + ys[iH3P] * ys[iO] * ratraw[iR032]
      + ys[iH3P] * ys[iCO] * ratraw[iR033] 
      - ys[iCP] * ys[iH2] * ratraw[iR035] 
      + ys[iH3P] * ys[iM] * ratraw[iR042]
      - ys[iH2P] * ys[iH2] * ratraw[iR043]
      + ys[iH3P] * ys[iH] * ratraw[iR044] 
      - ys[iH2] * ratraw[iP007]
      + ys[iH3P] * ratraw[iP008]
      - ys[iH2] * ratraw[iP015]
      - ys[iH2] * ratraw[iP016] ;
    
    dydt[iHP] = - ys[iH] * ys[iHP] * ratraw[iR003]
      + ys[iH] * ys[iH2P] * ratraw[iR004] 
      - ys[iHM] * ys[iHP] * ratraw[iR005]
      - ys[iH2] * ys[iHP] * ratraw[iR007]
      + ys[iH] * ys[iELEC] * ratraw[iR011] 
      - ys[iHP] * ys[iELEC] * ratraw[iR012]
      - ys[iHM] * ys[iHP] * ratraw[iR015]
      + ys[iHEP] * ys[iH] * ratraw[iR018]
      - ys[iHE] * ys[iHP] * ratraw[iR019]
      + ys[iH2] * ys[iHEP] * ratraw[iR022]
      - ys[iH2] * ys[iHP] * ratraw[iR026]
      + ys[iH2P] * ratraw[iP006]
      + ys[iH3P] * ratraw[iP008]
      + ys[iH] * ratraw[iP012];
    
    dydt[iH2P] = ys[iH] * ys[iHP] * ratraw[iR003]
      - ys[iH] * ys[iH2P] * ratraw[iR004]
      - ys[iH2P] * ys[iELEC] * ratraw[iR006]
      + ys[iH2] * ys[iHP] * ratraw[iR007]
      + ys[iHM] * ys[iHP] * ratraw[iR015]
      +ys[iH2] * ys[iHEP] * ratraw[iR021]
      - ys[iH2P] * ys[iH2] * ratraw[iR043]
      + ys[iH3P] * ys[iH] * ratraw[iR044]
      - ys[iH2P] * ratraw[iP006]
      + ys[iH3P] * ratraw[iP009]
      + ys[iH2] * ratraw[iP015] ;
    
    dydt[iH3P] = - ys[iH3P] * ys[iELEC] * ratraw[iR024]
      - ys[iH3P] * ys[iELEC] * ratraw[iR025]
      + ys[iH2] * ys[iHP] * ratraw[iR026] 
      - ys[iH3P] * ys[iC] * ratraw[iR031]
      - ys[iH3P] * ys[iO] * ratraw[iR032]
      - ys[iH3P] * ys[iCO] * ratraw[iR033]
      - ys[iH3P] * ys[iM] * ratraw[iR042]
      - ys[iH3P] * ratraw[iP008]
      - ys[iH3P] * ratraw[iP009]
      + ys[iH2P] * ys[iH2] * ratraw[iR043]
      - ys[iH3P] * ys[iH] * ratraw[iR044] ;
    
    dydt[iHE] = -ys[iHE] * ys[iELEC] * ratraw[iR016]
      + ys[iHEP] * ys[iELEC] * ratraw[iR017]
      + ys[iHEP] * ys[iH] * ratraw[iR018] 
      - ys[iHE] * ys[iHP] * ratraw[iR019]
      + ys[iH2] * ys[iHEP] * ratraw[iR021]
      +ys[iH2] * ys[iHEP] * ratraw[iR022]
      + ys[iHEP] * ys[iHM] * ratraw[iR023]
      + ys[iHEP] * ys[iCO] * ratraw[iR034]
      - ys[iHE] * ratraw[iP013] ;
    
    dydt[iHEP] = ys[iHE] * ys[iELEC] * ratraw[iR016] 
      - ys[iHEP] * ys[iELEC] * ratraw[iR017]
      - ys[iHEP] * ys[iH] * ratraw[iR018]
      + ys[iHE] * ys[iHP] * ratraw[iR019] 
      - ys[iH2] * ys[iHEP] * ratraw[iR021]
      - ys[iH2] * ys[iHEP] * ratraw[iR022] 
      - ys[iHEP] * ys[iHM] * ratraw[iR023] 
      - ys[iHEP] * ys[iCO] * ratraw[iR034]
      + ys[iHE] * ratraw[iP013] ;
    
    dydt[iC] = - ys[iH3P] * ys[iC] * ratraw[iR031]
      - ys[iC] * ys[iOHX] * ratraw[iR038]
      + ys[iCP] * ys[iELEC] * ratraw[iR039]
      + ys[iCHX] * ratraw[iP001]
      - ys[iC] * ratraw[iP010]
      + ys[iCO] * ratraw[iP011]
      - ys[iC] * ratraw[iP014] ;
    
    dydt[iCP] =  ys[iHEP] * ys[iCO] * ratraw[iR034] 
      - ys[iCP] * ys[iH2] * ratraw[iR035] 
      - ys[iCP] * ys[iOHX] * ratraw[iR036] 
      - ys[iCP] * ys[iELEC] * ratraw[iR039]
      + ys[iC] * ratraw[iP010]
      + ys[iC] * ratraw[iP014] ;
    
    dydt[iO] = - ys[iH3P] * ys[iO] * ratraw[iR032] 
      + ys[iHEP] * ys[iCO] * ratraw[iR034]
      - ys[iO] * ys[iCHX] * ratraw[iR037]
      + ys[iOHX] * ys[iO] * ratraw[iR045]
      + ys[iOHX] * ratraw[iP002]
      + ys[iCO] * ratraw[iP011] ;
    
    dydt[iCHX] = ys[iH3P] * ys[iC] * ratraw[iR031] 
      + ys[iCP] * ys[iH2] * ratraw[iR035] 
      - ys[iO] * ys[iCHX] * ratraw[iR037] 
      - ys[iCHX] * ratraw[iP001] ;
    
  

    dydt[iOHX] = ys[iH3P] * ys[iO] * ratraw[iR032]
      - ys[iCP] * ys[iOHX] * ratraw[iR036]
      - ys[iC] * ys[iOHX] * ratraw[iR038]
      - ys[iOHX] * ys[iO] * ratraw[iR045]
      - ys[iOHX] * ratraw[iP002] ; 
    
    dydt[iHCOP] = ys[iH3P] * ys[iCO] * ratraw[iR033]
      + ys[iCP] * ys[iOHX] * ratraw[iR036]
      - ys[iHCOP] * ys[iELEC] * ratraw[iR040]
      - ys[iHCOP] * ratraw[iP004] ;
    
    dydt[iCO] = - ys[iH3P] * ys[iCO] * ratraw[iR033] 
      - ys[iHEP] * ys[iCO] * ratraw[iR034]
      + ys[iO] * ys[iCHX] * ratraw[iR037]
      + ys[iC] * ys[iOHX] * ratraw[iR038]
      + ys[iHCOP] * ys[iELEC] * ratraw[iR040]
      + ys[iHCOP] * ratraw[iP004]
      - ys[iCO] * ratraw[iP011] ;
    
    dydt[iM] =  ys[iMP] * ys[iELEC] * ratraw[iR041]
      - ys[iH3P] * ys[iM] * ratraw[iR042] 
      - ys[iM] * ratraw[iP003] ;
    
    dydt[iMP] = - ys[iMP] * ys[iELEC] * ratraw[iR041]
      + ys[iH3P] * ys[iM] * ratraw[iR042]
      + ys[iM] * ratraw[iP003];


    //dydt[iELEC] = likely does not need filling in ;
    

#endif

  //cout << "dydt" << endl;
  //for(i=0;i<NSPECIES;i++]  cout << i << " " << dydt[i] << endl;
  //cout << "dydt end" <<endl;
    
    
};


 
void GetNetworkRates(const REAL Temp, const REAL Den, REAL* ys, REAL* RatRaw, double* rpar)
{
    int i;
    
    double lt1,lt2,lt3,lt4,lt5,lt6,lt7,lt8,lt9;
    double  lnt1,lnt2,lnt3,lnt4,lnt5,lnt6,lnt7,lnt8,lnt9;
    double  te,term,na,t3,naterm, capk, narate;
    double  ncrh, ncrh2, ncrhe, nt, ncr;
    double  ncr2, nv2, nlte2;
    double  nhh;
    double  xh,xh2,xhe,nlte,nv;
    double sqT,logT4,T4;
    double lnte1, lnte2, lnte3, lnte4, lnte5, lnte6, lnte7, lnte8, lnte9;
    
    
    // Radiation parameters
    double j21 = rpar[J21_par];
    // self-shielding factors should be determined here
    // (not external to chemistry) using passed in column densities
    // or passed in photodissociation rates (from ray tracing potentially)
    double fssh2 = 1.0;
    double fsshd = 1.0;
    //double j21 = RadParams[0];
    //double fssh2 = RadParams[1];
    //double fsshd = RadParams[2];
    
    for(i=0;i<NREACTION;i++) RatRaw[i] = 0.0;

    if(Temp < 0.1) {
      cout << "Temp is too low in GetNetworkRates... temp = " << Temp << endl;
      cout << "current ys = " << endl;
      for(int i=0;i<NINT;i++) cout << i << " " << ys[i] << endl;
      AbortChem();
    }
    
    sqT = pow(Temp,0.5);
    T4 = Temp / 1.0e4;
    logT4 = log10(T4);
    
    lt1 = log10(Temp);
    lt2 = pow(log10(Temp),2);
    lt3 = pow(log10(Temp),3);
    lt4 = pow(log10(Temp),4);
    lt5 = pow(log10(Temp),5);
    lt6 = pow(log10(Temp),6);
    lt7 = pow(log10(Temp),7);
    lt8 = pow(log10(Temp),8);
    lt9 = pow(log10(Temp),9);
    
    te  = Temp*0.00008617;
    t3  = Temp/300.0;
    
    lnte1 = log10(te);
    lnte2 = pow(log(te),2);
    lnte3 = pow(log(te),3);
    lnte4 = pow(log(te),4);
    lnte5 = pow(log(te),5);
    lnte6 = pow(log(te),6);
    lnte7 = pow(log(te),7);
    lnte8 = pow(log(te),8);
    lnte9 = pow(log(te),9);
 
    //double mh = 1.67e-24;
    narate = Den / mh / 1.3;
    naterm = Den / mh / 1.3;


    term = 0.0;

#if DO_CHEMISTRY == PRIMORDIAL
    
    int chem_cc_case = 1;

    ncrh  = 3.0-0.416*log10(Temp/10000.0)-0.327*pow((log10(Temp/10000.0)),2);
    ncrh  = pow(10.0,ncrh);
    ncrh2 = 4.845-1.3*log10(Temp/10000.0)+1.62*pow((log10(Temp/10000.0)),2);
    ncrh2 = pow(10.0,ncrh2);
    ncrhe = 5.0792*(1.0-1.23e-5*(Temp-2000.0));
    ncrhe = pow(10.0,ncrhe);
    ncr = 0.0; //FIX LATER
    nt  = 0.0; //FIX LATER
    
    nt = (ys[iHP]+ys[iH]+ys[iHM]+2.0*ys[iH2]+2.0*ys[iH2P]+ys[iHDP]+ys[iHD])*naterm;
    xh  = ys[iH];
    xh2 = ys[iH2];
    xhe = ys[iHE];
    
    
    //ncr = (1.0/(1.0+xhe))*(xh/ncrh+2.0*xh2/ncrh2+xhe/ncrhe);
    //ncr = 1.0/ncr;
    ncr = 0.0;
    ncr = xh/(ncrh) + xh2/ncrh2 + xhe/ncrhe;
    ncr = 1.0/ncr;
    
    
    ncr2 = ncr*100.0;
    
    
    nlte = (nt/ncr)/(1.0+(nt/ncr));
    nv   = 1.0/(1.0+(nt/ncr));
    
    nlte2 = (nt/ncr2)/(1.0+(nt/ncr2));
    nv2 = 1.0/(1.0+(nt/ncr2));
    
    capk = 1.05e-22*pow(Temp,-0.515)*exp(52000.0/Temp);
    
    //ARS special handling of dissociation rates, taken from Glover/Gadget
    double h2var0 = 0, h2var1=0, h2var2=0, yn, yncrh, yncrh2, h2_low_n, h2_high_n;
    double ykdhe, ykdh, ykdh2, ncrinv, abh2;
    double ch1, ch2, ch3, ch4, ch5, ch6, ch22, ch25, cl1 = 0, cl67 = 0;
    double f1, f2, alpha, t4log, temp_1eV;
    //double eV = 1.60219e-12;
    //double kboltz = 1.3806e-16;

    t4log = log10(Temp) - 4.e0;
    temp_1eV = eV / kboltz;

    abh2 = xh2;
    yn = narate;

     if (Temp <= 3.e2)
       ch1 = 1.1e-9 * pow(Temp,0.135e0) * exp(-5.2e4 / Temp);
     else
       ch1 = (3.7e-8 / pow(Temp,0.485e0)) * exp(-5.2e4 / Temp);

     if (Temp < 6.e4) {
       alpha = 1.e0 + 5.48e0 * temp_1eV / Temp;
       f1 = sqrt(Temp) * exp(-alpha);
       }
     else {
       alpha = 1.e0 + 5.48e0 * temp_1eV / 6e4;
       f1 = sqrt(6.e4) * exp(-alpha);
       }

      alpha = 1.e0 + 5.48e0 * temp_1eV / 4.5e3;
      f2 = sqrt(4.5e3) * exp(-alpha);

      h2_low_n  = 1.2e-16 * f1 / f2;
      h2_high_n = ch1;
    
    //Avoid division-by-zero. Note that this doesn't introduce any significant
    //naccuracy, since the test is only true when the dissociation rate is tiny.

    if (h2_high_n < 1.e-40)
      ch3 = 1.e-40;
    else
      ch3 = max(h2_low_n / h2_high_n, 1.e-40);

    ch2 =  6.5e-8 * exp(-5.2e4 / Temp) / pow(Temp,0.485);

    h2_high_n = ch2;

    if (Temp < 3.e4)
       h2_low_n  = 5.996e-30 * pow(Temp,4.1881) * exp(-54657.4e0 / Temp) / pow((1.e0 + 6.761e-06 * Temp),5.6881);
    else
       h2_low_n  = 5.996e-30 * pow(3.e4,4.1881) * exp(-54657.4e0 / 3.e4) / pow((1.e0 + 6.761e-06 * 3.e4),5.6881);

    if (h2_high_n < 1.e-40)
       ch4 = 1.e-40;
    else
       ch4 = max(h2_low_n / h2_high_n, 1.e-40);

    if (Temp < 3.e4) {
      yncrh  = pow(1.e1,( 3.e0 - 0.416e0 * t4log  - 0.327e0 * t4log * t4log ));
      yncrh2 = pow(1.e1,( 4.845e0 - 1.3e0 * t4log + 1.62e0 * t4log * t4log ));
      }
     else {
      yncrh = pow(1.e1,(3.e0 - 0.416e0 * (log10(Temp) - 4.e0) - 0.327e0 * (log10(Temp) - 4.e0) * (log10(Temp) - 4.e0)));
      yncrh2 = pow(1.e1,( 4.845e0 - 1.3e0 * (log10(Temp) - 4.e0) + 1.62e0 * (log10(Temp) - 4.e0) * (log10(Temp) - 4.e0) ));
        }

    ch5 = 1.e0 / yncrh;
    ch6 = 1.e0 / yncrh2;

    ncrinv   = (2.e0 * abh2 * (ch6 - ch5) + ch5);
    h2var0   = 1.e0 / ( 1.e0 + yn * ncrinv);
    h2var1   = pow(ch3,h2var0);
    h2var2   = pow(ch4,h2var0);
    ykdh     = ch1 * h2var1;
    ykdh2    = ch2 * h2var2;
    rpar[h2var0_par] = h2var0;

   //end ARS special handling
    
    // Loads rates for specific reactions
 
    lnt1 = log(te);
    lnt2 = pow(log(te),2);
    lnt3 = pow(log(te),3);
    lnt4 = pow(log(te),4);
    lnt5 = pow(log(te),5);
    lnt6 = pow(log(te),6);
    lnt7 = pow(log(te),7);
    lnt8 = pow(log(te),8);
    lnt9 = pow(log(te),9);
    

   //R001: H + e --> HM + gam
    if(Temp > 6000.0)
        term = pow(10,(-16.4199 + 0.1998*lt2 - 0.005447*lt4 + 0.000040415*lt6));
    else
        term = pow(10,(-17.8450 + 0.762*lt1 + 0.1523*lt2 - 0.03274*lt3));
    RatRaw[iR001]  = term * narate;
    
    //R002: H + HM --> H2 + e
    term = 1.3e-9;   //*t3,(-0.1)
    RatRaw[iR002]  = term * narate;
           
    //R003: H + HP --> H2P + gam
    term = pow(10,(-19.38 - 1.523*lt1 + 1.118*lt2 - 0.1269*lt3));
    RatRaw[iR003]  = term * narate;
    
    //R004: H + H2P --> H2 + HP
    term = 6.4e-10;
    RatRaw[iR004]  = term * narate;
           
    //R005: HP + HM --> H + H
    if(Temp < 2.e4)
    	term = 2.4e-6*pow(Temp,(-0.5))*(1.0+5.0e-5*Temp);
    else
    	term = 4.8e-6/sqrt(Temp);   //ARS added extra temp. constraint to this rate
    RatRaw[iR005]  = term * narate;
    
    
    
    //R006 : H2P + e --> H + H
    if(Temp > 617.0) term = 1.32e-6*pow(Temp,(-0.76));
    else term = 1.0e-8;
    RatRaw[iR006]  = term * narate;
    

    //R007: H2 + HP --> H + H2P
    if(Temp < 110.0)
        term = 0.0e0;
    else if(Temp > 3.0e4)
        term = (-3.3232183e-7 + 3.3735382e-7*log(3.0e4) - 1.4491368e-7*pow(log(3.0e4),2) + 3.4172805e-8*pow(log(3.0e4),3)- 4.7813720e-9*pow(log(3.0e4),4) + 3.9731542e-10*pow(log(3.0e4),5) -1.8171411e-11*pow(log(3.0e4),6) + 3.5311932e-13*pow(log(3.0e4),7))*exp(-21237.15/3.0e4);
    else
        term = (-3.3232183e-7 + 3.3735382e-7*log(Temp) - 1.4491368e-7*pow(log(Temp),2) + 3.4172805e-8*pow(log(Temp),3)- 4.7813720e-9*pow(log(Temp),4) + 3.9731542e-10*pow(log(Temp),5) -1.8171411e-11*pow(log(Temp),6) + 3.5311932e-13*pow(log(Temp),7))*exp(-21237.15/Temp);
                //    term = 2.4e-9*exp(-21200.0/Temp)
                //   term = exp(-24.2491469 + 3.4008244*lnt1 - 3.8980040*lnt2 + 2.0455878*lnt3 - 5.4161829*lnt4 &
                //           +8.4107750*lnt5 - 7.8790262*lnt6 + 4.1383984*lnt7 - 9.3634588*lnt8)
    RatRaw[iR007]  = term * narate;
                

    //R008: H2 + e --> H + H + e
    term = nv*log10(4.49e-9*pow(Temp,(0.11))*exp(-101858.0/Temp));
    term = term + nlte*log10(1.91e-9*pow(Temp,(0.136))*exp(-53407.1/Temp));
    term = pow(10,(term));
    //    term = 4.49e-9*Temp,(0.11)*exp(-101858.0/Temp)
    RatRaw[iR008]  = term * narate;
    
    //R009: H2 + H --> H + H + H
    term = nv*log10(6.67e-12*sqrt(Temp)*exp(-(1.0+63593.0/Temp)));
    term = term + nlte*log10(3.52e-9*exp(-43900.0/Temp)) ; //Simon's Rate
    term = pow(10,(term));
    //term = 6.67e-12*Temp,(0.5)*exp(-(1.0+63593.0/Temp))
    RatRaw[iR009]  = term * narate;
    
    //R010: H2 + H2 --> H + H + H2
    term = (5.996e-30*pow(Temp,(4.1881))/pow((1.0+6.761e-6*Temp),(5.6881)))*exp(-54657.4/Temp);
    term = nv*log10(term);
    term = term + nlte*log10(1.3e-9*exp(-53300.0/Temp)); //Simon's Rate
    term = pow(10,(term));
    RatRaw[iR010]  = term * narate;
    
    //R011: H2 + He --> H + H + He
    term = nv*(-27.029+3.801*lt1-29487.0/Temp);
    term = term + nlte*(-2.729-1.75*lt1-23474.0/Temp);
    //6.6e-10*Temp,(0.115)*exp(-52000.0/Temp))
    term = pow(10,(term));
    //   term = 10,(-27.029+3.801*lt1-29487.0/Temp);
    RatRaw[iR011]  = term * narate;
    
    // R012: H + e --> HP + e + e
    if(Temp < 2800.0)
        term = 0.0e0;
    else
        term = exp(-3.271396786e1 + 1.35365560e1*lnt1 - 5.73932875e0*lnt2 +1.56315498e0*lnt3 - 2.87705600e-1*lnt4 + 3.48255977e-2*lnt5 -2.63197617e-3*lnt6 + 1.11954395e-4*lnt7 -2.03914985e-6*lnt8);
    RatRaw[iR012]  = term * naterm;
    
    //R013: HP + e --> H + gamma
    if(chem_cc_case == 0)
        term = 1.269e-13*pow((315614.0/Temp),(1.503))*pow((1.0+pow((604625.0/Temp),(0.470))),(-1.923)); // Trying case A
    else if (chem_cc_case == 1)
        term = 2.753e-14*pow((315614.0/Temp),(1.500))*pow((1.0+pow((115188.0/Temp),(0.407))),(-2.242));  // This is case B
    else
        term = 2.753e-14*pow((315614.0/Temp),(1.500))*pow((1.0+pow((115188.0/Temp),(0.407))),(-2.242));  // This is case B
                                                          
    RatRaw[iR013]  = term * narate;
   
    //R014: HM + e --> H + e + e
    if(Temp < 100.0)
        term = 0.0e0;
    else
        term =exp(-1.801849334e1 + 2.36085220*lnt1  -2.82744300e-1*lnt2 + 1.62331664e-2*lnt3  -3.36501203e-2*lnt4 + 1.17832978e-2*lnt5 - 1.65619470e-3*lnt6   +1.06827520e-4*lnt7 - 2.63128581e-6*lnt8);
    RatRaw[iR014]  = term * narate;
   
    //R015: HM + H --> H + H + e
    if(te > 0.1)
        term = exp(-2.0372609e1 + 1.13944933e0*lnt1 - 1.4210135e-1*lnt2
                   +8.4644554e-3*lnt3 - 1.4327641e-3*lnt4 + 2.0122503e-4*lnt5
                   +8.6639632e-5*lnt6 - 2.5850097e-5*lnt7 + 2.4555012e-6*lnt8
                   -8.0683825e-8*lnt9);
    else
        term = 2.5634e-9*pow(te,(1.78186));
    RatRaw[iR015]  = term * narate;
    
    //R016: HP + HM --> H2P + e
    if(Temp > 8000.0)
        term = 9.6e-7*pow((Temp),(-0.90));
    else
        term = 6.9e-9*pow((Temp),(-0.35));
    RatRaw[iR016]  = term * narate;
    
    // R017: He + e --> HE+ + e + e
    if(Temp < 2800.0)
        term = 0.0e0;
    else
        term = exp(-4.409864886e1 + 2.391596563e1*lnt1 - 1.07532302e1*lnt2
                   +3.05803875e0*lnt3 - 5.68511890e-1*lnt4 + 6.79539123e-2*lnt5
                   -5.00905610e-3*lnt6 + 2.06723616e-4*lnt7 - 3.64916141e-6*lnt8);
    RatRaw[iR017]  = term * naterm;
    
    // R018: HEP + e --> HEPP + e + e
    if(Temp < 5500.0)
        term = 0.0e0;
    else
        term = exp(-6.87104099e1 + 4.393347633e1*lnt1 - 1.84806699e1*lnt2
                   +4.70162649e0*lnt3 - 7.6924663e-1*lnt4 + 8.113042e-2*lnt5
                   -5.32402063e-3*lnt6 + 1.97570531e-4*lnt7 - 3.16558106e-6*lnt8);
    
    RatRaw[iR018] = term * naterm;
    
    
    //R019: Hep + e --> He + gamma
    term = 0.32*(1.0e-11*pow((Temp),(-0.5))*(11.19 - 1.676*lt1 - 0.2852*lt2 + 0.04433*lt3)); // This is case B
    term = term + 0.68*(1.0e-11*pow((Temp),(-0.5))*(12.72 - 1.615*lt1 - 0.3162*lt2 + 0.0493*lt3)); // This is case A,
    term = term + 1.9e-3*pow(Temp,(-1.50))*exp(-473421.0/Temp)*(1.0+0.3*exp(-94684.0/Temp));
    //term = term + pow(Temp,(-1.5))*(5.966e-4*exp(-455600.0/Temp)+1.613e-4*exp(-555200.0/Temp)-2.223e-5*exp(-898200/Temp));
    RatRaw[iR019] = term * narate;

     //R020: Hepp + e --> Hep + gam
     if(chem_cc_case == 0 )
         term = 2.538e-13*pow((1262456.0/Temp),1.503)*pow((1.0+pow((2418500.0/Temp),0.407)),(-1.923)); // CASE A
     else if (chem_cc_case == 1)
         term = 5.506e-14*pow((1262456.0/Temp),1.500)*pow((1.0+pow((460752.0/Temp),0.407)),(-2.242)); // CASE B
     else
         term = 5.506e-14*pow((1262456.0/Temp),1.500)*pow((1.0+pow((460752.0/Temp),0.407)),(-2.242)); // CASE B
     
     RatRaw[iR020] = term * narate;
    
    
   
     //R021: H2P + HM --> H2 + H
     term = 1.4e-7*pow(t3,(-0.5));
     RatRaw[iR021] = term * narate;
     
     //R022: H2P + HM --> H + H + H
     term = 1.4e-7*pow(t3,(-0.5));
     RatRaw[iR022] = term * narate;
     
     //R023: H2 + e --> H + HM
     if(Temp < 500.0)
         term = 0.0e0;
     else
         term = 2.7e-8*pow(Temp,(-1.27))*exp(-43000.0/Temp);
     RatRaw[iR023] = term * narate;
     
     //R024: H2 + HeP --> He + H + HP
     term = 3.7e-14*exp(35.0/Temp);
     RatRaw[iR024] = term * narate;
     
     //R025: H2 + HeP --> HE + H2P
     term = 7.2e-15;
     RatRaw[iR025] = term * narate;
     
     //R026: H + HEP --> HE + HP + GAM
     term = 1.2e-15*pow(t3,(0.25));
     RatRaw[iR026] = term * narate;
     
     //R027: He + HP --> H + HeP
     if(Temp > 10000.0)
         term = 4.0e-37*pow(Temp,(4.74));
     else if(Temp < 1500.0)
         term = 0.0e0;
     else
         term = 1.26e-9*pow(Temp,(-0.75))*exp(-127500.0/Temp);
     
     RatRaw[iR027] = term * narate;
    
     //R028: HeP + HM --> He + H
     term = 2.32e-7*pow(t3,(-0.52))*exp(Temp/22400.0);
     RatRaw[iR028] = term * narate; 
     //ARS adding  high temp fit (see Stancil, Lepp & Dalgarno, 1998, ApJ, 509, 1)
     if(Temp > 2.24e4)
       RatRaw[iR028] = 2.32e-7 * pow(2.24e4/3e2,-0.52e0) * narate;
 
     //R029: HM + He --> H + He + e
     if(Temp < 250.0)
         term = 0.0e0;
     else
         term = 4.1e-17*pow(Temp,2.0)*exp(-19870.0/Temp);
     
    RatRaw[iR029] = term * narate;
    
     //R030: H + H + H --> H2 + H
     //  term = 1.44e-26*pow(Temp,(-1.54));  //7.7e-31*Temp,(-0.464)
     if(Temp < 300.0) 
         term = 1.14e-31*pow(Temp,(-0.38)); // -31   <-- ABN02 rate
     else
         term = 3.9e-31*pow(Temp,(-1.0)); // -31   <-- ABN02 rate
     
     //R030: 3H --> H2 + H (Forrey 2013)
     term = 6.0e-32*pow(Temp,-0.25) + 2.0e-31/sqrt(Temp);

     RatRaw[iR030] = term * narate * narate; //dens     

     //R031: H + H + H2 --> H2 + H2
     //   term = 1.44e-26*Temp,(-1.54)/8.0  ////9.625e-32*Temp,(-0.464)
     term = RatRaw[iR030]/8.0;
     RatRaw[iR031] = term; 
     //RatRaw[iR031] = term * narate * narate; //dens     

     //R032: H + H + He --> H2 + He
    term = 6.9e-32*pow(Temp,(-0.4));
     RatRaw[iR032] = term * narate * narate; //dens

   //ARS adding extra special handling of dissociation rates for high density:
   if(narate > 1.e8)
     {

     RatRaw[iR009] = RatRaw[iR030]/capk / narate;
     RatRaw[iR010] = RatRaw[iR031]/capk / narate;
     RatRaw[iR011] = RatRaw[iR032]/capk / narate;

     }
     
     //R033: DP + e --> D + gamma
     RatRaw[iR033] = RatRaw[iR013];
     
     //R034: D + HP --> H + DP
     if(Temp > 2.0e5)
         term = 3.44e-10*pow(Temp,(0.35));
     else
         term = 2.0e-10*pow(Temp,(0.402))*exp(-37.1/Temp)-3.31e-17*pow(Temp,(1.48));
     
     RatRaw[iR034] = term * narate;
     
     //R035: H + DP --> D + HP
     term = 2.06e-10*pow(Temp,(0.396))*exp(-33.0/Temp)+2.03e-9*pow(Temp,(-0.332));
     RatRaw[iR035] = term * narate;
     
     //R036: H + D --> HD + gam
     if(Temp > 200.0) 
         term = 1.0e-25*exp(507.207-370.889*log(Temp) + 104.854*pow(log(Temp),2.0)     -14.4192*pow(log(Temp),3.0) + 0.971469*pow(log(Temp),4.0) -   0.0258076*pow(log(Temp),5.0));
     else if (Temp < 10.0)
         term = 0.0;
     else
         term = 1.0e-25*(2.80202 - 6.63697*log(Temp) + 4.75619*pow(log(Temp),2.0)-1.39325*pow(log(Temp),3.0) + 0.178259*pow(log(Temp),4.0)-0.00817097*pow(log(Temp),5.0));
     
     RatRaw[iR036] = term * narate;
     
     //R037: H2 + D --> HD + H
     if(Temp > 2000.0) 
         term = 3.17e-10*exp(-5207.0/Temp);
     else
     {
         term = (-56.4737 + 5.888886*lt1 + 7.19692*lt2 + 2.25069*lt3 - 2.16903*lt4 + 0.317887*lt5);
         term = pow(10,term);
     }
     
     RatRaw[iR037] = term * narate;
     
     //R038: HDP + H --> HD + HP
     RatRaw[iR038] = RatRaw[iR004];
     
     //R039: H2 + DP --> HD + HP
     term = 4.17e-10 + 8.46e-10*lt1 - 1.37e-10*lt2;
     RatRaw[iR039] = term * narate;
     
     //R040: HD + H --> H2 + D
     if(Temp > 200.0) 
         term = 5.25e-11*exp(-4430.0/Temp + 173900.0/pow(Temp,2));
     else if(Temp < 50.0) 
         term = 0.0e0;
     else
         term = 5.25e-11*exp(-4430.0/Temp);
     
     RatRaw[iR040] = term * narate;
     
     //R041: HD + HP --> H2 + D+
     term = 1.1e-9*exp(-488.0/Temp);
     RatRaw[iR041] = term * narate;
     
     //R042: D + HP --> HDP + gam
     term = 3.9e-19*pow(t3,(1.8))*exp(20.0/Temp);
     RatRaw[iR042] = term * narate;
     
     //R043: H + DP --> HDP + gam
     term = 3.9e-19*pow(t3,(1.8))*exp(20.0/Temp);
     RatRaw[iR043] = term * narate;
     
     //R044: HDP + e --> H + D
     term = 7.2e-8*pow(Temp,(-0.5));
    RatRaw[iR044] = term * narate; //Had an 'l' at the end of narate? l;
     
     // R045: D + e --> DP + e + e
    RatRaw[iR045] = RatRaw[iR012];
    
     //R046: D + HEP --> HE + DP + gam
     term = 1.1e-15*pow(t3,(0.25));
     RatRaw[iR046] = term * narate;
     
     //R047: He + DP --> D + HeP
     if(Temp > 10000.0) 
         term = 5.9e-37*pow(Temp,(4.74));
     else if(Temp < 1500.0)
         term = 0.0e0;
     else
         term = 1.85e-9*pow(Temp,(-0.75))*exp(-127500.0/Temp);
     
     RatRaw[iR047] = term * narate;
     
     //R048: H2P + D --> HDP + H
     term = 1.07e-9*pow(t3,(0.062))*exp(-Temp/41400.0);
     RatRaw[iR048] = term * narate;
     
     //R049: D + HDP --> HD + DP
     term = 6.4e-10;
     RatRaw[iR049] = term * narate;
     
     //R050: HDP + H --> H2P + D
     term = 1.0e-9*exp(-154.0/Temp);
     RatRaw[iR050] = term * narate;
     
     //R051: D + e --> DM + gam
     RatRaw[iR051] = RatRaw[iR001];
     
     //R052: H + DM --> D + HM
     term = 6.4e-9*pow(t3,(0.41));
     RatRaw[iR052] = term * narate;
     
     //R053: D + HM --> D + DM
     term = 6.4e-9*pow(t3,(0.41));
     RatRaw[iR053] = term * narate;
     
     //R054: D + HM --> HD + e
     //term = 1.5e-9*pow(t3,(-0.1));
     RatRaw[iR054] = 0.5*RatRaw[iR002];
     
     //R055: H + DM --> HD + e
     //term = 1.5e-9*pow(t3,(-0.1));
     RatRaw[iR055] = 0.5*RatRaw[iR002];
     
     //R056: D + DM --> D2 + e
     //term = 1.6e-9*pow(t3,(-0.1));
     RatRaw[iR056] = RatRaw[iR002];
     
     //R057: HD + e --> DM + H
     if(Temp < 500.0) 
         term = 0.0e0;
     else
         term = 1.35e-9*pow(Temp,(-1.27))*exp(-43000.0/Temp);
     
     RatRaw[iR057] = term * narate;
     
     //R058: HD + e --> D + HM
     term = 1.35e-9*pow(Temp,(-1.27))*exp(-43000.0/Temp);
     RatRaw[iR058] = term * narate;
     
     //R059: D2 + e --> D + DM
     term = 6.7e-11*pow(Temp,(-1.27))*exp(-43000.0/Temp);
     RatRaw[iR059] = term * narate;
     
     //R060: HP + DM --> HDP + e
    term = 1.1e-9*pow(t3,(-0.4));
     RatRaw[iR060] = term * narate;
     
     //R061: DP + HM --> HDP + e
    term = 1.1e-9*pow(t3,(-0.4));
     RatRaw[iR061] = term * narate;
     
     //R062: DP + DM --> D2P + e
     term = 1.3e-9*pow(t3,(-0.4));
     RatRaw[iR062] = term * narate;
     
     //R063: DM + e --> D + e +e
     RatRaw[iR063] = RatRaw[iR014];
     
     //R064: DM + H --> D + H +e
     RatRaw[iR064]= RatRaw[iR015];
     
     //R065: DM + He --> D + He + e
     if(Temp < 250.0) 
         term = 0.0e0;
     else
         term = 1.5e-17*pow(Temp,2.0)*exp(-19870.0/Temp);
     
     RatRaw[iR065] = term * narate;
     
     //R066: DP + HM --> D + H
     RatRaw[iR066]=RatRaw[iR005];
     
     //R067: HP + DM --> D + H
     RatRaw[iR067]=RatRaw[iR005];
     
     //R068: DP + DM --> D + D
     RatRaw[iR068]=RatRaw[iR005];
     
     //R069: H2P + DM --> H2 + D
     term = 1.7e-7*pow(t3,(-0.5));
     RatRaw[iR069] = term * narate;
     
     //R070: H2P + DM --> H + H + D
     term = 1.7e-7*pow(t3,(-0.5));
     RatRaw[iR070] = term * narate;
     
     //R071: HDP + HM --> HD + H
     term = 1.5e-7*pow(t3,(-0.5));
     RatRaw[iR071] = term * narate;
     
     //R072: HDP + HM --> D + H + H
    term = 1.5e-7*pow(t3,(-0.5));
     RatRaw[iR072] = term * narate;
     
     //R073: HDP + DM --> HD + D
     term = 1.9e-7*pow(t3,(-0.5));
     RatRaw[iR073] = term * narate;
     
     //R074: HDP + DM --> D + H + D
     term = 1.9e-7*pow(t3,(-0.5));
     RatRaw[iR074] = term * narate;
     
     //R075: D2P + HM --> D2 + H
     term = 1.5e-7*pow(t3,(-0.5));
     RatRaw[iR075] = term * narate;
     
     //R076: D2P + HM --> D + D + H
    term = 1.5e-7*pow(t3,(-0.5));
     RatRaw[iR076] = term * narate;
     
     //R077: D2P + DM --> D2 + D
     term = 2.0e-7*pow(t3,(-0.5));
     RatRaw[iR077] = term * narate;
     
     //R078: D2P + DM --> D + D + D
     term = 2.0e-7*pow(t3,(-0.5));
     RatRaw[iR078] = term * narate;
     
     //R079: HEP + DM --> He + D
     term = 3.03e-7*pow(t3,(-0.52))*exp(Temp/22400.0);
     RatRaw[iR079] = term * narate;
     //ARS adding  high temp fit (see Stancil, Lepp & Dalgarno, 1998, ApJ, 509, 1)
     if(Temp > 2.24e4)
       RatRaw[iR079] = 2.32e-7 * pow(2.24e4/3e2,-0.52e0) * narate;
 
     //R080: D + DP --> D2P + gam
     term = 1.9e-19*pow(t3,(1.8))*exp(20.0/Temp);
     RatRaw[iR080] = term * narate;
     
     //R081: D + H2P --> H2 + DP
     term = 6.4e-10;
     RatRaw[iR081] = term * narate;
     
     //R082: H2P + D --> HD + HP
     term = 1.0e-9;
     RatRaw[iR082] = term * narate;
     
     //R083: HDP + H --> H2 + DP
     term = 1.0e-9;
     RatRaw[iR083] = term * narate;
     
     //R084: HDP + D --> D2P + H
     term = 1.0e-9;
     RatRaw[iR084] = term * narate;
     
     //R085: HDP + D --> D2 + HP
     term = 1.0e-9;
     RatRaw[iR085] = term * narate;
     
     //R086: D + D2P --> D2 + DP
     term = 6.4e-10;
     RatRaw[iR086] = term * narate;
     
     //R087: H + D2P --> D2 + HP
     term = 6.4e-10;
     RatRaw[iR087] = term * narate;
     
     //R088: D2P + H --> HDP + D
     term = 1.0e-9*exp(-472.0/Temp);
     RatRaw[iR088] = term * narate;
     
     //R089: D2P + H --> HD + DP
     term = 1.0e-9;
     RatRaw[iR089] = term * narate;
     
     //R090: H2 + DP --> D + H2P
     RatRaw[iR090] = RatRaw[iR007];
     
     
     //R091: H2 + DP --> HDP + H
     if(Temp < 230.0) 
         term = 0.0e0;
     else
         term = (1.04e-9 + 9.52e-9*(Temp/10000.0) - 1.81e-9*pow((Temp/10000.0),2))*exp(-21000.0/Temp);
     
     RatRaw[iR091] = term * narate;
     
     //R092: HD + HP --> H + HDP
     RatRaw[iR092] = RatRaw[iR007];
     
     //R093: HD + HP --> H2P + D
     if(Temp < 230.0) 
         term = 0.0e0;
     else
         term = 1.0e-9*exp(-21600.0/Temp);
     
    RatRaw[iR093] = term * narate;
    
     //R094: HD + DP --> D + HDP
     RatRaw[iR094] = RatRaw[iR007];
     
     //R095: HD + DP --> D2 + HP
     term = 1.0e-9;
     RatRaw[iR095] = term * narate;
     
     //R096: HD + DP --> D2P + H
     term = (3.54e-9 + 7.50e-10*(Temp/10000.0) - 2.92e-10*pow((Temp/10000.0),2))*exp(-21100.0/Temp);
     RatRaw[iR096] = term * narate;
     
     //R097: D2 + HP --> HD + DP
     term = 2.1e-9*exp(-491.0/Temp);
     RatRaw[iR097] = term * narate;
     
     //R098: D2 + HP --> HDP + D
     term = (5.18e-11 + 3.05e-9*(Temp/10000.0) - 5.42e-10*pow((Temp/10000.0),2))*exp(-20100.0/Temp);
     RatRaw[iR098] = term * narate;
     
     //R099: D2 + HP --> H + D2P
     RatRaw[iR099] = RatRaw[iR007];
     
     //R100: D2 + DP --> D2P + D
     RatRaw[iR100] = RatRaw[iR007];
     
     //R101: HD + HEP --> HE + HDP
     term = 7.2e-15;
     RatRaw[iR101] = term * narate;
     
     //R102: HD + HEP --> HE + HP + D
     term = 1.85e-14*exp(35.0/Temp);
     RatRaw[iR102] = term * narate;
     
     //R103: HD + HEP --> HE + HP + D
     term = 1.85e-14*exp(35.0/Temp);
     RatRaw[iR103] = term * narate;
     
     //R104: D2 + HeP --> HE + D2P
     term = 2.5e-14;
     RatRaw[iR104] = term * narate;
     
     //R105: D2 + HeP --> He + DP + D
     term = 1.1e-13*pow(t3,(-0.24));
     RatRaw[iR105] = term * narate;
     
     //R106: HD + D --> D2 + H
     term = 1.15e-11*exp(-3220.0/Temp);
     RatRaw[iR106] = term * narate;
     
     //R107: D2 + H --> HD + D
     if(Temp > 2200.0) 
         term = 2.67e-10*exp(-5945.0/Temp);
     else
     {
         term = -86.1558 + 4.53978*lt1 + 33.5707*lt2 - 13.0449*lt3 + 1.22017*lt4 + 0.0482453*lt5;
         term = pow(10,term);
     }
     
     RatRaw[iR107] = term * narate;
     
     //R108: HD + H --> H + D +H
     term = nv2*log10(6.67e-12*pow(Temp,(0.5))*exp(-(1.0+63593.0/Temp)));
     term = term + nlte2*log10(3.52e-9*exp(-43900.0/Temp));
     term = pow(10,(term));
     RatRaw[iR108 ] = term * narate; //RatRaw[iR009)
     
     //R109: HD + H2 --> H + D + H2
     term = (5.996e-30*pow(Temp,(4.1881))/pow((1.0+6.761e-6*Temp),(5.6881)))*exp(-54657.4/Temp);
     term = nv2*log10(term);
     term = term + nlte2*log10(1.3e-9*exp(-53300.0/Temp)); //Simon's Rate
     term = pow(10,(term));
     RatRaw[iR109] = term * narate;
     //   RatRaw[iR109RatRaw[iR010)
     
     //R110: HD + He --> H + D + He
     term = nv2*(-27.029+3.801*lt1-29487.0/Temp);
     term = term + nlte2*(-2.729-1.75*lt1-23474.0/Temp) ;  //6.6e-10*Temp,(0.115)*exp(-52000.0/Temp))
     term = pow(10,(term));
     //   term = pow(10,(-27.029+3.801*lt1-29487.0/Temp));
     RatRaw[iR110] = term * narate;
     //   RatRaw[iR110RatRaw[iR011)
     
     //R111: HD + e --> H + D + e
     term = 0.0e0;
     term = nv2*log10(5.09e-9*pow(Temp,(0.128))*exp(-103258.0/Temp));
     term = term + nlte2*log10(1.04e-9*pow(Temp,(0.218))*exp(-53070.7/Temp));
     term = pow(10,(term));
     //    term = 5.09e-9*pow(Temp,(0.128))*exp(-103258.0/Temp);
     RatRaw[iR111] = term * narate;
     
     //R112: D2 + H --> D + D + H
     RatRaw[iR112]=RatRaw[iR009];
     
     //R113: D2 + H2 --> D + D + H2
     RatRaw[iR113]=RatRaw[iR010];
     
     //R114: D2 + He --> D + D + He
     RatRaw[iR114]=RatRaw[iR011];
     
     //R115: D2 + e --> D + D + e
     term = 0.0e0;
     term = nv*log10(8.24e-9*pow(Temp,(0.126))*exp(-105388.0/Temp));
     term = term + nlte*log10(2.75e-9*pow(Temp,(0.163))*exp(-53339.7/Temp));
     term = pow(10,(term));
     //   term = 8.24e-9*pow(Temp,(0.126))*exp(-105388.0/Temp);
     RatRaw[iR115] = term * narate;
     
     //R116: HM + gam --> H + e
     term = 1.36e-11*j21;
     RatRaw[iR116]=term;
     
     //R117: DM + gam --> D + e
     term = 1.36e-11*j21;
     RatRaw[iR117]=term ;
     
     //R118: H2P + gam --> H + HP
     term = 4.11e-12*j21;
     RatRaw[iR118]=term;
     
     //R119: HDP + gam --> H + DP
     term = 2.05e-12*j21;
     RatRaw[iR119]=term;
     
     //R120: HDP + gam --> D + HP
     term = 2.05e-12*j21;
     RatRaw[iR120]=term;
     
     //R121: D2P + gam --> D + DP
     term = 4.11e-12*j21;
     RatRaw[iR121]=term;
     
     //R122: H2 + gam --> H + H
     term = 1.3e-12*fssh2*j21;
     RatRaw[iR122]=term;
     
     //R123: HD + gam --> H + D
     term = 1.45e-12*fsshd*j21;
     RatRaw[iR123]=term;
     
     //R124: D2 + gam --> D + D
     term = 1.3e-12*j21;
     RatRaw[iR124]=term;

#endif



#if DO_CHEMISTRY == CONTEMPORARY


    lnt1 = log(Temp);
    lnt2 = pow(log(Temp),2);
    lnt3 = pow(log(Temp),3);
    lnt4 = pow(log(Temp),4);
    lnt5 = pow(log(Temp),5);
    lnt6 = pow(log(Temp),6);
    lnt7 = pow(log(Temp),7);
    lnt8 = pow(log(Temp),8);
    lnt9 = pow(log(Temp),9);
    
    
    
    // GLOVER H/He network - what I need, numbering according to table
    // glover ordering           new ordering
    // 1-19                         1-19
    // 30                            20
    // 88-89                          21-22
    // 109                             23
    // 110-111 (may be overwritten by NL99 rates)   24,25
    // 141                              26
    // 154-156                          27, 28, 29
    // 165                              30
    // 54, 55                         
    
    // Glover photochemical rates
    // 166, 167, 168, 169, 170
    
    //NL99 photo rates
    // will only use
    // CHX + gamma -> CI + H
    // OHX + gamma -> OI + H
    // M + gamma -> M+ + e
    // and from glover:
    // 166, 167, 168, 169, 170, 171
    // 198
    // and CR rates
    // 199, 200, (H2 + cr from NL99)
    // 205, 
    // No induced cr ray rates from glover
    
    // from NL99
    // 1, cr + h2
    // 3, 4, 5, 
    // 7, 8, 9
    // both neut-neut: 10, 11
    // elec recomb
    // 14, c+ + e- -> cI +gamma
    // 15, 16
    // charge transfer, 17
    
    


   //R001: H + e --> HM + gam
    if(Temp > 6000.0)
        term = pow(10,(-16.4199 + 0.1998*lt2 - 0.005447*lt4 + 0.000040415*lt6));
    else
        term = pow(10,(-17.8450 + 0.762*lt1 + 0.1523*lt2 - 0.03274*lt3));
    RatRaw[iR001]  = term * narate;
  
  
    // H- + H -> H2 + e-
    if(Temp < 300.0) 
      term = 1.5e-9;
    else
      term = 4.0e-9 * pow(Temp,-0.17);
    RatRaw[iR002] = term * narate;
    

    // H + H+ -> H2+ + gamma
    term = pow(10.0,(-19.38 - 1.523 * lt1 + 1.118 * lt2 - 0.1269 * lt3));
    RatRaw[iR003] = term * narate;
  

    // H + H2+ -> H2 + H+
    RatRaw[iR004] = 6.4e-10  * narate;
  

    // H- + H+ -> H + H 
    term = 2.4e-6 / sqT * (1.0 + Temp / 2.0e4);
    RatRaw[iR005] = term * narate;
    

    // H2+ + e- -> H + H
    if(Temp < 617.0) 
      term = 1.0e-8;
    else
      term = 1.32e-6 * pow(Temp,-0.76);
    RatRaw[iR006] = term * narate;
  
    
    // H2 + H+ -> H2+ + H
    RatRaw[iR007] = (-3.3232183e-7 + 3.3735382e-7 * lnt1  
		     - 1.4491368e-7 * lnt2 + 3.4172805e-8 * lnt3 
		     - 4.7813720e-9 * lnt4 + 3.9731542e-10 * lnt5 
		     - 1.8171411e-11 * lnt6 + 3.5311932e-13 * lnt7)  
      * exp(-21237.15 / Temp) * narate;
    
  
    // H2 + e- -> H + H + e-
    RatRaw[iR008] = 3.73e-9 * pow(Temp,0.1121) * exp(-99430.0 / Temp) * narate;
  
    
  


 // Special handing of H2 dissociation rates:
  double noncr;
  double yh, yhe;
  double kl, kh;

  if(Temp  < 3.0e4) {
    // ncr, H
    ncrh = pow(10.0,3.0 - 0.416 * logT4 - 0.327 * pow(logT4,2));
  
    // ncr, H2
    ncrh2 = pow(10.0,4.845 - 1.3 * logT4 + 1.62 * pow(logT4,2));
  } else
    {
      ncrh = pow(10.0,3.0);
      ncrh2 = pow(10.0,4.845);
    }


  ncr = 1.0 / (ys[iH] / ncrh + ys[iH2] /  ncrh2);
  noncr = naterm / ncr;

  // k9, l: low rate for H2 + H -> 3H
  kl = 6.67e-12 * sqT * exp(-(1.0 + 63590.0 / Temp)); 
  
  // k9, h: high rate for H2 + H -> 3H
  kh = 3.52e-9 * exp(-43900.0 / Temp);
  
  // H2 + H -> H + H +H
  term  = pow(10.0,(( noncr / (1.0 + noncr)) * log10(kh) + (1.0 / (1.0 + noncr)) * log10(kl)));
  RatRaw[iR009]  = term * narate;

  // k10, l: low rate for H2 + H2 -> H2 + 2H
  kl = 5.996e-30 * pow(Temp,4.1881) / pow((1.0 + 6.761e-6 * Temp),5.6881) * exp(-54657.4 / Temp);
  
  // k10, h: high rate for H2 + H2 -> H2 + 2H
  kh = 1.3e-9 * exp(-53300.0 / Temp);
  
  // H2 + H2 -> H2 + H + H
  term = pow(10.0,( ( noncr / (1.0 + noncr)) * log10(kh) + (1.0 / (1.0 + noncr)) * log10(kl)));
  RatRaw[iR010]  = term * narate;

			  

    // R011: H + e --> H+ + e + e
    if(Temp < 2800.0)
      term = 0.0e0;
    else
      term = exp(-3.271396786e1 + 1.35365560e1*lnte1 - 5.73932875e0*lnte2 +1.56315498e0*lnte3 - 2.87705600e-1*lnte4 + 3.48255977e-2*lnte5 -2.63197617e-3*lnte6 + 1.11954395e-4*lnte7 -2.03914985e-6*lnte8);
    RatRaw[iR011]  = term * naterm;
  
  
  
    // case B electron recombination
    // H+ + e -> H + gamma
    term = 2.753e-14 * pow(315614.0 / Temp,1.5)  * pow(1.0 + pow(115188.0 / Temp,0.407),-2.242);
    RatRaw[iR012]  = term * naterm;
    

 

    // H- + e --> H + e + e
    if(Temp < 100.0)
        term = 0.0e0;
    else
      term =exp(-1.801849334e1 + 2.36085220*lnte1  -2.82744300e-1*lnte2 + 1.62331664e-2*lnte3  -3.36501203e-2*lnte4 + 1.17832978e-2*lnte5 - 1.65619470e-3*lnte6   +1.06827520e-4*lnte7 - 2.63128581e-6*lnte8);
    RatRaw[iR013]  = term * narate;
  
  
    
    // H- + H -> H + H + e-
  if(te < 0.1) 
    term = 2.5634e-9 * pow(te,1.78186);
  else
    term = exp(-2.0372609e1 
	       + 1.13944933 * lnte1 
	       - 1.4210135e-1 * lnte2
	       + 8.4644554e-3 * lnte3 
	       - 1.4327641e-3 * lnte4 
	       + 2.0122503e-4 * lnte5 
	       + 8.6639632e-5 * lnte6 
	       - 2.5850097e-5 * lnte7 
	       + 2.4555012e-6 * lnte8 
	       -8.0683825e-8 * lnte9);
  RatRaw[iR014]  = term * narate;
  
  
  // H- + H+ -> H2+ + e-
  if(Temp < 8000.0) 
    term = 6.9e-9 * pow(Temp,(-0.35));
  else
    term = 9.6e-7 * pow(Temp,(-0.90));
  RatRaw[iR015]  = term * narate;
  
  
  // He + e- -> He+ + e- + e-
  term = exp(-4.409864886e1  
	     + 2.391596563e1 * lnte1  
	     - 1.07532302e1 * lnte2  
	     + 3.05803875 * lnte3  
	     - 5.6851189e-1 * lnte4 
	     + 6.79539123e-2 * lnte5  
	     - 5.0090561e-3 * lnte6  
	     + 2.06723616e-4 * lnte7  
	     - 3.64916141e-6 * lnte8);
  RatRaw[iR016]  = term * narate;
  

  // case B He+ recombination
  // He+ + e- -> He + gamma
  term = 1.0e-11 / sqT * (11.19 - 1.676 * lt1 - 0.2852 * lt2 + 0.04433 * lt3);
  RatRaw[iR017]  = term * narate;
  
  // He+ + H -> He + H+
  term = 1.25e-15 * pow(t3,0.25);
  RatRaw[iR018]  = term * narate;
    
  // He + H+ -> He+ + H
  if(Temp < 1.0e4) 
    term = 1.26e-9 * pow(Temp,(-0.75)) * exp(-127500.0 / Temp);
  else
    term = 4.0e-37 * pow(Temp,(4.74));
  RatRaw[iR019]  = term * narate;
  
  
  // k30, l:
  kl = pow(10.0,(-27.029 + 3.801 * lt1 - 29487.0 / Temp))  ;
  
  // k30, h:
  kh = pow(10.0,(-2.729 - 1.75 * lt1 - 23474.0 / Temp))  ;
  
  // ncr, He (note: subtle typo in Glover et al. 2010)
  if(Temp < 1.0e4) 
    ncrhe = pow(10.0,(5.0792 * (1.0 - 1.23e-5 * (Temp - 2000.0))));
  else
    {
      double tfix = 1.0e4;
      ncrhe = pow(10.0,(5.0792 * (1.0 - 1.23e-5 * (tfix - 2000.0))));
    }

  // H2 + He -> H + H + He
  yhe = YHELIUM - ys[iHEP];
  ncr = 1.0 / (yhe / ncrhe + ys[iH2] /  ncrh2);
  noncr = naterm / ncr;
  term = pow(10.0,( ( noncr / (1.0 + noncr)) * log10(kh) + (1.0 / (1.0 + noncr)) * log10(kl)));
  RatRaw[iR020]  = term * narate;


  // H2 + He+ -> He + H2+
  term = 7.2e-15;
  RatRaw[iR021]  = term * narate;

  // H2 + He+ -> He + H + H+
  term = 3.7e-14 * exp(-35.0 / Temp);
  RatRaw[iR022]  = term * narate;


  // He+ + H- -> He + H
  if(Temp < 22400.0)
    term = 2.32e-7 * pow(t3,-0.52) * exp(Temp / 22400.0);
  else
    term = 2.32e-7 * pow(74.67,-0.52); 
  RatRaw[iR023]  = term * narate;
  
  // H3+ + e- -> H2 + H
  term = 2.34e-8 * pow(t3,-0.52);
  RatRaw[iR024]  = term * narate;
  
  // H3+ + e- -> H + H + H
  term = 4.36e-8 * pow(t3,-0.52);
  RatRaw[iR025]  = term * narate;


  // H2 + H+ -> H3+ + gamma
  term = 1.0e-16;
  RatRaw[iR026]  = term * narate;
  
  

  // 154 - 156: 3 body rates
  
  // H + H + H -> H2 + H
  if(Temp <  300.0)
    term = 1.32e-32 * pow(t3,-0.38); 
  else
    term = 1.32e-32 / t3;
  RatRaw[iR027]  = term * narate * narate;
  
  // H + H + H2 -> H2 + H2
  term = 2.8e-31 * pow(Temp,-0.6);
  RatRaw[iR028]  = term * narate * narate;

  // H + H + He -> H2 + He
  term = 6.9e-32 * pow(Temp,-0.4);
  RatRaw[iR029]  = term * narate * narate;


  // H2 formation on dust grains. 
  // H + H + dust -> H2 + dust
  double fa, sh, Tdust;  
  Tdust = rpar[tdust_par];
  fa = 1.0 / (1.0 + exp(7.5e2 * (1.0/75.0 - 1.0/Tdust))); 
  sh =  1.0 / (1.0 + 4.0e-2 * pow(Temp + Tdust,0.5) + 2.0e-3 * Temp + 8.0e-6 * Temp * Temp);
  term = 3.025e-17 * pow(Temp/100.0,0.5) * fa * sh;
  RatRaw[iR030]  = term * narate;


 
  
  
  // NL99 rates

    // from NL99
    // 1, cr + h2
    // 3, 4, 5, 
    // 7, 8, 9
    // both neut-neut: 10, 11
    // elec recomb
    // 14, c+ + e- -> cI +gamma
    // 15, 16
    // charge transfer, 17
    
  
  // H3+ + CI -> CHX + H2
  term = 2.0e-9;
  RatRaw[iR031] = term * narate;

  // H3+ + OI -> OHx + H2
  term = 8.0e-10;
  RatRaw[iR032] = term * narate;
  
  // H3+ + CO -> HCO+ + H2
  term = 1.7e-9;
  RatRaw[iR033] = 1.7e-9 * narate;
  
  // He+ + CO -> C+ + O + He
  term = 1.6e-9;
  RatRaw[iR034] = term * narate;
  
  // C+ + H2 -> CHX + H
  term = 4.0e-16;
  RatRaw[iR035] = term * narate;
  
  // C+ + OHX -> HCO+
  term = 1.0e-9;
  // hack munan's rate:
  //term = 9.15e-10 * (0.62 + 45.4 / pow(Temp,0.5));
  RatRaw[iR036] = term * narate;
  
  // OI + CHX -> CO + H 
  term = 2.0e-10;
  RatRaw[iR037] = term * narate;
  
  // CI + OHX -> CO + H
  term = 5.8e-12 * pow(Temp,0.5);
  // hack munan's rate:
  //term = 7.95e-10 * pow(Temp,-0.339) * exp(0.108 / Temp);
  RatRaw[iR038] = term * narate;

  // C+ + e- -> CI + gamma
  term = 1.4e-10 * pow(Temp,-0.61);
  RatRaw[iR039] = term * narate;
  
  // HCO+ + e- -> CO + H
  //term = 3.3e-5 / Temp;
  //More recent rate from Glover et al.
  term = 2.76e-7 * pow(Temp/300.0,-0.64);
  RatRaw[iR040] = term * narate;
  
  // M+ + e- -> M + gamma
  term = 3.8e-10 * pow(Temp,-0.65);
  RatRaw[iR041] = term * narate;
  
  // H3+ + M -> M+ + e + H2
  term = 2.0e-9;
  RatRaw[iR042] = term * narate;


 
  // H2+ + H2 -> H3+ + H, g54
  term = 2.24e-9 * pow(t3,0.042) * exp(-Temp / 46600.0);
  RatRaw[iR043] = term * narate;
  
  // H3+ + H -> H2+ + H2, g55
  term = 7.7e-9 * exp(-17560.0 / Temp);
  RatRaw[iR044] = term * narate;
  
  



  // munan's rate additions:
  
  // OHX + O -> O + O + H
  term = 3.5e-11;
  RatRaw[iR045] = term * narate;



  // Photo reactions
  double kb = 1.38e-16;
  double fdust, Av, fsh, fh2, fco, x, b5, b;
  double G0, NH, NH2, NCO, zeta;
  double Lshield;
  //fac = pi*kb / (mh * G)
  double ljfac = 3.892e15;
  
  
  G0 = rpar[g0_par];
  NH = rpar[NH_par];
  NH2 = rpar[NH2_par];
  NCO = rpar[NCO_par];
  zeta = rpar[zeta_par];
  
  // Jeans length shielding with temp cap:
  Lshield = pow(ljfac * min(Temp,40.0) / Den , 0.5);
 
  NH = Lshield * naterm;
  NH2 = Lshield * ys[iH2] * naterm;
  NCO = Lshield * ys[iCO] * naterm;
  
  Av = NH / 1.87e21;
  
    
    //NL99 photo rates
    // will only use
    // CHX + gamma -> CI + H
    // OHX + gamma -> OI + H
    // M + gamma -> M+ + e
  // HCO+ + gamma ->
  
  // CHX + gamma -> C + H
  fdust = exp(-1.5 * Av);
  RatRaw[iP001] = 1.0e-9 * G0 * fdust;
  
  // OHX + gamma -> O + H
  fdust = exp(-1.7 * Av);
  RatRaw[iP002] = 5.0e-10 * G0 * fdust;


  // M + gamma -> M+ + e
  fdust = exp(-1.9 * Av);
  RatRaw[iP003] = 2.0e-10 * G0 * fdust;
  
  // HCO+ + gamma -> CO + H
  fdust = exp(-2.5 * Av);
  RatRaw[iP004] = 1.5e-10 * G0 * fdust;

    // Glover photochemical rates
    // 166, 167, 168, 169, 170, 171, 198
  
  // H- + gamma -> H + e
  fdust = exp(-0.5 * Av);
  RatRaw[iP005] = 7.1e-7 * G0 * fdust;


  // H2+ + gamma -> H + H+
  fdust = exp(-1.9 * Av);
  RatRaw[iP006] = 1.1e-9 * G0 * fdust;

  // H2 + gamma -> H + H
  fdust = exp(-3.74 * Av);
  x = NH2 / 5.0e14;
  b = pow(kb * Temp / mh,0.5);
  b5 = b / 1.0e5;
  fsh = 0.965 / pow(1.0 + x / b5,2) + 0.035 / pow(1.0 + x,0.5) * exp(-8.5e-4 * pow(1.0 + x,0.5));
  RatRaw[iP007] = 5.6e-11 * fdust * fsh * G0;
  
  // H3+ + gamma -> H2 + H+
  fdust = exp(-1.8 * Av);
  RatRaw[iP008] = 4.9e-13 * G0 * fdust;
  
  // H3+ + gamma -> H2+ + H
  fdust = exp(-2.3 * Av);
  RatRaw[iP009] = 4.9e-13 * G0 * fdust;
  
  // C + gamma -> C+ + e-
  fdust = exp(-3.0 * Av);
  RatRaw[iP010] = 3.1e-10 * G0 * fdust;

  //CO + gamma -> C + O
  fdust = exp(-2.5 * Av);
  co_shielding_(&NH2, &NCO, &fh2, &fco);
  RatRaw[iP011] = 2.0e-10 * fco * fdust * fh2 * G0;

  

  // and CR rates
    // 199, 200, (H2 + cr from NL99)
    // 205, 
    // No induced cr ray rates from glover


  // H + cr -> H + + e-, g
  RatRaw[iP012] = zeta;
  
  // He + cr -> He+ + e-, g
  RatRaw[iP013] = 1.1 * zeta;
  
  // C + cr -> C+ + e-, g
  RatRaw[iP014] = 3.8 * zeta;

  // this reaction from NL99 makes no sense. 
  // H2 + cr -> H3+ + e + H, nl99
  //RatRaw[iP015] = 1.77 * zeta;

  // H2 + cr -> H2+ + e-
  RatRaw[iP015] = 2.0 * zeta;
    
  // H2 + cr -> H + H 
  RatRaw[iP016] = 0.22 * zeta;


#endif



#if DO_CHEMISTRY == SIMPLEMETAL



#endif



  if(Temp != Temp) {
    cout << "nan temp in GetNetworkRates!" << endl;
    cout << "ratraw. temp = " << Temp << endl;
    for(i=0;i<NREACTION;i++)  cout << i << " " << RatRaw[i] << endl;
    cout << "ratraw end" <<endl;
    AbortChem();
  }

};
