
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
using namespace std;

#include "ChemGlobals.H"

#include "ChemBurner.H"
#include "ChemNetwork.H"
#include "ChemIntegrate.H"

using namespace ChemNS;



double ysave[NINT];
double ymassB[NINT];
double ytotB[NSPECIES];
double RatRawB[NREACTION];
double ddyddx[NINT];
double advT = 0, firstDT, tempDT, endT, tmin, dtmax, ttot;
int jcounts = 0;

void Chem_Burner(double dt, double* eosIn, double* yIn, double* yOut, double* rpar)
{

    // Integration/Debug stuff
    int i;
    double smallx = 1.0e-20, bigx = 0.5, bigx2 = 1.0;
    double temp_before;
    ttot = dt; //ARS added "ttot" to save original value of dt
    dtmax = 0.2*ttot; 
    tmin  = dt*1.e-20; 


    //double ysave[NINT],ymass[NINT];
    for(i=0;i<NINT;i++)
      {
      ymassB[i] = yIn[i];
      ysave[i] = yIn[i];
    }
    //double ytot[NSPECIES];
    for(int i=0;i<NSPECIES;i++) {
      ytotB[i] = 0.0;
    }
  

    // Prepare the network rates (inputs: temp, density, array)
    //double RatRaw[NREACTION];
    GetNetworkRates(eosIn[1],eosIn[0],ymassB, RatRawB, rpar);

 
    // Get the values of derivatives for the ODE solver
    //double ddyddx[NINT];
    GetNetwork(dt,ddyddx,RatRawB,ymassB,rpar);

    // Given these, we can estimate the first timestep to try
    //double firstDT; 
    firstDT = dt;
    tempDT = 0.0;
    int isSmall = -1;
    for(i=0;i<NINT;i++)
    { 
	tempDT = fabs(max(ymassB[i],1.0e-4)/ddyddx[i]);
        if(tempDT < firstDT) if(ymassB[i] > smallx)
        {
            firstDT = tempDT;
            isSmall = i;
        }
    }
    if(firstDT > dtmax) firstDT = dtmax;

    
    int icounts = 0;//, jcounts=0;
    advT = 0.0; jcounts=0;
    while(dt>0) // Main integration loop
    {
        icounts = icounts + 1;
        
        // We have now set the time to integrate before updating photochem
        
        tmin = firstDT*1.0e-20;
        endT = firstDT;
        
	rpar[temp_par] = eosIn[1];
        ChemIntegrate(0,firstDT,tmin,endT,eosIn,ymassB,GetNetwork,rpar,RatRawB,jcounts);
	for(i=0;i<NINT;i++){
	  if(ymassB[i] != ymassB[i]){
	    cout << "ChemIntegrate produced nan!" << endl;
	    for(i=0;i<NINT;i++) cout << i << " " << ymassB[i] << endl;
	    cout << "eos[0] = " << eosIn[0] << endl;
	    cout << "eos[1] = " << eosIn[1] << endl;
	    cout << "eos[2] = " << eosIn[2] << endl;
	    AbortChem();
	  }
	}
	 
	
        // Prevent negative values
        for(i=0;i<NINT;i++) if(ymassB[i] < smallx ) ymassB[i] = smallx;

        //ARS also putting a check on the H2 and hP abundance 
        #if DO_CHEMISTRY == PRIMORDIAL
        if(ymassB[iH2] > bigx) ymassB[iH2] = bigx;
        if(ymassB[iHP] > bigx2) ymassB[iHP] = bigx2;        
        #endif

        // Cooling
	temp_before = eosIn[1];
	Cool_gas(eosIn,ymassB,RatRawB,rpar,endT); //updates e (=eosIn[2])
	fill_species(ymassB,ytotB);
	Chem_EOS(2,eosIn,ytotB,NSPECIES,6);
        rpar[temp_par] = eosIn[1];  //ARS updating temp value in rpar
	if(eosIn[1] < 0.1){
	  cout << "temp below 0.1 after cool species! not ok ... " << eosIn[1] << endl;
	  cout << "temp before cool_gas = " << temp_before << endl;
	  cout << "endT = " << endT << endl;
	  cout << "current abundances = " << endl;
	  for(int i=0;i<NINT;i++) cout << i << " " << ymassB[i] << endl;
	  AbortChem();
	}
        

    
        // Re-evaluate rates and network using new state values and species abundances
        GetNetworkRates(eosIn[1],eosIn[0],ymassB,RatRawB,rpar);
        GetNetwork(0.0,ddyddx,RatRawB,ymassB,rpar); // first entry time, not actually used
        
        // We have now successfully advanced 'endT' worth of time
        // Update time left
        advT = advT + endT;
        dt = dt - endT;
        firstDT = dt;
       
        // Re-evaluate firstDT
        for(i=0;i<NINT;i++)
        {
	    tempDT = fabs(max(ymassB[i],1.0e-4)/ddyddx[i]);
            if(tempDT < firstDT) if(ymassB[i] > smallx)
            {
                firstDT = tempDT;
                isSmall = i;
            }
        }
	if(tempDT < firstDT) 
	  {
	    firstDT = tempDT;
	    isSmall = iHP;
	  }
        if(firstDT > dtmax) firstDT = dtmax;

        //ARS adding a different temp cooling constraint
        double Tdot = 0.01*fabs(temp_before - eosIn[1])/firstDT;
        double tkDT = eosIn[1] / Tdot;
        if(tkDT < firstDT) firstDT = tkDT;

        if(firstDT > dt) firstDT = dt; // should never happen, but why not check.
        
    }

 
    // All done!
    for(i=0;i<NINT;i++)
    {
      yOut[i] = ymassB[i];
      if(yOut[i] < smallx) yOut[i] = smallx;
    }
 
  //ARS also putting a check on the H2 and HPabundance
  #if DO_CHEMISTRY == PRIMORDIAL
   if(yOut[iH2] > bigx) yOut[iH2] = bigx;
   if(yOut[iHP] > bigx2) yOut[iHP] = bigx2;
  #endif


};






// determine electron abunance from charge conservation:
// network dependent
void fill_species(double* ymass, double* ytot)
{

  for(int i=0;i<NSPECIES;i++) ytot[i]=0.0;
  for(int i=0;i<NINT;i++) ytot[i]=ymass[i];

#if DO_CHEMISTRY == CONTEMPORARY
  ytot[iELEC] = ytot[iHP] + ytot[iHEP] + ytot[iCP] + ytot[iHCOP] + ytot[iMP] - ytot[iHM] + ytot[iH2P] + ytot[iH3P];
  if(ytot[iELEC] < 1.0e-20) ytot[iELEC] = 1.0e-20;
#endif

#if DO_CHEMISTRY == PRIMORDIAL
  //ARS adding  extra lines here to make sure abundance totals add up correctly
  double DeutTot=0;

  ytot[iHE] = abhe - ytot[iHEP] - ytot[iHEPP];
  ytot[iH] = 1.- 2.e0 * ytot[iH2]  - ytot[iHP]  - ytot[iHD] - 2.*ytot[iH2P] - ytot[iHM];

  ytot[iELEC] = ytot[iHP] + ytot[iDP] + ytot[iHEP] + 2.0*ytot[iHEPP] + ytot[iHDP] + ytot[iH2P] + ytot[iD2P] - ytot[iHM] - ytot[iDM];

  for(int i=0;i<NSPECIES;i++) if(ytot[i] < 1.e-20) ytot[i] = 1.e-20;
#endif

};






void Cool_function(double ei, double* ytot, double* ratraw, double* rpar, double &cool_sum)  
{
    // What is returned:
    cool_sum = 0.0;
    
    double eosIn[6];
    eosIn[0] = rpar[rho_par];
    eosIn[2] = ei;
    // non ideal hack here to get gamma:
    eosIn[5] = rpar[gamma_par];
    Chem_EOS(2,eosIn,ytot,NSPECIES,3);
   
    double tmp = eosIn[1];
    
    double rho = rpar[rho_par];

    int i; 
    if(ei <= 0.0){
      cout << "entering Cool_function with negative ei =  " << ei << endl;
      cout << "current temp =  " << eosIn[1] << endl;
      cout << "current abundanes:  " << endl;
      for(i=0;i<NINT;i++) cout << i << " " << ytot[i] << endl;
      AbortChem();
      return;
    }
   
    double log10temp = log10(tmp);
    double sqt = pow(tmp,0.5);
    
    double T3 = tmp/1000.0;
    double lt1 = pow(log10(T3),1);
    double lt2 = pow(log10(T3),2);
    double lt3 = pow(log10(T3),3);
    double lt4 = pow(log10(T3),4);
    double lt5 = pow(log10(T3),5);
    

    double c1=0.0;
    double c2=0.0;


    double llow, lrlte, lvlte, llte;
    double t51, t52, t53, t54, t55, t56;
    double t12, t13, t14, t15, t16;

    // 4.606e23 = 1/(1.3 * mh)
    double naterm = rho * 4.606e23;
    double nin[NSPECIES];
    for(int i=0;i<NSPECIES;i++) nin[i] = ytot[i]*naterm;
    
    double lambda_h2 = 0.0;
    double lambda_co = 0.0;
    double lambda_gd = 0.0;
    double lambda_fs = 0.0;
    double lambda_atomic = 0.0;
    double gamma_grain_pe = 0.0;
    double gamma_cr = 0.0;
    
    //ARS updating netowrk rates again here, since they can change dramatically at high rho 
    GetNetworkRates(eosIn[1],eosIn[0],ytot, ratraw, rpar);

     // Go over all cooling rates for molecules

    // H-H2
    // above 10^4 K, all these collision rates take on constant values.
    
    double term=0.0;
    if(tmp < 10.0)  term = 0.0e0;
    else if(tmp >= 10.0 && tmp < 100.0)  {
        term = -16.818342 + 37.383713*lt1 + 58.145166*lt2 + 48.656103*lt3 + 20.159831*lt4 + 3.8479610*lt5;
        term = pow(10,(term));
    }
    else if(tmp >= 100.0 && tmp < 1000.0) {
        term = -24.311209 + 3.5692468*lt1 - 11.332860*lt2 - 27.850082*lt3 - 21.328264*lt4 - 4.2519023*lt5;
        term = pow(10,(term));
    }
    else if(tmp >= 1000.0 && tmp < 10000.0)  {
        term = -24.311209 + 4.6450521*lt1 - 3.7209846*lt2 + 5.9369081*lt3 - 5.5108047*lt4 + 1.5538288*lt5;
        term = pow(10,(term));
    }
    else {   
      term = -24.311209 + 4.6450521 - 3.7209846 + 5.9369081 - 5.5108047 + 1.5538288;
      term = pow(10,(term)); 
    }
        
    t51 = term;
 
    

  // H2 & H2
    term = 0.0e0;
    if(tmp < 10.0)  term = 0.0e0;
    else if(tmp >= 10.0 && tmp < 10000.0)  {
        term = -23.962112 + 2.09433740*lt1 - 0.77151436*lt2 + 0.43693353*lt3 - 0.14813216*lt4 - 0.033638326*lt5;
        term = pow(10,(term));
    } 
    else {
      term = -23.962112 + 2.09433740 - 0.77151436 + 0.43693353 - 0.14813216 - 0.033638326;
      term = pow(10,(term));
    }
    t52 = term;

    // H2 & He
    term = 0.0e0;
    if(tmp < 10.0)  term = 0.0e0;
    else if(tmp >= 10.0 && tmp < 10000.0)  {
        term = -23.689237 + 2.1892372*lt1 - 0.81520438*lt2 + 0.29036281*lt3  - 0.16596184*lt4 + 0.19191375*lt5;
        term = pow(10,(term));
    }
    else {
      term = -23.689237 + 2.1892372 - 0.81520438 + 0.29036281 - 0.16596184 + 0.19191375;
      term = pow(10,(term));
    }
    t53 = term;

    // H2 & HP
    term = 0.0e0;
    if(tmp < 10.0)  term = 0.0e0;
    else if(tmp >= 10.0 && tmp < 10000.0) {
        term = -21.716699 + 1.3865783*lt1 - 0.37915285*lt2 + 0.11453688*lt3 - 0.23214154*lt4 + 0.058538864*lt5;
        term = pow(10,(term)) ;
    }
    else {
      term = -21.716699 + 1.3865783 - 0.37915285 + 0.11453688 - 0.23214154 + 0.058538864;
      term = pow(10,(term)) ;
    }
    t54 = term;


    // H2 & e
    term = 0.0e0;
    if(tmp < 10.0)  term = 0.0e0;
    else if(tmp >= 10.0 && tmp < 200.0)  {
        term = -34.286155 - 48.537163*lt1 - 77.121176*lt2 - 51.352459*lt3 - 15.169160*lt4 - 0.98120322*lt5;
        term = pow(10,(term));
    }
    else if(tmp >= 200.0 && tmp < 10000.0)  {
        term = -22.190316 + 1.5728955*lt1 - 0.21335100*lt2 + 0.96149759*lt3 - 0.91023195*lt4 + 0.13749749*lt5;
        term = pow(10,(term));
    } 
    else {
      term = -22.190316 + 1.5728955 - 0.21335100 + 0.96149759 - 0.91023195 + 0.13749749;
      term = pow(10,(term));
    }
    t55 = term;

    
    if(tmp > 10.0){
      // factor_1
      t12 = 9.5e-22 * pow(T3,3.76);
      
      // factor_2 (dens dependance)
      t13 = exp(- pow(0.13 / T3,3));
      
      // factor_3
      t14 = (1.0 + 0.12 * pow(T3,2.1));
      
      // factor_4 (dens dependance)
      t15 = 3.0e-24 * exp(- 0.51 / T3);
      
      // LVLTE 
      t16 = 6.7e-19 * exp(-5.86 / T3) + 1.6e-18 * exp(-11.7 / T3);
      
      
      
      // H2 cooling rate based on Glover & Abel 2008:
      // LTE cooling rate still from Galli & Palla, which should still be good
      
      llow = t51 * nin[iH] + t52 * nin[iH2] + t53 * nin[iHE] + t54 * nin[iHP] + t55 * nin[iELEC];
      lrlte = t12 * t13 / t14 + t15;
      lvlte = t16;
      llte = lrlte + lvlte;
      lambda_h2 = llte * nin[iH2] / (1.0 + llte / llow);
    }
  
 
#if DO_CHEMISTRY == CONTEMPORARY
    
    // CO cooling
    double yco, yh2, ye, yc, ycp, yo, yh, yhp, yhep, yhe;
    yco = ytot[iCO];
    yh2 = ytot[iH2];
    ye = ytot[iELEC];
    yc = ytot[iC];
    ycp = ytot[iCP];
    yo = ytot[iO];
    yh = ytot[iH];
    yhp = ytot[iHP];
    yhe = ytot[iHE];
    yhep = ytot[iHEP];


    // using the alternative velocity divergence here, divv1
    double gradrho = rpar[gradrho_par];
    double vturb = rpar[vturb_par];
    double divv = rpar[divv1_par];
    co_cooling_(&yco,&yh2,&ye,&tmp,&naterm,&divv,&vturb,&gradrho,&lambda_co);
 

    
    // metal line fine-structure cooling - contains CMB temperature floor 
    // hardwired at TCMB = 2.725K
    double NH = rpar[NH_par];
    if(tmp > 1.0){
    fine_structure_cooling_(&yc, &ycp, &yo, &yh, &ye, &yhp, &yh2, &tmp, &naterm, &NH, &lambda_fs);
    }
 
    // atomic cooling
    double t1, t9;

    double temp_var = 5.5 - log10temp;
    double aux_var = exp(- pow(temp_var,2) / 3.0);
    double gff = 1.1 + 0.34 * aux_var;
 
    // free free 
    t1 =  1.42e-27 * gff * sqt;
    
    // Collisional excitation cooling of HI (ly-alpha)
    t9 = 0.0;
    if (tmp > 900.0){
      double XH2 = 1.0 / (1.0 + pow(tmp / 1.0e5,0.5));
      t9 = 7.5e-19 * XH2 * exp(-118348.0 / tmp);
    }
       
    lambda_atomic = t1 * (yhp + yhep) * ye 
			+ t9 * yh * ye;
    lambda_atomic = lambda_atomic * naterm * naterm;



    // dust-gas energy transfer
    // lambda_gd = positive for gas cooling (Tdust  < Tgas)
    //           = negative for gas heating (Tdust > Tgas)
    // if radiation it turned on, this term should be zero
    double alpha_gd = 3.2e-34; 
    double Tdust = rpar[tdust_par];
    lambda_gd = alpha_gd * pow(tmp,0.5) * (tmp - Tdust) * naterm * naterm;




    double Lshield;
    //fac = pi*kb / (mh * G)
    double ljfac = 3.892e15;
    
    // Jeans length shielding with temp cap:
    Lshield = pow(ljfac * min(tmp,40.0) / rho , 0.5);


    NH = Lshield * naterm;
    
    // dust grain photoelectric heating
    double G0 = rpar[g0_par];
    double G0_eff = G0 * exp(-0.5 * NH * 1.0e-21);
    double beta = 0.735 / pow(tmp,0.068);
    double fac = G0_eff * pow(tmp,0.5) / nin[iELEC];
    double eff = 4.9e-2 / (1.0 + 4.0e-3 * pow(fac,0.73))  
      + 3.7e-2 * pow(tmp / 1.0e4,0.7) / (1.0 + 2.0e-4 * fac);
    gamma_grain_pe = 1.0e-24 * eff * G0_eff * naterm
      - 4.65e-30 * pow(tmp,0.94) * pow(fac,beta) * nin[iELEC] * naterm;
    // simple photoelectric heating 
    //gamma_grain_pe = 4.0e-26 * G0_eff * naterm;
  
    
    // cosmic ray ionization heating rate:
    // see Krumholz 2014
    double zeta = rpar[zeta_par];
    double qion, qion_H2, qion_HI;
    double lognh = log10(naterm);
    double ev_in_erg = 1.602e-12;
    qion_HI = 1.04e-11 + 4.23e-11 * pow(ye / (ye + 0.07),0.5);
 
    if(lognh < 2.0) qion_H2 = 10.0 * ev_in_erg;
    else if(lognh >= 2.0 && lognh < 4.0) qion_H2 = (10.0 + 3.0 * (lognh - 2.0) / 2.0) * ev_in_erg;
    else if(lognh >= 4.0 && lognh < 7.0) qion_H2 = (13.0 + 4.0 * (lognh - 4.0) / 3.0) * ev_in_erg;
    else if(lognh >= 7.0 && lognh < 10.0) qion_H2 = (17.0 + (lognh - 7.0) / 3.0) * ev_in_erg;
    else qion_H2 = 18.0 * ev_in_erg;

    qion = yh * qion_HI + 2.0 * yh2 * qion_H2;
    gamma_cr = zeta * qion * naterm;
    

#endif
    
  
   
    cool_sum = lambda_h2 + lambda_gd + lambda_co + lambda_fs + lambda_atomic - gamma_grain_pe - gamma_cr;
   

    // simple optically thin cooling from Kim & Ostriker (for testing purposes):                            
    //cool_sum = (2.0e-19 * exp(-1.184e5 / (tmp + 1000.0))                                            
    //          + 2.8e-28 * pow(tmp,0.5) * exp(-92.0 / tmp))*naterm*naterm                             
    //          - 2.0e-26 * 0.2 * naterm;  


    // Need to make cooling negative and heating positive. And convert volumetric -> specific:
    cool_sum = -cool_sum / rho;

    if(cool_sum != cool_sum) {
      cout << "nan cool sum!" << endl;
      cout << "tmp = " << tmp << endl;
      cout << "rho = " << rho << endl;
      cout << "lambda h2 = " << lambda_h2 << endl;
      cout << "lambda gd = " << lambda_gd << endl;
      cout << "lambda co = " << lambda_co << endl;
      cout << "lambda fs = " << lambda_fs << endl;
      cout << "lambda atomic = " << lambda_atomic << endl;
      cout << "gamma pe = " << gamma_grain_pe << endl;
      cout << "gamma cr = " << gamma_cr << endl;
      
      cout << "***" << endl;
      cout << "ytot = " << endl;
      for(int i =0;i<NSPECIES;i++) cout << i << " " << ytot[i] << endl;
      AbortChem();
    }
    

};



void Cool_gas(double* eosIn, double* ymass, double* ratraw, double* rpar, double dt)
{
    int i;
    
    double tradmax = 1.0e9;
    double tradmin = 1.0e1;
    double dradmax = 1.0e-10;
    double dradmin = 1.0e-30;
    
    double rho = eosIn[0];
    double Temp = eosIn[1];
    double ei = eosIn[2];
    double mfrac = eosIn[4];
    
    
    double ytot[NSPECIES];
    fill_species(ymass,ytot);

    double dtHold = dt;
    double sdot;
    double dtmax;
    double escale;
    double eold;
    double sdotold;


    int nstep = 0;

    // integrate internal energy forward dt
    Cool_function(ei,ytot,ratraw,rpar,sdot);
    
    int jcounts = 0;
    double firstDt, endT, tempDt;
    double hdid, hnext, h;
    double x=0;
    int stepMax = 100000;
    double odescal = 1.0e-6;
    double yscal, ei0, eis;
    
    escale = ei;
    ei0 = ei;
    firstDt = dtHold;
    tempDt = 0.1 * ei / fabs(sdot);
    if(tempDt < firstDt) firstDt = tempDt;
    h = firstDt;
    eis = ei / escale;
    while(dtHold > 0.0){
      
    yscal = max(odescal,fabs(eis));
    hdid=0;
    hnext=0;
    eold = eis*escale;
    sdotold = sdot;
    for(int n=0;n<10;n++){
      StepMeSingle(eis,sdot,x,h,hdid,hnext,yscal,Cool_function,jcounts,ytot,ratraw,rpar,escale);

      if(eis>0.0) break;
      h=h/2.0;
      hdid=0.0;
      hnext=0.0;
      eis = eold/escale;
      sdot = sdotold;
      

      if(n > 5){
	cout << "more than 5 iterations of step me single. probably not right..." << " " << eis << endl;
	AbortChem();
      }
      
    }
    if(hdid > h){
      cout << "hdid greater than h in step me single! not OK: " << hdid << " " << h << endl;
      AbortChem();
    }

    
    if(eis <= 0.0 || hdid == 0.0) {
      cout << "negative ei produced by stepmesingle! or hdid=0, " << eis*escale << endl;
      cout << "input sdot = " << sdot << endl;
      cout << "eold = " << eold << endl;
      cout << "escale = " << escale << endl;
      cout << "jcounts = " << jcounts << endl;
      cout << "attempted timestep = " << h << endl;
      cout << "acahieved timestep = " << hdid << endl;
      cout << "temp = " << Temp << endl;
      cout << "current ytot = " << endl;
      for(int i =0;i<NSPECIES;i++) cout << i << " " << ytot[i] << endl;
      AbortChem();
    }

    
    dtHold = dtHold - hdid;
    if(dtHold < 0.0){
      cout << "dthold < 0.0 should never be here: " << dtHold << " " << hdid << " " << h << endl;
      cout << "original dt = " << dt << endl;
      cout << "x = " << x << endl;
      cout << "rpar density " << rpar[rho_par] << endl;
      AbortChem();
    }
    
    h = min(dtHold,hnext);
    Cool_function(eis*escale,ytot,ratraw,rpar,sdot);

  }

    eosIn[2] = eis*escale;// * escale;
    Chem_EOS(2,eosIn,ytot,NSPECIES,4);

    if(eosIn[1] < 0.1) {
      cout << "Temp is really low after cool_gas... temp = " << eosIn[1] << endl;
      eosIn[2] = ei0;// * escale;
      Chem_EOS(2,eosIn,ytot,NSPECIES,5);
      cout << "ei0 and Temp0 = " << ei0 << " " << eosIn[1] << endl; 

      cout << "current ytot = " << endl;
      for(int i =0;i<NSPECIES;i++) cout << i << " " << ytot[i] << endl;

      AbortChem();
    }
   
};


