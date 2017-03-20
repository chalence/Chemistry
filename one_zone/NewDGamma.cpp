/* -----------------------------------------------
 Determine Gamma Routine. Aaron Lee, 2015
 
 Optimized routine to determine either:
 (1) The adiabatic index and energy density of a cell
 (2) The adiabatic index and temperature of a cell
 
 The first procedure inputs the cell's composition 
 (mass fractions) from the tracers, density, and
 temperature of the cell. This allows us to analytically
 determine the adiabatic index and the thermal energy
 density. This is fast. 
 
 The second procedure inputs the cell's composition
 (mass fractions) from the tracers, density, and
 thermal energy density of the cell. The relation 
 connecting e to T and gamma is not a simple algebraic
 formula (a result of the adibatic index for molecular
 hydrogen also being a function of temperature). Therefore
 this procedure performs a Newton-Raphson routine to 
 simultaneously determine the temperature and gamma of
 the cell for the given value of thermal energy. This is
 not as fast as (1), but it has been optimized. 
 
 Both of these routines should return updated temperatures,
 energy densities, and gamma values (though only two of
 these quantities should change for a given procedure).
 
 In this routine, the adiabatic index for molecular hydrogen
 is a function of temperature, found by summing 800 terms in 
 its partition function and evluating the partition function
 as an array of 200 elements log-spaced from 10 K to 50,000 K
 
 ------------------------------------------------- */


#include "NewDGamma.H"

// Tolerance for Newton-Raphson routine
// NR iterations are performed on temperature until
// |T_old - T_new| < NRtol
// Testing shows this occurs in <4 iterations for NRtol = 1e-6
# define NRtol 0.000001
# define NRmax 100

// Some compile declarations
# define Navo 6.022140857e23 // Avogadro's Number
# define kBol 1.38064852e-16 // Boltz constant (cgs)
# define mpro 1.6726219e-24  // Mass of proton (cgs)
# define MRHE2H 0.24         // Mass ratio of H to He for contemporary environments


/* 
  The driver for determining gamma and either temperature or energy density
  Inputs:
  which_func = 1 or 2 (actually checks if which_func = 1, else chooses other branch)
            1 = Pass in tempearture, get out energy density and gamma
            2 = Pass in energy density, get out temperature and gamma
  rho =  mass density (per volume)
  TorEDEN = either thermal energy density (per volume) [could be a pointer to the cell's value]
            or temperature. You should pass in temperature here if which_func = 1, else pass in
            energy density.
  newgamma = adiabatic index of cell [could be a pointer to the cell's value]
  H2frac = mass fraction of H_2 gas [pass in from tracer field or pass 0 if not doing chemistry]
  HEfrac = mass fraction of Helium gas [pass in from tracer field or pass 0 if not doing chemistry]
      NOTE: If not doing chemistry, you can set the mass ratio of H to He above.
 
*/
double determineGamma(int which_func, double rho, double TorEDEN, double newgamma, double H2frac, double HEfrac)
{
    if(H2frac < 0)
    {
        //cout << "Inputs wrong, mass fractions < zero!! (values H2 = " << H2frac << ")" << endl;
        H2frac = 0;
    }
    if(HEfrac < 0)
    {
        //cout << "Inputs wrong, mass fractions < zero!! (values HE = " << HEfrac << ")" << endl;
        HEfrac = 0;
    }
    
    // Prepare mass fractions
    double X[3] = {0,H2frac,HEfrac};
    X[0] = 1.0 - H2frac - HEfrac;
    
#if DO_CHEMISTRY == NO
    X[1] = 1.0/(1.0+MRHE2H);
    X[2] = 0.0;
    X[3] = MRHE2H/(1.0+MRHE2H);
#endif
    
    // Checks, H fraction should never be less than zero!
    if(X[0]<0)
    {
        //cout << "H mass fraction less than zero!! (value = " << X[0] << ")" << endl;
        // error message goes here
        X[0] = 0;
    }
    
    
    if(which_func==1) // TorEDEN = Temperature
    {
        double Temp = TorEDEN;
        double e;
        
        calcGamma(which_func,Temp,e,newgamma,X); // newgamma and e have been updated
        
        // e is energy/g. Convert to eng/vol
        e = e*rho;
        
    }
    else // TorEDEN = energy volume density
    {
        
        // Convert to eng/g
        double e = TorEDEN/rho;
        
        // Composition mu
        double mu = X[0]/mpro + X[1]/(2.*mpro) + X[2]/(4.*mpro); // actually 1/mu
        mu = 1/mu;
        
        // Calculations, using old passed-in gamma, estimate Temperature for first guess
        double Temp = e*(newgamma-1.0)*mu/(Navo*kBol);
        
        calcGamma(which_func,Temp,e,newgamma,X); // newgamma and Temp have been updated
        
    }
    
    return newgamma;

}

/*
 Where all the calculations take place. 
 The energy / gram is a function of T and gamma
    e = 1/(gamma-1) * N_avg * k_Bol * T / mu
 where mu is the mean weight (function of composition). N_avo is Avogadro's num
 and k_Bol is Boltzmann's constant. 
 We have 1/(gamma-1) = mu * Sum ( X[i]/A[i] / (gamma[i]-1) )
    for species i, where A[i] is the mass of the specie
 Only H_2 has a gamma that depends on temperature, so 1/(gamma-1) can be written
 as a constant (the H and He terms) plus a term that is a function of temperature.
 Define Delta = 1/(gamma_H2 - 1). You can then define constants 0 and 1 so that
    e = [ Cst[0] + Cst[1]*Delta(T) ] * T
 For the part the involves a Newton-Raphson routine, the function we use is
    f(T) = e - [ Cst[0] + Cst[1]*Delta(T) ] * T
 Then the Newton-Raphson routine is 
    T_new = T_old - f(T_old)/f'(T_old)
 Where, by the chain rule
    f'(T) = -[Cst[0] + Cst[1]*Delta(T)] - Cst[1]*T*Delta'(T)
          = -[Cst[0] + Cst[1]*( Delta(T) + T*Delta'(T) )]
 (Note the overall minus sign if you think there's a typo below...)
 Again this routine alters the input values of Temp, e, and gamma
 */
void calcGamma(int which, double& Temp, double& e, double& newgamma, double X[])
{
    // Evaluate the constants here (they are functions of composition)
    double Cst[2];
    Cst[0] = Navo*kBol*( X[0]/mpro/(5./3.-1.) +  X[2]/4./mpro/(5./3.-1) );
    Cst[1] = X[1]*Navo*kBol/2./mpro;
    
    // Composition mu
    double mu = X[0]/mpro + X[1]/(2.*mpro) + X[2]/(4.*mpro); // actually 1/mu
    mu = 1.0/mu;
    
    if(which==1) // We know temperature, and want e and gamma. This one is easy
    {
        // New gamma
        double idx = getH2idx(Temp); // kept as double so we can linearly interpolate in a second.
        double curDel = LinInt(idx,DelGammaTable); // linear interpolation from table (in .H file)
        double H2gam = (curDel+1.0)/curDel; // gamma value for H_2
        newgamma = (mu/mpro)*( X[0]/1.0/(5./3.-1.) + X[1]/2.0/(H2gam-1.0) + X[2]/4.0/(5./3.-1.) ); // actually this is 1/(gamma-1)
        newgamma = (newgamma+1.0)/newgamma; // And this is gamma..
        
        // New e
        e = Navo*kBol*Temp/mu/(newgamma-1.0);
    }
    else if(which==2) // We know energy density, want T and gamma. More work needed
    {
        // We make an initial guess on T (or, initial guess on gamma_{H_2})
        // Either we could use the passed-in value or just make a guess
        // Rate of convergence is not a strong function of this initial guess
        double Told = Temp; //   e / (Cst[0] + 1.7*Cst[1]);
        double Tnew, curDel;
        
        int NRcur = 0;
        while(true)
        {
            double idx = getH2idx(Told); // kept as double so we can linearly interpolate in a second.
            curDel = LinInt(idx,DelGammaTable); // linear interpolation from table (in .H file)
            double curDelp = LinInt(idx,DelGammaPrimeTable); // linear interpolation from table (in .H file)
            
            // Newton-Raphson step. If you think there's a typo and the plus should be a minus,
            // remember that f'(T_old) had an overall minus sign (see above)
            Tnew = Told +
                ( e - (Cst[0]+Cst[1]*curDel)*Told )/(Cst[0]+Cst[1]*(curDel+curDelp*Told));
            //cout << "(" << NRcur << ") Told = " << Told << " and Tnew = " << Tnew << endl;
            
            // Check if we've converged enough
            if( fabs(Tnew-Told) <= NRtol ) break;
            
            // Else we try again
            Told = Tnew;
            NRcur++;
            
            if(NRcur==NRmax) // We aren't converging...
            {
                // Some error message should be put here.
                //cout << "We are not converging!!" << endl;
                break;
            }

        }
        
        // Out of the loop
        
        // Final temperature (technically using Tnew here but using the old value
        // of curDel, but since we're out of the Newton-Raphson routine, we're
        // convinced that these values are more or less the same)
        Temp = Tnew;
        
        // New total gamma
        double H2gam = (curDel+1.0)/curDel;
        newgamma = (mu/mpro)*( X[0]/1.0/(5./3.-1.) + X[1]/2.0/(H2gam-1.0) + X[2]/4.0/(5./3.-1.) ); // actually this is 1/(gamma-1)
        newgamma = (newgamma+1.0)/newgamma; // And this is gamma..
        
    }
}

// Linear interpolation in the table
double LinInt(double idx, double Table[])
{
    double value = 0;

    double y1 = Table[int(idx)];
    double y2 = Table[int(idx)+1];
    
    value = y1 + (y2-y1)*(idx-floor(idx));
    // ( idx - floor(idx) <= 1 )
    
    //cout << "Table looks at " << int(idx) << " and " << int(idx)+1 << endl;
    //cout << "Interpolated " << value << " between " << y1 << " and " << y2 << endl;

    return value;
}

// Linear interpolation to get the indices used for the table lookups
double getH2idx(double T)
{
    double idx=0;
    
    // Tables are log spaced from 10 K to 50,000 K
    // 201 entries, determined from 800 terms in the partition function
    idx = 200.0/(log10(50000)-log10(10)) * (log10(T)-log10(10)) ;
    //cout << "Index is " << idx << endl;
    
    if(idx<0) idx=0;
    if(idx>=200.0) idx=199.99999999;  // Can't be 200, since we will index int(idx) and int(idx)+1
    
    return idx;
}

// Undefine compile declarations
# undef NRtol
# undef NRmax
# undef Navo
# undef kBol
# undef mpro
# undef MRHE2H