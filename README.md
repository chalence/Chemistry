A chemistry (primordial and contemporary) routine for Orion2. Based off the package used in Gray and Scannapieco, 2010, Apj, 718, 417. Further updated by Chalence Safranek-Shrader for use with both primordial and contemporary chemistry, using the methodology of Nelson & Langer (1999). Primordial chemistry tracks the evolution of 16 species (including electrons), in particular molecular hydrogen and hydrogen deuteride. This scheme employs a robust 4th-order Runge Kutta scheme to evolve the species in time and allows for sub-cycling within both the chemistry and cooling update to ensure numerical stability. The species and network rates are taken from Glover and Abel 2008, MNRAS, 388, 1627, as well as Gray and Scannapieco (2010). 

Improvements for the Orion2 version include updated cooling tables from CLOUDY, a linear interpolation function when reading data from these tables, as well as general debugging and refactoring of the original Fortran code. Written in C++. The package is called for each cell. Each cell’s chemistry is independent of every other cell, giving excellent parallel scaling. 

Without ionization ray tracing, this package is called on the hydro time step. With ionization ray tracing particular to Pop III stars, this package is called after every ray trace (implementation still unfinished for use with ray-tracing). If there is no ray tracing, chemistry is called in AMRLevelOrion.cpp via the function ChemPrep(). 

As of now, the chemistry package requires 6 inputs, via ChemDriver:
- dt_in, the amount of time to evolve, 
- den_in, the mass density
- ei_in, the specific energy density (ei_in = erg/g)
- gamma_in, the value of the adiabatic index (gamma, determined external to chemistry, usually equal to 5/3),
- y_fromOrion, array of size NINT representing the species abundances (defined as the number density of that species divided by the number density of hydrogen nuclei
- rpar, array of size NRPAR containing additional user input needed to be passed around to the various chemistry functions (mainly used with contemporary chemistry, though easily extendable), labels of which can be found in ChemExternal.H

The tracer fields are used to track the abundances of each chemical species. For primordial, The ordering of the chemical species can be found in ChemGlobals.H. Changing the order has not been tested, so don't try it. The abundances should be passed into the chemistry routine in the exact same order as what is given in the ChemGlobals.H file, which is originally determined in the init.c file for your particular problem. The electron fraction is set by imposing charge conservation, and is not required to be stored as a tracer field. The chemistry package returns updated abudances and an updated specific energy. 


=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

DEFINITIONS.H flags: 
 - DO_CHEMISTRY  (default = NO) :: CONTEMPORARY = contemporary (CO) chemistry
   		 	    	   PRIMORDIAL = primordial, metal-free chemistry
 - CHEM_NSPECIES (default = 0) :: The number of species evolved and tracked by orion, manually set by the user depending on which network is used (should be changed to be automatic...) 
   		 	    	   PRIMORDIAL   = CHEM_NSPECIES = NINT = 15
  		 	    	   CONTEMPORARY = CHEM_NSPECIES = NINT = 17

Note: for both networks, the number of species integrated (NINT) is one less than the total number of considered species, with the one remaining species in both networks being electrons,whose abundance can always be calculated by enforcing charge conservation.

 The number of tracer fields NTRACER needs to be at least as big as CHEM_NSPECIES. At the moment, this package assumes the chemical species occupy the first CHEM_NSPECIES tracer field slots. 


=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Inside the chemistry package, there are a few knobs that can be set:

ChemDriver.cpp:
	double smallx : Floor chemical species abundance. Default = 1e-20.

ChemIntegrate.cpp — function ChemIntegrate()
	REAL odescal : When computing the errors in the Runge-Kutta routine, this value sets the lower limit for which abundances are considered. Default = 1e-6.
	
ChemIntegrate.H — function StepMe()
	REAL eps : Runge-Kutta error tolerance. Default = 1e-4


=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Brief description of each file and its contents
(this should really be cleaned up and better organized...)):

ChemGlobals.H     — Defines a namespace used for species and reaction indexing 
ChemDriver.H      — Holds EOS routine
ChemDriver.cpp    — Driver function that sets everything up and then calls the Burner routine
ChemBurner.H      — The real "meat" of chemistry. Contains reaction networks and cooling/heating processes. Also handles temperature/energy integration.
ChemBurner.cpp    — The main loop that re-evaluates abundance derivatives, calls the Runge-Kutta integrator, and calls the cooling routines
ChemIntegrate.H   — Includes some small matrix operation functions and StepMe, the function that advances by one Runge-Kutta step. 
ChemIntegrate.cpp — Houses the integration setup routine, which calls StepMe to evolve the abundances forward in time. 
Jacobian.cpp      — Calculates the jacobian matrix for the Runge-Kutta integration 


