#include "ChemGlobals.H"
#include "ChemIntegrate.H"
#include "ChemDriver.H"
#include <iostream>

using namespace std;

using namespace ChemNS;

void Jacobian(REAL x,REAL* yt,REAL (*dfdy)[NINT], REAL* ratraw, double* rpar);

void ChemIntegrate(REAL start,REAL firstDT,REAL tmin,REAL dt,REAL* eosIn, REAL* ymass, void (*derivs)(REAL,REAL*,REAL*,REAL*, double* ), double* rpar, REAL* RatRaw, int &jcounts)
{
    int i;
    
    //REAL dxsav = 1;
    //int kmax = 1;
    REAL odescal = 1.0e-6;
    int stepMax = 1000000, stepCur = 0;
    
    
    int nok = 0;
    int nbad = 0;
    
    // Store boundary conditions
    REAL y[NINT];
    REAL yscal[NINT];
    for(i=0;i<NINT;i++) y[i] = ymass[i];
    
    REAL dydx[NINT];
    
    REAL x = start;
    REAL h = firstDT, hdid = 0, hnext=0;
    
    for(stepCur=1;stepCur < stepMax ; stepCur++)
    {
     
        
        // Yet another check that everything is positive
        for(i=0;i<NINT;i++) y[i] = max(y[i],1.0e-20);
        
        // Update the derivatives
        //cout << "About to enter loop derivs" << endl;
        derivs(h,dydx,RatRaw,y,rpar);
        
        // Used in error checking 
        for(i=0;i<NINT;i++) yscal[i] = max(odescal,fabs(y[i]));
        
        // Check to make sure we're not overshooting the timestep
        if( (x+h-dt)*(x+h-start) > 0.0 ) h = dt-x;
        
        hdid = 0; hnext = 0;
        
        // Step away! (y,x,hdid,hnext,jcounts updated)
        StepMe(y,dydx,x,h,hdid,hnext,yscal,derivs,jcounts,RatRaw,rpar);
        
            
        // Tally up
        if( hdid == h ) nok = nok+1;
        else nbad = nbad + 1;
        
        
        // Checks
        if(fabs(h) < tmin)
        {
            // print an error message
            cout << "ERROR: with fabs(h) < tmin" << endl;
        }
        
        if(stepCur==stepMax-1)
        {
            // Exceed max number of steps...
            cout << "ERROR: Exceeded max number of steps" << endl;
        }
        
        // Are we done?
        if( stepCur==stepMax-1 || (x-dt)*(dt-start) >= 0)
        {
            for(i=0;i<NINT;i++) ymass[i] = y[i];
            break;
        }
        
        
        // If not, get ready to start over
        h = hnext;
        
        
    } // end of stepping loop
    
    //cout << "jcount = " << jcounts << " nok = " << nok << " and nbad = " << nbad << endl;
    return;
    
};



// Rosenbrock ODE semi-implicit method for a single variable 
// (specific internal energy in our case)
void StepMeSingle(REAL &y, REAL &dydx, REAL &x,REAL htry,REAL &hdid,REAL &hnext,
                  REAL &yscal, void (derivs)(REAL,REAL*,REAL*,REAL*,REAL&),
                  int &jcounts,REAL* ytot, REAL* ratraw, double* rpar, REAL escale)
{

    
    int i,j,jtry;
    
    int TryMax = 100;
    
    
    REAL a21 = 2.0e0, a31 = 48.0e0/25.0e0, a32 = 6.0e0/25.0e0, c21 = -8.0e0, c31 = 372.0e0/25.0e0, c32 = 12.e0/5.0e0, c41 = -112.0e0/125.0e0,   c42 = -54.0e0/125.0e0, c43 = -2.0e0/5.0e0,   c1x = 1.0e0/2.0e0, c2x = -3.0e0/2.0e0,  c3x = 121.0e0/50.0e0, c4x = 29.0e0/250.0e0, a2x = 1.0e0, a3x = 3.0e0/5.0e0;
    
    REAL eps = 1.0e-4 * escale;
    
    // Original values (needed if step fails)
    REAL xsav,ysav,dysav;
    ysav = y;
    dysav = dydx/escale;
    xsav = x;
    
    
    // Get the dense Jacobian
    REAL dfdy;
    // Jacobian(xsav,ysav,dfdy,ratraw,rpar);
    // Compute numerical jacobian of e here
    REAL yr, yl, dydxr, dydxl;
    REAL fr = 1.0 + 1.0e-2;
    REAL fl = 1.0 - 1.0e-2;
    yr = y * fr * escale;
    yl = y * fl * escale;
    derivs(yr,ytot,ratraw,rpar,dydxr);
    derivs(yl,ytot,ratraw,rpar,dydxl);
    dfdy = (dydxr - dydxl) / (yr - yl);
    ///
    
    REAL h = htry;
    
    
    // Needed for matrix solve
    REAL dmat;
    REAL av;
    REAL bb;
    
    
    REAL gamf = 0.5;
    // Solve equations loop
    for(jtry=1;jtry<TryMax;jtry++)
    {
        
        
        // Form the matrix
        REAL xx = 1.0/h/gamf;
        
	dmat = -dfdy;
	dmat = xx + dmat;
        
        
        // Set up and solve RHS for g1
        REAL g1;
        
	g1 =  dysav;
	av = dmat;
	bb = g1;
	
	g1 = bb/av;

	y = ysav + a21 * g1;
        
        // Set up and solve RHS for g2
        x = xsav + a2x * h ;
        REAL g2;
       
        derivs(y*escale,ytot,ratraw,rpar,dydx); 
	dydx = dydx / escale;

	g2 = dydx + c21 * g1 / h;
	av = dmat;
	bb = g2;
        
	g2 = bb/av;
        
	y = ysav + a31 * g1 + a32 * g2;
        
        
        // Set up and solve RHS for g3
        x = xsav + a3x*h;
        REAL g3;
        derivs(y*escale,ytot,ratraw,rpar,dydx);
	dydx = dydx / escale;

	g3 = dydx + (c31*g1 + c32*g2)/h;
	av = dmat;
	bb = g3;

	g3 = bb/av;
	
        // And finally, g4
        REAL g4;
	g4 = dydx + (c41*g1 + c42*g2 + c43*g3) / h;
	av = dmat;
	bb = g4;

	g4 = bb/av;
        
        
        // Compute 3rd and 4th order estimates for y
        REAL err;
        double b1 = 19.0e0/9.0e0, b2 = 1.0e0/2.0e0,
        b3 = 25.0e0/108.0e0, b4 = 125.0e0/108.0e0;
        double e1 = 17.0e0/54.0e0, e2 = 7.0e0/36.0e0,
        e3 = 0.0, e4 = 125.0e0/108.0e0;
        
	y = ysav + b1*g1 + b2*g2 + b3*g3 + b4*g4;
	err = e1*g1 + e2*g2 + e3*g3 + e4*g4;
        
        
        // Time advances
        x = xsav + h;
        
        if(x==xsav)
        {
            // Something has gone wrong
            cout << "ERROR: Uh oh, h=0 in temp ODE" << endl;
	    cout << "h = " << h << endl;
	    cout << "htry = " << htry << endl;
	    cout << "ysav = " << endl;
	    cout << i << " " << ysav << endl;
	    cout << "dysav = " << endl;
	    cout << i << " " << dysav << endl;
	    AbortChem();
        }
        
        
        // Determine the scaled accuracy
        REAL errmax = fabs(err/yscal) / eps;
       
        
        REAL safety = 0.9, grow= 1.5, pgrow= -0.25, shrnk= 0.5, pshrnk= -1.0/3.0;
        REAL errcon = 0.1296;
        if(errmax <= 1.0) // All went well
        {
            hdid = h;
            if(errmax > errcon) hnext = safety*h*pow(errmax,pgrow);
            else hnext = grow*h;
            jcounts = jcounts+jtry;
            
            if(jcounts>400000)
            {
                // Why is this taking so long?!
                cout << "Jcounts is too high in RS (jcounts="<<jcounts<<")" << endl;
            }
            
	   
            return;
            
        }
        else // Error too large, gots to try again
        {
	  
            hnext = safety*h*pow(errmax,pshrnk);
	    //if(y <= 0.0) hnext = shrnk*fabs(h);
	    //cout << "error too large in stem me single! " << h << " " << hnext << endl;
            double hsign = 1.0;
            if(h<0) hsign = -1.0;
            if(h<0) cout << "ERROR: hsign is < 0?! in RS" << endl; // should never happen
            h = hsign*max(fabs(hnext),shrnk*fabs(h));
            //cout << "Failed step " << htry << " " << h << " " << hnext << " " << hdid << " " << y << " " << errmax << endl;
            
            
        }
        
       
    } // End of solve equation loop

    // should never get here
    cout << "ERROR: Too many step tries in RS (jtry ="<<jtry<<")" <<endl;
    AbortChem();
};




void StepMe(REAL* y,REAL* dydx,REAL &x,REAL htry,REAL &hdid,REAL &hnext,REAL* yscal, void (*derivs)(REAL,REAL*,REAL*,REAL*,double*), int &jcounts, REAL* ratraw, double* rpar)
{

    
    int i,j,jtry;
    
    int TryMax = 30000000;
    
    
    REAL a21 = 2.0e0, a31 = 48.0e0/25.0e0, a32 = 6.0e0/25.0e0, c21 = -8.0e0, c31 = 372.0e0/25.0e0, c32 = 12.e0/5.0e0, c41 = -112.0e0/125.0e0,   c42 = -54.0e0/125.0e0, c43 = -2.0e0/5.0e0,   c1x = 1.0e0/2.0e0, c2x = -3.0e0/2.0e0,  c3x = 121.0e0/50.0e0, c4x = 29.0e0/250.0e0, a2x = 1.0e0, a3x = 3.0e0/5.0e0;
    
    REAL eps = 1.0e-4;
    
    
    // Original values (needed if step fails)
    REAL xsav,ysav[NINT], dysav[NINT];
    for(i=0;i<NINT;i++) ysav[i] = y[i];
    for(i=0;i<NINT;i++) dysav[i] = dydx[i];
    xsav = x;
    
    
    // Get the dense Jacobian
    REAL dfdy[NINT][NINT];
    Jacobian(xsav,ysav,dfdy,ratraw,rpar);
    REAL h = htry;
    
    
    // Needed for matrix solve
    REAL dmat[NINT][NINT];
    REAL av[NINT][NINT];
    REAL bb[NINT];
    
    
    REAL gamf = 0.5;
    // Solve equations loop
    for(jtry=1;jtry<TryMax;jtry++)
    {
        
        
        // Form the matrix
        REAL xx = 1.0/h/gamf;
        
        for(j=0;j<NINT;j++)
            for(i=0;i<NINT;i++)
                dmat[i][j] = -dfdy[i][j];
        for(i=0;i<NINT;i++) dmat[i][i] = xx + dmat[i][i];
        
        
        // Set up and solve RHS for g1
        REAL g1[NINT];
        for(i=0;i<NINT;i++) g1[i] = dysav[i];
        for(j=0;j<NINT;j++)
            for(i=0;i<NINT;i++)
                av[i][j] = dmat[i][j];
        for(i=0;i<NINT;i++) bb[i] = g1[i];
        
        
        leqs(av,bb);
        for(j=0;j<NINT;j++) g1[j] = bb[j];
        //g1[iELEC] = g1[iHP] + g1[iDP] + 2.0*g1[iHEPP] + g1[iHEP] + g1[iHDP] + g1[iH2P] + g1[iD2P] - g1[iHM] - g1[iDM];
        
        for(i=0;i<NINT;i++) y[i] = ysav[i] + a21 * g1[i];
        
        
        // Set up and solve RHS for g2
        x = xsav + a2x * h ;
        REAL g2[NINT];
        
        derivs(x,dydx,ratraw,y,rpar);
        for(i=0;i<NINT;i++) g2[i] = dydx[i] + c21*g1[i]/h;
        for(j=0;j<NINT;j++)
            for(i=0;i<NINT;i++)
                av[i][j]=dmat[i][j];
        for(i=0;i<NINT;i++) bb[i] = g2[i];
        leqs(av,bb);
        for(j=0;j<NINT;j++) g2[j] = bb[j];
        
        for(i=0;i<NINT;i++) y[i] = ysav[i] + a31*g1[i] + a32*g2[i];
        
        
        // Set up and solve RHS for g3
        x = xsav + a3x*h;
        REAL g3[NINT];
        derivs(x,dydx,ratraw,y,rpar);
        for(i=0;i<NINT;i++) g3[i] = dydx[i] + (c31*g1[i] + c32*g2[i])/h;
        for(j=0;j<NINT;j++)
            for(i=0;i<NINT;i++)
                av[i][j]=dmat[i][j];
        for(i=0;i<NINT;i++) bb[i] = g3[i];
        leqs(av,bb);
        for(j=0;j<NINT;j++) g3[j] = bb[j];

        // And finally, g4
        REAL g4[NINT];
        for(i=0;i<NINT;i++) g4[i]= dydx[i]+(c41*g1[i]+c42*g2[i]+c43*g3[i])/h;
        for(j=0;j<NINT;j++)
            for(i=0;i<NINT;i++)
                av[i][j]=dmat[i][j];
        for(i=0;i<NINT;i++) bb[i] = g4[i];
        leqs(av,bb);
        for(j=0;j<NINT;j++) g4[j] = bb[j];
        
        
        // Compute 3rd and 4th order estimates for y
        REAL err[NINT];
        double b1 = 19.0e0/9.0e0, b2 = 1.0e0/2.0e0,
        b3 = 25.0e0/108.0e0, b4 = 125.0e0/108.0e0;
        double e1 = 17.0e0/54.0e0, e2 = 7.0e0/36.0e0,
        e3 = 0.0, e4 = 125.0e0/108.0e0;
        
        for(i=0;i<NINT;i++)
        {
            y[i] = ysav[i] + b1*g1[i] + b2*g2[i] + b3*g3[i] + b4*g4[i];
            err[i] = e1*g1[i] + e2*g2[i] + e3*g3[i] + e4*g4[i];
        }
        
        
        
        // Time advances
        x = xsav + h;
        
        if(x==xsav)
        {
            // Something has gone wrong
            cout << "ERROR: Uh oh, h=0" << endl;
	    cout << "h = " << h << endl;
	    cout << "htry = " << htry << endl;
	    cout << "ysav = " << endl;
	    for(i=0;i<NINT;i++) cout << i << " " << ysav[i] << endl;
	    cout << "dysav = " << endl;
	    for(i=0;i<NINT;i++) cout << i << " " << dysav[i] << endl;

	    AbortChem();
        }
        
        
        // Determine the scaled accuracy
        REAL errmax = 0.0;
        for(i=0;i<NINT;i++){
	  errmax = max(errmax,fabs(err[i]/yscal[i]));
	  
	}
        errmax = errmax/eps;
       
        
        REAL safety = 0.9, grow= 1.5, pgrow= -0.25, shrnk= 0.5, pshrnk= -1.0/3.0;
        REAL errcon = 0.1296;
        if(errmax <= 1.0) // All went well
        {
            hdid = h;
            if(errmax > errcon) hnext = safety*h*pow(errmax,pgrow);
            else hnext = grow*h;
            jcounts = jcounts+jtry;
            
            if(jcounts>400000)
            {
                // Why is this taking so long?!
                cout << "Jcounts is too high (jcounts="<<jcounts<<")" << endl;
		cout << "h = " << h << endl;
		cout << "htry = " << htry << endl;
                cout << "temp = " << rpar[temp_par] << endl;  //ARS adding Temp print-out
		cout << "ysav = " << endl;
		for(i=0;i<NINT;i++) cout << i << " " << ysav[i] << endl;
		cout << "dysav = " << endl;
		for(i=0;i<NINT;i++) cout << i << " " << dysav[i] << endl;

		cout << "ratraw = " << endl;
		for(i=0;i<NREACTION;i++) cout << i << " " << ratraw[i] << endl;

		AbortChem();
            }
            
	   
            return;
            
        }
        else // Error too large, gots to try again
        {
            hnext = safety*h*pow(errmax,pshrnk);
            double hsign = 1.0;
            if(h<0) hsign = -1.0;
            if(h<0) cout << "ERROR: hsign is < 0?!" << endl; // should never happen
            h = hsign*max(fabs(hnext),shrnk*fabs(h));
            //cout << "Failed step " << h << " " << hnext << " " << hdid << endl;
        }
        
    } // End of solve equation loop

    // should never get here
    cout << "ERROR: Too many step tries (jtry ="<<jtry<<")" <<endl;
    AbortChem();
};



// Solves matrix equation:
//        a * x = b
// and puts solution in b
void leqs(REAL (*a)[NINT],REAL* b)
{
    int n1 = NINT-1;
    int i,j,k,l;
    REAL c,r;
    
    for(i=0;i<NINT;i++)
    {
        r = fabs(a[i][0]);
        for(j=1;j<NINT;j++)
        {
            c=fabs(a[i][j]);
            if(r<c) r=c;
        }
        
        for(j=0;j<NINT;j++) a[i][j] = a[i][j]/r;
        b[i] = b[i]/r;
    }
    
    for(j=0;j<n1;j++)
    {
        l = j+1;
        for(i=l;i<NINT;i++)
        {
            r = -a[i][j];
            if(r==0.0) continue; //goto onward;
            r=r/(a[j][j]);
            for(k=l;k<NINT;k++) a[i][k] = a[i][k] + r*a[j][k];
            b[i] = b[i] + r*b[j];
        //onward:
            continue;
        }
    }
    
    b[n1] = b[n1] / a[n1][n1]; //last entry
    
    for(l=0;l<n1;l++)
    {
        i= NINT-2-l;
        r=0;
        int imax = i+1;
        for(j=imax;j<NINT;j++)
        {
            int jj=i+NINT-j;
            r = r+(a[i][jj])*b[jj];
        }
        b[i]= (b[i]-r)/(a[i][i]);
    }
  
};
