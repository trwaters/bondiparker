#include <stdlib.h>
#include <iostream>											
#include <fstream>  										
#include "bondiparker.h"
#include "LambertW.h"
#include "LambertW.cc"

//#define DEBUG_SOLUTION  //dumps info for single solutions
//#define DEBUG_PSURVEY   //dumps info during the parameter survey

using namespace std;

void Solution::initialize(const int Nvals, double* x_user)
{

  if (x_user == NULL) // use default array length
    N = NMAX;
  else // determine the solution at user defined values
    N = Nvals; //ceil(sizeof(x_user)/sizeof(x_user[0]));
   //cout << "sizeof(x_user) = " << sizeof(x_user) << endl;
   //cout << "sizeof(x_user[0]) = " << sizeof(x_user[0]) << endl;

  // position array
  x_array = new double[N];

    // allocate memory
  for (int n = 0; n < MAX_ROOTS; n++)
  {
      M[n].subsonic = new double[N];
      M[n].supersonic = new double[N];
      M[n].transonic_inflow = new double[N];
      M[n].transonic_outflow = new double[N];

      v[n].subsonic = new double[N];
      v[n].supersonic = new double[N];
      v[n].transonic_inflow = new double[N];
      v[n].transonic_outflow = new double[N];

      d[n].subsonic = new double[N];
      d[n].supersonic = new double[N];
      d[n].transonic_inflow = new double[N];
      d[n].transonic_outflow = new double[N];

      P[n].subsonic = new double[N];
      P[n].supersonic = new double[N];
      P[n].transonic_inflow = new double[N];
      P[n].transonic_outflow = new double[N];

      econst[n].subsonic = new double[N];
      econst[n].supersonic = new double[N];
      econst[n].transonic_inflow = new double[N];
      econst[n].transonic_outflow = new double[N];
  }
    
}

/*----------------------------------------------------------------------------
 * Public function to be called from main to generate complete solutions
 *----------------------------------------------------------------------------*/
void Solution::run(const int Nvals, double* x_user=NULL, bool verbose=true, string path="do_not_print_output")
{
  clock_t start = clock();
  
  //Set necessary quantities
  cosi = cos(i_deg); //i_deg always shows up as cos(i_deg)
  

  //Allocate memory for solution Profiles
  initialize(Nvals,x_user);
   
  //Set the effective potential, effective force, and Keplerian speed
  setPotential(gravity);
  
  //Find critical points
  if (g==1.) pcritpt_func = &Solution::isothermal_critpt;
  else pcritpt_func = &Solution::polytropic_critpt;
  bracket_roots();
  root_bisection();

  //Tidy up isothermal solution
  if(g==1.) 
  {
    remove_iso_vel_minima_roots();  //for g==1, locations of v_min satisfy crit pt eqn 
    for(int i=0; i<nroots; i++) lcs[i] = lp;  //lc=lp for g==1 by definition
  }
  
  //Solve for Mach Number and other flow variables
  solutionProfiles(x_user); 

  //display terminal output
  if (verbose)
    cout_output(start); 
  
  //Print results to file:
  if (path != "do_not_print_output")
  { 
    filepath = path;
    if (L < 1. || L > 1.) {  //everywhere sub/supersonic curves occur for L<1
      print_Xtype_solns();   //everywhere sub/super-critical curves occur for L>1 
    }
    else if (L == 1.) {
      print_Ttype_solns();  //Ttype: transonic solutions satisfying the density BC
      // print_Xtype_solns();  //Xtype: the X-type transonic solution topology 
    }
  }

}

/*----------------------------------------------------------------------------
 * Setter for global parameter rho_bc
 *----------------------------------------------------------------------------
void setDensityBC(const double value);
{
    rho_bc = value;
}*/

/*----------------------------------------------------------------------------
 * Set the potential and force based on input from main
 *----------------------------------------------------------------------------*/
void Solution::setPotential(string choice)
{
  bool keplerian = (fabs(zeta-sqrt(lp))<TINY_NUMBER) ? true:false;
  if(choice == "Pseudo-Newtonian") {
      pUeff = &Solution::pseudoNewtonianUeff;
      pFeff = &Solution::pseudoNewtonianFeff;
      gravity = choice;
      if (keplerian) zeta = sqrt(8.*lp)/(1. + lp*bdy_());
  }
  else {
      pUeff = &Solution::newtonianUeff;
      pFeff = &Solution::newtonianFeff;
      gravity = "Newtonian";
      if (keplerian) zeta = sqrt(lp)/(1. + lp*bdy_());
  }
}
 
/*----------------------------------------------------------------------------
 * Build the remaining profiles from the Mach number and critical point
 *----------------------------------------------------------------------------*/
void Solution::solutionProfiles(double* x_user)
{
  double x,dx,xmin,xmax;
  
  xmin = 0.;
  xmax = outer_bdy; //(bdy_()==0.) ? OUTER_BDY:bdy_();
  // dx = (xmax-xmin)/(NMAX-1);
  dx = (xmax-xmin)/N;
  
#ifdef LOG_SPACING
  xmin = LOGXMIN;
  xmax = log10(xmax);
  // dx = (xmax-xmin)/(NMAX-1);
  dx = (xmax-xmin)/N;
#endif
  
  double Mo_sub = get_subsonic_MachNumber(0.);
  double Mo_sup = get_supersonic_MachNumber(0.);
  
  double Msub,vsub,dsub,Psub,ecsub;
  double Msup,vsup,dsup,Psup,ecsup;
  
  n_inflow = 0;
  n_outflow = 0;

  for (int n=0; n<nroots; n++)
  {
	  xc = xcs[n]; 
	  lc = lcs[n];
	  // cout << "\nNROOTS " << nroots << endl;
	  // cout << "ROOT " << n << endl;

    //set global variables tracking the inflow/outflow root
    if (ids[n] == 0)
      n_inflow = n;
    else
      n_outflow = n;

    // build the solution profiles
	  for (int i=0; i < N; i++)
		{
		
		if (x_user == NULL)
		{
		  x = xmin + i*dx;
	    #ifdef LOG_SPACING
	  	  x = pow(10.,x);
	    #endif
	    x_array[i] = x;
		}
		else 
		{
		  x = x_user[i];
		  x_array[i] = x;
		}
		  
		  //Build Mach Number Profiles
		  Msub = get_subsonic_MachNumber(x);
		  Msup = get_supersonic_MachNumber(x);
		  M[n].subsonic[i] = Msub; M[n].supersonic[i] = Msup;
		  
		  if (L == 1.) 
		  {
		    M[n].transonic_outflow[i] = (x<xc) ? Msub : Msup;
		    M[n].transonic_inflow[i]  = (x<xc) ? Msup : Msub;
		  }
		  
		  //Build profiles of remaining solution variables
		  getFlowVariables(Msub,x,&vsub,&dsub,&Psub,&ecsub); //subsonic values
		  getFlowVariables(Msup,x,&vsup,&dsup,&Psup,&ecsup); //supersonic values
		  if (g == 1.0) {
		    dsub = density_lambertWfct(x,Msub,Mo_sub);
		    dsup = density_lambertWfct(x,Msup,Mo_sup);
		  }
		  
		  v[n].subsonic[i] = vsub; v[n].supersonic[i] = vsup;
		  d[n].subsonic[i] = dsub; d[n].supersonic[i] = dsup;
		  P[n].subsonic[i] = Psub; P[n].supersonic[i] = Psup;
		  econst[n].subsonic[i] = ecsub; econst[n].supersonic[i] = ecsup;
		  
		  if (L == 1.) 
		  {
		    v[n].transonic_outflow[i] = (x<xc) ? vsub : vsup;
		    v[n].transonic_inflow[i]  = (x<xc) ? vsup : vsub;
		    
		    d[n].transonic_outflow[i] = (x<xc) ? dsub : dsup;
		    d[n].transonic_inflow[i]  = (x<xc) ? dsup : dsub;
		    
		    P[n].transonic_outflow[i] = (x<xc) ? Psub : Psup;
		    P[n].transonic_inflow[i]  = (x<xc) ? Psup : Psub;
		    
		    econst[n].transonic_outflow[i] = (x<xc) ? ecsub : ecsup;
		    econst[n].transonic_inflow[i]  = (x<xc) ? ecsup : ecsub;
		  }
		}
  }
}

/*----------------------------------------------------------------------------
 * Perform bisection in the subsonic regime
 *----------------------------------------------------------------------------*/
double Solution::get_subsonic_MachNumber(const double x)
{
  if(g==1.) return mach_lambertWfct(x,0); //principal real branch of Lambert-W fct
  else return mach_bisection(x,0.);
}

/*----------------------------------------------------------------------------
 * Perform bisection in the supersonic regime
 *----------------------------------------------------------------------------*/
double Solution::get_supersonic_MachNumber(const double x)
{
  if(g==1.) return mach_lambertWfct(x,-1); //-1 real branch of Lambert-W fct
  else return mach_bisection(x,1.);
}

/*------------------------------------------------------------------------------
 * Bracket all roots of the critical point equation
 *----------------------------------------------------------------------------*/
void Solution::bracket_roots()
{
  double x1,x2,dx,root,fprobe,xstep,fstep,fprod;
  int nfound,nbrak;
  double rc_bondi = 0.25*(5. - 3.*g);
  if (g <= 1.1 || g >= 1.5) nbrak = 1e2*UPPER_LIMIT; //use fine resolution
  else nbrak = 2e1*UPPER_LIMIT; //(2e1*OUTER_BDY was used in WP12)
  x1 = 0.;
  x2 = (bdy_() == 0.) ? UPPER_LIMIT:1.; //2.*rc_bondi; //For Bondi,all crit pts within Bondi radius 
  dx = (x2-x1)/nbrak;
  fprobe = (this->*pcritpt_func)(x1); //critical_point_func(x1);
  nfound = 0;
  nroots = MAX_ROOTS; //max number of roots sought
  xstep = x1;
  for (int j=1; j<nbrak; j++)
  {
      xstep = xstep+dx;
      fstep = (this->*pcritpt_func)(xstep); //critical_point_func(xstep);
      fprod = fprobe*fstep;
      if (fprod <=0.0)
      {
		xl[++nfound-1] = xstep-dx;
		xr[nfound - 1] = xstep;
		if (nfound == MAX_ROOTS) break; 
  	  }
      fprobe = fstep;
  }
  nroots = nfound; //number of roots actually found
  #ifdef DEBUG_SOLUTION
  cout << "nroots found " << nroots << "\tHEP " << lp << endl;
  #endif
}

/*------------------------------------------------------------------------------
 * Perform bisection on roots bracketed from above
 *----------------------------------------------------------------------------*/
void Solution::root_bisection()
{
  double x1,x2,dx,root,product;
  int iter = 0;
  root = 0.;
  for (int j=0; j<(nroots); j++) { //loop through all brackets and find the roots
    x1 = xl[j];
    x2 = xr[j];
    dx = x2-x1;
    while (fabs(dx) > TOL) 
    {
	  root = 0.5*(x1+x2);
	  product = (this->*pcritpt_func)(x1)*(this->*pcritpt_func)(root); 
	  if (product < 0.)
	  {
	    x2  = root;
	    dx = x2-x1;
	  }
	  else
	  {
	    x1  = root;
	    dx = x2-x1;
	  }
    }
    if (root>TINY_NUMBER) {xcs[j] = root; lcs[j] = lc; ids[j] = id;}
  }
}

/*------------------------------------------------------------------------------
 * Perform standard bisection to find the Mach Number
 *----------------------------------------------------------------------------*/
double Solution::mach_bisection(const double x, const double Mbound)
{
  double M1,M2,dM,root,product;
  M1 = (Mbound < 1.) ? TINY_NUMBER : 1e3;
  M2 = 1.; //M = 1 at sonic point by definition
  dM = M2-M1;
  while  (fabs(dM) > TOL) 
  {
    root = (M1+M2)/2.;
    product = mach_func(x,M1)*mach_func(x,root); 
    if (product < 0.) {
	  M2 = root;
	  dM = M2-M1;
	}
    else {
	  M1 = root;
	  dM = M2-M1;
	}
  }
  return root;
} 

/*------------------------------------------------------------------------------
 * Evaluate the polytropic critical point eqn, eqn (3.42) in WP12
 *----------------------------------------------------------------------------
double Solution::polytropic_critpt(const double x)
{
  double a,ao,aprime,so,wo;
  double Ueff_o,Ueff,Feff;
  double rhs,ec,gm1,gp1;
  double ro = 1./lp;
  gm1 = g-1.;
  gp1 = g+1.;
  Ueff_o = (this->*pUeff)(bdy_());
  Ueff = (this->*pUeff)(x);
  Feff = (this->*pFeff)(x);
  ao = area_(bdy_());
  a = area_(x);
  aprime = area_prime_(x);
  lc = lp*aprime/(a*Feff);
  so = lc/lp;
  wo = pow((a/ao),2.)*pow(lp/lc,gp1/gm1);
  rhs = .5*wo + Ueff_o + 1./gm1;
  ec = .5*gp1/gm1 + so*Ueff;
  
  //Classify roots: id==1 for outflow roots, 0 for inflow roots
  if(bdy_()==0.) id = (wo > 1.) ? 0 : 1; //inner BCs: wo < 1 is outflow
  else id = (wo > 1.) ? 1 : 0; //outer BCs: wo < 1 is inflow
  
  return ec/so - rhs;
} */

/*------------------------------------------------------------------------------
 * Evaluate the polytropic critical point eqn, eqn (3.42) in WP12,
 * generalized to include an arbitrary density boundary condition
 *----------------------------------------------------------------------------*/
double Solution::polytropic_critpt(const double x)
{
  double a,ao,aprime,so,wo;
  double Ueff_o,Ueff,Feff;
  double rhs,ec,gm1,gp1;
  double ro = 1./lp;
  gm1 = g-1.;
  gp1 = g+1.;
  Ueff_o = (this->*pUeff)(bdy_());
  Ueff = (this->*pUeff)(x);
  Feff = (this->*pFeff)(x);
  ao = area_(bdy_());
  a = area_(x);
  aprime = area_prime_(x);
  lc = lp*aprime/(a*Feff);
  so = lc/lp;
  wo = pow((a/ao),2.)*pow(lp/lc,gp1/gm1)*pow(rho_bc,-gp1);
  rhs = pow(rho_bc,gm1)*(0.5*wo + 1./gm1) + Ueff_o;
  ec = 0.5*gp1/gm1 + Ueff*lc/lp;
  
  //Classify roots: id==1 for outflow roots, 0 for inflow roots
  if(bdy_()==0.) id = (wo > 1.) ? 0 : 1; //inner BCs: wo < 1 is outflow
  else id = (wo > 1.) ? 1 : 0; //outer BCs: wo < 1 is inflow
  
  return ec*lp/lc - rhs;
}

/*------------------------------------------------------------------------------
 * Evaluate the isothermal critical point eqn, eqn (3.52) in WP12
 *----------------------------------------------------------------------------*/
double Solution::isothermal_critpt(const double x)
{
  double a,ao,aprime,wo,Feff;
  double ro = 1./lp;
  Feff = (this->*pFeff)(x);
  ao = area_(bdy_());
  a = area_(x);
  aprime = area_prime_(x);
  
  //Classify roots: id==1 for outflow roots, 0 for inflow roots
  if(bdy_()==0.) id = 1; //inner BCs: all roots should be outflow
  else id = 0; //outer BCs: all roots should be inflow
  
  return Feff - aprime/a;
}

/*------------------------------------------------------------------------------
 * Method to distinguish between velocity minima and critical points
 *----------------------------------------------------------------------------*/
void Solution::remove_iso_vel_minima_roots()
{
  double h = 1e-2;
  double machprime;

  for (int i=0; i<nroots; i++)
  { 
    xc = xcs[i]; //set xc for external functions
    machprime = ( mach_lambertWfct(xcs[i] + h, -1) - mach_lambertWfct(xcs[i] - h, 0) )/(2.*h);
    #ifdef DEBUG_SOLUTION
    cout << "root " << i << " has W_o(xc+h) = " << mach_lambertWfct(xcs[i],0) 
					<< "\t and W_-1(xc-h) = " << mach_lambertWfct(xcs[i],-1) << endl; 
	#endif
    if ( machprime!=machprime ) //if Mach derivative around xc is NaN
      ids[i] = 0; //then this crit pt marks a velocity minimum
    else ids[i] = 1;
  }
  
  int ncritpts = nroots;
  int j = 0;
  while (j<ncritpts) //remove any velocity minima roots
  {
    if(ids[j]==0) 
      for (int i=j; i<(ncritpts-1); i++)
      {
        ids[i] = ids[i+1]; 
    	xcs[i] = xcs[i+1]; 
    	nroots--;
      }
    j++;
  }
}

/*------------------------------------------------------------------------------
 * Evaluate the explicit solution, eqn (3.44) in WP12
 *----------------------------------------------------------------------------*/
double Solution::mach_func(const double x, const double M)
{
  double Ueff,Ueff_c;
  double Xfactor,FofM,Xofx;
  double so = lc/lp;
  double gexp = 2.*(g-1.)/(g+1.);
  double a = area_(x);
  double ac = area_(xc);
  double ro = 1./lp;
  Ueff = (this->*pUeff)(x);
  Ueff_c = (this->*pUeff)(xc);
  FofM = (1./(M*M*(g-1.)) + .5)*pow(M,4./(g+1.));
  Xfactor = 1./gexp + so*(Ueff_c - Ueff);
  Xofx = pow(a/ac,gexp)*Xfactor;
  return FofM - pow(L,-gexp)*Xofx;
}   

/*----------------------------------------------------------------------------
 * Build remaining flow variables from the Mach number and critical point
 *----------------------------------------------------------------------------*/
void Solution::getFlowVariables(const double M, const double x, double* v, double* d, double* P, double* ec)
{
  double gm1,gp1,units;
  double v2,w,s,T,Ueff;
  double a = area_(x);
  double ac = area_(xc);
  double ro = 1./lp;
  gm1 = g - 1.;
  gp1 = g + 1.;
  units=.5/lc;  	
  *v= pow(units,.5)*pow((ac/a),gm1/gp1)*pow(M,2./gp1); // escape speed units
  v2 = (*v)*(*v);
  w = M*M;
  s = 2.*lc*v2/w;
  T = (lp/lc)*s;
  *d = pow(T,1./gm1); 
  Ueff = (this->*pUeff)(x);
  *P = pow(*d,g)/(2.*lp*g);
  *ec = (lc/lp)*Ueff + s*(.5*w + 1./gm1);
}

/*------------------------------------------------------------------------------
 * Evaluate the isothermal solution, eqn (3.51) in WP12
 *----------------------------------------------------------------------------*/
double Solution::mach_lambertWfct(const double x, int branch)
{
  double Ueff,Ueff_c;
  double gammaB;
  double a = area_(x);
  double ac = area_(xc);
  double ao = area_(bdy_());
  double ro = 1./lp;
  Ueff = (this->*pUeff)(x);
  Ueff_c = (this->*pUeff)(xc);
  gammaB = (ac/ao)*exp(-0.5 - Ueff_c);
  double w = -LambertW(branch,-pow(L*gammaB*exp(Ueff)/(a/ao),2.)); //mach number squared
  return pow(w,0.5);
}  

/*------------------------------------------------------------------------------
 * Evaluate the isothermal density, eqn (3.56) in WP12
 *----------------------------------------------------------------------------*/
double Solution::density_lambertWfct(const double x, const double M, const double Mo)
{
  double a = area_(x);
  double ao = area_(bdy_());
  
  return Mo*ao/(M*a);
}  

/*------------------------------------------------------------------------------
 * Determine the value of the initial Mach Number given xc,lc
 *----------------------------------------------------------------------------*/	 
double Solution::getMo(const double xc, const double lc)
{
  double Mo;
  if(g==1.) {
    //Mo = (bdy_()==0.) ? mach_lambertWfct(0.,0) : mach_lambertWfct(0.,-1);
    Mo = mach_lambertWfct(0.,0);
  }
  else {
    double ac = area_(xc);
    double ao = area_(bdy_());
  	double gm1 = g-1.;
	double gp1 = g+1.;
	Mo = L*(ac/ao)*pow(lp/lc,.5*gp1/gm1)/pow(rho_bc,0.5*gp1);
  }
  return Mo;
}

/*------------------------------------------------------------------------------
 * Determine the velocity at chi=0
 *----------------------------------------------------------------------------*/   
double Solution::getv0(const double xc, const double lc)
{
  double v0,M0;
  if(g==1.) {
    M0 = mach_lambertWfct(0.,0);
  }
  else {
    double ac = area_(xc);
    double a0 = area_(0.);
    double gm1 = g-1.;
  double gp1 = g+1.;
  M0 = get_supersonic_MachNumber(0.);
  v0 = pow(L*(ac/a0),gm1/gp1)*sqrt(0.5/lc)*pow(M0,2./gp1);
  }
  return v0;
}

/*------------------------------------------------------------------------------
 * Determine the mass flux density (in units of rho_o c_o)
 *----------------------------------------------------------------------------*/   
double Solution::getmdot(const double xc, const double lc)
{
  double mfd,mfd_bondi;
  if(g==1.) {
    mfd = 1.;
  }
  else {
    double ac = area_(xc);
    // double a0 = area_(bdy_());
    double gm1 = g-1.;
  double gp1 = g+1.;
  mfd = ac*pow(lp/lc,0.5*gp1/gm1);
  // mfd = lc; //pow(lp/lc,0.5*gp1/gm1);
  // mfd_bondi = 0.25*pow(0.5*(5.-3.*g),-0.5*(5.-3.*g)/gm1);
  }
  return mfd;
}

/*------------------------------------------------------------------------------
 * Determine the value of the initial density given xc,lc
 *----------------------------------------------------------------------------*/	 
double Solution::getdo(const double xc, const double lc)
{
    double v,units,wo,s,T;
    double ac = area_(xc);
    double ao = area_(bdy_());
  	double gm1 = g-1.;
	double gp1 = g+1.;
	// double Mo = L*(ac/ao)*pow(lp/lc,.5*gp1/gm1);
	// wo = Mo*Mo;
	// units=.5/lc;
	// v= pow(units,.5)*pow((ac/ao),gm1/gp1)*pow(wo,1./gp1);
	// s = 2.*lc*v*v/wo;
 //    T = (lp/lc)*s;
 //    return pow(T,1./gm1); 

  double Mo = get_subsonic_MachNumber(bdy_());
  return pow(lp/lc,1./gm1)*pow(ac/ao/Mo,2./gp1);
}

/*------------------------------------------------------------------------------
 * Define the Newtonian Effective Potential: -GM/r + U_centrifugal
 *----------------------------------------------------------------------------*/
double Solution::newtonianUeff(const double x)
{
  double U,Urot;
  double ro = 1./lp;
  U  = -1./sqrt(x*x + 2.*x*ro*cosi + ro*ro);  
  Urot = 0.5*zeta*zeta*(ro+bdy_())*(ro+bdy_())/((ro+x*cosi)*(ro+x*cosi));
  return U + Urot; 
}

/*------------------------------------------------------------------------------
 * Define the (minus of the) Newtonian Effective Force: dUeff/dx
 *----------------------------------------------------------------------------*/
double Solution::newtonianFeff(const double x)
{
  double Fg,Fc;
  double ro = 1./lp;
  Fg = (x+ro*cosi)/pow(x*x + 2.*x*ro*cosi + ro*ro,1.5);
  Fc = zeta*zeta*(ro+bdy_())*(ro+bdy_())*cosi/pow((ro + x*cosi),3.); 
  return Fg - Fc;
}

/*------------------------------------------------------------------------------
 * Define the Paczynski-Wiita Effective Potential: -GM/(r-r_s) + U_centrifugal
 *----------------------------------------------------------------------------*/
double Solution::pseudoNewtonianUeff(const double x)
{
  double U,Urot;
  double ro = 1./lp;
  double r = sqrt(x*x + 2.*x*ro*cosi + ro*ro);
  U  = -1./(r-ro);  
  Urot = 0.5*zeta*zeta*(ro+bdy_())*(ro+bdy_())/((ro+x*cosi)*(ro+x*cosi));
  return U + Urot; 
}

/*------------------------------------------------------------------------------
 * Define the (minus of the) Paczynski-Wiita Effective Force: dUeff/dx
 *----------------------------------------------------------------------------*/
double Solution::pseudoNewtonianFeff(const double x)
{
  double Fg,Fc;
  double ro = 1./lp;
  double r = sqrt(x*x + 2.*x*ro*cosi + ro*ro);
  Fg = (x+ro*cosi)/((r-ro)*(r-ro)*r);
  Fc = zeta*zeta*(ro+bdy_())*(ro+bdy_())*cosi/pow((ro + x*cosi),3.); 
  return Fg - Fc;
}

/*------------------------------------------------------------------------------
 * Output summary to terminal
 *----------------------------------------------------------------------------*/
void Solution::cout_output(clock_t start)
{
  string model = identifyModelName();
  cout << "===============================================" << endl <<
  model.c_str() << endl;
  cout << "===============================================" << endl <<
  "Inputs:" << endl <<
  "HEP = " << lp << endl << 
  "Gamma = " << g << endl <<
  "Inclination angle = " << i_deg << " degrees" << endl <<
  "Rotation rate = ";
  bool keptest = (fabs(zeta-sqrt(lp))<TINY_NUMBER) ? true:false;
  if (keptest) { cout << "Keplerian" << endl;}
  else { cout << zeta << endl;}
  cout << "===============================================" << endl;
  
  cout<<"Calculations took " << 
  		(clock() - start) / (double)CLOCKS_PER_SEC << " seconds." << endl;
  if (nroots == 1)
  cout << endl << nroots <<
  		" critical point solution was found:" << endl;
  else
  cout << endl << nroots <<
  		" critical point solutions were found:" << endl;
  		
  if (nroots == 0) 
    cout << "Invalid input parameters.  Consult the parameter survey in WP12." << endl; 
   
  double Mo,v0,d0,Egas0,ec;
  string bdyspeed;
  
  for (int i=0; i<nroots; i++)
  { 
      Mo = getMo(xcs[i],lcs[i]);
      getFlowVariables(Mo, 0., &v0, &d0, &Egas0, &ec);
      if (ids[i] == 1) bdyspeed = (bdy_()==0.) ? "SUBSONIC":"SUPERSONIC";
      else if (ids[i] == 0) bdyspeed = (bdy_()==0.) ? "SUPERSONIC":"SUBSONIC";
      cout << "Solution " << i+1 << " is " 
    	<< bdyspeed << " at the boundary with:\n"
		<< "id= " << ids[i]  << "\t"
  		<< "xc= " << xcs[i]  << "\t"
  		<< "lc/ro= " << lp*xcs[i]  << "\t"
        //<< "lc= " << lcs[i] << "\t"
        << "d_o= " << d0 << "\t"
        << "M_o= " << Mo << "\t"
        << "v_o= " << v0 << "\t"
        << "P_o= " << Egas0 << "\t"
        << "ec = " << ec << endl;
  }
}

/*------------------------------------------------------------------------------
 * Print full solution topology to file
 *----------------------------------------------------------------------------*/
void Solution::print_Xtype_solns()
{
  ofstream results;
  for (int n=0; n<nroots; n++)
  {
      string tag = "X_";
      string filetowrite = default_filename(n+1,tag);
	  results.open(filetowrite.c_str(), ios::trunc);
	  results <<setw(10) << left << "x= " << 
				setw(10) << left << "\tM(sub)= " << 
				setw(10) << left << "\tM(super)= " <<
				setw(10) << left << "\tv(sub)= " <<
				setw(10) << left << "\tv(super)= " <<
				setw(10) << left << "\td(sub)= " <<
				setw(10) << left << "\td(super)= " << 
				setw(10) << left << "\tP(sub)= " <<
				setw(10) << left << "\tP(super)= " << 
				setw(10) << left << "\tec(sub)= " <<
				setw(10) << left << "\tec(super)= " << endl;
	  for (int i=0; i < N; i++)
	  {
		results << 	setprecision(10) << 
					setw(10) << left << x_array[i] << "\t" <<
					setw(10) << left << M[n].subsonic[i] << "\t" << 
					setw(10) << left << M[n].supersonic[i] << "\t" <<
					setw(10) << left << v[n].subsonic[i] << "\t" << 
					setw(10) << left << v[n].supersonic[i] << "\t" <<
					setw(10) << left << d[n].subsonic[i] << "\t" << 
					setw(10) << left << d[n].supersonic[i] << "\t" <<
					setw(10) << left << P[n].subsonic[i] << "\t" << 
					setw(10) << left << P[n].supersonic[i] << "\t" <<
					setw(10) << left << econst[n].subsonic[i] << "\t" << 
					setw(10) << left << econst[n].supersonic[i] << endl;
	  }
	  results.close();
	  filetowrite = filename;
  }
}

/*------------------------------------------------------------------------------
 * Print transonic solutions to file
 *----------------------------------------------------------------------------*/
void Solution::print_Ttype_solns()
{
  ofstream results;

  for (int n=0; n<nroots; n++)
  {
    string tag = "T_";
    string filetowrite = default_filename(n+1,tag);
	results.open(filetowrite.c_str(), ios::trunc);
	results <<setw(10) << left << "x= " << 
			  setw(10) << left << "\tMach Number= " << 
			  setw(10) << left << "\tvelocity= " << 
			  setw(10) << left << "\tdensity= " << 
			  setw(10) << left << "\tPressure= " << 
			  setw(10) << left << "\tCritPtConst= " << endl;
    if(ids[n]==1) {  //then solution satisfying BCs is subsonic at x=0
  	  for (int i=0; i < N; i++)
  	  {
  		results << 	setprecision(10) << 
  					setw(10) << left << x_array[i] << "\t" <<
  					setw(10) << left << M[n].transonic_outflow[i] << "\t" << 
  					setw(10) << left << v[n].transonic_outflow[i] << "\t" << 
  					setw(10) << left << d[n].transonic_outflow[i] << "\t" << 
  					setw(10) << left << P[n].transonic_outflow[i] << "\t" << 
  					setw(10) << left << econst[n].transonic_outflow[i] << endl;
  	  }
	  }
	  else if (g != 1.) 
    {  //solution satisfying BCs is supersonic at x=0
  	  for (int i=0; i < N; i++)
  	  {
  		results << 	setprecision(10) << 
  					setw(10) << left << x_array[i] << "\t" <<
  					setw(10) << left << M[n].transonic_inflow[i] << "\t" << 
  					setw(10) << left << v[n].transonic_inflow[i] << "\t" << 
  					setw(10) << left << d[n].transonic_inflow[i] << "\t" << 
  					setw(10) << left << P[n].transonic_inflow[i] << "\t" << 
  					setw(10) << left << econst[n].transonic_inflow[i] << endl;
  	  }
	  }
	results.close();
	filetowrite = filename;
  }
}

/*------------------------------------------------------------------------------
 * Determine the name of the model discussed in WP12 given the parameters
 *----------------------------------------------------------------------------*/
string Solution::identifyModelName()
{
  string modelName = classtag();
  if (i_deg == 0. && modelName == "CON") modelName = "Classic Parker Model";
  else if (i_deg == 0. && modelName == "CIA") modelName = "Cylindrical Parker Model";
  else if (i_deg == 0. && modelName == "BONDI") modelName = "Bondi Model";
  
  if (modelName == "CON") modelName = "Converging Model";
  else if (modelName == "CIA") modelName = "CIA Model";
  
  if (zeta != 0.) modelName += " with Rotation";
  
  if (gravity == "Pseudo-Newtonian") modelName += " (Paczynski-Witta Potential)";
  
  return modelName;
}

/*------------------------------------------------------------------------------
 * Construct the default filename for the solution file
 *----------------------------------------------------------------------------*/
string Solution::default_filename(int n, string label)
{
  string name = filepath + label;
  name += classtag() + "soln" + to_string<int>(n) + "of" + to_string<int>(nroots);
  if(label=="T_") name += (ids[n-1]==1) ? "(outflow)" : "(inflow)";
  name += "_HEP" + to_string<double>(lp);
  name += "_g" + to_string<double>(g);
  if (classtag() == "BONDI") name += "_at" + to_string<double>(outer_bdy);
  if (i_deg != 0.) name += "_i" + to_string<double>(i_deg/(DEG));
  int keptest = (fabs(zeta-sqrt(lp))<TINY_NUMBER) ? -1:1;
  if (keptest < 0) name += "Kep";
  else if (zeta != 0.) name += "_rot" + to_string<double>(zeta);
  if(gravity == "Pseudo-Newtonian") name += "(PW)";
  name += ".tab";
  return name;
}

/*----------------------------------------------------------------------------
 * Public function to be called from main to generate a parameter survey
 *----------------------------------------------------------------------------*/
void ParameterSurvey::run(string path="")
{
  clock_t start = clock(); 
  
  //set necessary values in Solution base class
  soln->cosi = cos(soln->i_deg);
  soln->filepath = path;
  
  //Set the effective potential and effective force
  soln->pUeff = &Solution::newtonianUeff;
  soln->pFeff = &Solution::newtonianFeff;
  
  //Determine the smallest possible HEP for the given gamma, i_deg, zeta
  double minHEP;
  if(useMin) minHEP = findMinHEP();
  else minHEP = soln->lp;
  
  //Output results
  string model = soln->identifyModelName();
  string startcondition = (useMin) ? "Minimum HEP" : "Custom HEP";
  cout << "=======================================================" << endl <<
  "Parameter Survey: " << model.c_str() << endl;
  cout << "=======================================================" << endl <<
  "Inputs:" << endl <<
  "Starting HEP = " << startcondition.c_str() << " = " << minHEP << endl <<
  "Gamma = " << soln->g << endl <<
  "Inclination angle = " << soln->i_deg << " degrees" << endl <<
  "Rotation rate = ";
  if (kep) { cout << "Keplerian" << endl;}
  else { cout << soln->zeta << endl;}
  cout << "=======================================================" << endl;
  
  //Run parameter survey
  doSurvey(minHEP); 
  
  if (useMin)
  cout <<"\nMINIMUM HEP FOUND = " << minHEP << endl;
  
  cout <<"\nCalculations took " << 
  		(clock() - start) / (double)CLOCKS_PER_SEC << " seconds." << endl;
  		
  print_psurvey();
}

/*----------------------------------------------------------------------------
 * Public function just to query if a solution exists for given inputs
 *----------------------------------------------------------------------------*/
void ParameterSurvey::query()
{
  clock_t start = clock();
  
  //set necessary values in Solution base class
  soln->cosi = cos(soln->i_deg); 
  if(kep) soln->zeta = sqrt(soln->lp)/(1. + soln->lp*soln->bdy_());
  
  //Set the effective potential and effective force
  soln->pUeff = &Solution::newtonianUeff;
  soln->pFeff = &Solution::newtonianFeff;
  
  //Find critical points
  if (soln->g==1.) soln->pcritpt_func = &Solution::isothermal_critpt;
  else soln->pcritpt_func = &Solution::polytropic_critpt;
  soln->bracket_roots();
  soln->root_bisection();
  
  double val1,val2;
  val1 = soln->polytropic_critpt(0.);
  val2 = soln->polytropic_critpt(1./soln->lp);
  cout << "FUNC(0) = " << val1 << endl;
  cout << "xc = " << soln->xcs[1] << endl;
  cout << "Ueff(xc) = " << soln->newtonianUeff(soln->xcs[1]) << endl;
  cout << "Mtest = " << soln->mach_func(soln->xcs[1],1.) << endl;
  //display terminal output
  soln->cout_output(start); 
}

/*------------------------------------------------------------------------------
 * Carry out the parameter survey of HEP values leading to transonic solutions
 *----------------------------------------------------------------------------*/
void ParameterSurvey::doSurvey(double minHEP)
{  
  minHEP = 0.01*ceil(minHEP*100.); //round HEP to nearest 100th
  
  int nroots_prior = 0;
  int k=0;
  while (k < HMAX) {
  	soln->lp = minHEP + HEP_STEP_SIZE*k;
  	hep[k] = soln->lp;
  
  	//Refresh crit pts
  	for (int n=0; n<MAX_ROOTS; n++) soln->xcs[n] = 0.;
	if(kep) soln->zeta = sqrt(soln->lp)/(1. + soln->lp*soln->bdy_());
	if (soln->g==1.) soln->pcritpt_func = &Solution::isothermal_critpt;
	else soln->pcritpt_func = &Solution::polytropic_critpt;
	soln->bracket_roots();
	/*
	while(soln->nroots < nroots_prior && soln->nbrackets < 1e3.*UPPER_LIMIT)
	{
	  soln->nbrackets*=2; 
	  soln->bracket_roots();
	}
	nroots_prior = soln->nroots;
	*/
	soln->root_bisection();
	if(soln->g==1.){
	  soln->remove_iso_vel_minima_roots();
	  for(int i=0; i<soln->nroots; i++) soln->lcs[i] = soln->lp;  //lc=lp for g==1 
	}
	
	for (int n=0; n<MAX_ROOTS; n++)
    {  
      //Note roots are found from closest to farthest. It is typically the close ones
      //that disappear so if there is only one root it is stored in the 2nd column
      //so that the roots line up when there are two
      if(soln->nroots>0){ 
		  if(soln->ids[n] == 0) {
		    psurvey[n].xc_inflow[k] = soln->xcs[n]; 
		    psurvey[n].lc_inflow[k] = soln->lcs[n]; 
		    psurvey[n].Mo_inflow[k] = soln->getMo(soln->xcs[n],soln->lcs[n]);
		    if(soln->g != 1.) 
		    psurvey[n].do_inflow[k] = soln->getdo(soln->xcs[n],soln->lcs[n]);
		  }
		  else if(soln->ids[n] == 1) {
		    psurvey[n].xc_outflow[k] = soln->xcs[n]; 
		    psurvey[n].lc_outflow[k] = soln->lcs[n]; 
		    psurvey[n].Mo_outflow[k] = soln->getMo(soln->xcs[n],soln->lcs[n]);
		    if(soln->g != 1.) 
		    psurvey[n].do_outflow[k] = soln->getdo(soln->xcs[n],soln->lcs[n]);
		  }  
	  } 
	}
	#ifdef DEBUG_PSURVEY
	if(k==0 || k%10==0) cout << "k = " << k << "\tHEP= " << soln->lp 
						<< "\txc1= " << soln->xcs[0] << " id1= " << soln->ids[0]
						<< "\txc2= " << soln->xcs[1] << " id2= " << soln->ids[1]
						<< "\txc3= " << soln->xcs[2] << " id3= " << soln->ids[2] << endl;
	#endif
	
	for (int n=0; n<MAX_ROOTS; n++) soln->ids[n] = -1; //reset ids to prevent storing previous roots
	++k;
  }
}

/*------------------------------------------------------------------------------
 * Determine the minimum value of the HEP yielding a transonic solution
 *----------------------------------------------------------------------------*/
double ParameterSurvey::findMinHEP()
{
  double lower,mid,upper;
  double xc,sqrt_lp;
  double g = soln->g;

  soln->nroots = 0;
  soln->lp = 0.0; //starting HEP
  
  //increase HEP until first root is found (HEP limit prevents infinite loop)
  while(soln->nroots == 0 && soln->lp < 200.) { 
      soln->lp += 0.1; 
      
	  //Refresh crit pts
	  if(kep) soln->zeta = sqrt(soln->lp)/(1. + soln->lp*soln->bdy_());
	  if (g==1.) soln->pcritpt_func = &Solution::isothermal_critpt;
	  else soln->pcritpt_func = &Solution::polytropic_critpt;
	  soln->bracket_roots();
	  soln->root_bisection();
	  if(g==1.) soln->remove_iso_vel_minima_roots(); 

	  xc = soln->xcs[0];
	  #ifdef DEBUG_PSURVEY
	  	cout << "HEP= " << soln->lp << "\txc= " << xc << endl;
	  #endif
    }  
  
  if(soln->lp>200.) {cout<< "No Solutions found for HEP<200!" << endl; return 0.;}
  upper = soln->lp;
  lower = upper - 0.1;
  mid = 0.5*(upper+lower);
  float diff = 0.5*(upper-lower);
  
  //now determine minimum HEP to high precision
  while (diff > 1e-6) { 
    soln->xcs[0] = 0.;
	soln->lp = mid;
	#ifdef DEBUG_PSURVEY
	cout << "HEP= " << soln->lp << "\tdiff= " << diff << endl;
	#endif
	//Refresh crit pts
	if(kep) soln->zeta = sqrt(soln->lp)/(1. + soln->lp*soln->bdy_());
	if (g==1.) soln->pcritpt_func = &Solution::isothermal_critpt;
	else soln->pcritpt_func = &Solution::polytropic_critpt;
	soln->bracket_roots();
	soln->root_bisection();
	if(g==1.) soln->remove_iso_vel_minima_roots(); 
	
	//recursively tighten the bounds 
	xc = soln->xcs[0];
    if (xc == 0.) lower = mid;
    else upper = mid;
    mid = 0.5*(upper+lower);
    diff = 0.5*(upper-lower);
    #ifdef DEBUG_PSURVEY
    cout << "Lower= " << lower << "\tUpper= " << upper << endl;
    #endif
  }
    
  if (xc == 0.) return upper;
  else return lower;
}

/*------------------------------------------------------------------------------
 * Print parameter survey file
 *----------------------------------------------------------------------------*/
void ParameterSurvey::print_psurvey()
{
  ofstream inflowcurves,outflowcurves;
  string name;
  
  name = soln->filepath + soln->classtag();
  name += "psurvey_inflow_g" + to_string<double>(soln->g);
  if (soln->classtag() == "BONDI") name += "_at" + to_string<double>(outer_bdy);
  if (soln->i_deg != 0.) name += "_i" + to_string<double>(soln->i_deg/(DEG));
  if (kep) name += "Kep";
  else if (soln->zeta != 0.) name += "_rot" + to_string<double>(soln->zeta);
  name += ".tab";
  inflowcurves.open(name.c_str(), ios::trunc);
      if (soln->g == 1.) 
      {
  				inflowcurves<<setw(10) << left << "HEP= " << 
				setw(10) << left << "\txc_1= " << 
				setw(10) << left << "\txc_2= " <<
				setw(10) << left << "\tMo_1= " <<
				setw(10) << left << "\tMo_2= " << endl;
	  }
	  else 
	  {
				inflowcurves<<setw(10) << left << "HEP= " << 
				setw(10) << left << "\txc_1= " << 
				setw(10) << left << "\txc_2= " <<
				setw(10) << left << "\tMo_1= " <<
				setw(10) << left << "\tMo_2= " << 
				setw(10) << left << "\tdo_1= " <<
				setw(10) << left << "\tdo_2= " << 
				setw(10) << left << "\tlc_1= " <<
				setw(10) << left << "\tlc_2= " << endl;
	  }
  for (int k=0; k < HMAX; k++) {
    if(psurvey[0].xc_inflow[k] != 0. || psurvey[1].xc_inflow[k] != 0.) {
      if (soln->g == 1.) 
      {
      			inflowcurves << setprecision(10) << 
				setw(10) << left << hep[k] << "\t" <<
				setw(10) << left << psurvey[0].xc_inflow[k] << "\t" << 
				setw(10) << left << psurvey[1].xc_inflow[k] << "\t" <<
				setw(10) << left << psurvey[0].Mo_inflow[k] << "\t" << 
				setw(10) << left << psurvey[1].Mo_inflow[k] << "\t" << endl;
	  }
	  else 
	  { 
	  			inflowcurves << setprecision(10) << 
				setw(10) << left << hep[k] << "\t" <<
				setw(10) << left << psurvey[0].xc_inflow[k] << "\t" << 
				setw(10) << left << psurvey[1].xc_inflow[k] << "\t" <<
				setw(10) << left << psurvey[0].Mo_inflow[k] << "\t" << 
				setw(10) << left << psurvey[1].Mo_inflow[k] << "\t" <<
				setw(10) << left << psurvey[0].do_inflow[k] << "\t" << 
				setw(10) << left << psurvey[1].do_inflow[k] << "\t" <<
				setw(10) << left << psurvey[0].lc_inflow[k] << "\t" << 
				setw(10) << left << psurvey[1].lc_inflow[k] << "\t" << endl;
	  }
	}
  }
  inflowcurves.close();
  
  name = soln->filepath + soln->classtag();
  name += "psurvey_outflow_g" + to_string<double>(soln->g);
  if(soln->classtag() == "BONDI") name += "_at" + to_string<double>(outer_bdy);
  if (soln->i_deg != 0.) name += "_i" + to_string<double>(soln->i_deg/(DEG));
  if (kep) name += "Kep";
  else if (soln->zeta != 0.) name += "_rot" + to_string<double>(soln->zeta);
  name += ".tab";
  outflowcurves.open(name.c_str(), ios::trunc);
      if (soln->g == 1.) 
      {
  				outflowcurves<<setw(10) << left << "HEP= " << 
				setw(10) << left << "\txc_1= " << 
				setw(10) << left << "\txc_2= " <<
				setw(10) << left << "\tMo_1= " <<
				setw(10) << left << "\tMo_2= " << endl;
	  }
	  else 
	  {
				outflowcurves<<setw(10) << left << "HEP= " << 
				setw(10) << left << "\txc_1= " << 
				setw(10) << left << "\txc_2= " <<
				setw(10) << left << "\tMo_1= " <<
				setw(10) << left << "\tMo_2= " << 
				setw(10) << left << "\tdo_1= " <<
				setw(10) << left << "\tdo_2= " << 
				setw(10) << left << "\tlc_1= " <<
				setw(10) << left << "\tlc_2= " << endl;
	  }
  for (int k=0; k < HMAX; k++) {
    if(psurvey[0].xc_outflow[k] != 0. || psurvey[1].xc_outflow[k] != 0.) {
      if (soln->g == 1.) 
      {
      			outflowcurves << setprecision(10) << 
				setw(10) << left << hep[k] << "\t" <<
				setw(10) << left << psurvey[0].xc_outflow[k] << "\t" << 
				setw(10) << left << psurvey[1].xc_outflow[k] << "\t" <<
				setw(10) << left << psurvey[0].Mo_outflow[k] << "\t" << 
				setw(10) << left << psurvey[1].Mo_outflow[k] << "\t" << endl;
	  }
	  else 
	  { 
	  			outflowcurves << setprecision(10) << 
				setw(10) << left << hep[k] << "\t" <<
				setw(10) << left << psurvey[0].xc_outflow[k] << "\t" << 
				setw(10) << left << psurvey[1].xc_outflow[k] << "\t" <<
				setw(10) << left << psurvey[0].Mo_outflow[k] << "\t" << 
				setw(10) << left << psurvey[1].Mo_outflow[k] << "\t" <<
				setw(10) << left << psurvey[0].do_outflow[k] << "\t" << 
				setw(10) << left << psurvey[1].do_outflow[k] << "\t" <<
				setw(10) << left << psurvey[0].lc_outflow[k] << "\t" << 
				setw(10) << left << psurvey[1].lc_outflow[k] << "\t" << endl;
	  }
	}
  }
  outflowcurves.close();
}

