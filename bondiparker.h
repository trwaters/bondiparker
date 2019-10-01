#include <cmath>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>
#include <time.h>                                        		

/* Parameters to adjust if running the ParameterSurvey class */
#define HEP_STEP_SIZE 1.	//HEP step size in ParameterSurvey
#define HMAX 1000			//number of HEP_STEP_SIZE steps in ParameterSurvey

/* Parameters controlling default position grid - ignored if a position array is inputted */
#define NMAX 1024			//number of points used for distance variable
#define LOG_SPACING			//define to use logarithmic spacing in distance variable
#define LOGXMIN -5.			//minimum distance is 10^LOGXMIN

/* Parameters used by root finders */	
#define UPPER_LIMIT	1e2		//max crit pt distance for winds (also max Mach number) 
#define MAX_ROOTS 4			//max number of critical points to bracket
#define TOL 1e-11  			//tolerance for root finding

/* Constants */
#define PI 3.14159265358979
#define DEG PI/180. 		//conversion factor: angles must be input in degrees! 
#define TINY_NUMBER 1e-13

typedef std::vector<double> array_t;  //dynamic array data type

 		
typedef struct
{

  double  *subsonic,
          *supersonic,
          *transonic_inflow,
          *transonic_outflow; 
} Profiles;
  

typedef struct 
{
  double	xc_inflow[HMAX],			//xc as a function of HEP for accretion solutions
  			xc_outflow[HMAX],			//xc as a function of HEP for wind solutions
  			lc_inflow[HMAX],			//corresponding critical point paramters
  			lc_outflow[HMAX],			 
			Mo_inflow[HMAX],			//corresponding initial Mach numbers
			Mo_outflow[HMAX],			
			do_inflow[HMAX],			//corresponding initial densities
			do_outflow[HMAX];
} HEPcurves;

class ParameterSurvey;



/* An abstract base class */
class Solution {
     			
  public: 
    Solution (double a, double b, double c=0., double d=0., double e=1.) :
     			lp(a),g(b),i_deg(c),zeta(d),L(e) {}
    void run(const int,double*,bool,std::string); 	//ties together all functions and outputs a complete solution
	std::string gravity;	//allows switching between various gravitational potentials
  double outer_bdy, rho_bc; 
	

	/* solution variables */ 
	Profiles M[MAX_ROOTS],					//struct of Mach Number profiles for each soln
  			 v[MAX_ROOTS],					//struct of velocity profiles for each soln
  			 d[MAX_ROOTS],					//struct of density profiles for each soln
  			 P[MAX_ROOTS],					//struct of energy profiles for each soln
  			 econst[MAX_ROOTS];				//struct of critical pt constant for each soln
  	double*	x_array; 					//distance along a streamline
  	int		ids[MAX_ROOTS],					//ID inflow vs. outflow transonic solutions
  			nroots,							//number of roots found
        n_inflow,           // index of the transonic inflow solution
        n_outflow,          // index of the transonic outflow solution
        N;                  // length of the position/solution arrays
    double xcs[MAX_ROOTS],         //critical points
           lcs[MAX_ROOTS];         //lamda_c's

    /* compute functions */
    double getMo(const double, const double);
    double getv0(const double, const double);
    double getmdot(const double, const double);
    double getdo(const double, const double);
  
  protected:
  
    /* input parameters */
	double	lp,								//HEP
			g,								//Gamma
			i_deg,							//inclination angle
			zeta,							//rotation parameter
			L;								//solution eigenvalue; L=1.0 for transonic
			
	/* derived parameters */
	double cosi;
	
	/* critical point parameters */
	double	lc,								//lambda_c
			xc,								//critical point location
			ec;								//critical point constant
	
	/* solution variables */ 				 
  	double	xl[MAX_ROOTS],					//left root brackets
  			xr[MAX_ROOTS];					//right root brackets

    
    
    int id,iter,nbrackets;
	

  /* Memory handling */
  void initialize(const int, double*);

	/* For output */
	std::string filename,filepath;
	virtual std::string classtag() = 0;
    std::string identifyModelName();
    std::string default_filename(int,std::string);
    void clearVariables();
    void cout_output(clock_t);
    void print_Xtype_solns();
    void print_Ttype_solns();
    
    /* Function pointers */
    typedef double (Solution::*p1ArgFunc)(const double);
    typedef void (Solution::*p0ArgFunc)();
    p1ArgFunc pcritpt_func,
    		  pUeff,
    		  pFeff;
    p0ArgFunc pmach_func; 
    
    /* Compute functions */
    void bracket_roots();
    void root_bisection();
    void solutionProfiles(double*);
    double get_subsonic_MachNumber(const double);
    double get_supersonic_MachNumber(const double);
    double mach_bisection(const double x, const double M);
    double mach_func(const double x, const double M);
    double mach_lambertWfct(const double, int);
    double density_lambertWfct(const double, const double, const double);
    double polytropic_critpt(const double x);
    void getFlowVariables(const double M, const double x, double* v, double* d, double* E, double* B);
    double isothermal_critpt(const double x);
    void remove_iso_vel_minima_roots();
    
    /* Allow for an arbitary density at the boundary */
    // void setDensityBC(const double value);

    /* Effective potentials and corresponding forces */
    void setPotential(std::string);
    double newtonianUeff(const double);
    double newtonianFeff(const double);
    double pseudoNewtonianUeff(const double);
    double pseudoNewtonianFeff(const double);
    
    /* Model-specific functions */
    virtual double bdy_() = 0;						//returns the location of the boundary
    virtual double area_(const double x) = 0; 		//area function of a given model
    virtual double area_prime_(const double x) = 0; //derivative of area function
    
    /* Allow function access to ... */
    friend class ParameterSurvey;
};

class CIAmodel : public Solution {
  public:
    CIAmodel (double a, double b, double c) : Solution(a,b,c*DEG,sqrt(a)){}
    CIAmodel (double a, double b, double c, double d) : Solution(a,b,c*DEG,d){}
    CIAmodel (double a, double b, double c, double d, double e) : Solution(a,b,c*DEG,d,e){}
  
  protected:
    double bdy_ () {return 0.;}
    double area_ (const double x) {return (1./lp + x*cosi); }
    double area_prime_ (const double x) { return cosi; }
    std::string classtag() {return "CIA";}
    
    friend class ParameterSurvey;
    
  private:
    double ro;
};

class CONmodel : public Solution {
  public:
    CONmodel (double a, double b, double c) : Solution(a,b,c*DEG,sqrt(a)){}
    CONmodel (double a, double b, double c, double d) : Solution(a,b,c*DEG,d){}
    CONmodel (double a, double b, double c, double d, double e) : Solution(a,b,c*DEG,d,e){}

  protected:
    double bdy_ () {return 0.;}
    double area_ (const double x) {ro = 1./lp; return (ro + x*cosi)*(ro + x*cosi); }
    double area_prime_ (const double x) {ro = 1./lp; return 2.*(ro + x*cosi)*cosi; }
    std::string classtag() {return "CON";}
    
    friend class ParameterSurvey;
  
  private:
    double ro;
};

class BONDImodel : public Solution {
  public:
    BONDImodel (double a, double b) : Solution(a,b,0.,0.){}
    BONDImodel (double a, double b, double c) : Solution(a,b,c*DEG,sqrt(a)){}
    BONDImodel (double a, double b, double c, double d) : Solution(a,b,c*DEG,d){}
    BONDImodel (double a, double b, double c, double d, double e) : Solution(a,b,c*DEG,d,e){}

  protected:
    double bdy_ () {return outer_bdy;}
    double area_ (const double x) {ro = 1./lp; return (ro + x*cosi)*(ro + x*cosi); }
    double area_prime_ (const double x) {ro = 1./lp; return 2.*(ro + x*cosi)*cosi; }
    std::string classtag() {return "BONDI";}
    
    friend class ParameterSurvey;
  
  private:
    double ro;
};

class ParameterSurvey {
  public:  
    ParameterSurvey (bool a, bool b, bool c) : kep(a), useMin(b), outer_bdy(c) {} 
    Solution *soln;
    void run(std::string);
    void query(); //use to query if a solution exists at a specific HEP value
    
  protected:
    
    /* parameter survey variables */
    double outer_bdy;
    double hep[HMAX];				//HEP (horizontal axis)
	HEPcurves psurvey[MAX_ROOTS];	//parameter survey variables 
    bool kep; 						//true -> Keplerian rotation on
    bool useMin;					//true -> use minimum HEP in survey
    
    /* compute functions */
    void doSurvey(double);
    double findMinHEP();
    
    /* print function */
    void print_psurvey();
    
  private:
    double testhep;
};

template <class T>
inline std::string to_string (const T& t)
{
   std::stringstream ss;
   ss << t;
   return ss.str();
}
