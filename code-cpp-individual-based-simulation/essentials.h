#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <stdlib.h>
#include <ctime>
#include <sys/time.h>

#include <gsl/gsl_rng.h> // random number generators from Gnu Scientific Library
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include <set>
#include <map>
#include <vector>
#include <algorithm>


#define POPSIZE 500000

// every 2 years, there's a new flu strain (antigen) that's just a little different from the
// one two years ago.  This is enough for 36 years of flu antigens
#define NUMANTIGENS 18


//#define HAVE_INLINE 	// this will make GSL vector operations faster


#ifndef ESSENTIALS
#define ESSENTIALS

using namespace std;

//
//
// this function contains the ode system
//
int func(double t, const double y[], double f[], void *params);


//void* jac;	// do this for C-compilation
//
// for C++ compilation we are replacing the void* declaration above with
// the inline dummy declaration below
inline int jac(double a1, const double* a2, double* a3, double* a4, void* a5)
{
    return 0;	
};


inline void write_to_file( const char* szFilename, vector< vector<double> >& vvDATA )
{
    FILE* fp = fopen( szFilename, "w" );
    int nr = vvDATA.size();	// number of rows
    int nc = vvDATA[0].size();
    
   for(int rr=0;rr<nr;rr++)
    {
	for(int cc=0;cc<nc;cc++)
	{
	    fprintf(fp, "%1.3f \t", vvDATA[rr][cc] );	
	}
	fprintf(fp, "\n");
    }

    fclose(fp);
    return;
}


inline int get_day_number( int d, int m, int y )
{
    int basedays;
    switch(y)
    {
      case 2009:
	basedays = 0; break;
      case 2010:
	basedays = 365; break;
      case 2011:
	basedays = 730; break;
      case 2012:
	basedays = 1095; break;
      case 2013:
	basedays = 1461; break;
      case 2014:
	basedays = 1826; break;
      case 2015:
	basedays = 2191; break;
      case 2016:
	basedays = 2556; break;
      default:
	assert(false);
    }
    
    int numleapdays = 0;
    if( y==2012 || y==2016 || y==2020 ) numleapdays = 1;
    
    int monthdays = 0;
    switch(m)
    {
      case 1:
	monthdays = 0; break;
      case 2:
	monthdays = 31; break;
      case 3:
	monthdays = 59+numleapdays; break;
      case 4:
	monthdays = 90+numleapdays; break;
      case 5:
	monthdays = 120+numleapdays; break;
      case 6:
	monthdays = 151+numleapdays; break;
      case 7:
	monthdays = 181+numleapdays; break;
      case 8:
	monthdays = 212+numleapdays; break;
      case 9:
	monthdays = 243+numleapdays; break;
      case 10:
	monthdays = 273+numleapdays; break;
      case 11:
	monthdays = 304+numleapdays; break;
      case 12:
	monthdays = 334+numleapdays; break;
      default:
	assert(false);
    }
  
    return basedays + monthdays + d;
}

inline double besttiter( double ti_4f, double ti_2f, double tiv_4f, double tiv_2f )
{
    bool bFourFoldFailed = false;
    double bt = -99.0;
    
    if( tiv_4f > 0.5 )
    {
        if( ti_4f > 0.0 )
	{
	    bt = ti_4f;
	}
	else
	{
	    bFourFoldFailed = true; 
	}
    }
    
    if( tiv_4f <= 0.5 || bFourFoldFailed )
    {
        if( tiv_2f > 0.5 && ti_2f > 0.0 )
	{
	    bt = ti_2f; 
	}
    }
    
    return bt;
}







#endif // ESSENTIALS
