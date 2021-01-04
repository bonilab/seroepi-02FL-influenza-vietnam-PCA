#ifndef PERSON
#define PERSON

#include "essentials.h"

#include <vector>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h> // random number generators from Gnu Scientific Library
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>  // this includes the gsl_ran_gaussian function


using namespace std;

extern gsl_rng *G_RNG;	

class person
{   
public:    
    explicit person();          	// constructor
    ~person();         	                        // destructor


    double age; // in years
    int last_infection_week;
    bool is_in_waning_period;
    bool oas_option;            // when this is true, 'original antigenic sin' is set to 'on', and your first 
                                // acquired titer is multiplied by the value 'oas_factor'
                                
                                
    
    static double beta;         // coefficient in logistic regression for probability of infection as a function of titer
    static bool boosting_is_on;
    static double boost_factor;
    
    static double log2_annual_waning;   // this is the 1-year expected antibody waning on a log2 scale
                                        // take 2^log2_annual_waning, which is 7.26 and this means that after
                                        // one your your titer is reduced about 7-fold
    
    
    // this is the antibody titer, on a scale of 10 to 1280 (or maybe 1810)
    // for 18 different influenza antigens (in order) from 0 to 17 (NUMANTIGENS = 18)
    vector<double> v_titer;

    // this is the antibody titer FLOOR, i.e. the lowest allowable titer for this particular antigen
    // for 18 different influenza antigens (in order) from 0 to 17 (NUMANTIGENS = 18)
    vector<double> v_titer_floor;
    
    // here, you pass the week of infection in, i.e. w=0, w=52, w=500, up to, w=1871
    // and the funtion gives you the probability that this person is infected if challenged
    double get_infection_probability( int week );

    
    // here, a person gets "selected for infection", or really, they get selected for having been
    // infected three weeks ago.  And, in this function below, we raise their antibody titer to their
    // specific infection; and we also raise their antibody titers to neighboring infections
    // according to the cross-reaction parameer sigma, which has to be between [0,1]
    //
    // also: if you pass in a value for the oas_factor *and* if oas_option is set to true, then the host's
    // infection will have a raised antibody titer compared to normal titers.
    void generate_abtiter_after_infection( int week, double sigma, double oas_factor=1.0 );

    
    // call this function once a week to update the person's antibody titers
    // so that they wane a little bit
    // first check the above 'is_in_waning_period' boolean; if it is false, there is 
    // no need to call this function
    void wane_all_antibody_titers( int week );
    

    
    // randomly determines if this individual dies; if he does, he is reincarnated as a 0-year old
    void randomly_die_or_not( void );
    
    
    
    
};

#endif // PERSON
