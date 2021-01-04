#include "person.h"

// Initialize static class members
double person::beta = 1000000.0;
bool person::boosting_is_on = false;
double person::boost_factor = 1.0;
double person::log2_annual_waning = 2.86;

// constructor
person::person( )
{
    age = 0.0;
    for(int i=0; i < NUMANTIGENS; i++)
    {
        v_titer.push_back( 10.0 );
        v_titer_floor.push_back( 10.0 );
    }
    is_in_waning_period = false;
    oas_option = false;
    last_infection_week = -1;
}


// destructor
person::~person()
{

    
}

double person::get_infection_probability( int week )
{
    // you pass in the global week, from 0 to 1871
    // this will give the year index from 0 to 35
    int yr = week/52;
    
    int ag_index = yr/2; // and this will give the antigen index from 0 to 17
    
    assert( ag_index >= 0 && ag_index <= 17 );
    
    // the antibody titer to this antigen (the currently circulating virus in this week)
    double ttr = v_titer[ ag_index ];
    double logttr = log2(ttr/10.0);
    
    // function from logistic regression showing infection prob as a function of titer
    double c = 11.8698; // from 10FL inference, using 2009 H3 titer against 2009 H3 infection
    assert( beta < 0.0 && beta > -2.0 );
    double prob = c*exp(beta*logttr) / (1.0 +  c*exp(beta*logttr));
    
    return prob;
}


// third argument defaults to 1.0 if nothing is passed in
void person::generate_abtiter_after_infection( int week, double sigma, double oas_factor )
{
    /*if( v_titer.size() != 18 )
    {
        fprintf( stderr, "\n\tv_titer length is %d \n\n", (int)v_titer.size() ); fflush(stderr);   
    }*/
    assert( v_titer.size() == 18 );
    
    bool is_their_first_infection = last_infection_week < 0 ? true : false; 
    
    last_infection_week = week;
    is_in_waning_period = true;
    
    // you pass in the global week, from 0 to 1871
    // this will give the year index from 0 to 35
    int yr = week/52;
    
    int ag_index = yr/2; // and this will give the antigen index from 0 to 17
    
    assert( ag_index >= 0 && ag_index <= 17 );
    
    // draw a random number to get the antibody titer for the corresponding antigen
    
    // from Stacy's 10FL study, we know that the mean antibody generated after infection
    // is 7.72 (on a log2 scale) with a standard deviation of 1.311
    double logttr = 7.72 + gsl_ran_gaussian( G_RNG, 1.311 );
    if( logttr<0.0 ) logttr = 0.0;
    double ttr = 10.0 * pow(2.0, logttr);
    
    // if it's their first infection, then via original antigenic sin, make
    // this titer larger so that even the waned titer many years later looks large
    if( is_their_first_infection ) ttr *= oas_factor;
    
    if( boosting_is_on )
    {
        // loop through all past antigens whose titers are above 50 (arbitrarily) so that they're definitely non-zero
        // and boost them up by a factor 'boost_factor' .. note that if boost_factor > 2.86 then these boosted
        // antibodies will not wane back to their previous level, and there will be some permanent boosting effect
        for( int ag=0; ag < ag_index; ag++ )
        {
            if( v_titer[ag] > 50.0 ) v_titer[ag] *= boost_factor;
        }
    }

    // one-year waning drop should be 7.26x
    double oyw = pow(2.0, log2_annual_waning);
    // assert( oyw > 7.26 && oyw < 7.27 );
    
    //NOTE  - in the loop below, we may overwrite some of the boosted titers above
    //      - but only if we are making them higher, not lower
    
    // now assign titers to the primary antigen and to neighboring antigens
    // loop through all neighboring antigens from 3-to-the-left to 3-to-the-right
    for( int k=-3; k<=3; k++ )
    {
        double T = ttr;
        if( k==-3 || k==3 ) T = ttr*sigma*sigma*sigma >= 10.0 ? ttr*sigma*sigma*sigma : 10.0;
        if( k==-2 || k==2 ) T = ttr*sigma*sigma >= 10.0 ? ttr*sigma*sigma : 10.0;
        if( k==-1 || k==1 ) T = ttr*sigma >= 10.0 ? ttr*sigma : 10.0;
        
        // now we have the new candidate titer T that we want to update the person's profile with
        // 1. we need to check that T is above the current titer (if it's not, no need to update anything)
        // 2. then, we need to check if we need to change (increase) the titer floor
        
        // check to make sure the index is in-bounds
        if( ag_index+k >= 0 && ag_index+k < NUMANTIGENS) 
        {
            // check that new titer (resulting from infection) is higher than pre-infection titer
            if( T > v_titer[ag_index+k] )
            {
                v_titer[ag_index+k] = T;
                
                // now check if we need to increase the floor
                double candidate_floor = T/oyw;
                if( candidate_floor > v_titer_floor[ag_index+k] ) v_titer_floor[ag_index+k] = candidate_floor;
                
            }
        }
    }
    
}






void person::wane_all_antibody_titers( int week )
{
    assert( v_titer.size() == 18 );
    
    // this is here so that a newly generated antibody titer is not waned
    if (week==last_infection_week) 
    {
        return;
    }
    
    // antibodies only wane for the first year
    // then they stay constant after that
    if( week - last_infection_week > 52 )
    {
        is_in_waning_period = false;
        return;
    }
    
    // from 10FL waning is 2.86 log2 units per year, or 0.054849315 per week
    // that means that we multiply the titer by the factor 0.96269498834 every week
    // TODO: do this calculation on the fly based on person::log2_annual_waning 
    
    for(int ag_index=0; ag_index < NUMANTIGENS; ag_index++)
    {
        v_titer[ag_index] = v_titer[ag_index] * 0.96269498834;
        if( v_titer[ag_index] < 10.0 ) v_titer[ag_index] = 10.0;
    }
    
}






void person::randomly_die_or_not( void )
{
    
    // probability of dying in one week
    /*double weekly_mortality[18] = {    
    0.0000254396,               // 0-4 age group
    0.0000092676,               // 5-9 age group
    0.0000188466,
    0.0000188466,
    0.0000287374,
    0.0000452329,
    0.0000485341,
    0.0000452329,
    0.0000644307,
    0.0001128153,
    0.0001420276,
    0.0002063076,
    0.0003035913,
    0.0004450196,
    0.0007499847,
    0.0010840304,
    0.0017896769,
    0.0035064367 };*/
    
    
    
    
    double four_weekly_mortality[18] = { 
    0.0001017544,           // 0-4 age group
    0.0000370698,
    0.0000753841,
    0.0000753841,
    0.0001149446,           // 20-24 age group
    0.0001809193,
    0.0001941221,
    0.0001809193,
    0.0002576981,           // 40-44 age group
    0.0004511848,
    0.0005679893,
    0.0008249750,
    0.0012138124,           // 60-64 age group
    0.0017788903,
    0.0029965658,
    0.0043290758,
    0.0071395130,           // 80-84 age group
    0.0139521487 };

    int ind = (int) (age/5.0);
    if( ind > 17) ind=17;
    
    // flip a coin to see if this person dies in this 4-week period
    double rv = gsl_ran_flat(G_RNG, 0.0, 1.0);
    
    if( rv < four_weekly_mortality[ind] ) // the person dies, so you need to re-incarnate them
    {
        // update their age
        age = 0.0;
        
        // zero-out their antibody titers
        for(int i=0; i < NUMANTIGENS; i++)
        {
            v_titer[i] = 10.0;
        }
        
        // reset their last infection week
        last_infection_week = -1;
        
        // flag them as not in a waning period
        is_in_waning_period = false;
        
    }

    
}













