#include "essentials.h"
#include "readdata.h"
#include "person.h"



// GLOBAL RANDOM NUMBER GENERATOR
gsl_rng *G_RNG;		

// GLOBAL VARIABLES
//int G_CLO_NUMITER = 5000; 	     // number of iterations in the Nelder-Mead minimizer routine
//bool G_CLO_PROFILE = false;     // determines whether a likelihood profile should be done (in addition to the MLE)
//int G_CLO_NUM_IC = 10;          // number of initial conditions to try in NM-searches; i.e. number of Nelder-Mead
                                // searches that will be done

                                


// cross protection parameter; you generate a titer of 1.0 to the strain you were infected with
// and sigma to the two neighboring strains, and sigma^2 to the next two neighboring strains, and sigma^3 after that
double G_CLO_SIGMA = 0.88;      // if you just do the simple averaging of the four Pearson correlation coefficients
                                // of the neighboring pairs of H3 strains, you get 0.88345

int G_CLO_AR = 25;              // this is the annual attack rate, 25 means 25% .... restricted to be an interger because
                                // we're not goint to try and differentiate between a 19% and 19.5% attack rate
                                
int G_CLO_REP = 1;           

double G_CLO_NEGBETA = 0.8;     // negative of coefficient in logistic regression
                                // that detemines shallowness of prob-infection-by-titer curve
                                // but remember that this curve is only valid for SYMPTOMATIC infection (beta = - 0.7797)
                                //
                                // for any infection (symptomatic or asymptomatic) neg_beta should be set to something
                                // between 0.1 and 0.3 (probably ... as this is not measured in the field typically)
                                // 
                                
                                
bool G_CLO_OVERSAMPLE_CHILDREN = false;
bool G_CLO_OAS = false;                 // original antigenic sin
bool G_CLO_BOOSTING = false;            //
double G_CLO_BOOSTINGFACTOR = 1.0;      //


                                

// GLOBAL STRUCTURES THAT HOLDS DATA
// vector< vector<double> > G_02FL_DATA;
// vector< vector<double> > G_10FL_DATA;
// ssedata G_SSEDATA;

// GLOBAL STRUCTURE THAT HOLDS THE PARTICULAR 10FL DATA WE ARE ANALYZING NOW
// map< int,vector<double> > G_10FL_SUBSET;
// int n; // this is the number of patients in the above subset

person GPOP[POPSIZE];
int GINDICES[POPSIZE];
int GINDICES_SHUFFLED[POPSIZE];

int GWEEKNUM = 0; // this is the global week number starts at 0, and it will go up to 36*52 - 1 = 1871

// globalpointer to the output file
FILE* GFOUT;


vector< vector<double> > GPERCENTINFECTED;  // this matrix has 9*52 = 468 weeks of infection


// FUNCTION DECLARATIONS
void ParseArgs(int argc, char **argv);
void initialize_population( void );
void infect_individuals_in_population( void );
void write_sample_to_outfile( int n );
void write_sample_to_outfile2( int n );


int main(int argc, char* argv[])
{
    ParseArgs(argc,argv);
    //assert( DIM==131 );
    
    // get random number seed from current time
    struct timeval pTV[1];
    gettimeofday( pTV, NULL );
    int seed = ((int)pTV->tv_sec) + ((int)pTV->tv_usec);  // this adds seconds to microseconds

    // make random number generator (RNG) the Mersenne Twister which has period 2^19937 - 1
    const gsl_rng_type *TT_RAND = gsl_rng_mt19937;
    G_RNG = gsl_rng_alloc(TT_RAND);
    gsl_rng_set( G_RNG, seed ); // seed the RNG    
    
    // read in the data that show what percent of the population each week is infected
    string fn("./vn_sentinel_data/ScaledWeeklyH3Cases_SouthernVN_25percentAR.tdl");
    readdata( fn , GPERCENTINFECTED, 3, false );
    
    // construct the output filename based on the parameters
    char filename[100];
    int b1 = (int)G_CLO_NEGBETA;
    int b2 = ((int)(10.0*G_CLO_NEGBETA)) % 10;
    sprintf(filename, "../code-matlab/B_Simulated/BMatrix.AR%02d.Sigma0p%d.nb%dp%d.rep%02d.csv", G_CLO_AR, (int)(G_CLO_SIGMA*100.0), b1, b2, G_CLO_REP);

    // open the file and write the header
    GFOUT = fopen( filename, "w" );  // fopen ( const char * filename, const char * mode );
    fprintf(GFOUT, "TTR_H1_1918 , TTR_H1_1977 , TTR_H1_1999 , TTR_H1_2007 ,  TTR_H1_2009 , TTR_H3_1968 , TTR_H3_2003 , TTR_H3_2005 , TTR_H3_2007 , TTR_H3_2009 , TTR_H3_2011, AGE , BIRTHYEAR , LOCATION_CODE , YEAR , 02FL_SAMPLE_ ID \n");


    // populate the indices array, one can be shuffled, the other is fixed
    for(int i=0; i<POPSIZE; i++ ) 
    {
        GINDICES[i]=i;
        GINDICES_SHUFFLED[i]=i;
    }
    
    // initialize the age in the population
    initialize_population();
    
    // initialize the beta value (a static member of person) that determines probability of infection
    person::beta = - G_CLO_NEGBETA;
    
    // initialize the back-boosting status in the population if this was selected on the command line
    person::boosting_is_on = G_CLO_BOOSTING;
    person::boost_factor = G_CLO_BOOSTINGFACTOR;
    
    
    // initialize the vector that will hold the sampling weeks
    vector<int> v_sampling_weeks;
    v_sampling_weeks.push_back(1150); // 1150 is in February 2010 (it is the 6th week of 2010)
    v_sampling_weeks.push_back(1159); // April 2010
    v_sampling_weeks.push_back(1167); // June 2010
    v_sampling_weeks.push_back(1176); // August 2010
    v_sampling_weeks.push_back(1185); // October 2010
    v_sampling_weeks.push_back(1193); // December 2010
    
    v_sampling_weeks.push_back(1202); // February 2011
    v_sampling_weeks.push_back(1211);
    v_sampling_weeks.push_back(1219);
    v_sampling_weeks.push_back(1228);
    v_sampling_weeks.push_back(1237);
    v_sampling_weeks.push_back(1245);
    
    v_sampling_weeks.push_back(1254); // February 2012
    v_sampling_weeks.push_back(1263);
    v_sampling_weeks.push_back(1271);
    v_sampling_weeks.push_back(1280);
    v_sampling_weeks.push_back(1289);
    v_sampling_weeks.push_back(1297);

    v_sampling_weeks.push_back(1306); // February 2013
    v_sampling_weeks.push_back(1315);
    v_sampling_weeks.push_back(1323);
    v_sampling_weeks.push_back(1332);
    v_sampling_weeks.push_back(1341);
    v_sampling_weeks.push_back(1349);

    v_sampling_weeks.push_back(1358); // February 2014
    v_sampling_weeks.push_back(1367);
    v_sampling_weeks.push_back(1375);
    v_sampling_weeks.push_back(1384);
    v_sampling_weeks.push_back(1393);
    v_sampling_weeks.push_back(1401);

    v_sampling_weeks.push_back(1410); // February 2015
    v_sampling_weeks.push_back(1419);
    v_sampling_weeks.push_back(1427);
    v_sampling_weeks.push_back(1436);
    v_sampling_weeks.push_back(1445);
    // v_sampling_weeks.push_back(1453); // no sampling was done in December 2015
    
    
    vector<int>::iterator it = v_sampling_weeks.begin();
    int num_sampling_weeks = v_sampling_weeks.size();
    
    printf("\n\t Num sampling weeks is %d \n", num_sampling_weeks );
    
    /*printf("\n\t Currently pointing to week %d \n", *it );
    it++;
    printf("\n\t Currently pointing to week %d \n", *it );*/

    
    // make sure that the difference in sampling weeks is always 8 or 9
    for( int i=1; i<num_sampling_weeks; i++ )
    {
        assert( v_sampling_weeks[i]-v_sampling_weeks[i-1] == 8 || v_sampling_weeks[i]-v_sampling_weeks[i-1]==9 );
    }

    
//     
//     ### TESTED ### THIS WORKS ###
//     for(int r=0; r<400; r++)
//     {
//         for(int c=0; c<3; c++)
//         {
//             printf("\t %1.4f", GPERCENTINFECTED[r][c] );   
//         }
//         printf("\n");
//     }
//


    
    //BEGIN MAIN LOOP
    for(GWEEKNUM=0; GWEEKNUM<=1871; GWEEKNUM++) // 1872 weeks is 36 years
    {
        
        // ### 1 ### select people to be infected -- just increase their antibody titers! -- no actual infection or transmission
        infect_individuals_in_population();
        
        
        
        // ### 2 ###  wane antibodies for everyone *except* the individuals just infected above
        for(int i=0; i<POPSIZE; i++)
        {
            if( GPOP[i].is_in_waning_period && GPOP[i].last_infection_week != GWEEKNUM ) 
            {
                GPOP[i].wane_all_antibody_titers( GWEEKNUM );
            }
        }
        
        
        
        // ### 3 ###  update everyone's age by one week
        for(int i=0; i<POPSIZE; i++)
        {
            GPOP[i].age += ( 7.0 / 365.0 );
        }
        
        
        
        // ### 4 ###  if it's a sampling week, add to the sampling file
        if( *it == GWEEKNUM )
        {
            write_sample_to_outfile2( 192 ); // this will sample 192 individuals, and then discard those younger than 6 months
            it++;
        }
        
        
        
        // ### 5 ###  do births and deaths - once every 4 weeks is enough
        if( GWEEKNUM%4 == 0 )
        {
            for(int i=0; i<POPSIZE; i++)
            {
                GPOP[i].randomly_die_or_not();
            }
        }
        
        
        
    }
    //END MAIN LOOP

    
    fclose(GFOUT);
    
    printf("\n\n\tSimulation finished.\n\n");
    

    return 0;
}



void infect_individuals_in_population( void )
{
    // NOTE you have to do "mod 468" below because there are only nine years of flu data
    //      so once the simulation gets past nine years, you simply cycle through the infection
    //      data gain
    //
    
    // inflation factor f
    // if the input attack rate is higher or lowe than 25%, then use this factor
    // to adjust the number of individuals to infect in the calculation below.
    double f = ( (double) G_CLO_AR ) / 25.0;
    
    
    // the GPERCENTINFECTED matrix tells you exactly what % of the population will be infected this week
    int num_individuals_to_infect = (int)( ((double) POPSIZE) * f * GPERCENTINFECTED[ GWEEKNUM%468 ][ 2 ] + 0.5);
    
    if( num_individuals_to_infect==0 )
    {
        if( GWEEKNUM%100==0 ) printf("\n\t week=%d \t number of infections: \t %d ", GWEEKNUM , num_individuals_to_infect );
        return;
    }
    else
    {
        if( GWEEKNUM%100==0 ) printf("\n\t week=%d \t number of infections: \t %d ", GWEEKNUM , num_individuals_to_infect );
    }

    int* dest;
    dest = new int[ num_individuals_to_infect ];

    int numleft = num_individuals_to_infect;

    while( numleft > 0 )
    {
        // do a random sample of indices of the individuals we want to 
        // try to infect_individuals_in_population
        
        // int gsl_ran_choose(const gsl_rng * r, void * dest, size_t k, void * src, size_t n, size_t size)
        gsl_ran_choose(G_RNG, dest, num_individuals_to_infect, GINDICES, POPSIZE, sizeof (int));
        
        // then, shuffle these numbers
        // gsl_ran_shuffle(const gsl_rng * r, void * base, size_t n, size_t size)
        gsl_ran_shuffle(G_RNG, dest, num_individuals_to_infect, sizeof (int));
        
        // now loop through all of the indices in the temporary array "dest"
        // and TRY to infect these individuals
        for(int i=0; i<num_individuals_to_infect; i++)
        {
            // if they're already been infected, skip and continue
            if( GPOP[ dest[i] ].last_infection_week == GWEEKNUM ) continue;
            
            // otherwise, try to infect person at index dest[i]
            double p  = GPOP[ dest[i] ].get_infection_probability( GWEEKNUM );   
            double rv = gsl_ran_flat(G_RNG, 0.0, 1.0);
            if( rv < p ) // person is successfully infected
            {
                if( G_CLO_OAS )
                {
                    GPOP[ dest[i] ].generate_abtiter_after_infection( GWEEKNUM, G_CLO_SIGMA, 2.0 ); // third argument is the OAS titer inflation factor
                }
                else
                {
                    GPOP[ dest[i] ].generate_abtiter_after_infection( GWEEKNUM, G_CLO_SIGMA ); // third argument is missing; it is 1.0 by default
                }
                GPOP[ dest[i] ].last_infection_week = GWEEKNUM;
                GPOP[ dest[i] ].is_in_waning_period = true;
                
                numleft--;
                if( numleft==0 ) break;
            }
            
        }
        
    }
    
    
    delete[] dest;
    
}






// initialize population with age structure as in real Vietnamese population
// used age distribution in 2014 with 17 age classes 
void initialize_population()
{
    // indices of dist represent the age classes in 5-year age bands: 0 is 0-4, 1 is 5-9, etc 
    // range: index: index*5 : index*5 + 4
    double dist[17] = { 8.33,
                        7.75,
                        7.39,
                        7.87,
                        9.73,
                        9.17,
                        8.45,
                        7.42,
                        7.02,
                        6.39,
                        5.74,
                        4.58,
                        3.05,
                        2.11,
                        1.61,
                        1.36,
                        2.61};
    
    unsigned int population[17];
    gsl_ran_multinomial(G_RNG, 17, POPSIZE, dist, population);
    
    // NOTE checked - this works
    // for(int i=0; i<17; i++ ) printf("\n\t %d \t %d", i , population[i] );
    
    int count = 0;
    
    // loop through each age band
    for(int i=0; i<17;i++)
    {
        int num_in_ageband = population[i];
        
        // draw a random age for each individual in this age band
        // and assign the age to the main population array GPOP
        for(int j=0; j< num_in_ageband; j++)
        {
            double temp = gsl_ran_flat(G_RNG, 0.0, 4.999999999999);
            GPOP[count].age = ( ((double)i) * 5.0 ) + temp;
            assert( count < POPSIZE );
            count++;
        }
    }

    // gsl_ran_shuffle(const gsl_rng * r, void * base, size_t n, size_t size)
    gsl_ran_shuffle (G_RNG, GPOP, POPSIZE, sizeof(person) );
    
    
}



// generate one bimonthly sample, as in 02FL, of size n which is normally 200, but we will do 190 here
// to account for the samples that had missing titers
void write_sample_to_outfile( int n )
{
    int* dest;
    dest = new int[ n ]; // this will hold the n indices of the individuals that are sampled

    gsl_ran_choose(G_RNG, dest, n, GINDICES, POPSIZE, sizeof (int));
    
    double a; // age placeholder in the loops below
    
    // loop through all sampled indices in the array 'dest'
    for(int i=0; i<n; i++ )
    {
        // age of this individual
        a = GPOP[ dest[i] ].age;
        
        // exclude individuals younger than 6 months
        if( a <= 0.50 ) continue;
        
        // randomly pick a value for the 1968 titer based on age
        // 1. get the expected mean titer by age (normal centered at age 62.23, w sd = 39.43
        double meanttr_byage_1968 = sqrt(6.28) * 39.43 * 737.8 * gsl_ran_gaussian_pdf( (62.23 - a) , 39.43);
        
        // 2. now create some variation around this mean age-based titer, based on empirical data
        double sd = a < 40.0 ? 4.00 : 6.25;
        double rv_deviation_from_zero = gsl_ran_gaussian(G_RNG, sd);
        double ttr_1968 = meanttr_byage_1968 + rv_deviation_from_zero;
        if( ttr_1968 > 1280.0 ) ttr_1968 = 1810.0;
        if( ttr_1968 <   10.0 ) ttr_1968 = 10.0;
        
        
        // the first six below are the 6 H3 titers, and the 7th one is the person's age
        fprintf(GFOUT, " -99 , -99 , -99 , -99 ,  -99 , %1.3f , %1.3f , %1.3f , %1.3f , %1.3f , %1.3f , %1.2f , -99 , -99 , -99 , -99\n",
                ttr_1968 , 
                GPOP[ dest[i] ].v_titer[12] <= 1280.0 ? GPOP[ dest[i] ].v_titer[12] : 1810.0, 
                GPOP[ dest[i] ].v_titer[13] <= 1280.0 ? GPOP[ dest[i] ].v_titer[13] : 1810.0,
                GPOP[ dest[i] ].v_titer[14] <= 1280.0 ? GPOP[ dest[i] ].v_titer[14] : 1810.0,
                GPOP[ dest[i] ].v_titer[15] <= 1280.0 ? GPOP[ dest[i] ].v_titer[15] : 1810.0,
                GPOP[ dest[i] ].v_titer[16] <= 1280.0 ? GPOP[ dest[i] ].v_titer[16] : 1810.0,
                a ); 
    }
    
    delete[] dest;
}






// generate one bimonthly sample, as in 02FL, of size n which is normally 200, but we will do ~190 here
// to account for the samples that had missing titers
void write_sample_to_outfile2( int n )
{
    
    // gsl_ran_shuffle(const gsl_rng * r, void * base, size_t n, size_t size)
    //
    // shuffle the indices, and then just walk through the indices left to right to 
    // select individuals of the appropriate ages
    gsl_ran_shuffle( G_RNG, GINDICES_SHUFFLED, POPSIZE, sizeof(int) );

    double a; // age placeholder in the loops below
    
    int i=0; // index to use in GINDICES_SHUFFLED array
    int j;
    
    int samplesize_under5 = (int) ( (((double)n) * (0.090)) + 0.5 );
    int samplesize_5to10  = (int) ( (((double)n) * (0.056)) + 0.5 );
    int samplesize_10plus = n - samplesize_5to10 - samplesize_under5;
    if( G_CLO_OVERSAMPLE_CHILDREN )
    {
        samplesize_under5 = (int) ( (((double)n) * (0.30)) + 0.5 );
        samplesize_5to10  = (int) ( (((double)n) * (0.30)) + 0.5 );
        samplesize_10plus = n - samplesize_5to10 - samplesize_under5;
    }
    
    //fprintf(stderr, "\n\tSample sizes are %d , %d , %d ", samplesize_under5 , samplesize_5to10 , samplesize_10plus );
    

    // while each age-group sampling need is still larger than zero, keep trying to sample from GINDICES_SHUFFLED
    while( samplesize_under5>0 || samplesize_5to10>0 || samplesize_10plus>0 )
    {
        j = GINDICES_SHUFFLED[i]; // j is the index to use in GPOP array
        
         // age of this individual
        a = GPOP[j].age;       
        
        // exclude individuals younger than 6 months
        if( a <= 0.50 ) { i++; assert( i < POPSIZE ); continue; }
        
        bool bWrite = false;
        
        if( a < 5.0  )
        {
            if( samplesize_under5 > 0 )
            {
                bWrite = true;
                samplesize_under5--;
            }
        }
        else if( a < 10.0 )
        {
            if( samplesize_5to10 > 0 )
            {
                bWrite = true;
                samplesize_5to10--;
            }
        }
        else
        {
            if( samplesize_10plus > 0 )
            {
                bWrite = true;
                samplesize_10plus--;
            }
        }
        
        // if the individual is in an age group that
        // still needs sampling bWrite will be set to true,
        // and the individual will be written to the file
        if( bWrite )
        {
            // randomly pick a value for the 1968 titer based on age
            // 1. get the expected mean titer by age (normal centered at age 62.23, w sd = 39.43
            double meanttr_byage_1968 = sqrt(6.28) * 39.43 * 737.8 * gsl_ran_gaussian_pdf( (62.23 - a) , 39.43);
            
            // 2. now create some variation around this mean age-based titer, based on empirical data
            double sd = a < 40.0 ? 4.00 : 6.25;
            double rv_deviation_from_zero = gsl_ran_gaussian(G_RNG, sd);
            double ttr_1968 = meanttr_byage_1968 + rv_deviation_from_zero;
            if( ttr_1968 > 1280.0 ) ttr_1968 = 1810.0;
            if( ttr_1968 <   10.0 ) ttr_1968 = 10.0;     
            
            // now draw a random variate for the standard deviation of this individual's titers
            double rv_sd = gsl_ran_gaussian(G_RNG, 0.6); // 0.6 is the standard deviation on a log2 titer scale for replicates of the same sample
            double rvexp = pow(2.0, rv_sd);
            
            ttr_1968 = ttr_1968*rvexp;
            if( ttr_1968 > 1280.0 ) ttr_1968 = 1280.0;
            if( ttr_1968 < 10.0 )   ttr_1968 = 10.0;
            
            double ttr_2003 = GPOP[j].v_titer[12] * rvexp;
            if( ttr_2003 > 1280.0 ) ttr_2003 = 1280.0;
            if( ttr_2003 < 10.0 )   ttr_2003 = 10.0;

            double ttr_2005 = GPOP[j].v_titer[13] * rvexp;
            if( ttr_2005 > 1280.0 ) ttr_2005 = 1280.0;
            if( ttr_2005 < 10.0 )   ttr_2005 = 10.0;

            double ttr_2007 = GPOP[j].v_titer[14] * rvexp;
            if( ttr_2007 > 1280.0 ) ttr_2007 = 1280.0;
            if( ttr_2007 < 10.0 )   ttr_2007 = 10.0;

            double ttr_2009 = GPOP[j].v_titer[15] * rvexp;
            if( ttr_2009 > 1280.0 ) ttr_2009 = 1280.0;
            if( ttr_2009 < 10.0 )   ttr_2009 = 10.0;

            double ttr_2011 = GPOP[j].v_titer[16] * rvexp;
            if( ttr_2011 > 1280.0 ) ttr_2011 = 1280.0;
            if( ttr_2011 < 10.0 )   ttr_2011 = 10.0;
            
    
    
            // the first six below are the 6 H3 titers, and the 7th one is the person's age
            fprintf(GFOUT, " -99 , -99 , -99 , -99 ,  -99 , %1.3f , %1.3f , %1.3f , %1.3f , %1.3f , %1.3f , %1.2f , -99 , -99 , -99 , -99\n",
                    ttr_1968 , 
                    ttr_2003 , 
                    ttr_2005 ,
                    ttr_2007 ,
                    ttr_2009 ,
                    ttr_2011 ,
                    a ); 
        }
        
        i++;
        assert( i < POPSIZE );
    }
    
}









// parses command line arguments
void ParseArgs(int argc, char **argv)
{
    int i, start;
    start=1;

    /*if( argc<start )
    { 
        PrintUsageModes(); 
        exit(-1);
    }
        
    if( argv[1][0]=='n' && argv[1][1]=='o' && argv[1][2]=='n' && argv[1][3]=='e' && argv[1][4]==0 )
    {
        //fprintf(stderr, "\n\n\tnot printing to Outfile\n\n");
    }
    else 
    {
        Outfile = fopen( argv[1], "w" );
    }
    
    prm_intro_day = atof( argv[2] );*/
    
    string str;
    for(i=start; i<argc; i++)
    {
        str = argv[i];
             if( str == "-ar" )           G_CLO_AR                      = atoi( argv[++i] );
        else if( str == "-rep" )          G_CLO_REP                     = atoi( argv[++i] );
        else if( str == "-sigma" )        G_CLO_SIGMA                   = atof( argv[++i] );
        else if( str == "-beta" )         G_CLO_NEGBETA                 = atof( argv[++i] );
        else if( str == "-osc" )          G_CLO_OVERSAMPLE_CHILDREN     = true;
        else if( str == "-oas" )          G_CLO_OAS                     = true;
        else if( str == "-boost" )        G_CLO_BOOSTING                = true;
        else if( str == "-boostfactor" )  G_CLO_BOOSTINGFACTOR          = atof( argv[++i] );
        else
        {
            fprintf(stderr, "\n\tUnknown option [%s] on command line.\n\n", argv[i]);
            exit(-1);
        }
    }

    return;
}
