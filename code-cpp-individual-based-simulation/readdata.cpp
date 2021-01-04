#include "readdata.h"

// struct tm string_to_tmstruct( string str, int format_type );
// int daysdiff( struct tm a, struct tm b );


bool isFloat( string myString ) 
{
    istringstream iss(myString);
    float f;
    iss >> noskipws >> f; // noskipws considers leading whitespace invalid
    // Check the entire string was consumed and if either failbit or badbit is set
    return iss.eof() && !iss.fail(); 
}

vector< vector<double> > column_filter( vector< vector<double> > &DATA, int col_index, double value )
{
    vector< vector<double> > vv;
    
    for( int r=0; r<DATA.size(); r++ )
    {
        if( DATA[r][col_index] == value )
        {
            vv.push_back( DATA[r] );   
        }
    }
    return vv;
}



void readdata( string filename, vector< vector<double> > &DATA, int numcol, bool header  )
{
    DATA.clear();
    assert( DATA.size() == 0 );

    ifstream infile( filename.c_str(), ifstream::in );
    
    //TODO need some error handling here if the file above is not found
	
    string strHeader;
    if( header ) getline(infile, strHeader);

    int col=0;
    int row=0;
    
    int counter=0;
    
    while( true )
    {
        counter++;
	  
        string str("");
        infile >> str;
	
        //fprintf( stderr, "\n\t%d \t%d \t%d \t %1.5f", (int)str.length(), row, col, atof( str.c_str() ) ); fflush(stderr);
	
        if( str.length()==0 && infile.eof() ) break;
	
        // quick & dirty check to make sure the string is a float
        //assert( isdigit( str[0] ) || str[0]=='-' || str[0]=='.' );
	
        // if this is the first column, then it's a new row and we need
        // to allocate a new vector
        if( col==0 )
        {
            vector<double> vd(numcol);
            DATA.push_back( vd );
        }

	if(isFloat(str))
        {
	    DATA[row][col] = atof( str.c_str() );
        }
        else
        {
            DATA[row][col] = -101.0;
        }
	
	
	if( col==numcol-1 ) // if the previous statement just read in the last column
	{
	    col=0;
	    row++;
	}
	else
	{
	    col++;
	}
	
    }
    
    infile.close();
    return;
}




void readdata_02FL( string filename, vector< vector<double> > &DATA, int numcol, bool header  )
{
    DATA.clear();
    assert( DATA.size() == 0 );

    int extracol = 5;
    
    ifstream infile( filename.c_str(), ifstream::in );
	
    string strHeader;
    if( header ) getline(infile, strHeader);

    int col=0;
    int row=0;
    
    while( true )
    {
	string str("");
        infile >> str;
        if( str.length()==0 && infile.eof() ) break;
	
	// if the string looks like an NA or NaN, change it to a double value of -99.0
	if( str.length()==2 || str.length()==3 )
	{
            if( (str[0]=='N' || str[0]=='n') && (str[1]=='A' || str[1]=='a') )
	    {
                str = "-99.0"; 
	    }
	}
	// quick-dirty sanity check here to make sure the string is a number
	assert( isdigit( str[0] ) || str[0]=='-' );
	
        // if this is the first column, then it's a new row and we need
        // to allocate a new vector
	if( col==0 )
	{
	    vector<double> vd(numcol+extracol);
	    DATA.push_back( vd );
	}

	DATA[row][col] = atof( str.c_str() );
	
	// double check here if you just read in col 17, it should be the same as column 1 ( redundant site codes in the data )
	if( col==17 )
	{
	    int site1 = (int)(DATA[row][1]+0.0001); 
	    int site2 = (int)(DATA[row][17]+0.0001); 
	    assert( site1==site2 );
	}
	
	
	if( col==numcol-1 ) // if the previous statement just read in the last column
	{
	    // precompute 4 special columns, reset to col to zero, and increment row
	  
	    // this is the integer day with Jan 1 2009 being day 1
	    int dn = get_day_number( ( (int) (DATA[row][7]+0.001) ), ( (int) (DATA[row][8]+0.001) ), ( (int) (DATA[row][9]+0.001) ) );
	    DATA[row][numcol] = (double) dn;
	    
	    // this is the floating point age
	    DATA[row][numcol+1] = 0.5*(DATA[row][4]+DATA[row][5]);   
	    
	    // this is the floating point date
            int basedays=0;
	    int y =  (int) (DATA[row][9]+0.001);
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

	    int dn_inyear = dn-basedays;
	    assert( dn_inyear > 0 && dn_inyear <= 366 );
	    double dDenom = (y==2012||y==2016||y==2020) ? 366.0 : 365.0; 
	    double dFloatingPointDate = DATA[row][9] + (  (((double)dn_inyear)-0.5) / dDenom  );
	    DATA[row][numcol+2] = dFloatingPointDate;   
	    
	    
	    // columns 46-49 are the H1-09-Titer
	    // columns 66-69 are the H3-09-Titer
	    // columns 70-73 are the H3-11-Titer
	    double h109 = besttiter( DATA[row][46], DATA[row][47], DATA[row][48], DATA[row][49] );
	    double h309 = besttiter( DATA[row][66], DATA[row][67], DATA[row][68], DATA[row][69] );
	    double h311 = besttiter( DATA[row][70], DATA[row][71], DATA[row][72], DATA[row][73] );
	    
	    // this is the best H1 titer
	    DATA[row][numcol+3] = h109;   
	    
	    // this is the best H3 titer
	    double gmt_h3;
	    if( h309 <= 0.0 )
	    {
	        gmt_h3 = h311; 
	    }
	    else
	    {
	        if( h311 <= 0.0 )
		{
		    gmt_h3 = h309; 
		}
		else
		{
		    gmt_h3 = sqrt( h309 * h311 );  
		}
	    }
	    DATA[row][numcol+4] = gmt_h3;   
	    
	    col=0;
	    row++;
	}
	else
	{
	    col++;
	}
	
    }
    
    infile.close();
    return;
}

