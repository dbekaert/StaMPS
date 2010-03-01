// *********************************************************************
// Calculate amplitude calibration constant for SLC files 
// ---------------------------------------------------------------------
// AUTHOR    : Andy Hooper
// ---------------------------------------------------------------------
// WRITTEN   : 04.08.2003
//
// Change History
// ==============================================
// 01/2009 MA Deprication Fix
// 03/2009 MA Fix for gcc 4.3.x
// ==============================================

#include <iostream>  
using namespace std;
     
#include <fstream>  
using namespace std;

#include <vector>  
using namespace std;
     
#include <cmath>  
using namespace std;
     
#include <cstdio>
using namespace std;     

#include <cstdlib>     
using namespace std;     

#include <complex>     
using namespace std;     

// =======================================================================
// Start of program 
// =======================================================================

//int main(long  argc, char *argv[] ) {    
int main(int  argc, char *argv[] ) {   // [MA]  long --> int for gcc 4.3.x

try {
 
  if (argc < 3)
  {	  
     cout << "Usage: calamp parmfile.in width parmfile.out" << "\n";
     cout << "  parmfile.in(input) SLC file names (complex float)" << endl;
     cout << "  width              width of SLCs" << endl;
     cout << "  parmfile.out(outut) SLC file names and calibration constants" << endl;
     throw "";
  }   
     
  const char *outfilename; // [MA] deprication fix
  if (argc < 4) 
     outfilename="parmfile.out";
  else outfilename = argv[3];   
     
  int width = atoi(argv[2]);

  ifstream ampfiles (argv[1], ios::in);
  if (! ampfiles.is_open()) 
  {	  
      cout << "Error opening file " << argv[1] << "\n"; 
      throw "";
  }   

  ofstream parmfile (outfilename, ios::out);
  if (! parmfile.is_open()) 
  {	  
      cout << "Error opening file " << outfilename << "\n"; 
      throw "";
  }   
      
  char ampfilename[256];
      
  complex<float>* buffer = new complex<float>[width];
  int linebytes = sizeof(complex<float>)*width;
  
  ampfiles >> ampfilename;
  
  while (! ampfiles.eof() ) 
  {
    float calib_factor=0;
    cout << "opening " << ampfilename << "...\n";
    ifstream ampfile (ampfilename, ios::in|ios::binary);

    if (! ampfile.is_open())
    {	    
      cout << "Error opening file " << ampfilename << "\n"; 
      throw "";
    }

    int i=0; 
    double sumamp=0; 
    ampfile.read (reinterpret_cast<char*>(buffer), linebytes);
    while (! ampfile.eof() ) 
    {
      i++;
      for (int j=0; j<width; j++) 
      {
         sumamp+=abs(buffer[j]);
      }
      ampfile.read (reinterpret_cast<char*>(buffer), linebytes);
    }		
    
    calib_factor = sumamp/i/width;

    ampfile.close(); 

    parmfile << ampfilename << " " << calib_factor << "\n";
    cout << "Mean amplitude = " << calib_factor << endl;
  
    ampfiles >> ampfilename;  
  }
  
  ampfiles.close();
  parmfile.close();
  
   
  
  }
  catch( ... ) {
    return(999);
  }

  return(0);
       
};

