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
// 01/2011 MCC neglecting pixel with  zero amplitude
// 12/2012 AH Add byteswap option
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
void cshortswap( complex<short>* f )
{
  char* b = reinterpret_cast<char*>(f);
  complex<short> f2;
  char* b2 = reinterpret_cast<char*>(&f2);
  b2[0] = b[1];
  b2[1] = b[0];
  b2[2] = b[3];
  b2[3] = b[2];
  f[0]=f2;
}

void cfloatswap( complex<float>* f )
{
  char* b = reinterpret_cast<char*>(f);
  complex<float> f2;
  char* b2 = reinterpret_cast<char*>(&f2);
  b2[0] = b[3];
  b2[1] = b[2];
  b2[2] = b[1];
  b2[3] = b[0];
  b2[4] = b[7];
  b2[5] = b[6];
  b2[6] = b[5];
  b2[7] = b[4];
  f[0]=f2;
}

int main(int  argc, char *argv[] ) {   // [MA]  long --> int for gcc 4.3.x

try {
 
  if (argc < 3)
  {	  
     cout << "Usage: calamp parmfile.in width parmfile.out precision byteswap maskfile" << "\n";
     cout << "  parmfile.in(input) SLC file names (complex float)" << endl;
     cout << "  width              width of SLCs" << endl;
     cout << "  parmfile.out(output) SLC file names and calibration constants" << endl;
     cout << "  precision(input) s or f (default)" << endl;
     cout << "  byteswap(input) 1 for to swap bytes, 0 otherwise (default)" << endl;
     cout << "  maskfile   (input)  mask rows and columns (optional)" << endl;
     throw "";
  }   
     
  const char *outfilename; // [MA] deprication fix
  if (argc < 4) 
     outfilename="parmfile.out";
  else outfilename = argv[3];   

  const char *prec;
  if (argc < 5) 
     prec="f";
  else prec = argv[4];   

  int byteswap;
  if (argc < 6) 
     byteswap=0;
  else byteswap = atoi(argv[5]);   

  const char *maskfilename; 
  if (argc < 7)
     maskfilename="";
  else maskfilename = argv[6];

     
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

  ifstream maskfile (maskfilename, ios::in);
  char mask_exists = 0;
  if (maskfile.is_open())
  {
      mask_exists=1;
  }

  char* maskline = new char[width];
  for (register int x=0; x<width; x++) // for each pixel in range
  {
      maskline[x] = 0;
  }

      
  char ampfilename[256];
      
  complex<float>* buffer = new complex<float>[width];
  complex<short>* buffers = reinterpret_cast<complex<short>*>(buffer);
  int linebytes;
  if (prec[0]=='s')
  {
     linebytes = sizeof(complex<short>)*width;
  }else linebytes = sizeof(complex<float>)*width;
  
  ampfiles >> ampfilename;
  
  while (! ampfiles.eof() ) // loop over SLC names
  {
    ifstream maskfile (maskfilename, ios::in);
    char mask_exists = 0;
    if (maskfile.is_open())
    {
        mask_exists=1;
        cout << "opening " << maskfilename << "...\n";
        maskfile.read (maskline, width);
    }


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
    double amp_pixel=0;
    long unsigned int nof_pixels=0;
    long unsigned int nof_zero_pixels=0;
    ampfile.read (reinterpret_cast<char*>(buffer), linebytes);
    while (! ampfile.eof() ) // loop to read all file using buffers
    {
      //i++;
      for (int j=0; j<width; j++) // loop over each read pixel pf the buffer
      { 
         complex<float> camp;
         if (prec[0]=='s')
         {
            if (byteswap == 1)
            {
               cshortswap(&buffers[j]);
            }
            camp=buffers[j];
         }
         else
         {
            camp=buffer[j];

            if (byteswap == 1)
            {
               cfloatswap(&camp);
            }
         }
         amp_pixel=abs(camp);
         if (amp_pixel >0.001 && maskline[j]==0)       //rejects pixels with low amplitude ~0
         {
          sumamp+=abs(camp);
          nof_pixels++;
         }else nof_zero_pixels++;

 
      }
      maskfile.read (maskline, width);
      ampfile.read (reinterpret_cast<char*>(buffer), linebytes);
    }		
    if ( nof_pixels != 0) 
    { 
     //calib_factor = sumamp/i/width;
      calib_factor= sumamp/nof_pixels;
    }
    else
    { 
     cout << "WARNING : SLC " << ampfilename << "has ZERO mean amplitude \n";
     calib_factor =0;
    }

    ampfile.close(); 

    parmfile << ampfilename << " " << calib_factor << "\n";
    cout << "Mean amplitude = " << calib_factor << endl;
    cout << "Number of pixels with zero amplitude = " <<  nof_zero_pixels   << "\n";
    cout << "Number of pixels with amplitude different than zero = " <<  nof_pixels   << "\n";
    
    ampfiles >> ampfilename; 
    if (mask_exists==1)
    {
        maskfile.close();
    }
     
  }
  
  ampfiles.close();
  parmfile.close();
  
   
  
  }
  catch( ... ) {
    return(999);
  }

  return(0);
       
};

