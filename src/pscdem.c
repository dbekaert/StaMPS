// *********************************************************************
// pscdem: Extract height for PS Candidates
// ---------------------------------------------------------------------
// AUTHOR    : Andy Hooper
// ---------------------------------------------------------------------
// WRITTEN   : 10.03.2005
//
// Change History
// ==============================================
// 03/2009 MA Fix for gcc 4.3.x
// 02/2010 AH Allow processing with no DEM
// 09/2015 AH Allow double precision
// ==============================================

#include <iostream>  
using namespace std;
     
#include <fstream>  
using namespace std;

#include <complex>  
using namespace std;
     
#include <string>  
using namespace std;
     
#include <cmath>  
using namespace std;
     
#include <cstdio>
using namespace std;     

#include <cstdlib>     
using namespace std;     

// =======================================================================
// Start of program 
// =======================================================================

//int main(long  argc, char *argv[] ) {    
int main(int  argc, char *argv[] ) {    // [MA]  long --> int for gcc 4.3.x


try {
 
  if (argc < 2)
  {	  
     cout << "Usage: pscdem pscands.1.ij pscands.1.hgt precision" << endl << endl;
     cout << "Input parameters:" << endl;
     cout << "  parmfile   (input)  width of dem files (range bins)" << endl;
     cout << "                      name of dem file (radar coords, float)" << endl;
     cout << "  pscands.1.ij (input)  location of PS candidiates" << endl;
     cout << "  pscands.1.hgt (output) height of PS candidiates" << endl << endl;
     cout << "  precision(input) d or f (default)" << endl;
     throw "";
  }   
     
//  char *ijname;
  const char *ijname;        // [MA] deprication fix
  if (argc < 3)
     ijname="pscands.1.ij";
  else ijname = argv[2];

//  char *outfilename;
  const char *outfilename;   // [MA] deprication fix
  if (argc < 4)
     outfilename="pscands.1.hgt";
  else outfilename = argv[3];
// precision
  const char *prec;
  if (argc < 5)
     prec="f";
  else prec = argv[4];

  ifstream parmfile (argv[1], ios::in);
  if (! parmfile.is_open()) 
  {	  
      cout << "pscdem: Error opening file " << argv[1] << endl; 
      throw "";
  }    

  ifstream psfile (ijname, ios::in|ios::binary);
  cout << "opening " << ijname << "...\n";

  if (! psfile.is_open())
  {	    
      cout << "pscdem: Error opening file " << ijname << endl; 
      throw "";
  }

  ofstream outfile(outfilename,ios::out);
  if (! outfile.is_open()) 
  {	  
      cout << "pscdem: Error opening file " << outfilename << endl; 
      throw "";
  }    

  char line[256];
  int num_files = 1;
  int width = 0;
  char ifgfilename[256];
  
  parmfile >> width;
  cout << "pscdem: width = " << width << "\n";	  
  parmfile.getline(ifgfilename,256);

  ifstream* ifgfile   = new ifstream[num_files];
  float* calib_factor = new float[num_files];
      
  parmfile >> ifgfilename;
  ifgfile[0].open (ifgfilename, ios::in|ios::binary);
  cout << "pscdem: opening " << ifgfilename << "...\n";

  int nodem_sw = 0;   
  if (! ifgfile[0].is_open())
  {	    
      cout << "pscdem: Warning: cannot open " << ifgfilename << " - assumed all zeros" << endl; 
      nodem_sw = 1;
  }

  char header[32];
  long magic=0x59a66a95;
  ifgfile[0].read(header,32);
  if (*reinterpret_cast<long*>(header) == magic)
      cout << "pscdem: sun raster file - skipping header\n";
  else ifgfile[0].seekg(ios::beg); 
  
  parmfile.close();
  
  char buffer[1000];

  int sizeofpixel; // 
  if (prec[0]=='d')
  {
      sizeofpixel=sizeof(double);
      cout << "pscdem: input file specified as double precision\n";
  }else 
  {
      sizeofpixel=sizeof(float);
      cout << "pscdem: input file specified as single precision\n";
  }

  char ifg_pixel[sizeof(double)]={0};
  int pscid=0;
  int x=0;
  int y=0;

  psfile >> pscid >> y >> x;
  psfile.getline(buffer,1000);

  long xyaddr_save = 0;

  while (! psfile.eof() ) 
  //for (int dummyi=0; dummyi<100; dummyi++) 
  {
    long xyaddr = (y*width+x)*sizeofpixel;

    if (nodem_sw == 0)
    {
      ifgfile[0].seekg(xyaddr, ios::beg);	    
      ifgfile[0].read (ifg_pixel, sizeofpixel);
    } 
    outfile.write(ifg_pixel, sizeofpixel);

    if (pscid/100000.0 == rint(pscid/100000.0))
      cout << "pscdem: " << pscid << " PS candidates processed\n";

    psfile >> pscid >> y >> x;
    psfile.getline(buffer,1000);
  }   
   
  
  }
  catch( char * str ) {
     cout << str << "\n";
     return(999);
  }   
  catch( ... ) {
    return(999);
  }

  return(0);
       
};

