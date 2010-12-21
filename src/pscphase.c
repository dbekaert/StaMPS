// *********************************************************************
// Extract phase for PS candidates from complex interferograms
// ---------------------------------------------------------------------
// AUTHOR    : Andy Hooper
// ---------------------------------------------------------------------
// WRITTEN   : 11.12.2004
//
// Change History
// ==============================================
// 03/2009 MA Fix for gcc 4.3.x
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
     cout << "Usage: pscphase parmfile pscands.1.ij pscands.1.ph" << endl << endl;
     cout << "Input parameters:" << endl;
     cout << "  parmfile   (input)  width of interferogram files (range bins)" << endl;
     cout << "                      list of interferogram files (complex float)" << endl;
     cout << "  pscands.1.ij (input)  location of permanent scatterer candidiates" << endl;
     cout << "  pscands.1.ph (output) phase of permanent scatterer candidates" << endl << endl;
     throw "";
  }   
     
//  char *ijname;
  const char *ijname;      // [MA] deprication fix
  if (argc < 3)
     ijname="pscands.1.ij";
  else ijname = argv[2];

//  char *outfilename;
  const char *outfilename; // [MA] deprication fix
  if (argc < 4)
     outfilename="pscands.1.ph";
  else outfilename = argv[3];

     
  ifstream parmfile (argv[1], ios::in);
  if (! parmfile.is_open()) 
  {	  
      cout << "Error opening file " << argv[1] << endl; 
      throw "";
  }    

  ifstream psfile (ijname, ios::in);
  cout << "opening " << ijname << "...\n";

  if (! psfile.is_open())
  {	    
      cout << "Error opening file " << ijname << endl; 
      throw "";
  }

  ofstream outfile(outfilename,ios::out);
  if (! outfile.is_open()) 
  {	  
      cout << "Error opening file " << outfilename << endl; 
      throw "";
  }    

  char line[256];
  int num_files = 0;
  int width = 0;
  char ifgfilename[256];
  
  parmfile >> width;
  cout << "width = " << width << "\n";	  
  parmfile.getline(ifgfilename,256);
  int savepos=parmfile.tellg();  
  parmfile.getline(ifgfilename,256);
  while (! parmfile.eof())
  {
      parmfile.getline(ifgfilename,256);
      num_files++;
  }    
  cout << "number of interferograms = " << num_files << "\n";
  parmfile.clear();
  parmfile.seekg(savepos);

  ifstream* ifgfile   = new ifstream[num_files];
  float* calib_factor = new float[num_files];
      
  for (int i=0; i<num_files; ++i) 
  {
    parmfile >> ifgfilename;
    ifgfile[i].open (ifgfilename, ios::in|ios::binary);
    cout << "opening " << ifgfilename << "...\n";

    if (! ifgfile[i].is_open())
    {	    
        cout << "Error opening file " << ifgfilename << endl; 
	throw "";
    }

    char header[32];
    long magic=0x59a66a95;
    ifgfile[i].read(header,32);
    if (*reinterpret_cast<long*>(header) == magic)
        cout << "sun raster file - skipping header\n";
    else ifgfile[i].seekg(ios::beg); 

  }
  
  parmfile.close();
  
  char buffer[1000];
  char ifg_pixel[sizeof(float)*2];;
  int pscid=0;
  int x=0;
  int y=0;

  psfile >> pscid >> y >> x;
  psfile.getline(buffer,1000);

  long xyaddr_save = 0;

  for ( int i=0; i<num_files; i++) 
  {
    while (! psfile.eof() ) 
    {
      long xyaddr = (y*width+x)*sizeof(float)*2;
      long xyaddr_diff = xyaddr - xyaddr_save;
      xyaddr_save = xyaddr;

      ifgfile[i].seekg(xyaddr, ios::beg);	    
      ifgfile[i].read (ifg_pixel, sizeof(float)*2);
      outfile.write(ifg_pixel, sizeof(float)*2); 

      //if (pscid/10000.0 == rint(pscid/10000.0))
      //  cout << pscid << " PS candidates processed\n";

      psfile >> pscid >> y >> x;
      psfile.getline(buffer,1000);
    }   
    cout << i+1 << " of " << num_files << " interferograms processed\n";
    psfile.clear();
    psfile.seekg(ios::beg);	    
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

