// *********************************************************************
// Select PS Candidates
// Input are SLC's
// ---------------------------------------------------------------------
// AUTHOR    : Andy Hooper
// ---------------------------------------------------------------------
// WRITTEN   : 04.08.2003
//  
// Change History
// ==============================================
// 03/2009 MA Fix for gcc 4.3.x
// ==============================================

#include <iostream>  
using namespace std;
     
#include <fstream>  
using namespace std;

#include <string>  
using namespace std;
     
#include <complex>  
using namespace std;
     
#include <vector>  
using namespace std;
     
#include <cmath>  
using namespace std;
     
#include <cstdio>
using namespace std;     

#include <cstdlib>     
using namespace std;     

#include <climits>     
using namespace std;     

// =======================================================================
// Start of program 
// =======================================================================

//int main(long  argc, char *argv[] ) {    
int main(int  argc, char *argv[] ) {   // [MA]  long --> int for gcc 4.3.x 

try {
 
  if (argc < 3)
  {	  
     cout << "Usage: selpsc parmfile patch.in pscands.1.ij pscands.1.da mean_amp.flt maskfile " << endl << endl;
     cout << "input parameters:" << endl;
     cout << "  parmfile (input) amplitude dispersion threshold" << endl;
     cout << "                   width of amplitude files (range bins)" << endl;
     cout << "                   SLC file names & calibration constants" << endl;
     cout << "  patch.in (input) location of patch in rg and az" << endl;
     cout << "  pscands.1.ij   (output) PS candidate locations" << endl;
     cout << "  pscands.1.da   (output) PS candidate amplitude dispersion" << endl << endl;
     cout << "  mean_amp.flt (output) mean amplitude of image" << endl << endl;
     cout << "  maskfile   (input)  mask rows and columns (optional)" << endl;
     throw "";
  }   
     
//  char *ijname;
  const char *ijname;          // [MA] deprication fix
  if (argc < 4) 
     ijname="pscands.1.ij";
  else ijname = argv[3];   
     
//  char *ampoutname;
//  if (argc < 4) 
//     ampoutname="pscands.1.amp";
//  else ampoutname = argv[3];   
     
//  char *daoutname;
  const char *daoutname;       // [MA] deprication fix
  if (argc < 5) 
     daoutname="pscands.1.da";
  else daoutname = argv[4];   
     
//  char *meanoutname;
  const char *meanoutname;    // [MA] deprication fix
  if (argc < 6) 
     meanoutname="mean_amp.flt";
  else meanoutname = argv[5];   
  
//  char *maskfilename;
  const char *maskfilename;   // [MA] deprication fix
  if (argc < 7) 
     maskfilename="";
  else maskfilename = argv[6];   
     
     
  ifstream parmfile (argv[1], ios::in);
  if (! parmfile.is_open()) 
  {	  
      cout << "Error opening file " << argv[1] << "\n"; 
      throw "";
  }    
  
      
  char line[256];
  int num_files = 0;
  
  int width = 0;
  float D_thresh = 0;
  int pick_higher = 0;

  parmfile >> D_thresh;
  cout << "dispersion threshold = " << D_thresh << "\n";
  float D_thresh_sq = D_thresh*D_thresh;
  if (D_thresh<0) { 
      pick_higher=1;
  }

  parmfile >> width;
  cout << "width = " << width << "\n";	  
  parmfile.getline(line,256);
  int savepos=parmfile.tellg();
  parmfile.getline(line,256);
  while (! parmfile.eof())
  {
      parmfile.getline(line,256);
      num_files++;
  }    
  //parmfile >> num_files;
  parmfile.clear();
  parmfile.seekg(savepos);
  char ampfilename[256];
  ifstream* ampfile   = new ifstream[num_files];
  float* calib_factor = new float[num_files];
      
  for (int i=0; i<num_files; ++i) 
  {
    parmfile >> ampfilename >> calib_factor[i];
    ampfile[i].open (ampfilename, ios::in|ios::binary);
    cout << "opening " << ampfilename << "...\n";

    if (! ampfile[i].is_open())
    {	    
        cout << "Error opening file " << ampfilename << "\n"; 
	throw "";
    }

    char header[32];
    long magic=0x59a66a95;
    ampfile[i].read(header,32);
    if (*reinterpret_cast<long*>(header) == magic)
        cout << "sun raster file - skipping header\n";
    else ampfile[i].seekg(ios::beg); 
  }
  
  parmfile.close();
  cout << "number of amplitude files = " << num_files << "\n";

  ifstream patchfile (argv[2], ios::in);
  if (! patchfile.is_open()) 
  {	  
      cout << "Error opening file " << argv[2] << "\n"; 
      throw "";
  }    

  int rg_start=0;
  int rg_end=INT_MAX;
  int az_start=0;
  int az_end=INT_MAX;
  patchfile >> rg_start;
  patchfile >> rg_end;
  patchfile >> az_start;
  patchfile >> az_end;
  patchfile.close();

  filebuf *pbuf;
  long size;
  long numlines;

  // get pointer to associated buffer object
  pbuf=ampfile[0].rdbuf();

  // get file size using buffer's members
  size=pbuf->pubseekoff (0,ios::end,ios::in);
  pbuf->pubseekpos (0,ios::in);
  numlines=size/width/sizeof(float)/2;

  cout << "nunber of lines per file = " << numlines << "\n";	  
  
  ifstream maskfile (maskfilename, ios::in);
  char mask_exists = 0;
  if (maskfile.is_open()) 
  {	  
      mask_exists=1;
  }    
  
  ofstream ijfile(ijname,ios::out);
  ofstream daoutfile(daoutname,ios::out);
  ofstream meanoutfile(meanoutname,ios::out);
 
  complex<float>* buffer = new complex<float>[num_files*width]; // used to store 1 line of all amp files

  char* maskline = new char[width];
  for (int x=0; x<width; x++) // for each pixel in range
  {
      maskline[x] = 0;
  }

  int linebytes = width*8;                      // bytes per line in amplitude files`
  int y=0;                                      // amplitude files line number
  int pscid=0;                                  // PS candidate ID number
  
  for ( int i=0; i<num_files; i++)              // read in first line from each amp file
  {
    ampfile[i].read (reinterpret_cast<char*>(&buffer[i*width]), linebytes);
  } 
     
  while (! ampfile[1].eof() && y < az_end) 
  {
     if (mask_exists==1) 
     {
         maskfile.read (maskline, width);
     }    
     if (y >= az_start-1)
       {
       for (int x=rg_start-1; x<rg_end; x++) // for each pixel in range
       {
     
        float sumamp = 0;
        float sumampdiffsq = 0;
        int i;

        for (i=0; i<num_files/2; i++)        // for each amp file
	{
           float amp1=abs(buffer[(i*2)*width+x])/calib_factor[i*2]; // get amp value
           float amp2=abs(buffer[(i*2+1)*width+x])/calib_factor[i*2+1]; // get amp value
           //cout << "amp1=" << amp1 << " amp2=" << amp2 << "\n";
           sumamp+=amp1;
           sumamp+=amp2;
           //cout << "sumamp=" << sumamp << "\n";
           sumampdiffsq+=(amp1-amp2)*(amp1-amp2);
           //cout << "sumsq=" << sumampdiffsq << "\n";
        }
        meanoutfile.write(reinterpret_cast<char*>(&sumamp),sizeof(float));	

        if (maskline[x]==0 && sumamp > 0)
        { 
	    float D_a=sqrt(sumampdiffsq/(num_files/2))/(sumamp/num_files); 
            if (pick_higher==0 && D_a<D_thresh ||                 \
                pick_higher==1 && D_a>=D_thresh) 
	    {
               ++pscid;

               ijfile << pscid << " " << y << " " << x << "\n"; 

               daoutfile << D_a << "\n";
               //cout << "Da=" << D_a << "\n";

            } // endif
         } // endif
       } // x++           
     } // endif

     for ( int i=0; i<num_files; i++)           // read in next line from each amp file
     {
        ampfile[i].read (reinterpret_cast<char*>(&buffer[i*width]), linebytes);
     } 
     
     y++;
     
     if (y/100.0 == rint(y/100.0))
        cout << y << " lines processed\n";
  }  
  ijfile.close();
  //ampoutfile.close();
  daoutfile.close();
  meanoutfile.close();
  if (mask_exists==1) 
  {	  
      maskfile.close();
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

