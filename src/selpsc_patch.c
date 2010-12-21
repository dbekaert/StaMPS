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
// 07/2010 A0 Update Need for Speed on patches
// 08/2010 MA Code optimization
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
  const char *ijname;  // [MA]
  if (argc < 4) 
     ijname="pscands.1.ij";
  else ijname = argv[3];   
     
//  char *ampoutname;
//  if (argc < 4) 
//     ampoutname="pscands.1.amp";
//  else ampoutname = argv[3];   
     
//  char *daoutname;
  const char *daoutname; // [MA]
  if (argc < 5) 
     daoutname="pscands.1.da";
  else daoutname = argv[4];   
     
//  char *meanoutname;
  const char *meanoutname; // [MA]
  if (argc < 6) 
     meanoutname="mean_amp.flt";
  else meanoutname = argv[5];   
  
//  char *maskfilename;
  const char *maskfilename; // [MA]
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
  register float* calib_factor = new float[num_files];
      
  for (register int i=0; i<num_files; ++i) 
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

  // [A0] determine size of a patch
  int patch_lines = az_end-az_start+1;
  int patch_width = rg_end-rg_start+1;
  const int sizeofelement = sizeof(float); // [MA] size of a pixel

  filebuf *pbuf;
  long size;
  long numlines;

  // get pointer to associated buffer object
  pbuf=ampfile[0].rdbuf();

  // get file size using buffer's members
  size=pbuf->pubseekoff (0,ios::end,ios::in);
  pbuf->pubseekpos (0,ios::in);
  numlines=size/width/sizeofelement/2;

  cout << "number of lines per file = " << numlines << "\n";	  
  
  cout << "patch lines = " << patch_lines  << endl;
  cout << "patch width = " << patch_width  <<  endl;

  ifstream maskfile (maskfilename, ios::in);
  char mask_exists = 0;
  if (maskfile.is_open()) 
  {	  
      mask_exists=1;
  }    
  
  ofstream ijfile(ijname,ios::out);
  ofstream daoutfile(daoutname,ios::out);
  ofstream meanoutfile(meanoutname,ios::out);
 
  //complex<float>* buffer = new complex<float>[num_files*width]; // used to store 1 line of all amp files
  register complex<float>* buffer = new complex<float>[num_files*patch_width]; // used to store 1 line of all amp files

  char* maskline = new char[patch_width];
  for (register int x=0; x<patch_width; x++) // for each pixel in range
  {
      maskline[x] = 0;
  }

  //int linebytes = width*8;                      // bytes per line in amplitude files
  const int linebytes = width*sizeofelement*2;  // bytes per line in amplitude files (SLCs)
  register int y=0;                                      // amplitude files line number
  register int pscid=0;                                  // PS candidate ID number

  const int patch_linebytes =  patch_width*sizeofelement*2;
  register long long pos_start;
  pos_start= (long long)(az_start-1+y)*linebytes+(rg_start-1)*sizeofelement*2; // define position of start of 1st line of patch
                                                                               // on SLC file.
  // Turn on following 3 lines for debuing pos_star [MA]
  //cout << "[debug ] " << pos_start << endl;
  //cout << "[debug ] " << (long long )(az_start-1+y)*linebytes << endl;
  //cout << "[debug ] " << (long long )(rg_start-1)*sizeofelement*2 << endl;
  
  for (register int i=0; i<num_files; i++)              // read in first line from each amp file
  {
    //ampfile[i].read (reinterpret_cast<char*>(&buffer[i*width]), linebytes);
    ampfile[i].seekg (pos_start, ios::beg);
    ampfile[i].read (reinterpret_cast<char*>(&buffer[i*patch_width]), patch_linebytes);
  } 
     
  while (! ampfile[1].eof() && y < patch_lines) 
  {
     if (mask_exists==1) 
     {
         //maskfile.read (maskline, width);
         maskfile.seekg (pos_start, ios::beg);      // [AO] set pointer to start of j-th line of mask file
         maskfile.read (maskline, patch_linebytes); // read from pointer nr_pixels*8 bytes
     }    
     if (y >=0) // was (y >= az_start-1)
       {
       for (register int x=0; x<patch_width; x++) // for each pixel in range (width of the patch)
       {
     
        register float sumamp = 0;
        register float sumampsq = 0;

        for (register int i=0; i<num_files; i++)        // for each amp file
	{
           //float amp=abs(buffer[i*width+x])/calib_factor[i]; // get amp value
           register float amp=abs(buffer[i*patch_width+x])/calib_factor[i]; // get amp value
           sumamp+=amp;
           sumampsq+=amp*amp;
        }
        
        //meanoutfile.write(reinterpret_cast<char*>(&sumamp),sizeofelement);	
        meanoutfile.write(reinterpret_cast<char*>(&sumamp),sizeof(float));  // datatypes should be always the same	

        if (maskline[x]==0 && sumamp > 0)
        { 
	    register float D_sq=num_files*sumampsq/(sumamp*sumamp) - 1; // var/mean^2
            if (pick_higher==0 && D_sq<D_thresh_sq ||                 \
                pick_higher==1 && D_sq>=D_thresh_sq) 
	    {
               ++pscid;

               ijfile << pscid << " " << (az_start-1)+y << " " << (rg_start-1)+x << "\n"; 

	       register float D_a = sqrt(D_sq);
               daoutfile << D_a << "\n";

            } // endif
         } // endif
       } // x++           
     } // endif y

     y++;

     for (register int i=0; i<num_files; i++)           // read in next line from each amp file
     {
        //pos_start=(long long)(az_start-1+y)*linebytes+(rg_start-1)*sizeofelement*2; // get pos of the next line of patch
        //ampfile[i].seekg (pos_start, ios::beg);
        ampfile[i].seekg (linebytes-patch_linebytes, ios::cur);  // [MA]
        ampfile[i].read (reinterpret_cast<char*>(&buffer[i*patch_width]), patch_linebytes);
     } 
     
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

