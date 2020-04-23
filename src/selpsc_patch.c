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
// 01/2010 MCC Drop low amplitudes
// 12/2012 AH Add byteswap and short options
// 12/2012 AH Correct mask processing
// 08/2017 AH add "else" in  mask processing
// ==============================================

#include <string.h> 
using namespace std;

#include <iostream>  
using namespace std;
     
#include <fstream>  
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

void shortswap( short* f )
{
  char* b = reinterpret_cast<char*>(f);
  short f2;
  char* b2 = reinterpret_cast<char*>(&f2);
  b2[0] = b[1];
  b2[1] = b[0];
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

void floatswap( float* f )
{
  char* b = reinterpret_cast<char*>(f);
  float f2;
  char* b2 = reinterpret_cast<char*>(&f2);
  b2[0] = b[3];
  b2[1] = b[2];
  b2[2] = b[1];
  b2[3] = b[0];
  f[0]=f2;
}

void longswap( int32_t* f )
{
  char* b = reinterpret_cast<char*>(f);
  int32_t f2;
  char* b2 = reinterpret_cast<char*>(&f2);
  b2[0] = b[3];
  b2[1] = b[2];
  b2[2] = b[1];
  b2[3] = b[0];
  f[0]=f2;
}

//int main(long  argc, char *argv[] ) {    
int main(int  argc, char *argv[] ) {   // [MA]  long --> int for gcc 4.3.x 

try {
 
  if (argc < 3)
  {	  
     cout << "Usage: selpsc parmfile patch.in pscands.1.ij pscands.1.da mean_amp.flt precision byteswap maskfile " << endl << endl;
     cout << "input parameters:" << endl;
     cout << "  parmfile (input) amplitude dispersion threshold" << endl;
     cout << "                   width of amplitude files (range bins)" << endl;
     cout << "                   SLC file names & calibration constants" << endl;
     cout << "  patch.in (input) location of patch in rg and az" << endl;
     cout << "  pscands.1.ij   (output) PS candidate locations" << endl;
     cout << "  pscands.1.da   (output) PS candidate amplitude dispersion" << endl << endl;
     cout << "  mean_amp.flt (output) mean amplitude of image" << endl << endl;
     cout << "  precision(input) s or f (default)" << endl;
     cout << "  byteswap   (input) 1 for to swap bytes, 0 otherwise (default)" << endl;
     cout << "  maskfile   (input)  mask rows and columns (optional)" << endl;
     cout << "  master amplitude (input) in case files in parmfile are ifgs not SLCs (optional)" << endl;
     throw "";
  }   
     
//  char *ijname;
  const char *ijname;  // [MA]
  if (argc < 4) 
     ijname="pscands.1.ij";
  else ijname = argv[3];   

  char jiname[256]; // float format big endian for gamma
  strcpy (jiname,ijname);
  strcat (jiname,".int");
  //MCC
  char ijname0[256]; //used to store PS with at least one amplitude =0. This can happen if frames do not fully overlap
                  
  strcpy (ijname0,ijname);
  strcat (ijname0,"0");
  cout << "file name for zero amplitude PS: " << ijname0 << "\n";
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
  
  const char *prec;
  if (argc < 7)
     prec="f";
  else prec = argv[6];

  int byteswap;
  if (argc < 8)
     byteswap=0;
  else byteswap = atoi(argv[7]);

//  char *maskfilename;
  const char *maskfilename; // [MA]
  if (argc < 9) 
     maskfilename="";
  else maskfilename = argv[8];   


  char masterampfilename[256]="0000";
  char masteramp_exists = 0;
  if (argc < 10) 
  {
  }
  else
  { 
      ifstream masterparmfile (argv[9], ios::in);
      if (! masterparmfile.is_open()) 
      {	  
          cout << "Error opening file " << argv[9] << "\n"; 
          throw "";
      }    
      masterparmfile >> masterampfilename;
  }

  ifstream masterampfile (masterampfilename, ios::in);
  if (masterampfile.is_open()) 
  {	  
      masteramp_exists=1;
      cout << "opening " << masterampfilename << "...\n";
  }    
     
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

  const int sizeoffloat=4; // [MA] size of a pixel
  int sizeofelement; // [MA] size of a pixel
  if (prec[0]=='s')
  {
      sizeofelement = sizeof(short);
  }else sizeofelement = sizeof(float);

  const int linebytes = width*sizeofelement*2;  // bytes per line in amplitude files (SLCs)
  const int patch_linebytes =  patch_width*sizeofelement*2;
  const int patch_amp_linebytes =  patch_width*sizeofelement;

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
      cout << "opening " << maskfilename << "...\n";
  }    

  
  
  ofstream ijfile(ijname,ios::out);
  ofstream jifile(jiname,ios::out);
  ofstream ijfile0(ijname0,ios::out);
  ofstream daoutfile(daoutname,ios::out);
  ofstream meanoutfile(meanoutname,ios::out);
 
  //complex<float>* buffer = new complex<float>[num_files*width]; // used to store 1 line of all amp files
  char* buffer = new char[num_files*patch_linebytes]; // used to store 1 line of all amp files
  complex<float>* bufferf = reinterpret_cast<complex<float>*>(buffer); 
  complex<short>* buffers = reinterpret_cast<complex<short>*>(buffer);


  char* maskline = new char[patch_width];
  for (register int x=0; x<patch_width; x++) // for each pixel in range
  {
      maskline[x] = 0;
  }

  char* masterampline = new char[patch_linebytes]; // used to store 1 line of all amp files
  complex<float>* masterlinef = reinterpret_cast<complex<float>*>(masterampline); 
  complex<short>* masterlines = reinterpret_cast<complex<short>*>(masterampline);
  for (register int x=0; x<patch_width; x++) // for each pixel in range
  {
      if (prec[0]=='s')
      {
          masterlines[x] = 1;
          if (byteswap == 1)
          {
             cshortswap(&masterlines[x]);
          }
      }
      else
      {
          masterlinef[x] = 1;
          if (byteswap == 1)
          {
             cfloatswap(&masterlinef[x]);
          }
      }
  }

  //int linebytes = width*8;                      // bytes per line in amplitude files
  register int y=0;                                      // amplitude files line number
  register int pscid=0;                                  // PS candidate ID number

  register long long pix_start;
  register long long pos_start;
  pix_start= (long long)(az_start-1)*width+(rg_start-1); // define pixel number of start of 1st line of patch
  pos_start= pix_start*sizeofelement*2; // define position of start of 1st line of patch
                                                                               // on SLC file.
  // Turn on following 3 lines for debuing pos_star [MA]
  //cout << "[debug ] " << pos_start << endl;
  //cout << "[debug ] " << (long long )(az_start-1+y)*linebytes << endl;
  //cout << "[debug ] " << (long long )(rg_start-1)*sizeofelement*2 << endl;
  
  for (register int i=0; i<num_files; i++)              // read in first line from each amp file
  {
    //ampfile[i].read (reinterpret_cast<char*>(&buffer[i*width]), linebytes);
    ampfile[i].seekg (pos_start, ios::beg);
    ampfile[i].read (&buffer[i*patch_linebytes], patch_linebytes);
  } 
  if (mask_exists==1) 
  {
      //maskfile.read (maskline, width);
      maskfile.seekg (pix_start, ios::beg);      // set pointer to start of patch in mask file
      maskfile.read (maskline, patch_width); // read from pointer nr_pixels
  }    
  if (masteramp_exists==1) 
  {
      masterampfile.seekg (pos_start, ios::beg);      // set pointer to start of patch in mask file
      masterampfile.read (&masterampline[0], patch_linebytes); // read from pointer nr_pixels
  }    
     
     
  while (! ampfile[1].eof() && y < patch_lines) 
  {
     if (y >=0) // was (y >= az_start-1)
       {
       for (register int x=0; x<patch_width; x++) // for each pixel in range (width of the patch)
       {
     
        register float sumamp = 0;
        register float sumampsq = 0;
        int amp_0 =0;

        complex<float> master_amp; //  master amp value
        if (prec[0]=='s')
        {
            if (byteswap == 1)
            {
               cshortswap(&masterlines[x]);
            }
            master_amp=masterlines[x];
        }
        else
        {
            master_amp=masterlinef[x]; // get amp value
            if (byteswap == 1)
            {
               cfloatswap(&master_amp);
            }
        }
        //cout << "master_amp: " << abs(master_amp) << endl;
        if (abs(master_amp)==0)
        {
            master_amp=1;
        }

        for (register int i=0; i<num_files; i++)        // for each amp file
	   {
           complex<float> camp; //  amp value
           //float amp=abs(buffer[i*width+x])/calib_factor[i]; // get amp value
           if (prec[0]=='s')
           {
               if (byteswap == 1)
               {
                  cshortswap(&buffers[i*patch_width+x]);
               }
               camp=buffers[i*patch_width+x];
           }
           else
           {
               camp=bufferf[i*patch_width+x]; // get amp value
               if (byteswap == 1)
               {
                  cfloatswap(&camp);
               }
           }

           //cout << "camp: " << abs(camp) << " calib " << calib_factor[i] << " master " << abs(master_amp) << endl ;
         
           register float amp=abs(camp)/calib_factor[i]/abs(master_amp); // get amp value
           //cout << "amp: " << amp << endl ;
           if (amp <=0.00005) // do not use amp = 0 values for calculating the AD and set flag to 1
           {
            amp_0=1;
            sumamp=0;
            //cout << "coord amp zero : az  " << (az_start-1)+y << ", rng  " << (rg_start-1)+x << "\n";  
            //cout << "  amp : " << amp << " \n" ;
            continue  ; 
           }else
           {
           sumamp+=amp;
           sumampsq+=amp*amp;
           }
         }
	
        meanoutfile.write(reinterpret_cast<char*>(&sumamp),sizeoffloat);	
 //       cout << "amp0: " << amp_0  << "sumamp " << sumamp <<" \n" ;
        if (maskline[x]==0 && sumamp > 0)
        {
      //Amplitude disperion^2 
	    register float D_sq=num_files*sumampsq/(sumamp*sumamp) - 1; // var/mean^2
        if (pick_higher==0 && D_sq<D_thresh_sq  ||                 \
            pick_higher==1 && D_sq>=D_thresh_sq) 
	      {
            if (amp_0 != 1)
             {
             // cout << "selected! \n";
             // cout << "amp0: " << amp_0  << ", sumamp: " << sumamp <<" \n" ;

               ++pscid;

               ijfile << pscid << " " << (az_start-1)+y << " " << (rg_start-1)+x << "\n";      
               int32_t J=(rg_start-1)+x;
               int32_t I=(az_start-1)+y;
               longswap(&J);
               longswap(&I);
               jifile.write(reinterpret_cast<char*>(&J), sizeof(int32_t));
               jifile.write(reinterpret_cast<char*>(&I), sizeof(int32_t));

	       register float D_a = sqrt(D_sq);
               daoutfile << D_a << "\n";
             }
            else //Keeping track of PS with zero amplitude values
             {
              ijfile0 << pscid << " " << (az_start-1)+y << " " << (rg_start-1)+x << "\n";
             } // end if amp_0
             
            
         } // endif pick_highe
         } // endif maskline
       } //for loop x++           
     } // endif y >=0

     y++;

     for (register int i=0; i<num_files; i++)           // read in next line from each amp file
     {
        //pos_start=(long long)(az_start-1+y)*linebytes+(rg_start-1)*sizeofelement*2; // get pos of the next line of patch
        //ampfile[i].seekg (pos_start, ios::beg);
        ampfile[i].seekg (linebytes-patch_linebytes, ios::cur);  // [MA]
        ampfile[i].read (&buffer[i*patch_linebytes], patch_linebytes);
     } 
     if (mask_exists==1) 
     {
     maskfile.seekg (width-patch_width, ios::cur);  // [MA]
     maskfile.read (maskline, patch_width);
     }
     if (masteramp_exists==1) 
     {
         masterampfile.seekg (linebytes-patch_linebytes, ios::cur);  // 
         masterampfile.read (&masterampline[0], patch_linebytes); // read from pointer nr_pixels
     }    
     
     
     if (y/100.0 == rint(y/100.0))
        cout << y << " lines processed\n";
        //cout << pscid  << " selected PS \n";
        //cout << D_thresh_sq << "D_thresh_sq \n";
  }  
  ijfile.close();
  jifile.close();
  ijfile0.close();
  //ampoutfile.close();
  daoutfile.close();
  meanoutfile.close();
  if (mask_exists==1) 
  {	  
      maskfile.close();
  }    
  if (masteramp_exists==1) 
  {	  
      masterampfile.close();
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

