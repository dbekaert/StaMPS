// *********************************************************************
// cpxsum - program to sum (or subtract) 2 complex phase files 
// ---------------------------------------------------------------------
// AUTHOR    : Andy Hooper
// ---------------------------------------------------------------------
// WRITTEN   : 04.07.2006
//
// Change History
// ==================================================
// 03/2009 MA Fix for gcc 4.3.x
// 03/2014 AH add int before main for latest compiler
// ==================================================

#include <iostream>
using namespace std;

#include <fstream>
using namespace std;

#include <complex> 
using namespace std;

#include <cstring> 
using namespace std;

#include <cstdlib>  // [MA] gcc 4.3.x for exit and atoi

int main(int  argc, char *argv[] )    
{
  cout << "\ncpxsum, Andy Hooper, July 2006\n\n";

  int i,j,nlf,width,sumflag,normflag;

  if (argc < 5 ) {
    printf("cpxsum multiplies a complex file by a 2nd file (or its conjugate)\n\n");
    printf("usage: cpxsum INFILE INFILE2 OUTFILE WIDTH FORMAT SUMFLAG NORMFLAG\n");
    printf("       FORMAT refers to INFILE2 and can be cr4 (default) or r4\n");
    printf("       SUMFLAG is 1 for sum (default) or -1 for conjugate\n");
    printf("       NORMFLAG is 1 to normalize the amplitude of 2nd file\n\n");
    exit(0);
  }
  
//  char* infile1;
  const char* infile1; // [MA] gcc 4.3.x 
  infile1=argv[1];
//  char* infile2;
  const char* infile2; // [MA] gcc 4.3.x 
  infile2=argv[2];
//  char* outfile;
  const char* outfile; // [MA] gcc 4.3.x 
  outfile=argv[3];
//  char* formind;
  const char* formind;  // [MA] deprication fix
  width=atoi(argv[4]);
  if (argc > 5 ) {
      formind=argv[5];
  }
  else {
      formind="cr4";
  } 
  if (argc > 6 ) {
      sumflag=atoi(argv[6]);
  }
  else {
      sumflag=1;
  } 
  if (argc > 7 ) {
      normflag=atoi(argv[7]);
  }
  else {
      normflag=0;
  } 

   /* open files */
    ifstream fpin1;
    //cout << "opening " << infile1 << "...\n";
    fpin1.open (infile1, ios::in|ios::binary);
    if (! fpin1.is_open())
    {
        cout << "Error opening file " << infile1 << "\n";
        exit(4);
    }


    ifstream fpin2;
    //cout << "opening " << infile2 << "...\n";
    fpin2.open (infile2, ios::in|ios::binary);
    if (! fpin2.is_open())
    {
        cout << "Error opening file " << infile2 << "\n";
        exit(4);
    }

    ofstream fpout;
    //cout << "opening " << outfile << "...\n";
    fpout.open (outfile, ios::out);
    if (! fpout.is_open())
    {
        cout << "Error opening file " << outfile << "\n";
        exit(4);
    }


 // get pointer to associated buffer object
 filebuf* pbuf=fpin1.rdbuf();

 // get file size using buffer's members
 long size=pbuf->pubseekoff (0,ios::end,ios::in);
 pbuf->pubseekpos (0,ios::in);
 nlf=size/width/sizeof(complex<float>);

 cout << "Number of lines in " << infile1 << " = " << nlf << "\n";


 /* allocate memory */
 complex<float>* inrec1 = new complex<float>[width];
 int file2size; 
 if (strncmp(formind,"cr4",3)==0) {
     file2size=sizeof(complex<float>);
     cout << "format of " << infile2 << " assumed complex single\n";
 }
 else {
     file2size=sizeof(float);
     cout << "format of " << infile2 << " assumed real single\n";
 }
 char* inrec2 = new char[width*file2size]; 
 fpin1.read (reinterpret_cast<char*>(inrec1), width*sizeof(complex<float>));
 fpin2.read (reinterpret_cast<char*>(inrec2), width*file2size);
 for (i=0;i<nlf;i++) {
   for (j=0;j<width;j++) {
      complex<float> cpxphase2;
      if (strncmp(formind,"cr4",3)==0) {
          cpxphase2 = *reinterpret_cast<complex<float>*>(&inrec2[j*file2size]); 
      }
      else {
          float phase2 = *reinterpret_cast<float*>(&inrec2[j*file2size]); 
          complex<float> J(0,1);
          cpxphase2 = exp(J*phase2); 
      }

      if (sumflag < 0)
          {cpxphase2 = conj(cpxphase2);
      }
      if ((normflag > 0) && (cpxphase2 != complex<float>(0,0)) )
          {cpxphase2 = cpxphase2/abs(cpxphase2);
      }
      inrec1[j] *= cpxphase2;
   }
   fpout.write(reinterpret_cast<char*>(inrec1),width*sizeof(complex<float>));
   fpin1.read (reinterpret_cast<char*>(inrec1),width*sizeof(complex<float>));
   fpin2.read (reinterpret_cast<char*>(inrec2), width*file2size);
 }
     
 /* clean up */
cout << "Output written to " << outfile << "\n";
fpin1.close();
fpin2.close();
fpout.close();
return(0);
}
