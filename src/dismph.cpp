#include "XInfo.h"
#include "CDisp.h"
#include "CDispComp.h"
#include "CGetData.h"
#include <string.h>
#include <stdlib.h>

using namespace std;

XInfo *xinfo ;
CDispComp *cdi ;
CGetData *cdg ;


int main (int argc, char *argv[]) 
{
	char     infile [240] ;
	int 	 status, ns=5120, nl=0, startl=0, flipflag=0, expflag=1;
	float 	 expval = 0.3, scalefac=1. ;

	if (argc < 3) {
		cout << "Usage : dismph infile nsamps <startl=0> <nlines> <scale=1.> <expval=0.3> <flipflag=0> " << endl ;
		return (1) ;
	}
	argc--;
	strcpy (infile, *++argv) ;
	argc-- ;
	if (argc) 
	{
		ns=atoi (*++argv) ;
		argc-- ;
	}
	if (argc) {
		startl = atoi (*++argv) ;
		argc-- ;
	}
	if (argc) {
		nl = atoi (*++argv) ;
		argc-- ;
	}
	if (argc) {
		scalefac = atof (*++argv) ;
		argc-- ;
	}
	if (argc) {
		expflag = 1 ;
		expval = atof (*++argv) ;
		argc-- ;
	}
	if (argc) {
		flipflag = atoi (*++argv) ;
		argc-- ;
	}

	xinfo = new XInfo (argc, argv) ;
	cdi = new CDispComp (xinfo->toplevel, infile, xinfo) ;
	cdg = new CGetData ;
	cdg->setparams (infile, ns, startl, nl) ;
	cdg->getarrayMPH (scalefac, expflag, expval, flipflag) ;
	cdi->LoadArrays (cdg->mag, cdg->phase, cdg->nsamps, cdg->nlines) ;
	cdg->DeleteMPH () ;
	cdi->StartDisp () ;
	status = xinfo->startloop() ;
	return (0) ;
}


 
void closeup_fb (Widget w, XtPointer client_data, XmAnyCallbackStruct *cbs)
{
	delete cdi ;
	delete cdg ;
   	delete xinfo ;
	exit (0) ;   
     
    
 
}      
