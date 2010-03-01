#include "CDisp.h"

#ifndef CDM
#define CDM

class CDispMag : public CDisp {
	public : 
		
		CDispMag (Widget w, char *infile, XInfo *xi) ;
		~CDispMag () ;
		int LoadArray (unsigned char *, int, int) ;
		int StartDisp () ;
		XImage *ximage_mag ;
		unsigned char *mag ;
	
		// overidden CDisp base functions
                void expose_zm (Widget w, XtPointer xinfoeq, void *cbs) ;
                static void expose_dazm (Widget w, XtPointer xinfoeq, void *cbs) ;  

} ;

#endif

