#include "CDisp.h"
#ifndef CDComp 
#define CDComp

class CDispComp : public CDisp {
	public :
		XImage *ximage_mag ;
		XImage *ximage_phase ;
		XImage *ximage_both ;
		int    DispArr ;
		unsigned char *phase, *mag, *both ;
		~CDispComp () ;
		CDispComp (Widget w, char *infile, XInfo *xi) ;
		int LoadArrays (unsigned char *r, unsigned char *g, int ns, int nl) ;
		int WriteArrays (unsigned char *r, unsigned char *g, int ns, int nl) ;
		int StartDisp () ;

		// overidden CDisp base functions
		void expose_zm (Widget w, XtPointer xinfoeq, void *cbs) ;  
		static void expose_dazm (Widget w, XtPointer xinfoeq, void *cbs) ;

		static void click_da (Widget w, XtPointer xinfoeq, void *cbs) ; 
		void click (Widget w, XtPointer xinfoeq, void *cbs) ;
} ;

#endif
