#include "XInfo.h"
#include <Xm/DrawingA.h>
#include <Xm/ScrolledW.h>
#include <Xm/RowColumn.h>
#include <Xm/Xm.h>
#include <X11/Core.h>
#ifndef CD
#define CD

class CDisp {
	public : 
		CDisp () ;
		CDisp (Widget toplevel, char *dispfile, XInfo *xinfo_ptr) ;
		virtual ~CDisp () ;

	Widget  sw ;
	Widget  zoomshell, dispshell ;
	Widget  da ;
	Widget  da_zm ;	
	XImage  *ximage ;
	XInfo   *xinfo_ptr ;
	Pixmap  pmap, pmap_zm ;
	int	ns, nl ;
	int     xstart, ystart ;
	int	zm_fac ;
	//int     zm_samps, zm_lines ;
	Dimension     zm_samps, zm_lines ;

	unsigned char *rgb ;

	// functions for displaying 
	int DispMono (int ns, int nl) ;

	// callback functions 
	static void expose_da (Widget w, XtPointer xinfoeq, void *cbs) ;
	void expose (Widget w, XtPointer xinfoeq, void *cbs) ;

	static void resize_da (Widget w, XtPointer xinfoeq, void *cbs) ;
	void resize (Widget w, XtPointer xinfoeq, void *cbs) ;

	static void click_da (Widget w, XtPointer xinfoeq, void *cbs) ;
	virtual void click (Widget w, XtPointer xinfoeq, void *cbs) ;

	static void click_dazm (Widget w, XtPointer xinfoeq, void *cbs) ;
	virtual void click_zm (Widget w, XtPointer xinfoeq, void *cbs) ;

	static void expose_dazm (Widget w, XtPointer xinfoeq, void *cbs) ;
	virtual void expose_zm (Widget w, XtPointer xinfoeq, void *cbs) ;



}
 
;
#endif 
