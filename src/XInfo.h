#include <Xm/RowColumn.h>
#include <Xm/Text.h>
#include <iostream>

#ifndef XIdef
#define XIdef

class XInfo {
	public :
		XtAppContext 	app ;
		XVisualInfo  	vinfo ;
		GC	     	imgGC ;
		Visual 		*visual ;
		Display		*dpy ;
		Screen 		*screen ;

		Window 		rootwin ;
		Widget		toplevel ;
		Widget		curswin ;
		Widget		rowcol ;
		int		scr_num ;
		Cardinal	depth ;
	XInfo (int argc, char *argv[]) ;
	~XInfo () ;
	int startloop () ;
}


;
#endif
