#include <Xm/Protocols.h>
#include "XInfo.h"

using namespace std;

void closeup_fb (Widget w, XtPointer client_data, XmAnyCallbackStruct *cbs) ; 

XInfo::XInfo (int argc, char *argv[]) 
{
	Atom  WM_DELETE_WINDOW ;
	toplevel = XtAppInitialize (&app, "Programs", NULL,
                0, &argc, argv, NULL, NULL, 0) ; 
/*
	XtVaSetValues (toplevel, 
		XmNtitle, "Blipppy",
		NULL) ;
*/

	rowcol = XtVaCreateManagedWidget ("rowcol",
		xmRowColumnWidgetClass, toplevel,
		NULL ) ;

	curswin = XtVaCreateManagedWidget ("textf",
                xmTextWidgetClass, rowcol,
                XmNwidth, 300,
                XmNeditMode, XmMULTI_LINE_EDIT,
                XmNeditable, False,
                XmNrows, 10,
                NULL) ;   

	// get the depth and default visual
	dpy = XtDisplay (toplevel) ;
        screen = XtScreen (toplevel) ;
        scr_num = DefaultScreen (dpy) ;
        rootwin = RootWindowOfScreen (screen) ;
        depth = DefaultDepthOfScreen (screen) ;   
	imgGC = XCreateGC (dpy, rootwin, 0, NULL) ; 
	// now need to get a visual, grab the default visual
        visual = DefaultVisual (dpy, scr_num) ;
        int status = XMatchVisualInfo (dpy, scr_num, depth, TrueColor,
                &vinfo) ;
	if (status == 0) {
		cout << "Count not get True Color Visual " << endl ;
		return ;
	}
 
        cout << vinfo.c_class << endl ;
        if (vinfo.c_class == TrueColor)
        {
                cout << "True Color visual is grabbed " << endl ;
                cout << "Depth is : " << depth << endl ;
        }
	WM_DELETE_WINDOW = XmInternAtom (XtDisplay (toplevel),
                "WM_DELETE_WINDOW", False) ;  
	XmAddWMProtocolCallback (toplevel, WM_DELETE_WINDOW,
                XtCallbackProc (closeup_fb),
                this) ; 
 
}

XInfo::~XInfo () 
{
	XtDestroyWidget (toplevel) ;
}

int XInfo::startloop () {
        XtRealizeWidget (toplevel) ;
        XtAppMainLoop (app) ;
	cout << "Status is 0" << endl ;
	return (0) ;
}  

