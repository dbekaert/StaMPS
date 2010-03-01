//  CDisp class will display an array to a frame containing a drawing area in
//  a scroll window
#include "CDisp.h"
#include "CDispMag.h"
#include <stdio.h>
#include <Xm/Protocols.h>
#include <Xm/Xm.h>
#include <X11/Core.h>

using namespace std;

FILE *cdispinputfile;

void closeup_fb (Widget w, XtPointer client_data, XmAnyCallbackStruct *cbs) ;



CDisp::CDisp()
{
        ximage = NULL ;
}

CDisp::~CDisp ()
{
	if (ximage) XDestroyImage (ximage) ;
	if (dispshell) XtPopdown (dispshell) ;
	if (zoomshell) XtPopdown (zoomshell) ;
	if (pmap) XFreePixmap (xinfo_ptr->dpy, pmap) ;
}

CDisp::CDisp (Widget toplevel, char *dispfile, XInfo  *xiptr)
{
	char temptitle [240] ;
	Atom WM_DELETE_WINDOW ;

        cdispinputfile = fopen(dispfile,"r");

        rgb =NULL ;
        pmap=0 ;
        pmap_zm = 0 ;
        xstart = 0 ;
        ystart = 0 ;
        zm_fac = 4 ;
        zm_samps=0 ;
        zm_lines=0 ;
        ximage = NULL ;
        // given the top level window, generate the top shell containing
        // the drawing area, etc.

	xinfo_ptr = xiptr ;

	XtVaSetValues (xinfo_ptr->toplevel,
		XmNtitle, dispfile,
		NULL) ;

        dispshell = XtVaCreatePopupShell ("Topshell",
                topLevelShellWidgetClass, xinfo_ptr->toplevel,
                XmNtitle, dispfile,
                NULL) ;


        // make scroll window
        sw = XtVaCreateManagedWidget ("scrollwin",
                xmScrolledWindowWidgetClass, dispshell,
                XmNwidth, 600,
                XmNheight, 400,
                XmNscrollingPolicy, XmAUTOMATIC,
                NULL) ;

        // put drawing area in scroll window
	da = XmCreateDrawingArea (sw, "test", NULL, 0) ;
	XtVaSetValues (da, XmNtraversalOn, False,
                XmNheight, 512,
                XmNwidth, 512,
                NULL) ; 
	XtAddCallback (da, XmNinputCallback, XtCallbackProc (click_da), this ) ; 
	XtAddCallback (da, XmNexposeCallback, XtCallbackProc (expose_da), this ) ; 

	sprintf (temptitle, "Zoom : %s", dispfile) ;
	zoomshell = XtVaCreatePopupShell ("topshell",
                topLevelShellWidgetClass, xinfo_ptr->toplevel,
                XmNtitle, "ZOOOM",
                NULL) ;
	
	XtVaSetValues (zoomshell,
		XmNtitle, temptitle,
		NULL) ;
	da_zm = XmCreateDrawingArea (zoomshell, "test", NULL, 0) ;
	XtVaSetValues (da_zm, XmNtraversalOn, False,
                XmNheight, 256,
                XmNwidth, 256,
                NULL) ; 
	XtAddCallback (da_zm, XmNinputCallback, XtCallbackProc (click_dazm), this ) ; 
	XtAddCallback (da_zm, XmNexposeCallback, XtCallbackProc (expose_dazm), this ) ; 
	
	WM_DELETE_WINDOW = XmInternAtom (XtDisplay (toplevel),
                "WM_DELETE_WINDOW", False) ;
        XmAddWMProtocolCallback (dispshell, WM_DELETE_WINDOW,
                XtCallbackProc (closeup_fb),
                this) ;
        XmAddWMProtocolCallback (zoomshell, WM_DELETE_WINDOW,
                XtCallbackProc (closeup_fb),
                this) ;


	XtManageChild (da) ;
	XtManageChild (da_zm) ;
	XtPopup (dispshell, XtGrabNone) ;
	XtPopup (zoomshell, XtGrabNone) ;
}

void CDisp::expose (Widget w, XtPointer  xinfoeq, void *cbs)
{
        //int samps=0 ;
        Dimension samps=0 ;
	//int lines=0 ;
	Dimension lines=0 ;
        XtVaGetValues (da, XmNwidth, &samps, XmNheight, &lines, NULL) ;
	XCopyArea (xinfo_ptr->dpy, pmap, XtWindow (w), xinfo_ptr->imgGC,
               0, 0, samps, lines, 0, 0) ;
}

	

void  CDisp::click_zm (Widget w, XtPointer xinfoeq, void *cbs)
{  
	char cursbuf [80] ;
	int    xloc, yloc, newx, newy ;
	int    xbuffer, ybuffer ;
	XmDrawingAreaCallbackStruct *dcbs ;
	dcbs = (XmDrawingAreaCallbackStruct *) cbs ;
	XEvent  *event =  dcbs->event ;

        if (event->xany.type==ButtonPress) return ;
	xloc = event->xbutton.x ;
	yloc = event->xbutton.y ;
	newx = xstart + xloc / zm_fac ;
	newy = ystart + yloc / zm_fac ;
	xloc = newx ;
	yloc = newy ;
        cout << "xloc  yloc : " << xloc << "  " << yloc << endl ;

        sprintf (cursbuf, "X :   %5d      Y :    %5d", xloc, yloc) ;

        XmTextSetString (xinfo_ptr->curswin, cursbuf) ;
        XtVaGetValues (da_zm, XmNwidth, &zm_samps, XmNheight, &zm_lines, NULL) ;
            xstart = xloc - zm_samps / zm_fac / 2 ;
            ystart = yloc - zm_lines / zm_fac / 2 ;
	xbuffer = zm_samps / zm_fac  ; 
	ybuffer = zm_lines / zm_fac  ; 
	xstart = (xstart<0) ? 0:xstart ;
	xstart = (xstart>ns - xbuffer-1) ? ns - xbuffer-1 : xstart ;
	ystart = (ystart<0) ? 0:ystart ;
	ystart = (ystart>nl - ybuffer-1) ? nl - ybuffer-1 : ystart ;
        XClearArea (xinfo_ptr->dpy, XtWindow (da_zm), 0, 0, 0, 0, True) ;

}

void  CDisp::click (Widget w, XtPointer xinfoeq, void *cbs)
{
        char cursbuf [80] ;
        int    xloc, yloc ;
        int    samps, lines ;
	int    xbuffer, ybuffer ;
        XmDrawingAreaCallbackStruct *dcbs ;

	float value[1];
	unsigned char  bvalue[1];

        dcbs = (XmDrawingAreaCallbackStruct *) cbs ;
        XEvent  *event =  dcbs->event ;
        xloc = event->xbutton.x ;
        yloc = event->xbutton.y ;

        samps = ns ;
        lines = nl ;

        if (event->xany.type==ButtonPress) return ;

        // cout << "xloc  yloc : " << xloc << "  " << yloc << endl ;
        fseek ( cdispinputfile, yloc*samps*4+xloc*4, SEEK_SET);
        fread (value,4,1,cdispinputfile);
        fseek ( cdispinputfile, yloc*samps+xloc, SEEK_SET);
        fread (bvalue,1,1,cdispinputfile);

	printf ("xloc  yloc : %d %d, Value: %f (%d)\n",xloc,yloc, value[0],bvalue[0]);
	sprintf (cursbuf, "X, Y :  %5d  %5d   Value %f (%d)", xloc, yloc, value[0], bvalue[0]) ;
        XmTextSetString (xinfo_ptr->curswin, cursbuf) ;
        XtVaGetValues (da_zm, XmNwidth, &zm_samps, XmNheight, &zm_lines, NULL) ;
        xstart = xloc - zm_samps / zm_fac / 2 ;
        ystart = yloc - zm_samps / zm_fac / 2 ;
	xbuffer = zm_samps / zm_fac  ; 
	ybuffer = zm_lines / zm_fac  ; 
	xstart = (xstart<0) ? 0:xstart ;
        xstart = (xstart>ns - xbuffer-1) ? ns - xbuffer-1 : xstart ;
        ystart = (ystart<0) ? 0:ystart ;
        ystart = (ystart>nl - ybuffer-1) ? nl - ybuffer-1 : ystart ; 
	XClearArea (xinfo_ptr->dpy, XtWindow (da_zm), 0, 0, 0, 0, True) ;



}

void  CDisp::expose_zm (Widget w, XtPointer xinfoeq, void *cbs)
{
        unsigned char *rgbzm, *iptr, *optr, *inarr ;
        int     ib, i, j, izm, jzm ;
        long    isamploc, osamploc ;
        XImage  *ximage_zm ;

        // using the zoom factor, the start coords, and the zoom window size
        XtVaGetValues (da_zm, XmNwidth, &zm_samps, XmNheight, &zm_lines, NULL) ;
        int fr_samps = zm_samps / zm_fac ;
        int fr_lines = zm_lines / zm_fac ;


        // allocate the memory for the data array to hold the zoomed memory
        inarr = rgb ;
        rgbzm = new unsigned char [zm_samps * zm_lines * 3] ;
        for (i=0; i<fr_lines; i++) {
                for (j=0; j<fr_samps; j++) {
                        isamploc = (ystart + i) * ns * 3 + (xstart + j) * 3  ;
                        for (izm= 0; izm<zm_fac; izm++) {
                                for (jzm=0; jzm<zm_fac; jzm++) {
					osamploc = (i * zm_fac +  izm) * zm_samps * 3 +
                                                (j * 3* zm_fac + jzm * 3) ;
                                        for (ib=0; ib<3; ib++) {
                                                iptr = inarr + isamploc + ib ;
                                                optr = rgbzm + osamploc + ib ;
                                                *optr = *iptr ;
                                        }
                                }
                        }
                }
        }
        ximage_zm = XCreateImage (xinfo_ptr->dpy, xinfo_ptr->visual,
                xinfo_ptr->depth,
                ZPixmap, 0,
                (char *) rgbzm, zm_samps, zm_lines, 8, zm_samps * 3) ;
        ximage_zm->bits_per_pixel= 24 ;
        if (pmap_zm)
                XFreePixmap (xinfo_ptr->dpy, pmap_zm) ;
        pmap_zm = XCreatePixmap (xinfo_ptr->dpy, xinfo_ptr->rootwin, zm_samps, zm_lines, 24) ;
        XPutImage (xinfo_ptr->dpy, pmap_zm, xinfo_ptr->imgGC, ximage_zm, 0, 0, 0, 0,
                zm_samps, zm_lines) ;
        XCopyArea (xinfo_ptr->dpy, pmap_zm, XtWindow (w), xinfo_ptr->imgGC,
               0, 0, zm_samps, zm_lines, 0, 0) ;
        XDestroyImage (ximage_zm) ;

}

void CDisp::click_dazm (Widget w, XtPointer xinfoeq, void *cbs) 
{
	CDisp *cditmp = (CDisp *) xinfoeq ;
	cditmp->click_zm (w, xinfoeq, cbs) ;
}

void CDisp::click_da (Widget w, XtPointer xinfoeq, void *cbs)
{
        CDisp *thisptr = (CDisp *) xinfoeq ;
        thisptr->click (w, xinfoeq, cbs) ;
}

void CDisp::expose_dazm (Widget w, XtPointer xinfoeq, void *cbs)
{
        CDisp *thisptr = (CDisp *) xinfoeq ;
        thisptr->expose_zm (w, xinfoeq, cbs) ;
}

void CDisp::expose_da (Widget w, XtPointer xinfoeq, void *cbs)
{
        CDisp *thisptr = (CDisp *) xinfoeq ;
        thisptr->expose (w, xinfoeq, cbs) ;
}
