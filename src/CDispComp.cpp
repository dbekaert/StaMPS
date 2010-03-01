#include "CDispComp.h"
#include "CDisp.h"
#include <math.h>
#include <stdio.h>

using namespace std;

CDispComp::CDispComp(Widget w, char *ifile, XInfo *xi) : CDisp(w, ifile, xi) 
{
	// the rg file is the famous Zebker amp file
	phase = NULL ;
	mag = NULL ;
	both = NULL ;
	DispArr = 2 ;
}

int CDispComp::LoadArrays (unsigned char *mg, unsigned char *ph, int nsamps, int nlines) 
{
	// phase runs from -PI to PI, we map this from 0 to 255 
	// this function creates the arrays which are used to create the ximage arrays
	int 	i, j, phval  ;
	float 	ival ;
	float 	red [360], grn[360], blue[360] ; // arrays for the phase lut
	float   mgval ;
	double PI, PI2 ;
	ns = nsamps ;
	nl = nlines ;

	PI = atan (1.) * 4. ;
	PI2 = atan (1.) * 2. ;
	
	// assign values to the phase:
	for (i=0; i<120; i++)
	{
		red [i] = float(i) * 2.13 * 155./255. + 100 ;
		grn [i] = float (119. -i) * 2.13 * 155./255. + 100. ;
		blue [i] = 255. ;
	}
	for (i=120; i<240; i++)
	{
	 	ival = i - 120. ;
		red [i] = 255. ;
		grn [i] = float (ival) *  2.13 * 155./255. + 100. ;
		blue [i] = float (239. - i) * 2.13 * 155./ 255. + 100. ; 
	}
	for (i=240; i<360; i++)
	{
	 	ival = i - 240. ;
		red [i] = float (359. - i)  * 2.13 * 155. / 255. + 100. ;
		grn [i] = 255. ;
		blue [i] = float (ival) * 2.13 * 155. / 255. + 100. ;
	}

	phase = new unsigned char [ns * nl * 3] ;
	mag = new unsigned char [ns * nl * 3] ;
	both = new unsigned char [ns * nl * 3] ;
	ximage_phase = XCreateImage (xinfo_ptr->dpy, xinfo_ptr->visual,
                xinfo_ptr->depth,
                ZPixmap, 0,
                (char *) phase, ns, nl, 8, ns * 3) ; 
	ximage_mag = XCreateImage (xinfo_ptr->dpy, xinfo_ptr->visual,
                xinfo_ptr->depth,
                ZPixmap, 0,
                (char *) mag, ns, nl, 8, ns * 3) ; 
	ximage_both = XCreateImage (xinfo_ptr->dpy, xinfo_ptr->visual,
                xinfo_ptr->depth,
                ZPixmap, 0,
                (char *) both, ns, nl, 8, ns * 3) ; 
	ximage_phase->bits_per_pixel=24 ;
	ximage_mag->bits_per_pixel=24 ;
	ximage_both->bits_per_pixel=24 ;

	if(ImageByteOrder(xinfo_ptr->dpy) == LSBFirst){
	   for (i=0; i<nl; i++) {
		for (j=0; j<ns; j++) {
			phval = int(float (*(ph + i * ns + j)) /255. * 360. +0.5) ;
			mgval = *(mg + i * ns + j) ;
			*(phase+i * ns * 3+j*3 +2) = (unsigned char) (red [phval]) ;
			*(phase+i * ns * 3+j*3 +1) = (unsigned char) grn [phval] ;
			*(phase+i * ns * 3+j*3 +0) = (unsigned char) blue [phval] ;
			*(mag+i * ns * 3+j*3 +2) = (unsigned char)mgval ;
			*(mag+i * ns * 3+j*3 +1) = (unsigned char)mgval ;
			*(mag+i * ns * 3+j*3 +0) = (unsigned char)mgval ;
			*(both+i * ns * 3+j*3+2) = (unsigned char)(red [phval] * (float) mgval / 255.) ;
			*(both+i * ns * 3+j*3 +1) = (unsigned char) (grn [phval] * (float) mgval / 255.) ;
			*(both+i * ns * 3+j*3 +0) = (unsigned char) (blue [phval] * (float)mgval / 255.) ;
		}
	   }
	} else {
	   for (i=0; i<nl; i++) {
		for (j=0; j<ns; j++) {
			phval = int(float (*(ph + i * ns + j)) /255. * 360. +0.5) ;
			mgval = *(mg + i * ns + j) ;
			*(phase+i * ns * 3+j*3 +0) = (unsigned char) (red [phval]) ;
			*(phase+i * ns * 3+j*3 +1) = (unsigned char) grn [phval] ;
			*(phase+i * ns * 3+j*3 +2) = (unsigned char) blue [phval] ;
			*(mag+i * ns * 3+j*3 +0) = (unsigned char)mgval ;
			*(mag+i * ns * 3+j*3 +1) = (unsigned char)mgval ;
			*(mag+i * ns * 3+j*3 +2) = (unsigned char)mgval ;
			*(both+i * ns * 3+j*3 +0) = (unsigned char)(red [phval] * (float) mgval / 255.) ;
			*(both+i * ns * 3+j*3 +1) = (unsigned char) (grn [phval] * (float) mgval / 255.) ;
			*(both+i * ns * 3+j*3 +2) = (unsigned char) (blue [phval] * (float)mgval / 255.) ;
		}
	   }
	}
	return (1) ;
}
			
int CDispComp::WriteArrays (unsigned char *mg, unsigned char *ph, int nsamps, int nlines) 
{
	// phase runs from -PI to PI, we map this from 0 to 255 
	// this function creates the arrays which are used to create the ximage arrays
	int 	i, j, phval  ;
	float 	ival ;
	float 	red [360], grn[360], blue[360] ; // arrays for the phase lut
	float   mgval ;
	double PI, PI2 ;
	FILE *Of;
	ns = nsamps ;
	nl = nlines ;

	PI = atan (1.) * 4. ;
	PI2 = atan (1.) * 2. ;
	
	// assign values to the phase:
	for (i=0; i<120; i++)
	{
		red [i] = float(i) * 2.13 * 155./255. + 100 ;
		grn [i] = float (119. -i) * 2.13 * 155./255. + 100. ;
		blue [i] = 255. ;
	}
	for (i=120; i<240; i++)
	{
	 	ival = i - 120. ;
		red [i] = 255. ;
		grn [i] = float (ival) *  2.13 * 155./255. + 100. ;
		blue [i] = float (239. - i) * 2.13 * 155./ 255. + 100. ; 
	}
	for (i=240; i<360; i++)
	{
	 	ival = i - 240. ;
		red [i] = float (359. - i)  * 2.13 * 155. / 255. + 100. ;
		grn [i] = 255. ;
		blue [i] = float (ival) * 2.13 * 155. / 255. + 100. ;
	}

	//test arrays of color
        //Of = fopen("rgbtable.txt","w") ; 
        //if (Of == NULL){
        //    fprintf(stderr, "file open problem!\n");
        //    return(1);
        //}
        //for (i=0; i<360;i++)
        //{
        //    fprintf(Of,"%d %f %f %f\n",i,(double)red[i],grn[i],blue[i]) ;
	//}
	//fclose(Of);


	phase = new unsigned char [ns * nl * 3] ;
	mag = new unsigned char [ns * nl * 3] ;
	both = new unsigned char [ns * nl * 3] ;

	if(ImageByteOrder(xinfo_ptr->dpy) == LSBFirst){
	   for (i=0; i<nl; i++) {
		for (j=0; j<ns; j++) {
			phval = int(float (*(ph + i * ns + j)) /255. * 360. +0.5) ;
			mgval = *(mg + i * ns + j) ;
			*(phase+i * ns * 3+j*3 +2) = (unsigned char) (red [phval]) ;
			*(phase+i * ns * 3+j*3 +1) = (unsigned char) grn [phval] ;
			*(phase+i * ns * 3+j*3 +0) = (unsigned char) blue [phval] ;
			*(mag+i * ns * 3+j*3 +2) = (unsigned char)mgval ;
			*(mag+i * ns * 3+j*3 +1) = (unsigned char)mgval ;
			*(mag+i * ns * 3+j*3 +0) = (unsigned char)mgval ;
			*(both+i * ns * 3+j*3+2) = (unsigned char)(red [phval] * (float) mgval / 255.) ;
			*(both+i * ns * 3+j*3 +1) = (unsigned char) (grn [phval] * (float) mgval / 255.) ;
			*(both+i * ns * 3+j*3 +0) = (unsigned char) (blue [phval] * (float)mgval / 255.) ;
		}
	   }
	} else {
	   for (i=0; i<nl; i++) {
		for (j=0; j<ns; j++) {
			phval = int(float (*(ph + i * ns + j)) /255. * 360. +0.5) ;
			mgval = *(mg + i * ns + j) ;
			*(phase+i * ns * 3+j*3 +0) = (unsigned char) (red [phval]) ;
			*(phase+i * ns * 3+j*3 +1) = (unsigned char) grn [phval] ;
			*(phase+i * ns * 3+j*3 +2) = (unsigned char) blue [phval] ;
			*(mag+i * ns * 3+j*3 +0) = (unsigned char)mgval ;
			*(mag+i * ns * 3+j*3 +1) = (unsigned char)mgval ;
			*(mag+i * ns * 3+j*3 +2) = (unsigned char)mgval ;
			*(both+i * ns * 3+j*3 +0) = (unsigned char)(red [phval] * (float) mgval / 255.) ;
			*(both+i * ns * 3+j*3 +1) = (unsigned char) (grn [phval] * (float) mgval / 255.) ;
			*(both+i * ns * 3+j*3 +2) = (unsigned char) (blue [phval] * (float)mgval / 255.) ;
		}
	   }
	}

	//write to disk
        Of = fopen("dismph.dat","w") ; 
        if (Of == NULL){
            fprintf(stderr, "file open problem!\n");
            return(1);
        }
        fwrite((char *)both, sizeof(char), nl*ns*3, Of) ;
	fclose(Of);

        Of = fopen("dismph.ppm","w") ; 
        if (Of == NULL){
            fprintf(stderr, "file open problem!\n");
            return(1);
        }
        fprintf(Of,"P6 %d %d 255\n",ns,nl);
        fwrite((char *)both, sizeof(char), nl*ns*3, Of) ;
	fclose(Of);



	return (1) ;
}
			
 
int CDispComp::StartDisp() 
{
	int samps ;
	int lines ;
 

 
        samps = ns ;
        lines = nl ;
        XtVaSetValues (da, XmNwidth, samps, XmNheight, lines, NULL) ;
	if (pmap) XFreePixmap (xinfo_ptr->dpy, pmap) ;
        pmap = XCreatePixmap (xinfo_ptr->dpy, xinfo_ptr->rootwin, samps, lines, 24) ;
        XPutImage (xinfo_ptr->dpy, pmap, xinfo_ptr->imgGC, ximage_both, 0, 0, 0, 0, samps, lines) ;
 
        return (1) ;
}   
		

CDispComp::~CDispComp () {
	XDestroyImage (ximage_phase) ;
	XDestroyImage (ximage_mag) ;
	XDestroyImage (ximage_both) ;
/*
	delete [] phase ;
	delete [] mag ;
	delete [] both ;
*/
	if (pmap) XFreePixmap (xinfo_ptr->dpy, pmap) ;
	if (pmap_zm) XFreePixmap (xinfo_ptr->dpy, pmap_zm) ;
	pmap = 0 ;
	pmap_zm = 0 ;
}


void  CDispComp::expose_zm (Widget w, XtPointer xinfoeq, void *cbs)
{
        unsigned char *rgbzm, *iptr, *optr, *inarr=NULL ;
        int     ib, i, j, izm, jzm ;
        long    isamploc, osamploc ;
        XImage  *ximage_zm ;
 
        // using the zoom factor, the start coords, and the zoom window size
        XtVaGetValues (da_zm, XmNwidth, &zm_samps, XmNheight, &zm_lines, NULL) ;
        int fr_samps = zm_samps / zm_fac ;
        int fr_lines = zm_lines / zm_fac ;
 
 
        // allocate the memory for the data array to hold the zoomed memory
	if (DispArr == 0) inarr = phase ;
	if (DispArr == 1) inarr = mag ;
	if (DispArr == 2) inarr = both ;
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

void  CDispComp::click (Widget w, XtPointer xinfoeq, void *cbs)
{
        char cursbuf [80] ;
        int    xloc, yloc ;
	int    samps, lines ;
	int	xbuffer, ybuffer ;
        XmDrawingAreaCallbackStruct *dcbs ;
        dcbs = (XmDrawingAreaCallbackStruct *) cbs ;
        XEvent  *event =  dcbs->event ;
        xloc = event->xbutton.x ;
        yloc = event->xbutton.y ;

	samps = ns ;
	lines = nl ;

	if (event->xany.type==ButtonPress) return ;
	if (event->xbutton.button==3) {
		DispArr++ ; 
		if (DispArr > 2) DispArr = 0 ; 
		switch (DispArr) {
		    case 0 :
			XPutImage (xinfo_ptr->dpy, pmap, xinfo_ptr->imgGC, ximage_phase, 0, 0, 0, 0,
                		samps, lines) ; 
			break ;
		    case 1 :
			XPutImage (xinfo_ptr->dpy, pmap, xinfo_ptr->imgGC, ximage_mag, 0, 0, 0, 0,
                		samps, lines) ; 
			break ;
		    case 2 :
			XPutImage (xinfo_ptr->dpy, pmap, xinfo_ptr->imgGC, ximage_both, 0, 0, 0, 0,
                		samps, lines) ; 
			break ;
		}
        	XClearArea (xinfo_ptr->dpy, XtWindow (da), 0, 0, 0, 0, True) ;
	}
 
	cout << "xloc  yloc : " << xloc << "  " << yloc << endl ;
        sprintf (cursbuf, "X :   %5d      Y :    %5d", xloc, yloc) ;
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
 
void CDispComp::click_da (Widget w, XtPointer xinfoeq, void *cbs)
{
        CDispComp *thisptr = (CDispComp *) xinfoeq ;
        thisptr->click (w, xinfoeq, cbs) ;
}             

void CDispComp::expose_dazm (Widget w, XtPointer xinfoeq, void *cbs)
{
        CDispComp *thisptr = (CDispComp *) xinfoeq ;
        thisptr->expose_zm (w, xinfoeq, cbs) ;
}
 
/*
void CDisp::expose_da (Widget w, XtPointer xinfoeq, void *cbs)
{
        CDispComp *thisptr = (CDispComp *) xinfoeq ;
        thisptr->expose (w, xinfoeq, cbs) ;
}
*/
/*
void CDispComp::resize_da (Widget w, XtPointer xinfoeq, void *cbs)
{
        CDisp *thisptr = (CDisp *) xinfoeq ;
        thisptr->resize (w, xinfoeq, cbs) ;
}
*/
