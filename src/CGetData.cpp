#include "CGetData.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>

using namespace std;

void bytescale (unsigned char *indat, unsigned char *outdat, int ns, int nl, float
        scale, float offset, int dtype, int flipflag) ;

int CGetData::bpp [4] = {1, 2, 4, 4} ;

CGetData::CGetData () {
	red = NULL ;
	grn = NULL ;
	mag = NULL ;
	phase = NULL ;
	hgt = NULL ;
	raw = NULL ;
	} ;

int CGetData::setparams (char *ifile, int ns, int startl, int nl) 
{
	strcpy (infile, ifile) ;
	nsamps = ns ;
	startline = startl ;
	nlines   = nl ;

	ifstream ifil (infile, ios::in) ;
	if (ifil.bad()) {
		cout << "Problem with : " << infile << endl ;
		return (-1) ;
	}
	
	ifil.close() ;
	return (1) ;


	
	
}

CGetData::~CGetData () {
	if (red) delete [] red ;
	if (grn) delete [] grn ;
	if (mag) delete [] mag ;
	if (phase) delete [] phase ;
	if (hgt) delete [] hgt ;
	if (raw) delete [] raw ;
}

int CGetData::getarrayMag (int dtype, float minval, float maxval, int flipflag) 
{
	// get a byte array based upon the min and max value
        unsigned char 	*bytearr ; 
        unsigned char  	*temparr ; 
        int  		i, j, i_out, j_out, hflip=0, vflip=0 ;
	int  		bytesperpixel ;
        long 		lpos ;
	float 		scalemag, offmag ;

        ifstream ifil (infile, ios::in) ;
 
	bytesperpixel = bpp [dtype] ;

	// scale and offsets for going from orig data to byte array
	// first subtract offmag then mult by scalemag
	scalemag = 255. / (maxval - minval) ;
	offmag = -minval * scalemag ;
	cout << "Scale value  : " << scalemag << endl ;
	cout << "Offset value : " << offmag << endl ;
	//
        // exponent apply flag and exponent value are class members

	if (flipflag)
	{
		hflip = flipflag %2 ;
		vflip = flipflag /2 ;
	}

	// allocate the input array and the output arrays
	temparr = new unsigned char [nsamps * bytesperpixel] ;
	bytearr = new unsigned char [nsamps] ;

        ifil.seekg (long(0), ios::end ) ;
        lpos = ifil.tellg () ;
        totlines = lpos / (nsamps * bytesperpixel) ;
 
	if (nlines > totlines-startline) {
                nlines = totlines - startline ;
	}
        if (nlines ==0) {
                nlines = totlines - startline ;
        }
	mag = new unsigned char [nsamps * nlines] ;
        ifil.seekg (long(startline) * nsamps * bytesperpixel, ios::beg) ;
 
        for (i =startline; i<nlines; i++) {
		// read in a line of data then convert
                ifil.read ((char*)temparr, nsamps * bytesperpixel) ;
		i_out = (vflip)? nlines -i - 1: i ;
		bytescale (temparr, bytearr, nsamps, 1, scalemag, offmag, dtype, 0) ;
		for (j=0; j<nsamps; j++) {
			j_out = (hflip)? nsamps -j - 1: j ;
			*(mag + i_out * nsamps + j_out) = *(bytearr + j) ;
		}
	}
	delete [] temparr ;
	delete [] bytearr ;

	ifil.close() ;
	return (1) ;
}
		


	
int CGetData::getarrayMag (float scalefac, int expfl, float expv, int flipflag, int fifthflag)
{
	// get the magnitude from the standard complex array
        unsigned char *optr0 ;
        int  i, j, i_out, j_out, count=0, hflip=0, vflip=0 ;
        long lpos ;
        float *temparr,  *iptr0, *iptr1 ;
	float scalemag ;
        float b1 ;
        double totmag, magval, redval ;
        ifstream ifil (infile, ios::in) ;
 
        // exponent apply flag and exponent value are class members
        expflag = expfl ;
        expval = expv ;

	if (flipflag)
	{
		hflip = flipflag %2 ;
		vflip = flipflag /2 ;
	}
 
        ifil.seekg (long(0), ios::end ) ;
        lpos = ifil.tellg () ;
        totlines = lpos / (nsamps * 8) ;
 
	if (nlines > totlines-startline) {
                nlines = totlines - startline ;
	}
        if (nlines ==0) {
                nlines = totlines - startline ;
        }
 
        mag = new unsigned char [nsamps * nlines] ;
        temparr = new float [nsamps * 2] ;
 
        // get the scaling parameters
        if (expflag) cout << "Applying exponent : " << expval << endl ;
        if (!expflag) cout << "Not Applying exponent  " << endl  ;
        totmag = 0. ;
        for (i =startline; i<nlines; i+=32) {
                ifil.seekg (long(i) * nsamps * 8, ios::beg) ;
                ifil.read ((char*)temparr, nsamps * 8) ;
                for (j=32; j<nsamps-32; j+= 32) {
                        count += 1 ;
			magval = getmag (*(temparr+ j*2), *(temparr + j*2+1));

                        if (expflag)
                        {
                        	totmag += pow ( double(magval), double(expval)) ;
                        }
                        else {
                        	totmag += magval ;
                        }
                }
        }
 
        // get the avereage of the red and grn
        totmag /= count ;
	scalemag = scalefac * 127. / totmag ;
 
        cout << "Scaling parameter   :   "   << scalemag << endl ;
 
        ifil.seekg (startline * nsamps * 8, ios::beg) ;
        for (i=0; i<nlines; i++) {
                ifil.read ((char*)temparr, nsamps * 8) ;
                for (j=0; j<nsamps; j++)
                {
			i_out = (vflip) ? nlines -i  -1 : i ;
			j_out = (hflip) ? nsamps -j  -1 : j ;
                        optr0 = mag + i_out * nsamps + j_out ;
                        iptr0 = temparr + j * 2 ;
                        iptr1 = temparr + j * 2 + 1;
			magval = getmag (*iptr0, *iptr1) ;
                        if (expflag) {
                                b1 = pow (double(magval), double(expval)) ;
                        }
                        else {
                                b1 = magval ;
                        }
                        magval = b1 * scalemag ;
			redval = magval ;
                        redval = (redval < 0) ? 0. : redval ;
                        redval = (redval > 255) ? 255. : redval ;
                        *optr0 = (unsigned char) redval ;
                }
        }
        delete [] temparr ;
        ifil.close() ;
        return (1) ;
}          

int CGetData::getarrayMPH (float scalefac, int expfl, float expv, int flipflag)
{
        unsigned char *optr0, *optr1 ;
        int  i, j, count=0, hflip=0, vflip=0, i_out, j_out;
        long lpos ;
        float *temparr,  *iptr0, *iptr1 ;
	float scalemag ;
        float b1 ;
        double totmag, magval=0, phaseval=0, redval, TWOPI ;
        ifstream ifil (infile, ios::in) ;
 
        // exponent apply flag and exponent value are class members
        expflag = expfl ;
        expval = expv ;

	if (flipflag) {
		hflip = flipflag % 2 ;	
		vflip = flipflag / 2 ;	
	}
	TWOPI = 8. * atan (1.) ;
 
        ifil.seekg (long(0), ios::end ) ;
        lpos = ifil.tellg () ;
        totlines = lpos / (nsamps * 8) ;
 
	if (nlines > totlines-startline) {
                nlines = totlines - startline ;
	}
        if (nlines ==0) {
                nlines = totlines - startline ;
        }
 
        mag = new unsigned char [nsamps * nlines] ;
	phase = new unsigned char [nsamps * nlines] ;
        temparr = new float [nsamps * 2] ;
 
        // get the scaling parameters for magnitude
        if (expflag) cout << "Applying exponent : " << expval << endl ;
        if (!expflag) cout << "Not Applying exponent  " << endl  ;
        totmag = 0. ;
        for (i =startline; i<startline+nlines; i+=32) {
                ifil.seekg (long(i) * nsamps * 8, ios::beg) ;
                ifil.read ((char*)temparr, nsamps * 8) ;
                for (j=32; j<nsamps-32; j+= 32) {
                        count += 1 ;
			magval = getmag (*(temparr+ j*2), *(temparr + j*2+1));

                        if (expflag)
                        {
                        	totmag += pow ( double(magval), double(expval)) ;
                        }
                        else {
                        	totmag += magval ;
                        }
                }
        }
 
        // get the avereage of the red and grn
        totmag /= count ;
	scalemag = scalefac * 150. / totmag ;
 
        cout << "Scaling parameter   :   "   << scalemag << endl ;
	double phasescale = 255. / TWOPI ; 
 
        ifil.seekg (startline * nsamps * 8, ios::beg) ;
        for (i=0; i<nlines; i++) {
                ifil.read ((char*)temparr, nsamps * 8) ;
                for (j=0; j<nsamps; j++)
                {
			i_out = (vflip) ? nlines -i  -1 : i ;
			j_out = (hflip) ? nsamps -j  -1 : j ;
                        optr0 = mag + i_out * nsamps + j_out ;
                        optr1 = phase + i_out * nsamps + j_out ;
                        iptr0 = temparr + j * 2 ;
                        iptr1 = temparr + j * 2 + 1;
			getmph (*iptr0, *iptr1, &magval, &phaseval) ;
			if (phaseval < 0.) phaseval += TWOPI ;
                        if (expflag) {
                                b1 = pow (double(magval), double(expval)) ;
                        }
                        else {
                                b1 = magval ;
                        }
                        magval = b1 * scalemag ;
			phaseval *= phasescale ; 
			redval = magval ;
                        redval = (redval < 0) ? 0. : redval ;
                        redval = (redval > 255) ? 255. : redval ;
                        *optr0 = (unsigned char) redval ;
			*optr1 = (unsigned char) (phaseval) ; 
                }
        }
        delete [] temparr ;
        ifil.close() ;
        return (1) ;
}          

int CGetData::getarrayHgt (float scalefac, int expfl, float expv, float cont_interval, int flipflag)
{
	// this reads in the bil float file with the back scatter as band
	// 1 and the height as band 2 
	// the arrays created are the height and the bs arrays as byte
	// a contour interval is applied to the height so that it can 
	// be treated in the same fashion as the cyclical phase
        unsigned char *optr1 ;
        int  i, j, count=0, hgtval, hflip=0, vflip=0, i_out, j_out ;
        long lpos ;
        float *temparr,  *iptr0, *iptr1 ;
	float scalemag, scalehgt ;
        float b1 ;
        double totmag, magval=0, redval, contours ;
        ifstream ifil (infile, ios::in) ;
 
        // exponent apply flag and exponent value are class members
        expflag = expfl ;
        expval = expv ;

	if (flipflag) {
		hflip = flipflag % 2 ;	
		vflip = flipflag / 2 ;	
	}
 
        ifil.seekg (long(0), ios::end ) ;
        lpos = ifil.tellg () ;
        totlines = lpos / (nsamps * 8) ;
 
	if (nlines > totlines-startline) {
                nlines = totlines - startline ;
	}
        if (nlines ==0) {
                nlines = totlines - startline ;
        }
 
        mag = new unsigned char [nsamps * nlines] ;
	hgt = new unsigned char [nsamps * nlines] ;
        temparr = new float [nsamps * 2] ;
 
        // get the scaling parameters for magnitude
        if (expflag) cout << "Applying exponent : " << expval << endl ;
        if (!expflag) cout << "Not Applying exponent  " << endl  ;
        totmag = 0. ;
        for (i =startline; i<startline+nlines; i+=32) {
                ifil.seekg (long(i) * nsamps * 4 * 2, ios::beg) ;
		// just read in the magnitude portion of the line (first band)
                ifil.read ((char*)temparr, nsamps * 4) ;
                for (j=32; j<nsamps-32; j+= 32) {
                        count += 1 ;
			magval = *(temparr+ j) ;

                        if (expflag)
                        {
                        	totmag += pow ( double(magval), double(expval)) ;
                        }
                        else {
                        	totmag += magval ;
                        }
                }
        }
 
        // get the avereage of the red and grn
        totmag /= count ;
	scalemag = scalefac * 150. / totmag ;
	scalehgt = 255. / (float) cont_interval ;
 
        cout << "Magnitude scaling parameter   :   "   << scalemag << endl ;
        cout << "Height scaling parameter      :   "   << scalehgt << endl ;

	double min=1.E18, max=-1.E18; 
 
        ifil.seekg (startline * nsamps * 4 * 2, ios::beg) ;
        for (i=0; i<nlines; i++) {
		// first read in the magnitude
                ifil.read ((char*)temparr, nsamps * 4) ;
                for (j=0; j<nsamps; j++)
                {
			i_out = (vflip) ? nlines -i  -1 : i ;
			j_out = (hflip) ? nsamps -j  -1 : j ;
                        optr1 = mag + i_out * nsamps + j_out ;
                        iptr0 = temparr + j  ;
			magval = *iptr0 ;
			//cout << magval << endl;

                        if (expflag) {
                                b1 = pow (double(magval), double(expval)) ;
                        }
                        else {
                                b1 = magval ;
                        }
                        magval = b1 * scalemag ;
			redval = magval ;
                        redval = (redval < 0) ? 0. : redval ;
                        redval = (redval > 255) ? 255. : redval ;
                        *optr1 = (unsigned char) redval ;
                }
		// then the height
                ifil.read ((char*)temparr, nsamps * 4) ;
                for (j=0; j<nsamps; j++)
                {
			i_out = (vflip) ? nlines -i  -1 : i ;
			j_out = (hflip) ? nsamps -j  -1 : j ;
                        optr1 = hgt + i_out * nsamps + j_out ;
                        iptr1 = temparr + j ;
			if (*iptr1 > max) max = *iptr1 ;
			if (*iptr1 < min) min = *iptr1 ;
			// added 18mar02 hz for real contour interval
			contours = (*iptr1) / cont_interval;
			contours = (contours - floor(contours)) * cont_interval;
			hgtval = int (contours * scalehgt);
			//  test lines
			//cout << (*iptr1) << "  " << hgtval << "  " << contours << endl;

			// ** hgtval = int (float(int(*iptr1) % cont_interval) * scalehgt) ; 
			hgtval = (hgtval < 0) ? 0 : hgtval ;
			hgtval = (hgtval > 255) ? 255 : hgtval ;
			*optr1 = (unsigned char) (hgtval) ; 
                }
        }
        delete [] temparr ;
	cout << "Min hgt is " << min << endl ;
	cout << "Max hgt is " << max << endl ;
        ifil.close() ;
        return (1) ;
}          

int CGetData::getarrayByte (int min, int max, int flipflag)
{
        unsigned char *optr0 ;
        int  i, j, i_out, j_out, vflip=0, hflip=0 ;
        long lpos ;
        unsigned char *temparr ;
	float scalebyte, offset=0. ;
        double redval ;
        ifstream ifil (infile, ios::in) ;
 
        // exponent apply flag and exponent value are class members
 
        ifil.seekg (long(0), ios::end ) ;
        lpos = ifil.tellg () ;
        totlines = lpos / (nsamps) ;
 
	if (flipflag) {
		hflip = flipflag % 2 ;	
		vflip = flipflag / 2 ;	
	}
	if (nlines > totlines-startline) {
                nlines = totlines - startline ;
	}
        if (nlines ==0) {
                nlines = totlines - startline ;
        }
 
        mag = new unsigned char [nsamps * nlines] ;
        temparr = new unsigned char [nsamps] ;
	scalebyte = 255 / (max -min)  ;
	offset = min ;
 
        // get the scaling parameters
        cout << "Applying scalefactor  : " << scalebyte << endl ;
        cout << "Applying offset       : " << offset  ;
 
 
 
        ifil.seekg (startline * nsamps, ios::beg) ;
        for (i=0; i<nlines; i++) {
                ifil.read ((char*)temparr, nsamps) ;
                for (j=0; j<nsamps; j++)
                {
			i_out = (vflip) ? nlines -i  -1 : i ;
			j_out = (hflip) ? nsamps -j  -1 : j ;
			optr0 = mag + i_out * nsamps + j_out ;
			redval = (*(temparr + j) - offset) * scalebyte ;
			redval = (redval <0)? 0 : redval ;
			redval = (redval >255)? 255 : redval ;
                        *optr0 = (unsigned char) redval ;
                }
        }
        delete [] temparr ;
        ifil.close() ;
        return (1) ;
}          

int CGetData::Raw2Mag (int dtype, float min, float max, int flipflag) {
	float scale, offset ;
	scale = 255. / (max - min) ; 
	offset = -min * scale ;
	
	cout << "Raw2Mag scaling " << endl ;
	cout << "Scale  : " << scale << endl ;
	cout << "Offset : " << offset << endl ;

        mag = new unsigned char [nsamps * nlines] ;
	bytescale (raw, mag, nsamps, nlines, scale, offset, dtype, flipflag) ;
	return (1) ;
}
	

int CGetData::getarrayRaw (int dtype) 
{
	int bytespersample ;
	long lpos, npix = 0 ;

        ifstream ifil (infile, ios::in) ;
 
	bytespersample = bpp [dtype] ;
 
        ifil.seekg (long(0), ios::end ) ;
        lpos = ifil.tellg () ;
        totlines = lpos / (nsamps * bytespersample) ;
	if (nlines > totlines-startline) {
                nlines = totlines - startline ;
	}
	if (nlines ==0) {
		nlines = totlines - startline ;
	}
 
	npix = nsamps * (totlines - startline) * bytespersample ;

	raw = new unsigned char [npix] ;
	if (raw == NULL) {
		cout << "Could not allocate memory for the raw data array " << endl ;
		return (-1) ;
	}

	ifil.seekg (startline * nsamps * bytespersample, ios::beg) ;
	ifil.read ((char*)raw, npix) ;
	ifil.close () ;
	return (1) ;
}
 

int CGetData::getarrayRG (float scalefac, int expfl, float expv, int flipflag) 
{
	unsigned char *optr0, *optr1 ;
	int  i, j, count=0, hflip=0, vflip=0, i_out, j_out;
	long lpos ; 
	float *temparr, redval, grnval, *iptr0, *iptr1 ;
	float b1, b2 ;
	double totred, totgrn ;
	ifstream ifil (infile, ios::in) ;

	// exponent apply flag and exponent value are class members
	expflag = expfl ;
	expval = expv ;

	if (flipflag) {
		hflip = flipflag % 2 ;	
		vflip = flipflag / 2 ;	
	}
	

	ifil.seekg (long(0), ios::end ) ;
	lpos = ifil.tellg () ;
	totlines = lpos / (nsamps * 8) ;
	
	if (nlines > totlines-startline) {
                nlines = totlines - startline ;
	}
	if (nlines ==0) {
		nlines = totlines - startline ;
	}
 
	red = new unsigned char [nsamps * nlines] ;
	grn = new unsigned char [nsamps * nlines] ;
	temparr = new float [nsamps * 2] ;

	// get the scaling parameters 
	if (expflag) cout << "Applying exponent : " << expval << endl ;
	if (!expflag) cout << "Not Applying exponent  " << endl  ;
	totred = 0. ;
	totgrn = 0. ;
	for (i =startline; i<startline+nlines; i+=32) {
		ifil.seekg (long(i) * nsamps * 8, ios::beg) ;
		ifil.read ((char*)temparr, nsamps * 8) ;
		for (j=32; j<nsamps-32; j+= 32) {
			count += 1 ;
			if (expflag) 
			{
			totred += pow (double (*(temparr + j*2)), double(expval)) ;
			totgrn += pow (double (*(temparr + j*2+1)), double(expval)) ;
			}
			else {
			totred += *(temparr + j*2) ;
			totgrn += *(temparr + j*2+1) ;
			}
		}
	}

	// get the avereage of the red and grn 
	totred /= count ;
	totgrn /= count ;
	scalered = scalefac * 150./totred ;
	scalegrn = scalefac * 150./totgrn ;

	// 
	cout << "Scaling parameters   Red :   "   << scalered << "    Green :   " << scalegrn << endl ; 

	ifil.seekg (long(startline)*nsamps*8, ios::beg) ;
	for (i=0; i<nlines; i++) {
		ifil.read ((char*)temparr, nsamps * 8) ;
		for (j=0; j<nsamps; j++)	
		{
			i_out = (vflip) ? nlines -i  -1 : i ;
			j_out = (hflip) ? nsamps -j  -1 : j ;
			optr0 = red + i_out * nsamps + j_out ;
			optr1 = grn + i_out * nsamps + j_out ;
			
			iptr0 = temparr + j * 2 ; 
			iptr1 = temparr + j * 2 + 1; 
			if (expflag) {
				b1 = pow (*iptr0, expval) ;
				b2 = pow (*iptr1, expval) ;
			}
			else {
				b1 = *iptr0 ;
				b2 = *iptr1 ;
			}
			redval = b1 * scalered ;
			redval = (redval < 0) ? 0. : redval ; 
			redval = (redval > 255) ? 255. : redval ; 
			*optr0 = (unsigned char) redval ;
			grnval = b2 * scalegrn ;
			grnval = (grnval < 0) ? 0. : grnval ; 
			grnval = (grnval > 255) ? 255. : grnval ; 
			*optr1 = (unsigned char) grnval ;
		}
	}
	delete [] temparr ;
	ifil.close() ;
	return (1) ;
}


int CGetData::DeleteRG () 
{
	if (red) delete [] red ;
	if (grn) delete [] grn ;
	red = NULL ;
	grn = NULL ;
	return (1) ;
}
int CGetData::DeleteMag () 
{
	if (mag) delete [] mag ;
	mag = NULL ;
	return (1) ;
}
int CGetData::DeleteRaw () 
{
	if (raw) delete [] raw ;
	raw = NULL ;
	return (1) ;
}
int CGetData::DeleteHgt () 
{
	if (hgt) delete [] hgt ;
	hgt = NULL ;
	if (mag) delete [] mag ;
	mag = NULL ;
	return (1) ;
}
	
int CGetData::DeleteMPH () 
{
	if (mag) delete [] mag ;
	if (phase) delete [] phase ;
	mag = NULL ;
	phase = NULL ;
	return (1) ;
}
	
			
	

	
	
	
	
	
