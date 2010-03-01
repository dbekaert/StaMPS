#include <iostream>
#include <fstream>
#include <math.h>

#ifndef CGD 
#define CGD

class CGetData {
	public :
	CGetData () ;
	~CGetData () ;
	int setparams (char *ifile, int ns, int startl, int nl) ;
	int getarrayRG (float scalefac, int expfl, float exval, int flipflag) ;
	int getarrayMag (float scalefac, int expfl, float exval, int flipflag, int fifthflag) ;
	int getarrayMag (int dtype, float minval, float maxval, int flipflag) ;
	int getarrayHgt (float scalefac, int expfl, float exval, float cont_interval, int flipflag) ;
	int getarrayMPH (float scalefac, int expfl, float exval, int flipflag) ;
	int getarrayByte (int min, int maxxval, int flipflag) ;
	int Raw2Mag (int dtype, float min, float max, int flipflag) ;

	int getarrayRaw (int dtype) ;
	int DeleteRG () ;
	int DeleteMag () ;
	int DeleteMPH () ;
	int DeleteHgt () ;
	int DeleteRaw () ;
	

	// file parameters
	char infile [240] ;
	int nlines, nsamps, totlines, startline ;
	static int bpp [] ;

	unsigned char *red ;
	unsigned char *grn ;
	unsigned char *mag ;
	unsigned char *phase ;
	unsigned char *hgt ;
	unsigned char *raw ;
	int expflag ;
	int dtype ;
	float expval ;
	
	float scalered, scalegrn ;
	inline double getmag (float re, float im) ;
	inline void getmph (float re, float im, double * mg, double *ph) ;
} ;

inline double CGetData::getmag (float re, float im) {
	return (sqrt(re * re + im * im)) ;
}

inline void CGetData::getmph (float re, float im, double *mag, double *ph) {
	*mag =  (sqrt(re * re + im * im)) ;
	*ph  =  atan2 (double(im), double (re)) ;
}


#endif 


