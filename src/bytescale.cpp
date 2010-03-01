#include <iostream>
#include <fstream>


template <class k>
void tbytescale (k *indat, unsigned char *outdat, int ns, int nl, float scale, float offset, int flipflag) 
{

	unsigned char *outptr ;
	int 	i, j, hflip=0, vflip=0, i_out, j_out ;
	k	*datptr ;
	float	fval ;

	if (flipflag %2) hflip=1 ;
	if (flipflag /2) vflip=1 ;

	for (i=0; i<nl; i++) 
	{
		i_out = (vflip)? nl -i - 1 : i ;
	for (j=0; j<ns; j++) 
	{
		j_out = (hflip)? ns -j - 1 : j ;
		datptr = indat + i * ns + j ;
		outptr = outdat +i_out * ns + j_out ;
		fval = (float) *datptr * scale + offset ;
		fval = (fval>255.) ? 255. : fval ;
		fval = (fval<0.) ? 0. : fval ;
		*outptr = (unsigned char) fval ;
	}
	}


	return ;
}

void bytescale (unsigned char *indat, unsigned char *outdat, int ns, int nl, float
	scale, float offset, int dtype, int flipflag) 
{
	if (dtype==0)
		tbytescale (&indat[0], outdat, ns, nl, scale, offset, flipflag) ;
	if (dtype==1)
		tbytescale ((short *) &indat[0], outdat, ns, nl, scale, offset, flipflag) ;
	if (dtype==2)
		tbytescale ((int *) &indat[0], outdat, ns, nl, scale, offset, flipflag) ;
	if (dtype==3)
		tbytescale ((float *)&indat[0], outdat, ns, nl, scale, offset, flipflag) ;
	return ;
}
