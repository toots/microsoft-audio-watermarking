/*

(c) Microsoft Corporation. All rights reserved. 

*/
//
//		HTHRES.H	-	Computation of hearing thresholds
//
//		(c) 1999 Microsoft Corp.
//
//		History:
//
//		10/Mar/98 - Rico Malvar, Matlab version
//		14/Jun/99 - Wenqing Jiang, C version

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define	MAXN			2048 		// max length of MCLT block
#define	C_EPS       	1e-16
#define	C_NCBANDS   	25			// No. of Bark subbands
#define	Dabs        	70.
#define	Rfac		  	10.			// dB masking @ same frequency
#define	Thmin		  	80.

float thrabs( float f);
void hthres(float *Ht, float *X, int Nbands, long fs);
