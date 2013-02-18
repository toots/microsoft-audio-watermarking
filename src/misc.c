/*

(c) Microsoft Corporation. All rights reserved. 

*/
//
//		MISC.C	-	Watermarking of audio - miscellaneous functions
//
//		(c) 1999 Microsoft Corp.
//
//		History:
//
//		14/Jun/99 - Darko Kirovski, first version, based on his Matlab code
//		29/Jun/99 - Rico Malvar, using analog warping function
//		30/Jun/99 - Rico Malvar, miscellaneous codes to read/write key/ID bits
//		02/Jul/99 - Darko & Rico - new width control for nonlinear subbands

#include "wmark.h"

// ------------------------------------------------------------
// Routine to warp the frequency scale, to support frequency scale changes

float warpfreq(float x, float a)
{
	return((exp(a*x) - 1) / (exp(a) - 1) );
}

// ------------------------------------------------------------
// Routine to generate subband limits

void gensubindex(SSBANDS *ssbands, long fs, float scalefactor) {
	int		bit;
	int		fc;
	float	freq, x, y, a;

	// First define bits for SS subbands, on which the chips carrying
	// the CCI bits will be inserted
	//
	// The band from FWCCIMIN to FWMAX is divided in CHIPSPERBLOCK subbands
	//
	//	 FWCCIMIN                                     FWMAX
	//	    |---------|---------|--- . . . ---|---------|
	//          :         :                       :
	//        bit 0     bit 1             bit CHIPSPERBLOCK-1
	//
	// The bands are not uniformly distributed. The parameter SPEEDUP controls
	// the subband with distortion: higgher subbands should be wider as we
	// increase speedup.
	// The main trick: figure out how much the last subband width needs to
	// be amplified; that's the parameter g

    // first determine a from max gradient of warp
	a = 1.7 * log(1 + SPEEDUP * CHIPSPERBLOCK) + 0.5 * (SPEEDUP * CHIPSPERBLOCK); 
	fc = (int) ((FWMIN * scalefactor * 2.0 / fs) * (float) NFREQ);
	bit = 0;
    // start at index corresponding to FWCCIMIN
	ssbands->fstart[0] = fc; 
	while (bit < CHIPSPERBLOCK) {
		x = (bit + 1.0) / CHIPSPERBLOCK;
		y = warpfreq(x, a);
        // analog frequency, in Hz
		freq = FWMIN + (FWMAX-FWMIN) * y; 
        // convert to index
		fc = (int) (freq * scalefactor * (2.0 / fs) * NFREQ); 
		ssbands->fend[bit] = fc;
		ssbands->fmiddle[bit] = (ssbands->fstart[bit] + fc) / 2;
		// ??? ssbands->width[bit] = ssbands->fend[bit] - ssbands->fstart[bit] + 1;
		bit++;
		if (bit == CHIPSPERBLOCK) break;
		ssbands->fstart[bit] = fc + 1;
	}
#ifdef PRINT_SUBBAND_INDEX
    {
		FILE *f;
		int i;
		f = fopen("subband_index.m", "w");
		fprintf(f, "sbi=[\n");
		for (bit = 0; bit < CHIPSPERBLOCK; bit++) {
			for (i = 0; i < SRFREQ; i++) 
				fprintf(f, "%d %d %d\t", ssbands->fstart[bit], ssbands->fmiddle[bit], ssbands->fend[bit]);
			fprintf(f, "\n");
		}
		fprintf(f, "];\n");
		fclose(f);
	}
#endif
}

// ------------------------------------------------------------
// Routines to convert seconds to MM:SS:MSMS string

static char tmpstr[200];

char *ftimestr(float secs)
{
	int		i;
	int		min, sec, ms;

	// Get minutes, seconds, and milliseconds

	min = (int) (secs/60);
	sec = (int) (secs-60*min);
	ms  = (int) (0.5 + (secs-60*min-sec) * 1000.0);

	for (i = 0; i < 31; i++) tmpstr[i] = 0;
	sprintf(tmpstr,"%d:%02d:%03d", min, sec, ms);

	return(tmpstr);
}






// ---------------------------------------------------------------
// Routines not used currently


/*
char *ftimestr2(float secs)
{
	int		i;
	int		min, sec, cs;

	// Get minutes, seconds, and milliseconds

	min = (int) (secs/60);
	sec = (int) (secs-60*min);
	cs  = (int) (0.5 + (secs-60*min-sec) * 100.0);

	for (i = 0; i < 31; i++) tmpstr[i] = 0;
	sprintf(tmpstr,"%d:%02d:%02d", min, sec, cs);

	return(tmpstr);
}


// ------------------------------------------------------------
// Routine to convert CCI to binary string

char *fbinstr(int cci)
{
	char	b3, b2, b1, b0;

	b3 = (cci & 0x08)? '1' : '0';
	b2 = (cci & 0x04)? '1' : '0';
	b1 = (cci & 0x02)? '1' : '0';
	b0 = (cci & 0x01)? '1' : '0';
	sprintf(tmpstr, "%c%c%c%c%c", b3, b2, b1, b0, '\0');

	return(tmpstr);
}

char *fbinstr2(int cc)
{
	char	b1, b0;

	b1 = (cc & 0x02)? '1' : '0';
	b0 = (cc & 0x01)? '1' : '0';
	sprintf(tmpstr, "%c%c%c", b1, b0, '\0');

	return(tmpstr);
}


// ------------------------------------------------------------
// Routine to convert ISRC to hex string

char *fl64str(LONG64 isrc)
{
	sprintf(tmpstr, "%07X%08X", 
		(DWORD) (isrc >> 32), (DWORD) (isrc & 0xFFFFFFFF));

	return(tmpstr);
}

// approximation of log10(erfc())

float	log10erfc(float x) {
	float 	t, u;
	u = fabs(x);
	t = 1.0 / (1.0 + 0.5 * u);
	u = log10(t) + 0.43429448190325 *
		 (-u * u - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * 
		 (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * 
		 (-1.13520398 + t * (1.48851587 + t * 
		 (-0.82215223 + t * 0.17087277)))))))));
	return u;
}

*/
/*
void gensubindex(SSBANDS *ssbands, PLBANDS *plbands, 
		long fs, int NFREQ, float scalefactor)
*/
	// Now define bits for payload subbands: "Usage" and "ISRC"
	//
	// The band from FWMIN to FWMAX is divided in NMODBITS subbands
	//
	//	  FWMIN                                       FWMAX
	//	    |---------|---------|--- . . . ---|---------|
	//          :         :                       :
	//        bit 0     bit 1               bit NMODBITS-1
	//
	// The bands are not uniformly distributed. The parameter SPEEDUP controls
	// the subband with distortion: higgher subbands should be wider as we
	// increase speedup.
	/*
	N = NMODBITS;
	// The main trick: figure out how much the last subband width needs to
	// be amplified; that's the parameter g
	g = 1 + SPEEDUP * N;
	a = 1.7 * log(g) + 0.5 * (g-1); // determine a from max gradient of warp
	fc = (int) ((FWMIN * scf * 2.0 / fs) * (float) NFREQ);
	bit = 0;
	plbands->fstart[0] = fc;
	while (bit < NMODBITS) {
		x = (bit + 1.0) / N;
		y = warpfreq(x, a);
		freq = FWMIN + (FWMAX-FWMIN) * y; // analog frequency, in Hz
		fc = (int) (freq * scf * (2.0 / fs) * NFREQ); // convert to index
		plbands->fend[bit] = fc;
		plbands->fmiddle[bit] = (plbands->fstart[bit] + fc) / 2;
		plbands->width[bit] = plbands->fend[bit] - plbands->fstart[bit] + 1;
		bit++;
		if (bit == NMODBITS) break;
		plbands->fstart[bit] = fc + 1;
	}
	*/
/*
	printf("FWCCIMIN = %d    FWMIN = %d    FWMAX = %d\n", FWCCIMIN, FWMIN, FWMAX);
	printf ("Total bits: %d  a = %g\n", NMODBITS, a);
	for (bit = 0; bit < NMODBITS; bit++) {
	 	printf("BIT %3d    FSTART = %4d  FMIDDLE = %4d  FEND = %4d  width = %4d\n", 
		bit, plbands->fstart[bit], plbands->fmiddle[bit], plbands->fend[bit],
		plbands->fend[bit]-plbands->fstart[bit]+1);
	}
	exit(0);
	*/
/*
	printf("FWMIN = %d  FWMAX = %d\n", FWMIN, FWMAX);
	printf ("Total bits: %d  a = %g\n", BITSPERBLOCK, a);
	for (bit = 0; bit < BITSPERBLOCK; bit++) {
	 	printf("BIT %3d    FSTART = %4d  FMIDDLE = %4d  FEND = %4d  width = %4d\n", 
		bit, ssbands->fstart[bit], ssbands->fmiddle[bit], ssbands->fend[bit],
		ssbands->width[bit]);
	}
	*/
