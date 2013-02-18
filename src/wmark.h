/*

(c) Microsoft Corporation. All rights reserved. 

*/
//
//		WMADD.H	-	Adding watermark to wave files
//
//		(c) 1999 Microsoft Corp.
//
//		History:
//
//		14/Jun/99 - Darko Kirovski, frist version, based on Matlab code
//		21/Jun/99 - Rico, add log10erfc() and support for multichannel files
//		24/Jun/99 - Rico, support for DVD-RAW files
//		27/Jun/99 - Rico, support for SDMI 15-second window and insertion of
//						payload bits (CCI + Usg + ISRC)
//		29/Jun/99 - Darko, new pointer structures
//		02/Jul/99 - Darko & Rico - new width control for nonlinear subbands

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <io.h>
#include <limits.h>
#include "waveio.h"
#include "fxform.h"
//#include "fdct.h"

// ignore double to float conversion warning
#pragma warning( disable : 4996 4244 4305)

// x(t) - mclt block of time domain samples
// x(w) - mclt block of frequency domain samples

//-------------WMADD DEFINITIONS----------------
//#define PRINT_SUBBAND_INDEX
//#define PRINT_CHESS_CHIPS
//#define PRINT_PERMUTATIONS
//#define PRINT_REPEAT
//----------------------------------------------

//-------------WMTEST DEFINITIONS---------------
//#define PRINT_TIME_POINTERS
#define TO_CEPSTRUM_WITH_FDISTV
#define SPECTRUM_NORMALIZATION
//#define SMALL_MCLT_WINDOW
//----------------------------------------------


// Constants
#define	MS_WAV				31
#define	WRONG_FORMAT		32
#define	NONE				(-1)		// value for no detection
#define	MAX_CHAR			1024		// max Hex numbers to be embedded in 
										// an audio clip
#define BUFFSIZE			327680   	// file I/O buffer
#define	INVMAXX				(1.0/((float) LONG_MAX)) 
										// scale factor for .wav samples
// General watermarking parameters

#define	NCHMAX				2			// max # of input channels
#ifdef SMALL_MCLT_WINDOW
	#define	NFREQ			1024		// max |x(w)| 
	#define LOG2_NFREQ		10			// log2 size of |x(w)| 
	#define CHIPSPERBLOCK	120			// chips added in a block [60 to 100]
	#define CF				5			// lowpass filter on the cepstrum
	#define DECISIONBAR		7.75   		// if (nc>) then watermark detected
	#define PREDECISIONBAR	2.25		// if (nc>) then watermark detected
	#define PM				20			// cepstrum amplitude clip 
#else
	#define	NFREQ			2048		// max |x(w)| 
	#define LOG2_NFREQ		11			// log2 size of |x(w)| 
	#define CHIPSPERBLOCK	120			// chips added in a block [60 to 100]
	#define CF				10			// lowpass filter on the cepstrum
	#define DECISIONBAR		7.75   		// if (nc>) then watermark detected
	#define PREDECISIONBAR	0.50		// if (nc>) then watermark detected
	#define PM				3			// cepstrum amplitude clip 
#endif
#define FWMIN				2000		// min freq of a chip subband in Hz
#define	FWMAX				7200		// max freq of a chip subband in Hz
#define	WINDOWSIZE			11.14558	// window size in seconds
#define	BITSPERWINDOW		4			// bits exored in a window
#define	FRAMESPERWINDOW		24			// window size in frames 
// COMMENT: number of blocks per window proportional to |x(w)|
#define	MAXBPF				32			// max number of blocks in a frame:
// COMMENT: in this scheme 10 for 4096 and 20 for 2048 mclt timeblocks
#define NWATERMARKS			16			// #CCIbits = log2(NWATERMARKS)
#define	OFFSET				4			// chip magnitude change in dB
#define SPEEDUP				0.01		// resilience to dynamic freq stretch

// Pre-echo control parameters

#define NSEC				8			// subblocks per block
#define ETHR				1E-8		// to stabilize energy ratio 
#define	ERLIM				150.		// max energy ratio of two subbblocks

// Watermark detection parameters

#define DBCUT				-4.0		// how many dBs to cut
#define SRTIME				7			// how many time scales to search 
#define SRFREQ				7			// how many freq scales to search
#define TIME_RESILIENCE		10.			// resilience to time scale in %
#define FREQ_RESILIENCE		5.			// resilience to freq shift in %
#define FDNOISEFLOOR		1e-4		// min subband energy to correlate

// Search steps

#define BASIC_STEP			0.3			// basic search step in frames 
// COMMENT: load_buffer() is written such that integrationarea=2*BASIC_STEP 
// #define SEARCHSTEP			0.13931973	// search step in time
#define	iSTUBBORNESS		WINDOWSIZE/2
#define	sSTUBBORNESS		WINDOWSIZE*3/4
#define	STUBBORNESS			WINDOWSIZE/2
#define xSTUBBORNESS		WINDOWSIZE*3/8
#define sPROGRESSSTEP		WINDOWSIZE*3/4
#define uPROGRESSSTEP		WINDOWSIZE*2/4

// General macros

#define	DB(x)			pow(10.0,(x)/20.0)
#define error_exit		exit(1)  /* error return with DOS exit code */
/*
#define DBLTOSHRT(x)	((x) > (float) SHRT_MAX ? SHRT_MAX : \
                      	((x) < (float) SHRT_MIN ? SHRT_MIN : \
                      	((x) < 0 ? (int) ((x) - 0.5) : \
                      	(int) ((x) + 0.5))))
*/

// Set the state of internal RN generator
// A better set of seeds should be defined later

#define	INITIALIZE_SEED(seed, i)	srand((71912317*(i*(i+2)*(seed+1))+(i+4779)*317*(seed))%15991)

// 64-bit integers
// typedef	__int64		LONG64;

// Payload structure

typedef struct {
	int				bitsloaded;
	unsigned int	xCCI[MAX_CHAR];
	unsigned int	xLOAD[MAX_CHAR];
} WMBITS;

// Structure for SS bits insertion

typedef struct {
   int fstart[CHIPSPERBLOCK];	// subband lower limit for SS bit
   int fend[CHIPSPERBLOCK];		// subband upper limit for SS bit
   int fmiddle[CHIPSPERBLOCK];	// pointers to center of subbands
} SSBANDS;

// Structure for SS bits detection

typedef struct {
   int fmiddle[CHIPSPERBLOCK][SRFREQ];	// pointers to center of subbands
   int cbe[CHIPSPERBLOCK][SRFREQ];
   int cbs[CHIPSPERBLOCK][SRFREQ];
} SSDBANDS;

// Structures for detection

typedef struct {				// circular buffer for MCLT data
	float		*buffer;
	short int	*ht;
	long		length;
	long 		pointer;
} CIRCULAR_BUFFER;

// Structure that contains the correlation values

typedef struct {
	int		card[BITSPERWINDOW][2];
	float   corr[BITSPERWINDOW][2];
	// float	std[PAYLOADBITS][2];
	float 	sumsquares;
} PAYLOAD;

typedef struct {
	float	dtime;
	float 	nc;
	int		cci;
	int		load;
	float 	payload[BITSPERWINDOW];
} RESULT;

// Function prototypes

void add_watermark(int (* fn_getblock) (WSAMPLE *wp, int npts), 
				   void (* fn_putblock) (WSAMPLE *wp, int npts),
				   WMBITS *wmbits, 
				   long nsamples, 
				   int nch, 
				   long fs);

void gensubindex (SSBANDS *ssbands, long fs, float scalefactor);

float log10erfc(float x);

char *ftimestr(float secs);

// char *ftimestr2(float secs);
// char *fbinstr(int cci);
// char *fbinstr2(int cc);
// char *fl64str(LONG64 isrc);
