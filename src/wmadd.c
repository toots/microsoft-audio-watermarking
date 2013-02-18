/*

(c) Microsoft Corporation. All rights reserved. 

*/
//
//		WMADD.C	-	Adding watermark to wave files
//
//		(c) 1999 Microsoft Corp.
//
//		History:
//
//		14/Jun/99 - Darko Kirovski, first version, based on his Matlab code
//		16/Jun/99 - Rico Malvar, cleanup, remove Ymag, Xmag (speedup)
//		21/Jun/99 - Rico, add log10erfc() and support for multichannel files
//		24/Jun/99 - Rico, getblock and putblock function pointers
//		27/Jun/99 - Rico, support for SDMI 15-second window and insertion of
//					payload bits (CCI + Usage + ISRC)

#include "wmark.h"

// Long vectors, should be declared outside the routines
WSAMPLE	xwav[NCHMAX*NFREQ];		// I/O vector
float	xx[NCHMAX][NFREQ]; 		// time-domain samples
float	Xc[NCHMAX][NFREQ];		// MCLT coefficients, cosine part
float	Xs[NCHMAX][NFREQ];		// MCLT coefficients, sine part
float	h[NFREQ];				// MCLT window
float	ya[NCHMAX][3*NFREQ/2];	// MCLT internal buffers - analysis
float	ys[NCHMAX][3*NFREQ/2];	// MCLT internal buffers - synthesis
float	chip[NFREQ][FRAMESPERWINDOW];	// watermarking chips, one per MCLT component
								// (used for insertion of CCI btis)
SSBANDS	ssbands;				// pointers to subbands for SS bits
int		permutation[CHIPSPERBLOCK][BITSPERWINDOW];

// ------------------------------------------------------------
// Routine to generate SS bits
void getchips(int seed, int cci, int usage) {
	int		bit, i, j, k, freq, rows;
	int		ssbit, frame;
	int		toembed[BITSPERWINDOW];
	int		framesperbit;

	// how the watermark is created?
	// we have hardwired two lookup tables of rows for
	// chess watermarks of length 2 and four bits.
	// we create the larger 6 and 8 bit wide tables  
	// by crossconcatenating the smaller ones.

	// look up tables of size two and four blocks in a row
	int		lut2[2][2] = {   1, 0, 
							 0, 1};
	int		lut4[6][4] = {	1, 1, 0, 0, 
							1, 0, 1, 0, 
							1, 0, 0, 1, 
							0, 0, 1, 1, 
							0, 1, 0, 1, 
							0, 1, 1, 0}; 
	// look up tables of size 8 and 6 blocks in a row
	int		lut8[36][8];
	int		lut6[12][6];
	int		p_lut;

	// create the lookup tables for watermark creation
	// lookup tables of size 6 and 8 are made from 
	// lookup tables fo size 2 and four
	for (i = 0; i < 6; i++)
		for (j = 0; j < 6; j++) 
			for (k = 0; k < 8; k++) 
				if (k < 4) lut8[i*6+j][k] = lut4[i][k];
				else lut8[i*6+j][k] = lut4[j][k-4];
	for (i = 0; i < 6; i++)
		for (j = 0; j < 2; j++) 
			for (k = 0; k < 6; k++) 
				if (k < 4) lut6[i*2+j][k] = lut4[i][k];
				else lut6[i*2+j][k] = lut2[j][k-4];

	// Prepare seed for the random number generator; it
	// should be a function of the CCI bits and the input seed.
	// For now, we just add the 4-bit value of CCI to the seed.

	// Initialize RN generator with the highest bit of cci
	// there are 16 different watermarks 0 <= i <= 15.
	i = cci & 0xF;  
	INITIALIZE_SEED(seed, i);

	// encoding the usage and cci bits
	for (i = 0; i < BITSPERWINDOW; i++) {
		toembed[i] = (usage & (0x1 << i)) >> i;
	}

	// period equals the number of chip blocks per PAYLOADBIT frame
	framesperbit = (int) FRAMESPERWINDOW / BITSPERWINDOW;
	// depending on period the appropriate lookup table is selected
	// first we select its dimension (number of different rows)
	if (framesperbit == 6) rows = 12;
	else if (framesperbit == 8) rows = 36;
	else if (framesperbit == 4) rows = 6;
	// else exit(-1); // this case is an error
	for (bit = 0; bit < CHIPSPERBLOCK; bit++) {
		for (frame = 0; frame < FRAMESPERWINDOW; frame++) {
			// point to different row in the lookup table
			// after each period.
			if (frame % framesperbit == 0) {
				p_lut = (int) ((rand()*rows)/RAND_MAX);
			}
			// get the actual bit to be embedded at [bit,frame]
			if (framesperbit == 6) ssbit = lut6[p_lut][frame%6];
			else if (framesperbit == 8) ssbit = lut8[p_lut][frame%8];
			else if (framesperbit == 4) ssbit = lut4[p_lut][frame%4];
			// define chip from ssbit
			for (freq = ssbands.fstart[bit]; freq <= ssbands.fend[bit]; freq++) {
				// --- embedding the actual chip  or toembed = 0x0.
				// chip[freq][frame] = ssbit;
				// --- xoring with the communication channel (toembed) without permutations.
				// chip[freq][frame] = ssbit ^ toembed[(int) frame/period];
				// --- xoring + permutations of the communication channel
				// for each subband (bit from CHIPSPERBLOCK).
				chip[freq][frame] = ssbit ^ toembed[permutation[bit][(int) frame/framesperbit]];
			}
		}
	}
#ifdef PRINT_CHESS_CHIPS
    {
		int i,j;
		FILE *f;
		f = fopen("chess_chips.m", "w");
		fprintf(f, "cc=[\n");
		for (i = 0; i < NFREQ; i++) {
			for (j = 0; j < FRAMESPERWINDOW; j++) {
				fprintf(f, "%d ", (int)chip[i][j]);
			}
			fprintf(f, "\n");
		}
		fprintf(f, "];\n");
		fclose(f);
	}
#endif
}

// ------------------------------------------------------------
// Routine to check for pre-echo problems
int getwmflag(float *xx)
{
	int		i, j, scale;
	float	er, max_energy, min_energy, temp;

	// The case for not watermarking a block. If abrupt changes in energy are
	// detected in a block, watermark is not inserted to avoid pre-echoes

	scale = NFREQ / NSEC;
	max_energy = 0;
	min_energy = (float) 1e+20;
	for (i = 0; i < NSEC; i++) {
		temp = 0;
		for (j = 0; j < scale; j++) {
			temp += xx[i*scale + j] * xx[i*scale + j];
		}
		temp /= scale;
		if (temp > max_energy) max_energy = temp;
		if (temp < min_energy) min_energy = temp;
	}
	er = max_energy/(min_energy + ETHR);
	if (er < ERLIM) return(1); else return(0);
}

// creating the secret permutation
void create_permutation (int seed) {
	int bit, i, j, count, temp;
	srand(seed);
	for (bit = 0; bit < CHIPSPERBLOCK; bit++) {
		for (i = 0; i < BITSPERWINDOW; i++) 
			permutation[bit][i] = -1;
		for (i = 0; i < BITSPERWINDOW; i++) {
			temp = rand() % (BITSPERWINDOW - i);
			count = 0;
			for (j = 0; j < BITSPERWINDOW; j++) {
				if (permutation[bit][j] == -1) {
					if (count == temp) {
						permutation[bit][j] = i;
						break;
					} else count++;
				}
			}
		}
	}
#ifdef PRINT_PERMUTATIONS
    {
		int i,j;
		FILE *f;
		f = fopen("permutations.m", "w");
		fprintf(f, "p=[\n");
		for (i = 0; i < CHIPSPERBLOCK; i++) {
			for (j = 0; j < BITSPERWINDOW; j++) {
				fprintf(f, "%d ", (int)chip[i][j]);
			}
			fprintf(f, "\n");
		}
		fprintf(f, "];\n");
		fclose(f);
	}
#endif
}

// --------------------------------------------------------
// Routine to add watermak
void add_watermark(
			int  (* fn_getblock) (WSAMPLE *wp, int npts),
			void (* fn_putblock) (WSAMPLE *wp, int npts),
			WMBITS *wmbits, long nsamples, int nch, long fs) {

	int		i,j,freq,block,NBlocks,bsamples,ch;
	float	blocktlength;		// how many seconds per block
	float	frametlength;		// how many seconds per frame
	float	frametime;			// how many seconds within frame
	float	windowtime;			// how many seconds within the window
	long	window;				// window counter
	long	frame;				// frame counter
	float	wmgain;				// gain factor for each watermark
	int     freqmin; 			// index of min frequency to watermark
	int     freqmax; 			// index of max frequency to watermark
	int		nblknw[NCHMAX];		// how many blocks were not watermarked per channel
	int		nwindows;			// how many 15-sec windows
	float   mul, invmul;		// signal gain due to watermarking

	mul = DB(OFFSET);
	invmul = 1/mul;

	// Setting parameters
	NBlocks = (nsamples + NFREQ -1) / NFREQ;
	nwindows = (int) ceil(nsamples / (fs * WINDOWSIZE));
	freqmin = (int) ((FWMIN*2.0/fs) * (float) NFREQ);
	freqmax = (int) ((FWMAX*2.0/fs) * (float) NFREQ);

	// Initialize counters
	blocktlength = NFREQ / (float) fs;
	frametlength = (float) WINDOWSIZE / FRAMESPERWINDOW; 
	frametime = 0.0;
	windowtime = 0.0;
	window = 0;
	frame = 0;
	for (ch = 0; ch < nch; ch++) 
		nblknw[ch] = 0;

	//	Initialize MCLT window
	mlt_sine_window(h, NFREQ);

	// Initialize subband pointers; for insertion, no frequency shift
	gensubindex(&ssbands, fs, 1); 

	// Get first set of watermarking chips
	create_permutation(window);
	getchips(window, wmbits->xCCI[window], wmbits->xLOAD[window]);

	printf("\nEmbedding watermarks...\n");

	// Watermark insertion, block by block
#ifdef PRINT_REPEAT
	printf("Processing window  1 out of %2d. [Embedding=%X %X] \r", nwindows, wmbits->xCCI[window], wmbits->xLOAD[window]);
#else
    printf("T=%3.3fsec Processing window  1 out of %2d. [Embedding=%X%X] \n", WINDOWSIZE * window + windowtime, nwindows, wmbits->xCCI[window], wmbits->xLOAD[window]);
#endif

	for (block = 0; block < NBlocks; block++) {

		// Read one signal block (all channels)
		bsamples = fn_getblock(xwav, NFREQ * nch);

		// Copy sample data to signal vector xx and compute MCLT
		// Multiple channels are interlevd in xwav vector
		for (ch = 0; ch < nch; ch++) {
            
            // signal MCLT computation
			for (i = 0, j = ch; i < NFREQ; i++, j += nch) {
				xx[ch][i] = xwav[j] * INVMAXX;
			}
			fmclt(Xc[ch], Xs[ch], xx[ch], ya[ch], h, NFREQ);

			// Compute watermark flag; set to 1 if block is
			// "good enough " (no pre-echo) to be watermarked
			// Modify frequencies only up to FWMAX.
			// All channels are watermarked with the same chips
			if (getwmflag(xx[ch]) == 1) {
				for (freq = freqmin; freq <= freqmax; freq++) {
					if (chip[freq][frame] == 1) wmgain = mul;
					else wmgain = invmul;
					Xc[ch][freq] *= wmgain;
					Xs[ch][freq] *= wmgain;
				}
			} else nblknw[ch] += 1;

    		// Inverse MCLT and convert signal vector xx to sample data,
	    	// with channel interleaving
			fimclt(Xc[ch], Xs[ch], xx[ch], ys[ch], h, NFREQ);
			for (i = 0, j = ch; i < NFREQ; i++, j += nch) {
				xwav[j] = DBLTOWSAMPLE(xx[ch][i] * WSAMPLE_MAX);
			}
		}
	
		// Write the block of watermarked samples
		if (block == 0) {
			// do nothing - don't write first output block, to align output
		} else if (block == (NBlocks - 1)) {
			fn_putblock(xwav, bsamples);
			// write an additional blank output block
			for (i = 0; i < (NFREQ *nch); i++) xwav[i] = 0;
			fn_putblock(xwav, NFREQ * nch);
		} else {
			fn_putblock(xwav, bsamples);
		}

		// Update frame & window pointers
		frametime  += blocktlength;
		windowtime += blocktlength;

		// Check if current block should finish current frame
		if (frametime > frametlength) {

			// Yes, this was the last block of a frame
			frametime -= frametlength; // update counting for next frame
			frame++; // increment frame counter
		}

		// Check if current frame finishes a window
		if (windowtime > WINDOWSIZE) {
			windowtime -= (float) WINDOWSIZE;
			window++;	// increment window counter
			frame = 0;	// reset frame counter
#ifdef PRINT_REPEAT
            printf("Processing window %2d out of %2d. [Embedding=%X%X] \r", window+1, nwindows, wmbits->xCCI[window], wmbits->xLOAD[window]);
#else
			printf("T=%3.3fsec Processing window %2d out of %2d. [Embedding=%X%X] \n", WINDOWSIZE * window + windowtime, window+1, nwindows, wmbits->xCCI[window], wmbits->xLOAD[window]);
#endif
			// Get next set of spread-spectrum chips
			// and mext ramp pattern
			// getchips(window, wmbits->cci);
			create_permutation(window);
			getchips(window, wmbits->xCCI[window], wmbits->xLOAD[window]);
		}
	}
	printf("<< size of last window is %3.3f%% of watermark length\n", 100.*windowtime/(float)WINDOWSIZE);
	printf("\n");
	
	for (ch = 0; ch < nch; ch++) {
		if (nblknw[ch] > 0) {
			printf("In channel %2d, %4.1f%% of blocks were not watermarked.\n",
				ch+1, nblknw[ch] * 100.0 / NBlocks);
		}
	}
	
}


