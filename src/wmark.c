/*

(c) Microsoft Corporation. All rights reserved. 

*/

//
//		WMARK.C	-	Adding watermark to wave files - main program
//
//		(c) 1999 Microsoft Corp.
//
//		History:
//
//		14/Jun/99 - Darko Kirovski, first version, based on his Matlab code
//		16/Jun/99 - Rico Malvar, cleanup
//		21/Jun/99 - Rico, add support for multichannel files
//		24/Jun/99 - Rico, support for DVD-RAW files
//		27/Jun/99 - Rico, support for SDMI 15-second window and insertion of
//						payload bits (CCI + Usg + ISRC)

#include <string.h>
#include "wmark.h"

// Global variables - needed for block R/W functions
FILE		*rawin;			// pointer to input raw file
FILE		*rawout;		// pointer to output raw file
WAVEFILE    *win;			// pointer to input .wav file
WAVEFILE    *wout;			// pointer to output .wav file

// Function to read a block of samples from a wave file
int getblock_wav(WSAMPLE *xwav, int npts)
{
	return (wavegetblock(xwav, npts, win));
}

// Function to write a block of samples to a wave file
void putblock_wav(WSAMPLE *xwav, int npts)
{
	waveputblock(xwav, npts, wout);
}

// --------------------------------
// Routine to read the payload bits
void read_wmbits(char *arg, WMBITS *wmbits)
{
	int				i, result;
	unsigned char	temp;

	wmbits->bitsloaded = (int) strlen(arg);
	for (i = 0; i < MAX_CHAR; i++) {
		wmbits->xCCI[i] =  0;
		wmbits->xLOAD[i] = 0;
	}
	for (i = 0; i < wmbits->bitsloaded; i++) {
		temp = arg[i];
		if (temp >= 48 && temp <= 57)
			result = temp - 48;
		else if (temp >= 65 && temp <= 70)
			result = temp - 55;
		else if (temp >= 97 && temp <= 102)
			result = temp - 87;
		if (i % 2 == 0) wmbits->xCCI[i>>1] =  result;
		else wmbits->xLOAD[i>>1] = result;
	}
		
	// Report values
	printf("WM bits: CCI  = ");
	for (i = 0; i < (wmbits->bitsloaded + 1) >> 1; i++)
		printf("%X", wmbits->xCCI[i]);
	printf("\nWM bits: LOAD = ");
	for (i = 0; i < wmbits->bitsloaded >> 1; i++)
		printf("%X", wmbits->xLOAD[i]);
	printf("\n");
}


// --------------------------------
// The main program
int  main(int argc, char *argv[])
{
	clock_t     start, finish;	// needed for time measurement
	float       duration;		// needed for time measurement
	long		nsamples;     	// no. of samples in .wav file
	int			ssize;			// sample size, in bits
	int			nchannels;		// no. of channels in .wav file
	long        Nblocks;		// number of blocks
	int			ok;				// parameter check flag
	long		fs;				// sampling frequency
	int			input_format;	// set to DVD_RAW or MS_WAV
	int			output_format;	// set to DVD_RAW or MS_WAV
	WMBITS		wmbits;			// watermarking bits
	int			nwindows;		// how many 15-sec windows

	// pointers to functions that read and write blocks of samples
	int  		(* fn_getblock) (WSAMPLE *wp, int npts);
	void 		(* fn_putblock) (WSAMPLE *wp, int npts);


	printf("(c) Microsoft Corporation. All rights reserved.\n");

	// Help message
	printf("\n"
		"Audio watermarking (c) 1999 Microsoft Corp.\n"
		"** CONFIDENTIAL AND PROPRIETARY **\n");

	if (argc != 4) {
		printf(
        "\n"
 	    "Usage:  watermark  input_file  output_file  HEX2embed \n"
		"\n"
        "Arguments:\n"
        "\n"
        "   input_file   : should be in .wav format\n"
        "   output_file  : will be in .wav format\n"
        "   HEX2embed    : hex string to be embedded in the audio clip\n"
		"                  (two digits per 11-second window)\n"
        "\n"
		"Example: watermark  test.wav  testw.wav  A087CD7\n\n"
		);
		exit(0);
	}

	// Read watermarking bits
	// now  this function reads the message from the command line.
	read_wmbits(argv[3], &wmbits);

	// Check input file type
	if (strstr(argv[1], "wav") || strstr(argv[1], "WAV"))
		input_format = MS_WAV;
	else
		input_format = WRONG_FORMAT;

	// Check output file type
	if (strstr(argv[2], "wav") || strstr(argv[2], "WAV"))
		output_format = MS_WAV;
	else
		output_format = WRONG_FORMAT;

	// Open input file
	if (input_format == MS_WAV) {
		win = wavefopen(argv[1], "r", BUFFSIZE);
   		if (win == NULL) waveperror(1);
		fn_getblock = getblock_wav;
		nchannels = wavenchannels(win);
		fs = wavefs(win);
		nsamples = wavensamples(win);
		ssize = waveslength(win);

	} else { // Wrong format
		fprintf(stderr, "Error: WrongFormat>%s\n", argv[1]);
      	error_exit;
	}

	// Compute number of blocks & define parameters
	if (nchannels > NCHMAX) {
		fprintf(stderr, 
			"\nSorry, this program handles at most %d channels.\n", NCHMAX);
		error_exit;
	}
	switch (fs) {
		// the only allowed sampling frequencies
		case 44100:
			ok = 1; break;
		default: 
			ok = 0;
	}
	if (!ok) {
		fprintf(stderr,
			"\nSorry, a sampling frequency of %.1f kHz is not supported.\n", fs/1e3);
		error_exit;
	}

	Nblocks = (nsamples + NFREQ - 1) / NFREQ;
	nwindows = (int) (nsamples / (fs * WINDOWSIZE));

	// Open output file & define .wav header parameters
	if (output_format == MS_WAV) {
		wout = wavefopen(argv[2], "w", BUFFSIZE);
		if (wout == NULL) waveperror(1);
		wavesetparms(wout, fs, ssize, nchannels);
		fn_putblock = putblock_wav;

	} else { // Wrong Format
		fprintf(stderr, "Error: WrongFormat>%s\n", argv[2]);
		error_exit;
	}

	printf(
		"Input  file: %s\n"
		"Output file: %s\n"
		"Number of channels = %d\n"
		"Sampling frequency = %d Hz\n"
		"Sound clip length = %s (M:S:D) = %d samples\n"
		"Number of %.0f-second windows = %d \n"
		"Block length = %d samples = %.1f ms\n"        
		"Number of blocks = %d\n", 
		argv[1], argv[2], nchannels, fs, 
		ftimestr(nsamples / (float) fs), nsamples, 
		WINDOWSIZE, nwindows, NFREQ, NFREQ*1000.0/fs,	Nblocks);

#ifdef PRINT_SUBBAND_INDEX
    printf("Subband indices located in <subband_index.m> as variable <sbi>.\n");
#endif
#ifdef PRINT_CHESS_CHIPS
    printf("Chess chips located in <chess_chips.m> as variable <cc>.\n");
#endif
#ifdef PRINT_PERMUTATIONS
    printf("Permutations located in <permutations.m> as variable <p>.\n");
#endif

	// Reset processing time counter
	start = clock();

	// Embed watermark
	add_watermark(fn_getblock, fn_putblock, &wmbits, nsamples, nchannels, fs);

	// Measure elapsed time
	finish = clock();
	duration = (float) (finish - start) / CLOCKS_PER_SEC;
	printf("Elapsed time: %s (M:S:D) = %3.0f%% of real time.\n",
		ftimestr(duration), 100.0 * (duration * (float) fs/nsamples));

	// Close files and bye!

	if (input_format == MS_WAV)
		if (wavefclose(win)) waveperror(0);
	if (output_format == MS_WAV)
		if (wavefclose(wout)) waveperror(0);
	fcloseall();
}
