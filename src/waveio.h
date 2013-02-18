/*

(c) Microsoft Corporation. All rights reserved. 

*/
/*--------------------------------------------------------------------*
 *    WAVEIO.H  -  Definitions for module for accessing Microsoft     *
 *                 .wav files                                         *
 *                                                                    *
 *    (c) 1996 H. S. Malvar & 1999 Microsoft Corp.                    *
 *                                                                    *
 *    History:                                                        *
 *                                                                    *
 *    16/Jan/96 - Henrique Malvar - first version, limitations:       *
 *                + within RIFF wave subchunks, only 'data' and 'fmt' *
 *                  are considered, others are ignored.               *
 *                + only PCM (uncompressed) file types supported.     *   
 *		20/Jun/99 - Rico Malvar, added support for 24-bit and 32-bit    *
 *                files, change waveform sample type to long integer  *
 *--------------------------------------------------------------------*/
 
#include <limits.h>

typedef  unsigned char  BYTE;    /*  8-bit value, unsigned */
typedef  signed short   WORD;    /* 16-bit value, signed */
typedef  unsigned short UWORD;   /* 16-bit value, unsigned */
typedef  unsigned long  DWORD;   /* 32-bit value, unsigned */
typedef  signed long	WSAMPLE; /* 32-bit waveform sample value */

#define	WSAMPLE_MAX		LONG_MAX
#define	WSAMPLE_MIN		LONG_MIN

#define  FMT_PCM  		1

#define  W_READ   		0
#define  W_WRITE  		1
#define  W_RDWR   		2


/*--------------------------------------------------------------------*
 *    Brief info on RIFF WAVE files                                   *
 *                                                                    *
 *    The Resource Interchange File Format (RIFF) specifies that data *
 *    must be written to the file in chunks of data.  A chunk is a    *
 *    structure of the form                                           *
 *                                                                    *
 *    typedef struct{                                                 *
 *      char    chunk_ID[4];  <- four-character chunk ID              *
 *      ulong   csize;        <- size of chunk in bytes, unsigned     *
 *      uchar   data[csize];  <- the actual data                      *
 *    } chunk;                                                        *
 *                                                                    *
 *    Chunks with ID "RIFF" and "LIST" have an additional field:      *
 *                                                                    *
 *    typedef struct{                                                 *
 *      char    chunk_ID[4];  <- four-character chunk ID              *
 *      ulong   csize;        <- size of chunk in bytes, unsigned     *
 *      union {                                                       *
 *        char   form_ID[4];  <- four-character form identification   *
 *        uchar  data[csize]; <- the actual data                      *
 *      }                                                             *
 *    } chunk;                                                        *
 *                                                                    *
 *    The first chunk of a RIFF file must be a RIFF chunk. Within the *
 *    RIFF chunks, mandatory subchunks define file properties and     *
 *    carry signal data.  LIST chunks contain additional related      *
 *    information, such as copyright notices, text, icons, etc.       *
 *                                                                    *
 *    A wave file (.wav) contains one RIFF chunk with many subchunks, *
 *    in the form                                                     *
 *                                                                    *
 *    RIFF (WAVE) {                                                   *
 *       fmt    chunk;   <- mandatory, specifies sampling info        *
 *       fact   chunk;   <- optional chunk, e.g. for compressed types *
 *       other  chunks;  <- additional optional chunks                *
 *       data   chunk;   <- the data chunk containing waveform data   *
 *       fact   chunk;   <- optional chunk, e.g. for compressed types *
 *       DISP   chunk;   <- used for icons, for example               *
 *       LIST (INFO) chunk; <- optional, with many kinds of info      *
 *    }                                                               *
 *--------------------------------------------------------------------*/
 
/*--------------------------------------------------------------------*
 *    Format structure for the data portion of 'fmt ' subchunks       *
 *--------------------------------------------------------------------*/
 
typedef  struct {    /* format structure */
   short fmttag;     /* data format, should be = 1 for PCM */
							/* for other tags there may be additional fields */
   WORD  nchannels;  /* number of channels */
							/* 1 = mono, 2 = stereo (Left,Right) */
							/* 3 = (L,R,Center), 4 = (FrontL,FR,RL,RR) */
							/* 4 = (L,C,R,Surround), 6 = (LC,L,C,RC,R,S) */
   DWORD sampfreq;   /* sampling frequency in Hertz; most common
                        values are 11025, 22050, and 44100 */
   DWORD bytefreq;   /* sampling frequency * bytesps */
   WORD  bytesps;    /* bytes per set of nchannels samples */
   WORD  bitsps;     /* bits per sample - 8, 12, 16, 24, or 32 */
} wfmt;
 
 
/*--------------------------------------------------------------------*
 *    WAVEFILE pointer structure                                      *
 *--------------------------------------------------------------------*/
 
typedef  struct {    /* structure for wave file control */
   FILE     *fp;     /* pointer for corresponding DOS file */
   wfmt     *fmt;    /* pointer for corresponding format stucture */
   char     *name;   /* pointer for file name string */
   char     mode;    /* W_READ or W_WRITE  */
   DWORD    fsize;   /* format chunk size */
   DWORD    dpos;    /* file offset of data chunk */
   DWORD    dsize;   /* size of data chunk in bytes */
   DWORD    nsamp;   /* number of data samples per channel */
} WAVEFILE;
 
 
/*--------------------------------------------------------------------*
 *    wavefname - File name                                           *
 *                                                                    *
 *    Return a pointer to a string containing the wave file name      *
 *                                                                    *
 *    wp is a pointer to the opened wave file control structure.      *
 *--------------------------------------------------------------------*/
 
#define wavefname(wp) ((wp)->name)
 

/*--------------------------------------------------------------------*
 *    wavefs - sampling frequency                                     *
 *                                                                    *
 *    Usage: fs = wavefs(wp);                                         *
 *                                                                    *
 *    fs is a long that will contain the sampling frequency.          *
 *                                                                    *
 *    wp is a pointer to the opened wave file control structure.      *
 *--------------------------------------------------------------------*/
 
#define wavefs(wp) ((wp)->fmt->sampfreq)
 
 
/*--------------------------------------------------------------------*
 *    waveslength - length in bits of each sample                     *
 *                                                                    *
 *    Usage: length = waveslength(wp);                                *
 *                                                                    *
 *    length has the number of bits per sample, typically 8 to 32     *
 *                                                                    *
 *    wp is a pointer to the opened wave file control structure.      *
 *--------------------------------------------------------------------*/
 
#define waveslength(wp) ((wp)->fmt->bitsps)
 
 
/*--------------------------------------------------------------------*
 *    wavenchannels - number of channels, e.g. 1 if mono, 2 if stereo *
 *                    more channels are possible.                     *
 *                                                                    *
 *    Usage: st = wavenchannels(wp);                                  *
 *                                                                    *
 *    wp is a pointer to the opened wave file control structure.      *
 *--------------------------------------------------------------------*/
 
#define wavenchannels(wp) ((wp)->fmt->nchannels)
 
 
/*--------------------------------------------------------------------*
 *    wavensamples - no. of waveform samples in each channel          *
 *                                                                    *
 *    Usage: n = wavensamples(wp);                                    *
 *                                                                    *
 *    n has the number of waveform samples available in each channel  *
 *                                                                    *
 *    wp is a pointer to the opened wave file control structure.      *
 *--------------------------------------------------------------------*/
 
#define wavensamples(wp) ((wp)->nsamp)

 
/*--------------------------------------------------------------------*
 *    waveget8  - get  8-bit sample value from wave file              *
 *    waveget12 - get 12-bit sample value from wave file              *
 *    waveget16 - get 16-bit sample value from wave file              *
 *    waveget24 - get 24-bit sample value from wave file              *
 *    waveget32 - get 32-bit sample value from wave file              *
 *                                                                    *
 *    Usage: value = waveget8(wp);                                    *
 *           value = waveget12(wp);                                   *
 *           value = waveget16(wp);                                   *
 *           value = waveget24(wp);                                   *
 *           value = waveget32(wp);                                   *
 *                                                                    *
 *    wp is a pointer to the opened wave file control structure.      *
 *                                                                    *
 *    The returned value has type WSAMPLE.                            *
 *                                                                    *
 *    Note that 8-bit wave data is unsigned, others are signed, in    *
 *    little-endian order (LS byte first).                            *
 *                                                                    *
 *    These are implemented as macros that call getc.  For most       *
 *    compilers, getc is also a macro, and this leads to very fast    *
 *    sample reading, with full buffering.                            *
 *--------------------------------------------------------------------*/

#define  waveget8(wp) ((getc((wp)->fp) - 0x0080) << 24)

static   FILE     *ffftmp;
static	union {
				WSAMPLE	sample;
				short		sampint;
				BYTE		byte[4];
			} xxxtmp;

#define  waveget12(wp) ((ffftmp  = (wp)->fp, \
							    xxxtmp.sampint = 0, \
								 xxxtmp.byte[2] = getc(ffftmp), \
                         xxxtmp.byte[3] = getc(ffftmp)), \
                         (xxxtmp.sample))

#define  waveget16(wp) ((ffftmp  = (wp)->fp, \
							    xxxtmp.sampint = 0, \
								 xxxtmp.byte[2] = getc(ffftmp), \
                         xxxtmp.byte[3] = getc(ffftmp)), \
                         (xxxtmp.sample))

#define  waveget24(wp) ((ffftmp  = (wp)->fp, \
							    xxxtmp.sampint = 0, \
								 xxxtmp.byte[1] = getc(ffftmp), \
								 xxxtmp.byte[2] = getc(ffftmp), \
                         xxxtmp.byte[3] = getc(ffftmp)), \
                         (xxxtmp.sample))

#define  waveget32(wp) ((ffftmp  = (wp)->fp, \
								 xxxtmp.byte[0] = getc(ffftmp), \
								 xxxtmp.byte[1] = getc(ffftmp), \
								 xxxtmp.byte[2] = getc(ffftmp), \
                         xxxtmp.byte[3] = getc(ffftmp)), \
                         (xxxtmp.sample))


/*--------------------------------------------------------------------*
 *    waveput8  - write  8-bit sample value to wave file              *
 *    waveput12 - write 12-bit sample value to wave file              *
 *    waveput16 - write 16-bit sample value to wave file              *
 *    waveput24 - write 24-bit sample value to wave file              *
 *    waveput32 - write 32-bit sample value to wave file              *
 *                                                                    *
 *    Usage: waveput8(value, wp);                                     *
 *           waveput12(value, wp);                                    *
 *           waveput16(value, wp);                                    *
 *           waveput24(value, wp);                                    *
 *           waveput32(value, wp);                                    *
 *                                                                    *
 *    value: value to be written, of type WSAMPLE.                    *
 *    wp:    pointer to the opened wave file control structure.       *
 *                                                                    *
 *    For 8-bit wave files, value is written as unsigned bytes.       *
 *                                                                    *
 *    This is implemented as a macro that calls  putc.  For most      *
 *    compilers, putc is also a macro, and this leads to very fast    *
 *    sample writing, with full buffering.                            *
 *--------------------------------------------------------------------*/

#define  waveput8(v, wp)  {putc((((v) >> 24) + 0x80), (wp)->fp);}

#define  waveput12(v, wp)	{xxxtmp.sample = (v); ffftmp = (wp)->fp; \
									 xxxtmp.byte[2] &= 0xF0; \
               				 putc(xxxtmp.byte[2], ffftmp); \
									 putc(xxxtmp.byte[3], ffftmp);}

#define  waveput16(v, wp)	{xxxtmp.sample = (v); ffftmp = (wp)->fp; \
               				 putc(xxxtmp.byte[2], ffftmp); \
									 putc(xxxtmp.byte[3], ffftmp);}

#define  waveput24(v, wp)	{xxxtmp.sample = (v); ffftmp = (wp)->fp; \
               				 putc(xxxtmp.byte[1], ffftmp); \
               				 putc(xxxtmp.byte[2], ffftmp); \
									 putc(xxxtmp.byte[3], ffftmp);}

#define  waveput32(v, wp)	{xxxtmp.sample = (v); ffftmp = (wp)->fp; \
               				 putc(xxxtmp.byte[0], ffftmp); \
               				 putc(xxxtmp.byte[1], ffftmp); \
               				 putc(xxxtmp.byte[2], ffftmp); \
									 putc(xxxtmp.byte[3], ffftmp);}


/*--------------------------------------------------------------------*
 *    DBLTOWSAMPLE - macro to convert float values to waveform       *
 *    sample values, with rounding and clipping                       *
 *--------------------------------------------------------------------*/

#define DBLTOWSAMPLE(x) ((x) > (float) WSAMPLE_MAX ? WSAMPLE_MAX : \
                        ((x) < (float) WSAMPLE_MIN ? WSAMPLE_MIN : \
                        ((x) < 0 ? (WSAMPLE) ((x) - 0.5) : \
                         (WSAMPLE) ((x) + 0.5))))


/*--------------------------------------------------------------------*
 *    Function prototypes for functions in waveio.c                   *
 *--------------------------------------------------------------------*/

char     *wavestrerror(void);
void     waveperror(int);
WAVEFILE *wavefopen(char *, char *, size_t);
int      wavefclose(WAVEFILE *);
void     wavesetparms(WAVEFILE *, long, int, int);
int      wavegetblock(WSAMPLE *, int, WAVEFILE *);
int      waveputblock(WSAMPLE *, int, WAVEFILE *);
WSAMPLE	wavegetvalue(WAVEFILE *);
void 		waveputvalue(WSAMPLE, WAVEFILE *);
int      wavegotosample(WAVEFILE *, long);
