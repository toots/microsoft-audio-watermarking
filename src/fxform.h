/*

(c) Microsoft Corporation. All rights reserved. 

*/
/*--------------------------------------------------------------------*
 *    FXFORM.H  -  Fast Sinusoidal Transforms                         *
 *                                                                    *
 *    (c) 1991 Henrique Malvar & (c) 1998 Microsoft Corp.             *
 *                                                                    *
 *    History:                                                        *
 *                                                                    *
 *     1/Aug/91 - Rico Malvar, first version.                         *
 *     5/Oct/98 - RM, first MS version.                               *
 *     6/Oct/98 - RM, added MCLT.                                     *
 *--------------------------------------------------------------------*/
 
void fft(float *u, int nfft);
void fdctiv(float *x, int n);
void fdstiv(float *x, int n);
void mlt_sine_window(float *h, int n);
void fmlt(float *x, float *ya, float * ha, int n);
void fimlt(float *x, float *ys, float * hs, int n);
void fmclt(float *Xc, float *Xs, float *x, float *ya, float *ha, int n);
void fimclt(float *Xc, float *Xs, float *x, float *ys, float *hs, int n);
