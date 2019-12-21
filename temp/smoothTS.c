//
//  main.c
//  PhenoRice
//
//  Created by Fabio Pavesi on 03/04/13.
//  Copyright (c) 2013 Fabio Pavesi. All rights reserved.
//

#include <R.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <memory.h>
#include <time.h>
#include <sys/time.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include "smoothTS.h"

#define TAPPO -99999.9
#define VECTOR_SIZE 1000

#define CLOUDY 		1e4
#define CLEAR		1e-4
#define A_BIT_CLOUDY	1e-2

// 2nd degree poly coefficients
#define C(i) (gsl_vector_get(c,(i)))
// covariance matrix elements
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

// standard value for unweighted smoothing
#define NO_WEIGHT 1.0

// all of the following need be odd numbers
#define SMOOTHING_WINDOW_SIZE 7
#define DESPIKE_WINDOW_SIZE 5
#define DESPIKE_SMOOTHING_WINDOW_SIZE 7

void useC(int *i) {
    i[0] = 11;
}

static double *adattaPesi(double *in, int quanti) {
    double *out = malloc(77 * sizeof(double));
    
    for ( int i = 0; i < quanti; i++ ) {
        switch( (int)in[i] ) {
            case 3:
                out[i] = CLEAR;
                break;
            case 1:
                out[i] = CLOUDY;
                break;
            case 2:
                out[i] = A_BIT_CLOUDY;
        }
    }
    return out;
}

static void dumpVector(double vector[], int length) {
    for ( int i = 0; i < length; i++ ) {
        printf("%d. %g\n", i, vector[i]);
    }
}

/* Return 1 if the difference is negative, otherwise 0.  */
static int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1)
{
    long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);
    result->tv_sec = diff / 1000000;
    result->tv_usec = diff % 1000000;
    
    return (diff<0);
}

static void timeval_print(struct timeval *tv)
{
    char buffer[30];
    time_t curtime;
    
    printf("%ld.%06ld", tv->tv_sec, tv->tv_usec);
    curtime = tv->tv_sec;
    strftime(buffer, 30, "%m-%d-%Y  %H:%M:%S", localtime(&curtime));
    printf(" = %s.%06ld\n", buffer, tv->tv_usec);
}

static double smooth(double *yi, double *ei, int startIndex, int numberOfSamples) {
    gsl_matrix *X, *cov;
    gsl_vector *y, *w, *c;
    double xn, yn, en;
    double chisq;
    double result;
    
    X = gsl_matrix_alloc (numberOfSamples, 3);
    y = gsl_vector_alloc (numberOfSamples);
    w = gsl_vector_alloc (numberOfSamples);
    
    c = gsl_vector_alloc (3);
    cov = gsl_matrix_alloc (3, 3);
    
    for ( int i = 0; i < numberOfSamples; i++ ) {
        // xn = xi[i+startIndex]; //*(xi + ( i * sizeof(double) + startIndex ));
        xn = -(numberOfSamples / 2) + i;
        yn = yi[i+startIndex]; //i*(yi + ( i * sizeof(double) + startIndex ));
        en = ei[i+startIndex]; // *(ei + ( i * sizeof(double) + startIndex ));
        gsl_matrix_set (X, i, 0, 1.0);
        gsl_matrix_set (X, i, 1, xn);
        gsl_matrix_set (X, i, 2, xn * xn);
        
        gsl_vector_set (y, i, yn);
        // gsl_vector_set (w, i, 1.0/(en * en));
        gsl_vector_set (w, i, en);
        
    }
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc (numberOfSamples, 3);
    gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);
    gsl_multifit_linear_free (work);

    result = C(0);

    gsl_vector_free(c);
    gsl_vector_free(y);
    gsl_vector_free(w);
    gsl_matrix_free(X);
    gsl_matrix_free(cov);
     
    return result;
}

static gsl_vector *smooth_unweighted(double *xi, double *yi, int startIndex, int numberOfSamples) {
    gsl_matrix *X, *cov;
    gsl_vector *y, *w, *c;
    double xn, yn, en;
    double chisq;
    
    X = gsl_matrix_alloc (numberOfSamples, 3);
    y = gsl_vector_alloc (numberOfSamples);
    w = gsl_vector_alloc (numberOfSamples);
    
    c = gsl_vector_alloc (3);
    cov = gsl_matrix_alloc (3, 3);
    
    for ( int i = 0; i < numberOfSamples; i++ ) {
        // xn = xi[i+startIndex]; //*(xi + ( i * sizeof(double) + startIndex ));
        xn = xi[i+startIndex];
        yn = yi[i+startIndex]; //i*(yi + ( i * sizeof(double) + startIndex ));
        
        en = NO_WEIGHT;
        
        gsl_matrix_set (X, i, 0, 1.0);
        gsl_matrix_set (X, i, 1, xn);
        gsl_matrix_set (X, i, 2, xn * xn);
        
        gsl_vector_set (y, i, yn);
        gsl_vector_set (w, i, 1.0/(en * en));
        
    }
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc (numberOfSamples, 3);
    gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);
    gsl_multifit_linear_free (work);
    
    
     gsl_vector_free(y);
     gsl_vector_free(w);
     gsl_matrix_free(X);
     gsl_matrix_free(cov);
     
    return c;
}

static double despike(double *yi, double stdDevMultiplier, int startIndex, int numberOfSamples, int numberOfSmoothingSamples) {
    gsl_vector *c = NULL;
    double x[VECTOR_SIZE];
    double smoothingSample[VECTOR_SIZE];
    
    double sample[VECTOR_SIZE];
    double out[VECTOR_SIZE];
    
    
    int middleIndex = startIndex + (numberOfSmoothingSamples) / 2;
    int sampleStart = startIndex + (numberOfSmoothingSamples - numberOfSamples) / 2;
    int sampleEnd = startIndex + numberOfSamples;
    
    int smoothingSampleStart = startIndex;
    int smoothingSampleEnd = startIndex + numberOfSmoothingSamples;
    
    //    printf("Smoothing sample start: %d\n", smoothingSampleStart);
    //    printf("Sample start: %d, end: %d\n", sampleStart, sampleEnd);
    //    printf("Middle index: %d -> %g\n", middleIndex, yi[middleIndex]);
    
    for ( int i = smoothingSampleStart; i < smoothingSampleEnd; i++ ) {
        if ( i >= sampleStart && i <= sampleEnd ) {
            sample[i - sampleStart] = yi[i];
            //            printf("sample[%d] = %g\n", i - sampleStart, yi[i]);
        } else {
            sample[i - sampleStart] = -999999.0;
        }
        smoothingSample[i - smoothingSampleStart] = yi[i];
        x[i - smoothingSampleStart] = i - middleIndex;
    }
    //    dumpVector(sample, numberOfSamples);
    //    dumpVector(smoothingSample, numberOfSmoothingSamples);
    for ( int i = middleIndex - startIndex; i < numberOfSmoothingSamples - 1; i++ ) {
        smoothingSample[i] = smoothingSample[i+1];
        //        printf("smoothing sample[%d] = %g\n", i , smoothingSample[i+1]);
        x[i] = x[i + 1];
    }
    for ( int i = (numberOfSamples) / 2; i < numberOfSamples - 1; i++ ) {
        sample[i] = sample[i+1];
        //        printf("sample[%d] = %g\n", i, sample[i+1]);
    }
    //    dumpVector(sample, numberOfSamples);
    //    dumpVector(smoothingSample, numberOfSmoothingSamples);
    /*
     dumpVector(x, numberOfSmoothingSamples - 1);
     dumpVector(smoothingSample, numberOfSmoothingSamples - 1);
     */
    double mean = gsl_stats_mean(sample, 1, numberOfSamples - 1);
    double stdDeviation = gsl_stats_sd(sample, 1, numberOfSamples - 1);
    
    //    printf("Valore centrale: %g, media: %g, deviazione standard: %g\n", yi[middleIndex], mean, stdDeviation);
    
    if ( fabs( yi[middleIndex] - mean ) > stdDevMultiplier * stdDeviation ) {
        //        printf("spike: \n");
        c = smooth_unweighted(x, smoothingSample, 0, numberOfSmoothingSamples - 1);
        // out[middleIndex] = C(0);
        //        printf("Sostituisco il valore centrale con: ", C(0));
	double result = C(0);
	
	gsl_vector_free(c);

        return result;
    }
    
    return yi[middleIndex];
}

static double *removeSmallDrops(double *y, int length) {
    double *res = (double *)malloc(length * sizeof(double));
    for ( int i = 1; i < length - 1; i++ ) {
        if ( y[i] < y[i+1] && y[i] < y[i-1] ) {
            res[i] = ( y[i+1] + y[i-1] ) / 2.0;
        } else {
            res[i] = y[i];
        }
    }
    return res;
}

void smoothTS(double *y, double *w) {
    static int runCounter = 0;
    double *wAdattato;
    double despiked[77];
    double *deSmallDropped1;
    double *deSmallDropped2;
    double smoothed[77];
    
    // Rprintf("Smoother v0.1\n");
    gsl_vector *c;
    int numSamples = 7;
    
    struct timeval tvBegin, tvEnd, tvDiff;
    gettimeofday(&tvBegin, NULL);
    //    timeval_print(&tvBegin);
    
    // wAdattato = adattaPesi(w, 77);
    wAdattato = w;

    // Rprintf("1");
    double res;
    for ( int i = 0; i < 3; i++ ) {
        despiked[i] = y[i];
        smoothed[i] = y[i];
    }
    for ( int i = 0; i < 77 - 3; i++ ) {
        // c = smooth(y, wAdattato, i, numSamples);
        despiked[i+3] = despike(y, 2.0, i, 5, 7);
    }
    
    // Rprintf("2");
    deSmallDropped1 = removeSmallDrops(despiked, 77);
    for ( int i = 0; i < 3; i++ ) {
        deSmallDropped1[i] = y[i];
	deSmallDropped1[76 - i] = y[76 - i];
    }
    // Rprintf("3");
    // deSmallDropped1 = removeSmallDrops(despiked, 77);
    // Rprintf("4");
    deSmallDropped2 = removeSmallDrops(deSmallDropped1, 77);
    for ( int i = 0; i < 3; i++ ) {
        deSmallDropped2[i] = y[i];
	deSmallDropped2[76 - i] = y[76 - i];
    }
    
    for ( int i = 0; i < 3; i++ ) {
        smoothed[76 - i] = y[76 - i];
    }
    // Rprintf("5");
    for ( int i = 0; i < 77 - 3; i++ ) {
        y[i+3] = smoothed[i+3] = smooth(deSmallDropped2, wAdattato, i, 7);
    }
    
    free(deSmallDropped1);
    free(deSmallDropped2);

    gettimeofday(&tvEnd, NULL);
    timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
    // Rprintf("\n\n%ld.%06d\n", tvDiff.tv_sec, tvDiff.tv_usec);
    
    // Rprintf("6");
    /*
    char filename[255];
    sprintf(filename, "elaborazone_%.07d.csv", runCounter++);
    FILE *f = fopen(filename, "w");
    fprintf(f, "original, despiked, smalldrop1, smalldrop2, smoothed\n");
    for ( int i = 0; i < 77; i++ ) {
        fprintf(f, "%g, %g, %g, %g, %g, %g\n", y[i], w[i], despiked[i], deSmallDropped1[i], deSmallDropped2[i], smoothed[i]);
    }
    fclose(f);
     */
    // return 0;
}


#define BUFSIZE 4096

int main(int argc, char **argv) {
    struct timeval tvBegin, tvEnd, tvDiff;
    double dato[77];
    double peso[77];
    enum states { IDLE, NEW_RECORD, WAIT_FIELD1, READING_FIELD1, WAIT_FIELD2, READING_FIELD2 };
    enum states state;
    char buffer[BUFSIZE];
    int bufferCount = 0;
    FILE *f, *timing;
    // f = fopen("ProfilesEviNoise_res300_2010h26v06_gm.txt", "r");
    f = stdin;
    char c;
    int row = 0;
    timing = fopen("timing.log", "w");
    
    state = WAIT_FIELD1;
    while ( (c = fgetc(f)) != EOF ) {
        switch (c) {
            case ' ':
            case '\t':
            case '\n':
            case '\r':
                if ( state == READING_FIELD2 ) {
                    buffer[bufferCount++] = '\0';
                    sscanf(buffer, "%lg", &peso[row++]);
                    bufferCount = 0;
                    state = WAIT_FIELD1;
                    
                }
                if ( state == READING_FIELD1 ) {
                    buffer[bufferCount++] = '\0';
                    sscanf(buffer, "%lg", &dato[row]);
                    bufferCount = 0;
                    state = WAIT_FIELD2;
                }
                break;
            default:
                if ( state == WAIT_FIELD1 ) {
                    state = READING_FIELD1;
                }
                if ( state == WAIT_FIELD2 ) {
                    state = READING_FIELD2;
                }
                buffer[bufferCount++] = c;
                break;
        }
        if ( state == WAIT_FIELD1 && row == 77 ) {
            gettimeofday(&tvBegin, NULL);
            //            timeval_print(&tvBegin);
            
            smoothTS(dato, peso);
            
            gettimeofday(&tvEnd, NULL);
            timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
            // fprintf(timing, "%ld.%06d\n", tvDiff.tv_sec, tvDiff.tv_usec);
            
            row = 0;
        }
    }
    fclose(timing);
}


