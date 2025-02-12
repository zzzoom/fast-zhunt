#include "delta_linking.h"

#include <math.h>
#include <stdlib.h>

static double *bztwist; // Read-only after init
static const double _k_rt = -0.2521201; /* -1100/4363 */
static const double sigma = 16.94800353; /* 10/RT */
static const double explimit = -600.0;

void delta_linking_init(int max_dinucleotides)
{
    static double a = 0.357, b = 0.4; /* a = 2 * (1/10.5 + 1/12) */
    double ab;

    bztwist = (double*)calloc(max_dinucleotides, sizeof(double));
    ab = b + b;
    for (int i = 0; i < max_dinucleotides; i++) {
        ab += a;
        bztwist[i] = ab;
    }
}

void delta_linking_destroy(void)
{
    free(bztwist);
}

static void delta_linking_exponent(double dl, int terms, const double* logcoef, double* expmini_out, double* exponent) {
    double expmini = 0.0;
    #pragma omp simd reduction(min:expmini)
    for (int i = 0; i < terms; i++) {
        double z = dl - bztwist[i];
        exponent[i] = z = logcoef[i] + _k_rt * z * z;
        if (z < expmini) {
            expmini = z;
        }
    }
    *expmini_out = (expmini < explimit) ? explimit - expmini : 0.0;
}

static double delta_linking(double dl, double deltatwist, const double* logcoef, int terms)
{
    double expmini;
    double exponent[terms];
    delta_linking_exponent(dl, terms, logcoef, &expmini, exponent);

    double sump = 0.0;
    double sumq = 0.0;
    #pragma omp simd reduction(+:sumq,sump)
    for (int i = 0; i < terms; i++) {
        double z = exp(exponent[i] + expmini);
        sumq += z;
        sump += bztwist[i] * z;
    }
    sumq += exp(_k_rt * dl * dl + sigma + expmini);
    return deltatwist - sump / sumq;
}

static double linear_search_dl(double x1, double x2, double tole, double deltatwist, const double* logcoef, int terms)
{
    double f = delta_linking(x1, deltatwist, logcoef, terms);
    double fmid = delta_linking(x2, deltatwist, logcoef, terms);
    if (f * fmid >= 0.0) {
        return x2;
    }

    double dx, xmid;
    double x = (f < 0.0) ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
    do {
        dx *= 0.5;
        xmid = x + dx;
        fmid = delta_linking(xmid, deltatwist, logcoef, terms);
        if (fmid <= 0.0) {
            x = xmid;
        }
    } while (fabs(dx) > tole);
    return x;
}

void delta_linking_logcoef(int dinucleotides, const double* best_bzenergy, double* logcoef)
{
    double bzenergy_scratch[dinucleotides];

    for (int i = 0; i < dinucleotides; i++) {
        bzenergy_scratch[i] = 1.0;
    }
    for (int i = 0; i < dinucleotides; i++) {
        double sum = 0.0;
        for (int j = 0; j < dinucleotides - i; j++) {
            bzenergy_scratch[j] *= best_bzenergy[i + j];
            sum += bzenergy_scratch[j];
        }
        logcoef[i] = log(sum);
    }
}

double find_delta_linking(int dinucleotides, double dtwist, const double* logcoef)
{
    return linear_search_dl(10.0, 50.0, 0.001, dtwist, logcoef, dinucleotides);
}

double delta_linking_slope(double dl, const double* logcoef, int terms)
{
    double sump, sump1, sumq, sumq1, x, y, z;

    double expmini;
    double exponent[terms];
    delta_linking_exponent(dl, terms, logcoef, &expmini, exponent);

    sump = sump1 = sumq = sumq1 = 0.0;
    x = 2.0 * _k_rt;
    for (int i = 0; i < terms; i++) {
        z = dl - bztwist[i];
        y = exp(exponent[i] + expmini);
        sumq += y;
        sump += bztwist[i] * y;
        y *= z * x;
        sumq1 += y;
        sump1 += bztwist[i] * y;
    }
    y = exp(_k_rt * dl * dl + sigma + expmini);
    sumq += y;
    sumq1 += x * dl * y;
    return (sump1 - sump * sumq1 / sumq) / sumq;
} /* slope at delta linking = dl */
