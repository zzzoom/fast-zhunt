#include "delta_linking.h"

#include <stdlib.h>
#include <math.h>

int terms;
static double *bztwist, *logcoef, *exponent;
static double *bzenergy_scratch;
static const double _k_rt = -0.2521201; /* -1100/4363 */
static const double sigma = 16.94800353; /* 10/RT */
static double deltatwist;
static const double explimit = -600.0;

void delta_linking_init(int dinucleotides)
{
    static double a = 0.357, b = 0.4; /* a = 2 * (1/10.5 + 1/12) */
    double ab;

    bztwist = (double*)calloc(dinucleotides, sizeof(double));
    ab = b + b;
    for (int i = 0; i < dinucleotides; i++) {
        ab += a;
        bztwist[i] = ab;
    }
    logcoef = (double*)calloc(dinucleotides, sizeof(double));
    exponent = (double*)calloc(dinucleotides, sizeof(double));
    bzenergy_scratch = (double*)calloc(dinucleotides, sizeof(double));
}

void delta_linking_destroy(void)
{
    free(exponent);
    free(logcoef);
    free(bztwist);
    free(bzenergy_scratch);
}

static double linear_search(double x1, double x2, double tole, double (*func)())
{
    double f = func(x1);
    double fmid = func(x2);
    if (f * fmid >= 0.0) {
        return x2;
    }

    double dx, xmid;
    double x = (f < 0.0) ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
    do {
        dx *= 0.5;
        xmid = x + dx;
        fmid = func(xmid);
        if (fmid <= 0.0) {
            x = xmid;
        }
    } while (fabs(dx) > tole);
    return x;
}

static double delta_linking(double dl)
{
    double expmini = 0.0;
    for (int i = 0; i < terms; i++) {
        double z = dl - bztwist[i];
        exponent[i] = z = logcoef[i] + _k_rt * z * z;
        if (z < expmini) {
            expmini = z;
        }
    }
    expmini = (expmini < explimit) ? explimit - expmini : 0.0;
    double sump = 0.0;
    double sumq = 0.0;
    for (int i = 0; i < terms; i++) {
        double z = exp(exponent[i] + expmini);
        sumq += z;
        sump += bztwist[i] * z;
    }
    sumq += exp(_k_rt * dl * dl + sigma + expmini);
    return deltatwist - sump / sumq;
}

double delta_linking_slope(double dl)
{
    double sump, sump1, sumq, sumq1, x, y, z;

    double expmini = 0.0;
    for (int i = 0; i < terms; i++) {
        z = dl - bztwist[i];
        exponent[i] = z = logcoef[i] + _k_rt * z * z;
        if (z < expmini) {
            expmini = z;
        }
    }
    expmini = (expmini < explimit) ? explimit - expmini : 0.0;
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

double find_delta_linking(int dinucleotides, double dtwist, const double* best_bzenergy)
{
    double sum;
    int i, j;

    for (i = 0; i < dinucleotides; i++)
        bzenergy_scratch[i] = 1.0;
    for (i = 0; i < dinucleotides; i++) {
        sum = 0.0;
        for (j = 0; j < dinucleotides - i; j++) {
            bzenergy_scratch[j] *= best_bzenergy[i + j];
            sum += bzenergy_scratch[j];
        }
        logcoef[i] = log(sum);
    }
    terms = dinucleotides;
    deltatwist = dtwist;
    return linear_search(10.0, 50.0, 0.001, delta_linking);
}
