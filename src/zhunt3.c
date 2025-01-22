/*
Z-HUNT-2 computer program, Sept.19, 1990
*/

/*
serialized i/o to allow to run against large datasets
Updated i/o permissions to match posix. campt 1/10/2000
*/

/*
Turbo C compiler Ver. 1.5
Compact model, 8086/80286 instruction set, math emulator/8087/80287
Speed option, unsigned char
Run on IBM PC XT/AT or compatibles
*/

/*
Written by Ping-jung Chou, under the instruction of Pui S. Ho, according to the
paper "A computer aided thermodynamic approach for predicting the formation of
Z-DNA in naturally occurring sequences",
The EMBO Journal, Vol.5, No.10, pp2737-2744, 1986
With 0.22 kcal/mol/dinuc for mCG (Zacharias et al, Biochemistry, 1988, 2970)
*/

#include "delta_linking.h"

#include <error.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* creates a temporary file and mmaps the sequence into it */
#ifdef USE_MMAP
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#endif

/* Delta BZ Energy of Dinucleotide */
static const double dbzed[4][16] = {
    /* AS-AS */
    { 4.40, 6.20, 3.40, 5.20, 2.50, 4.40, 1.40, 3.30, 3.30, 5.20, 2.40, 4.20, 1.40, 3.40, 0.66, 2.40 },
    /* AS-SA */
    { 6.20, 6.20, 5.20, 5.20, 6.20, 6.20, 5.20, 5.20, 5.20, 5.20, 4.00, 4.00, 5.20, 5.20, 4.00, 4.00 },
    /* SA-AS */
    { 6.20, 6.20, 5.20, 5.20, 6.20, 6.20, 5.20, 5.20, 5.20, 5.20, 4.00, 4.00, 5.20, 5.20, 4.00, 4.00 },
    /* SA-SA */
    { 4.40, 2.50, 3.30, 1.40, 6.20, 4.40, 5.20, 3.40, 3.40, 1.40, 2.40, 0.66, 5.20, 3.30, 4.20, 2.40 }
};
static double expdbzed[4][16]; /* exp(-dbzed/rt) */
static int* bzindex; /* dinucleotides */

static double *bzenergy, *best_bzenergy; /* dinucleotides */
static double best_esum; /* assigned before call to anti_syn_energy() */
static char* best_antisyn; /* nucleotides */

char *tempstr, *sequence;
#ifdef USE_MMAP
int sequencefile;
#endif

static double assign_probability(double dl);

static void analyze_zscore(char* filename);
static void calculate_zscore(double a, int maxdinucleotides, int min, int max, char* filename);

static void best_anti_syn(char* antisyn, int dinucleotides, double esum);
static void anti_syn_energy(char* antisyn, int din, int dinucleotides, double esum);

static FILE* open_file(int mode, char* filename, char* typestr);
static void assign_bzenergy_index(int nucleotides, char seq[]);
static unsigned input_sequence(FILE* file, int nucleotides, int showfile);

static void assign_bzenergy_index(int nucleotides, char seq[])
{
    int i, j, idx;
    char c1, c2;

    i = j = 0;
    do {
        c1 = seq[i++];
        c2 = seq[i++];
        switch (c1) {
        case 'a':
            switch (c2) {
            case 'a':
                idx = 0;
                break;
            case 't':
                idx = 1;
                break;
            case 'g':
                idx = 2;
                break;
            case 'c':
                idx = 3;
            }
            break;
        case 't':
            switch (c2) {
            case 'a':
                idx = 4;
                break;
            case 't':
                idx = 5;
                break;
            case 'g':
                idx = 6;
                break;
            case 'c':
                idx = 7;
            }
            break;
        case 'g':
            switch (c2) {
            case 'a':
                idx = 8;
                break;
            case 't':
                idx = 9;
                break;
            case 'g':
                idx = 10;
                break;
            case 'c':
                idx = 11;
            }
            break;
        case 'c':
            switch (c2) {
            case 'a':
                idx = 12;
                break;
            case 't':
                idx = 13;
                break;
            case 'g':
                idx = 14;
                break;
            case 'c':
                idx = 15;
            }
        }
        bzindex[j++] = idx;
    } while (i < nucleotides);
}

static void antisyn_string_bzenergy(const char* antisyn_string, int dinucleotides, double* bzenergy)
{
    if (dinucleotides == 0) {
        return;
    }

    int i = antisyn_string[0] == 'A' ? 0 : 3;
    bzenergy[0] = expdbzed[i][bzindex[0]];
    for (int din = 1; din < dinucleotides; ++din) {
        if (antisyn_string[2*din] == 'A') {
            i = (antisyn_string[2*din - 2] == 'A') ? 0 : 2;
        } else if (antisyn_string[2*din] == 'S') {
            i = (antisyn_string[2*din - 2] == 'S') ? 3 : 1;
        } else {
            error(1, 1, "antisyn_string_bzenergy: shouldn't be here");
        }
        bzenergy[din] = expdbzed[i][bzindex[din]];
    }
}

static void calculate_bzenergy(const char* antisyn, int dinucleotides, double* bzenergy)
{
    if (dinucleotides == 0) {
        return;
    }

    int i = antisyn[0] == 0 ? 0 : 3;
    bzenergy[0] = expdbzed[i][bzindex[0]];
    for (int din = 1; din < dinucleotides; ++din) {
        if (antisyn[din] == 0) {
            i = (antisyn[din - 1] == 0) ? 0 : 2;
        } else if (antisyn[din] == 1) {
            i = (antisyn[din - 1] == 1) ? 3 : 1;
        } else {
            error(1, 1, "calculate_bzenergy: shouldn't be here");
        }
        bzenergy[din] = expdbzed[i][bzindex[din]];
    }
}

static void antisyn_string(const char* antisyn, int dinucleotides, char* dest)
{
    for (int din = 0; din < dinucleotides; ++din) {
        if (antisyn[din] == 0) {
            dest[2 * din] = 'A';
            dest[2 * din + 1] = 'S';
        } else {
            dest[2 * din] = 'S';
            dest[2 * din + 1] = 'A';
        }
    }
    dest[2 * dinucleotides] = '\0';
}

static void best_anti_syn(char* antisyn, int dinucleotides, double esum)
{
    if (esum < best_esum) {
        best_esum = esum;
        antisyn_string(antisyn, dinucleotides, best_antisyn);
    }
}

static void anti_syn_energy(char* antisyn, int din, int dinucleotides, double esum)
{
    antisyn[din] = 0;
    int i1 = (din == 0) ? 0 : ((antisyn[din - 1] == 0) ? 0 : 2);
    /*
    if (din > 0) {
      printf("%c%c-%c%c %d\n", antisyn[nucleotides-2], antisyn[nucleotides-1], antisyn[nucleotides], antisyn[nucleotides+1], i1);
    }
    */
    double e1 = dbzed[i1][bzindex[din]];

    if (din + 1 == dinucleotides) {
        best_anti_syn(antisyn, dinucleotides, esum + e1);
    } else {
        anti_syn_energy(antisyn, din + 1, dinucleotides, esum + e1);
    }

    antisyn[din] = 1;
    int i2 = (din == 0) ? 3 : ((antisyn[din - 1] == 1) ? 3 : 1);
    /*
    if (din > 0) {
      printf("%c%c-%c%c %d\n", antisyn[nucleotides-2], antisyn[nucleotides-1], antisyn[nucleotides], antisyn[nucleotides+1], i2);
    }
    */
    double e2 = dbzed[i2][bzindex[din]];

    if (din + 1 == dinucleotides) {
        best_anti_syn(antisyn, dinucleotides, esum + e2);
    } else {
        anti_syn_energy(antisyn, din + 1, dinucleotides, esum + e2);
    }
}

static FILE* open_file(int mode, char* filename, char* typestr)
{
    static char* rwstr[] = { "w", "r" };
    char* fullfile;

    FILE* file;
    file = NULL;
    fullfile = (char*)malloc(sizeof(char) * (strlen(filename) + strlen(typestr) + 2));
    strcpy(fullfile, filename);
    if (strlen(typestr) != 0) {
        strcat(fullfile, ".");
        strcat(fullfile, typestr);
    }
    printf("opening %s\n", fullfile);

    file = fopen(fullfile, rwstr[mode]);
    free(fullfile);
    return file;
}

static unsigned input_sequence(FILE* file, int nucleotides, int showfile)
{
    unsigned length, i, j;
    char c;
#ifdef USE_MMAP
    FILE* OUTPUT;
#endif

    printf("inputting sequence\n");

    length = 0; /* count how many bases */
    while (j = 0, fgets(tempstr, 128, file) != NULL) {
        while ((c = tempstr[j++]) != 0) {
            if (c == 'a' || c == 't' || c == 'g' || c == 'c' || c == 'A' || c == 'T' || c == 'G' || c == 'C') {
                length++;
            }
        }
    }

#ifndef USE_MMAP
    sequence = (char*)malloc(length + nucleotides);
#else
    sequence = (char*)malloc(nucleotides);
    OUTPUT = (FILE*)fopen("/usr/local/apache/htdocs/zhunt/temp", "w");
#endif

    rewind(file);

    if (showfile) {
        printf("\n");
    }
    i = 0;
    while (j = 0, fgets(tempstr, 128, file) != NULL) {
        while ((c = tempstr[j++]) != 0) {
            if (c == 'a' || c == 't' || c == 'g' || c == 'c' || c == 'A' || c == 'T' || c == 'G' || c == 'C') {
                if (c == 'A') {
                    c = 'a';
                }
                if (c == 'T') {
                    c = 't';
                }
                if (c == 'G') {
                    c = 'g';
                };
                if (c == 'C') {
                    c = 'c';
                }
#ifndef USE_MMAP
                sequence[i++] = c;
#else
                fprintf(OUTPUT, "%c", c);
                if (i < nucleotides) {
                    sequence[i++] = c;
                }
#endif
            }
        }
    }
    for (j = 0; j < nucleotides; j++) /* assume circular nucleotides */
    {
#ifndef USE_MMAP
        sequence[i++] = sequence[j];
#else
        fprintf(OUTPUT, "%c", sequence[j]);
#endif
    }
#ifdef USE_MMAP
    close(OUTPUT);
    sequencefile = open("/usr/local/apache/htdocs/zhunt/temp", O_RDWR, NULL);
    free(sequence);
    sequence = (char*)mmap(0, length, PROT_READ | PROT_WRITE, MAP_SHARED, sequencefile, 0);
#endif

    return length;
}

/* calculate the probability of the value 'dl' in a Gaussian distribution */
/* from "Data Reduction and Error Analysis for the Physical Science" */
/* Philip R. Bevington, 1969, McGraw-Hill, Inc */

static double assign_probability(double dl)
{
    static double average = 29.6537135;
    static double stdv = 2.71997;
    static double _sqrt2 = 0.70710678118654752440; /* 1/sqrt(2) */
    static double _sqrtpi = 0.564189583546; /* 1/sqrt(pi) */

    double x, y, z, k, sum;

    z = fabs(dl - average) / stdv;
    x = z * _sqrt2;
    y = _sqrtpi * exp(-x * x);
    z *= z;
    k = 1.0;
    sum = 0.0;
    do {
        sum += x;
        k += 2.0;
        x *= z / k;
    } while (sum + x > sum);
    z = 0.5 - y * sum; /* probability of each tail */
    return (dl > average) ? z : 1.0 / z;
}

int main(int argc, char* argv[])
{
    static double a = 0.357;
    static double rt = 0.59004; /* 0.00198*298 */

    int i, j, nucleotides, dinucleotides;
    int min, max;

    if (argc < 5) {
        printf("usage: zhunt windowsize minsize maxsize datafile\n");
        exit(1);
    }
    tempstr = (char*)malloc(128);
    dinucleotides = atoi((char*)argv[1]);
    min = atoi((char*)argv[2]);
    max = atoi((char*)argv[3]);

    printf("dinucleotides %d\n", dinucleotides);
    printf("min/max %d %d\n", min, max);
    printf("operating on %s\n", (char*)argv[4]);

    nucleotides = 2 * dinucleotides;

    delta_linking_init(dinucleotides);

    best_antisyn = (char*)malloc(nucleotides + 1);

    bzindex = (int*)calloc(dinucleotides, sizeof(int));
    best_bzenergy = (double*)calloc(dinucleotides, sizeof(double));

    for (i = 0; i < 4; i++)
        for (j = 0; j < 16; j++)
            expdbzed[i][j] = exp(-dbzed[i][j] / rt);

    calculate_zscore(a, dinucleotides, min, max, (char*)argv[4]);
#ifndef PROB_ONLY
    analyze_zscore((char*)argv[4]);
#endif

    free(best_bzenergy);
    free(bzindex);
    free(best_antisyn);
    free(tempstr);
    delta_linking_destroy();
    return 0;
}

static void calculate_zscore(double a, int maxdinucleotides, int min, int max, char* filename)
{
    static double pideg = 57.29577951; /* 180/pi */
    char* bestantisyn;
    FILE* file;
    unsigned seqlength, i;
    int fromdin, todin, din, nucleotides;
    long begintime, endtime;
    double dl, slope, probability, bestdl;
    float initesum;

    fromdin = min;
    todin = max;
    printf("calculating zscore\n");

    file = open_file(1, filename, "");
    if (file == NULL) {
        printf("couldn't open %s!\n", filename);
        return;
    }
    seqlength = input_sequence(file, 2 * maxdinucleotides, 0);
    fclose(file);

    file = open_file(0, filename, "Z-SCORE");
    if (file == NULL) {
#ifndef USE_MMAP
        free(sequence);
#else
        munmap(sequence, seqlength);
        close(sequencefile);
#endif
        return;
    }

    if (todin > maxdinucleotides) {
        todin = maxdinucleotides;
    }
    if (fromdin > todin) {
        fromdin = todin;
    }
    nucleotides = 2 * todin;

#ifndef PROB_ONLY
    fprintf(file, "%s %u %d %d\n", filename, seqlength, fromdin, todin);
#endif

    a /= 2.0;
    initesum = 10.0 * todin;
    bestantisyn = (char*)malloc(nucleotides + 1);
    char* antisyn = (char*)malloc(nucleotides + 1);

    time(&begintime);
    for (i = 0; i < seqlength; i++) {
        assign_bzenergy_index(nucleotides, sequence + i);
        bestdl = 50.0;
        for (din = fromdin; din <= todin; din++) {
            best_esum = initesum;
            antisyn[2 * din] = 0;
            anti_syn_energy(antisyn, 0, din, 0.0); /* esum = 0.0 */
            antisyn_string_bzenergy(best_antisyn, din, best_bzenergy);
            dl = find_delta_linking(din, a * (double)din, best_bzenergy);
            if (dl < bestdl) {
                bestdl = dl;
                strcpy(bestantisyn, best_antisyn);
            }
        }
#ifndef PROB_ONLY
        slope = atan(delta_linking_slope(bestdl)) * pideg;
#endif
        probability = assign_probability(bestdl);
#ifndef PROB_ONLY
        fprintf(file, " %7.3lf %7.3lf %le %s\n", bestdl, slope, probability, bestantisyn);
#else
        fprintf(file, " %7.3lf %le\n", bestdl, probability);
#endif
    }
    time(&endtime);

    free(bestantisyn);
    free(antisyn);
    fclose(file);
    printf("\n run time=%ld sec\n", endtime - begintime);
#ifndef USE_MMAP
    free(sequence);
#else
    munmap(sequence, seqlength);
    close(sequencefile);
#endif
}

static void analyze_zscore(char* filename)
{
    float *dl, *slope, *probability;
    unsigned seqlength, i;
    char** antisyn;
    FILE* file;
    int fromdin, todin, nucleotides;

    printf("analyzing_zscore\n");

    file = open_file(1, filename, "Z-SCORE");
    if (file == NULL) {
        printf("couldn't open %s.Z-SCORE!\n", filename);
        return;
    }
    fscanf(file, "%s %u %d %d", tempstr, &seqlength, &fromdin, &todin);
    nucleotides = 2 * todin;
    dl = (float*)calloc(seqlength, sizeof(float));
    slope = (float*)calloc(seqlength, sizeof(float));
    probability = (float*)calloc(seqlength, sizeof(float));
    antisyn = (char**)calloc(seqlength, sizeof(char*));

    for (i = 0; i < seqlength; i++) {
        fscanf(file, "%f %f %f %s", dl + i, slope + i, probability + i, tempstr);
        antisyn[i] = strdup(tempstr);
    }
    fclose(file);

    file = open_file(1, filename, "");
    input_sequence(file, nucleotides, 0);
    fclose(file);

    free(dl);
    free(slope);
    free(probability);
    for (i = 0; i < seqlength; i++)
        free(antisyn[i]);
    free(antisyn);
#ifndef USE_MMAP
    free(sequence);
#else
    munmap(sequence, seqlength);
    close(sequencefile);
#endif
}