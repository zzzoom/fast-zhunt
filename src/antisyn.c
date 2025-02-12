#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int* bzindex; /* dinucleotides */

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
/* Integer version of the above for exact sums */
static const int int_dbzed[4][16] = {
    /* AS-AS */
    { 440, 620, 340, 520, 250, 440, 140, 330, 330, 520, 240, 420, 140, 340, 66, 240 },
    /* AS-SA */
    { 620, 620, 520, 520, 620, 620, 520, 520, 520, 520, 400, 400, 520, 520, 400, 400 },
    /* SA-AS */
    { 620, 620, 520, 520, 620, 620, 520, 520, 520, 520, 400, 400, 520, 520, 400, 400 },
    /* SA-SA */
    { 440, 250, 330, 140, 620, 440, 520, 340, 340, 140, 240, 66, 520, 330, 420, 240 }
};
static double expdbzed[4][16]; /* exp(-dbzed/rt) */

typedef int64_t esum_t;
typedef struct {
    esum_t esum;
    char* antisyn;
} Candidate;

Candidate g_best0;
Candidate g_best0_prev;
Candidate g_best1;

void antisyn_init(int dinucleotides)
{
    static double rt = 0.59004; /* 0.00198*298 */

    g_best0.antisyn = (char*)malloc(2 * dinucleotides + 1);
    g_best0_prev.antisyn = (char*)malloc(2 * dinucleotides + 1);
    g_best1.antisyn = (char*)malloc(2 * dinucleotides + 1);

    bzindex = (int*)calloc(dinucleotides, sizeof(int));

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 16; j++) {
            expdbzed[i][j] = exp(-dbzed[i][j] / rt);
        }
    }
}

void antisyn_destroy(void)
{
    free(g_best0.antisyn);
    free(g_best0_prev.antisyn);
    free(g_best1.antisyn);
    free(bzindex);
}

void assign_bzenergy_index(int nucleotides, char seq[])
{
    int i = 0;
    int j = 0;
    int idx = 0;
    do {
        char c1 = seq[i++];
        char c2 = seq[i++];
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

static void antisyn_01_bzenergy(const char* antisyn, int dinucleotides, double* bzenergy)
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
            fprintf(stderr, "antisyn_bzenergy: wrong value %d in antisyn\n", antisyn[din]);
            exit(EXIT_FAILURE);
        }
        bzenergy[din] = expdbzed[i][bzindex[din]];
    }
}

void antisyn_bzenergy(const char* antisyn_string, int dinucleotides, double* bzenergy)
{
    if (dinucleotides == 0) {
        return;
    }

    int i = antisyn_string[0] == 'A' ? 0 : 3;
    bzenergy[0] = expdbzed[i][bzindex[0]];
    for (int din = 1; din < dinucleotides; ++din) {
        if (antisyn_string[2 * din] == 'A') {
            i = (antisyn_string[2 * din - 2] == 'A') ? 0 : 2;
        } else if (antisyn_string[2 * din] == 'S') {
            i = (antisyn_string[2 * din - 2] == 'S') ? 3 : 1;
        } else {
            fprintf(stderr, "antisyn_bzenergy: wrong value %c in antisyn\n", antisyn_string[2 * din]);
            exit(EXIT_FAILURE);
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
        } else if (antisyn[din] == 1) {
            dest[2 * din] = 'S';
            dest[2 * din + 1] = 'A';
        } else {
            fprintf(stderr, "antisyn_string: wrong value %d in antisyn\n", antisyn[din]);
            exit(EXIT_FAILURE);
        }
    }
    dest[2 * dinucleotides] = '\0';
}

void find_best_antisyn(int dinucleotides, char* antisyn_out)
{
    if (dinucleotides < 1) {
        return;
    }

    g_best0.esum = int_dbzed[0][bzindex[0]];
    g_best0.antisyn[0] = 0;

    g_best1.esum = int_dbzed[3][bzindex[0]];
    g_best1.antisyn[0] = 1;

    for (int din = 1; din < dinucleotides; ++din) {
        const esum_t dbzed00 = int_dbzed[0][bzindex[din]];
        const esum_t dbzed01 = int_dbzed[1][bzindex[din]];
        const esum_t dbzed10 = int_dbzed[2][bzindex[din]];
        const esum_t dbzed11 = int_dbzed[3][bzindex[din]];

        const esum_t prev_best0 = g_best0.esum;
        const esum_t prev_best1 = g_best1.esum;
        // g_best0 is dirty when processing g_best1
        memcpy(g_best0_prev.antisyn, g_best0.antisyn, din);

        esum_t esum00 = prev_best0 + dbzed00;
        esum_t esum10 = prev_best1 + dbzed10;
        // relatively expensive comparison to preserve original order
        if ((esum00 < esum10) || ((esum00 == esum10) && (strncmp(g_best0.antisyn, g_best1.antisyn, din) <= 0))) {
            g_best0.esum = esum00;
        } else {
            g_best0.esum = esum10;
            memcpy(g_best0.antisyn, g_best1.antisyn, din);
        }
        g_best0.antisyn[din] = 0;

        esum_t esum01 = prev_best0 + dbzed01;
        esum_t esum11 = prev_best1 + dbzed11;
        if ((esum11 < esum01) || ((esum11 == esum01) && (strncmp(g_best1.antisyn, g_best0.antisyn, din) < 0))) {
            g_best1.esum = esum11;
        } else {
            g_best1.esum = esum01;
            memcpy(g_best1.antisyn, g_best0_prev.antisyn, din);
        }
        g_best1.antisyn[din] = 1;
    }

    Candidate* best = g_best0.esum <= g_best1.esum ? &g_best0 : &g_best1;
    antisyn_string(best->antisyn, dinucleotides, antisyn_out);
}
