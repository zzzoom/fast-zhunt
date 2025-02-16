#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t esum_t;

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

void antisyn_init()
{
    static double rt = 0.59004; /* 0.00198*298 */
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 16; j++) {
            expdbzed[i][j] = exp(-dbzed[i][j] / rt);
        }
    }
}

void antisyn_destroy(void) {}

void assign_bzenergy_index(int nucleotides, const char* seq, int* bzindex)
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

static void antisyn_01_bzenergy(int dinucleotides, const char* antisyn, const int* bzindex, double* bzenergy)
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

void antisyn_bzenergy(int dinucleotides, const char* antisyn_string, const int* bzindex, double* bzenergy)
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

void find_best_antisyn(int dinucleotides, const int* bzindex, char* antisyn_out)
{
    if (dinucleotides < 1) {
        return;
    }

    esum_t best0_esum = int_dbzed[0][bzindex[0]];
    char best0_antisyn[dinucleotides];
    char best0_prev_antisyn[dinucleotides];
    best0_antisyn[0] = 0;

    esum_t best1_esum = int_dbzed[3][bzindex[0]];
    char best1_antisyn[dinucleotides];
    best1_antisyn[0] = 1;

    for (int din = 1; din < dinucleotides; ++din) {
        const esum_t dbzed00 = int_dbzed[0][bzindex[din]];
        const esum_t dbzed01 = int_dbzed[1][bzindex[din]];
        const esum_t dbzed10 = int_dbzed[2][bzindex[din]];
        const esum_t dbzed11 = int_dbzed[3][bzindex[din]];

        const esum_t prev_best0 = best0_esum;
        const esum_t prev_best1 = best1_esum;
        // best0 is dirty when processing best1, copy the original value
        memcpy(best0_prev_antisyn, best0_antisyn, din);

        esum_t esum00 = prev_best0 + dbzed00;
        esum_t esum10 = prev_best1 + dbzed10;
        // relatively expensive comparison to preserve original order
        if ((esum00 < esum10) || ((esum00 == esum10) && (strncmp(best0_antisyn, best1_antisyn, din) <= 0))) {
            best0_esum = esum00;
        } else {
            best0_esum = esum10;
            memcpy(best0_antisyn, best1_antisyn, din);
        }
        best0_antisyn[din] = 0;

        esum_t esum01 = prev_best0 + dbzed01;
        esum_t esum11 = prev_best1 + dbzed11;
        if ((esum11 < esum01) || ((esum11 == esum01) && (strncmp(best1_antisyn, best0_antisyn, din) < 0))) {
            best1_esum = esum11;
        } else {
            best1_esum = esum01;
            memcpy(best1_antisyn, best0_prev_antisyn, din);
        }
        best1_antisyn[din] = 1;
    }

    char* best_antisyn = best0_esum <= best1_esum ? best0_antisyn : best1_antisyn;
    antisyn_string(best_antisyn, dinucleotides, antisyn_out);
}
