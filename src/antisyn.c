#include <error.h>
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

Candidate g_best;
Candidate g_scratch;
Candidate g_best0;
Candidate g_best0_prev;
Candidate g_best1;

void antisyn_init(int dinucleotides)
{
    static double rt = 0.59004; /* 0.00198*298 */

    g_best.antisyn = (char*)malloc(2 * dinucleotides + 1);
    g_scratch.antisyn = (char*)malloc(2 * dinucleotides + 1);

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
    free(g_best.antisyn);
    free(g_best0.antisyn);
    free(g_best0_prev.antisyn);
    free(g_best1.antisyn);
    free(g_scratch.antisyn);
    free(bzindex);
    // free(best_bzenergy);
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

static void antisyn_bzenergy(const char* antisyn, int dinucleotides, double* bzenergy)
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
            fprintf(stderr, "din: %d, total: %d, as: %d\n", din, dinucleotides, antisyn[din]);
            fflush(stderr);
            error(1, 1, "antisyn_bzenergy: shouldn't be here");
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

static void best_anti_syn(const char *antisyn, int dinucleotides, esum_t esum)
{
    if (esum < g_best.esum) {
        g_best.esum = esum;
        memcpy(g_best.antisyn, antisyn, dinucleotides);
    }
}

/*
#define ANTISYN_BATCH_COMPARABLE_RESULTS
#define ANTISYN_BATCH_ROUNDS 4
#define ANTISYN_BATCH_SIZE (1 << ANTISYN_BATCH_ROUNDS)
#define ANTISYN_BATCH_MASK (ANTISYN_BATCH_SIZE - 1)
static double batch_dbzed[ANTISYN_BATCH_ROUNDS][ANTISYN_BATCH_SIZE];

// Run the last 2^ANTISYN_BATCH_SIZE elements together instead of recursively
// in order to vectorize operations
// Compile with aggresive vectorization to get gains (-O3 -march=native)
static void anti_syn_energy_batch(char* antisyn, double esum, int din, int dinucleotides)
{
    if (din + ANTISYN_BATCH_ROUNDS != dinucleotides) {
        error(1, 1, "anti_syn_energy_batch: wrong depth");
    }

    // Calculate indices and gather dbzeds into batch_dbzed
    int antisyn_prev = antisyn[din - 1];
    int antisyn_prev_bits = antisyn_prev << ANTISYN_BATCH_ROUNDS;
    for (int round = 0; round < ANTISYN_BATCH_ROUNDS; ++round) {
        for (int antisyn_index = 0; antisyn_index < ANTISYN_BATCH_SIZE; ++antisyn_index)
        {
            int antisyn_bits = antisyn_prev_bits + antisyn_index;
            int bzidx = (antisyn_bits >> (ANTISYN_BATCH_ROUNDS - round - 1)) & 0x3;
            batch_dbzed[round][antisyn_index] = dbzed[bzidx][bzindex[din+round]];
        }
    }

    #ifdef ANTISYN_BATCH_COMPARABLE_RESULTS
    // Accumulate (((esum + batch[0]) + ...) + batch[n]) as original
    for (int antisyn_index = 0; antisyn_index < ANTISYN_BATCH_SIZE; ++antisyn_index) {
        batch_dbzed[0][antisyn_index] += esum;
    }
    #endif

    // Easy to vectorize per-lane esum reduction
    for (int round = 1; round < ANTISYN_BATCH_ROUNDS; ++round) {
        for (int antisyn_index = 0; antisyn_index < ANTISYN_BATCH_SIZE; ++antisyn_index) {
            batch_dbzed[0][antisyn_index] += batch_dbzed[round][antisyn_index];
        }
    }

    // Find minimum esum in batch
    double min_esum = batch_dbzed[0][0];
    int min_index = 0;
    for (int antisyn_index = 0; antisyn_index < ANTISYN_BATCH_SIZE; ++antisyn_index) {
        if (batch_dbzed[0][antisyn_index] < min_esum) {
            min_esum = batch_dbzed[0][antisyn_index];
            min_index = antisyn_index;
        }
    }

    #ifdef ANTISYN_BATCH_COMPARABLE_RESULTS
    double best_batch_esum = min_esum;
    #else
    // Accumulate esum + (((batch[0] + batch[1]) + ...) + batch[n])
    double best_batch_esum = esum + min_esum;
    #endif

    if (best_batch_esum < g_best.esum) {
        g_best.esum = best_batch_esum;
        for (int i = 0; i < din; ++i) {
            g_best.antisyn[i] = antisyn[i];
        }
        for (int i = 0; i < ANTISYN_BATCH_ROUNDS; ++i) {
            g_best.antisyn[din + i] = (min_index >> (ANTISYN_BATCH_ROUNDS - i - 1)) & 0x1;
        }
    }
}

static void anti_syn_energy_rec(char* antisyn, double esum, int din, int dinucleotides)
{
    antisyn[din] = 0;
    int i0 = (din == 0) ? 0 : ((antisyn[din - 1] == 0) ? 0 : 2);
    double e0 = dbzed[i0][bzindex[din]];
    double esum0 = esum + e0;

    if (din + 1 == dinucleotides) {
        best_anti_syn(antisyn, dinucleotides, esum0);
    } else if (din + 1 + ANTISYN_BATCH_ROUNDS == dinucleotides) {
        anti_syn_energy_batch(antisyn, esum0, din + 1, dinucleotides);
    } else {
        anti_syn_energy_rec(antisyn, esum0, din + 1, dinucleotides);
    }

    antisyn[din] = 1;
    int i1 = (din == 0) ? 3 : ((antisyn[din - 1] == 1) ? 3 : 1);
    double e1 = dbzed[i1][bzindex[din]];
    double esum1 = esum + e1;

    if (din + 1 == dinucleotides) {
        best_anti_syn(antisyn, dinucleotides, esum1);
    } else if (din + 1 + ANTISYN_BATCH_ROUNDS == dinucleotides) {
        anti_syn_energy_batch(antisyn, esum1, din + 1, dinucleotides);
    } else {
        anti_syn_energy_rec(antisyn, esum1, din + 1, dinucleotides);
    }
}
*/

/*
void anti_syn_energy(int dinucleotides, double max_esum, char* antisyn_out, double* bzenergy_out)
{
    g_best.esum = max_esum;
    anti_syn_energy_rec(g_scratch.antisyn, 0.0, 0, dinucleotides);

    antisyn_string(g_best.antisyn, dinucleotides, antisyn_out);
    antisyn_bzenergy(g_best.antisyn, dinucleotides, bzenergy_out);
}
*/

void anti_syn_energy(int dinucleotides, double max_esum, char* antisyn_out, double* bzenergy_out)
{
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
        memcpy(g_best0_prev.antisyn, g_best0.antisyn, din);

        esum_t esum00 = prev_best0 + dbzed00;
        esum_t esum10 = prev_best1 + dbzed10;
        // relatively expensive comparison to preserve comparable results
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

    if (g_best0.esum <= g_best1.esum) {
        g_best.esum = g_best0.esum;
        antisyn_string(g_best0.antisyn, dinucleotides, antisyn_out);
        antisyn_bzenergy(g_best0.antisyn, dinucleotides, bzenergy_out);
    } else {
        g_best.esum = g_best1.esum;
        antisyn_string(g_best1.antisyn, dinucleotides, antisyn_out);
        antisyn_bzenergy(g_best1.antisyn, dinucleotides, bzenergy_out);
    }
}
