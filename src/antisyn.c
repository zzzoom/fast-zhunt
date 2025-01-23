#include <error.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

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

// static double *bzenergy, *best_bzenergy; /* dinucleotides */
typedef struct {
    double esum;
    char* antisyn;
} Candidate;

Candidate g_best;
Candidate g_scratch;

void antisyn_init(int dinucleotides)
{
    static double rt = 0.59004; /* 0.00198*298 */

    g_best.antisyn = (char*)malloc(2 * dinucleotides + 1);
    g_scratch.antisyn = (char*)malloc(2 * dinucleotides + 1);

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
    free(g_scratch.antisyn);
    free(bzindex);
    // free(best_bzenergy);
}

void assign_bzenergy_index(int nucleotides, char seq[])
{
    int i, j, idx;
    i = j = 0;
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

static void best_anti_syn(const char *antisyn, int dinucleotides, double esum)
{
    if (esum < g_best.esum) {
        g_best.esum = esum;
        memcpy(g_best.antisyn, antisyn, dinucleotides);
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
    } else {
        anti_syn_energy_rec(antisyn, esum0, din + 1, dinucleotides);
    }

    antisyn[din] = 1;
    int i1 = (din == 0) ? 3 : ((antisyn[din - 1] == 1) ? 3 : 1);
    double e1 = dbzed[i1][bzindex[din]];
    double esum1 = esum + e1;

    
    if (din + 1 == dinucleotides) {
        best_anti_syn(antisyn, dinucleotides, esum1);
    } else {
        anti_syn_energy_rec(antisyn, esum1, din + 1, dinucleotides);
    }
}

void anti_syn_energy(int dinucleotides, double max_esum, char* antisyn_out, double* bzenergy_out)
{
    g_best.esum = max_esum;
    anti_syn_energy_rec(g_scratch.antisyn, 0.0, 0, dinucleotides);

    antisyn_string(g_best.antisyn, dinucleotides, antisyn_out);
    antisyn_bzenergy(g_best.antisyn, dinucleotides, bzenergy_out);
}