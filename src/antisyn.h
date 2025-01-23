#pragma once

void antisyn_init(int dinucleotides);
void antisyn_destroy(void);

void assign_bzenergy_index(int nucleotides, char seq[]);
void best_anti_syn(char* antisyn, int dinucleotides, double esum);

void anti_syn_energy(int dinucleotides, double max_esum, char* antisyn_out, double* bzenergy_out);
