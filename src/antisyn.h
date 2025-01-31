#pragma once

void antisyn_init(int dinucleotides);
void antisyn_destroy(void);

void assign_bzenergy_index(int nucleotides, char seq[]);
void anti_syn_energy(int dinucleotides, char* antisyn_out, double* bzenergy_out);