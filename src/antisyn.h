#pragma once

void antisyn_init(int dinucleotides);
void antisyn_destroy(void);

void assign_bzenergy_index(int nucleotides, char seq[]);
void find_best_antisyn(int dinucleotides, char* antisyn_out);
void antisyn_bzenergy(const char* antisyn, int dinucleotides, double* bzenergy);