#pragma once

void antisyn_init(void);
void antisyn_destroy(void);

void assign_bzenergy_index(int nucleotides, const char* seq, int* bzindex);
void find_best_antisyn(int dinucleotides, const int* bzindex, char* antisyn_out);
void antisyn_bzenergy(int dinucleotides, const char* antisyn, const int* bzindex, double* bzenergy);