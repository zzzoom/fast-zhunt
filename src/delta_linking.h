#pragma once

void delta_linking_init(int dinucleotides);
void delta_linking_destroy(void);
double delta_linking_slope(double dl, int terms);
double find_delta_linking(int dinucleotides, double deltatwist, const double* best_bzenergy);