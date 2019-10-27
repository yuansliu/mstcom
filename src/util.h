#ifndef MSTCOM_UTIL_H_H
#define MSTCOM_UTIL_H_H

#include <stdio.h>
#include <stdint.h>

void gen_eight_level_bin_table(char *eight_level_bin_tab);
void gen_binary_table(char *binary_tab, unsigned int thr = 20, unsigned int low = 6, unsigned int high = 40);
// void gen_binary_table(char *binary_tab, unsigned int thr, unsigned int low, unsigned int high);

#endif