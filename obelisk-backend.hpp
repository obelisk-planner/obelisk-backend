#ifndef _OBELISK_BACKEND_
#define _OBELISK_BACKEND_
#include "obelisk-compressed-rrmat.hpp"

  using namespace std;
    
    int swap_pointers(double** x, double** y);
    int get_rr_inv_pair(double** input_mat, double** inv, compressed_rrmat* lhs, int dim);
    int print_nz(int size, double* array);
    int active_set_alg(compressed_rrmat* lhs, int num_rows, double** inv, double* solution, double* kkt_rhs);
    int update_inv(compressed_rrmat* lhs, int* num_rows, double** inv, int inv_cols, int new_active);
    int main_algorithm(double* solution);
#endif
