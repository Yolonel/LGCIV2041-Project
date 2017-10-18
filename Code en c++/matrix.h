#ifndef _MATRIX_H_
#define _MATRIX_H_
struct matrix_t;
matrix_t * create_matrix (int n, int m);
void delete_matrix (matrix_t *K);
void add_to_matrix (matrix_t *K, int n, int m, double val);
double get_from_matrix (matrix_t *K, int n, int m);
bool solve_linear_system (matrix_t *K, double *f, double *x);
#endif
