#ifndef matrix_t_h
#define matrix_t_h


typedef struct{
	//nbre ligne
	int n;
	//nbre colonne
	int m;
	//elment of matrix_t
	double **mData;
}matrix_t;



matrix_t * create_matrix(int n, int m);

void delete_matrix(matrix_t *K);

void add_to_matrix(matrix_t *K, int n, int m, double val);

int solve_linear_system(matrix_t *K, double *f, double *x);

void print_m(matrix_t *K);

void print(double *K, int n);


#endif