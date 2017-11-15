#include "matrix.h"

extern "C" {
  void dgesv(int *N, int *nrhs, double *A, int *lda, int *ipiv,
	     double *b, int *ldb, int *info);
}

struct matrix_t {
  int _n,_m;
  double *_data;
  matrix_t (int n, int m) : _n(n), _m(m){
    _data = new double [n*m];
    for (int i=0;i<n*m;i++)_data[i] = 0.0;
  }
};

matrix_t * create_matrix (int n, int m){
  return new matrix_t (n,m);
}
void delete_matrix (matrix_t *K) {
  delete K;
}

void add_to_matrix (matrix_t *K, int n, int m, double val){
  K->_data[n+K->_n*m] += val;
}

double get_from_matrix (matrix_t *K, int n, int m){
  return K->_data[n+K->_n*m];
}

bool solve_linear_system(matrix_t *K, double *f, double *x)
{
  int N = K->_n, nrhs = 1, lda = N, ldb = N, info;
  int *ipiv = new int[N];
  for(int i = 0; i < N; i++) x[i] = f[i];
  dgesv(&N, &nrhs, K->_data, &lda, ipiv, x, &ldb, &info);
  delete [] ipiv;
  if(info == 0) return true;
  return false;
}
