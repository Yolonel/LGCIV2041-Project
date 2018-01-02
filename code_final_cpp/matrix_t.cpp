#include "matrix_t.hpp"
#include <math.h>
#include <cstdio>

//
// Matrix_t Structure
//


// Matrix_t contains :
//	int n 		 : number of rows
//  int m 		 : number of columns
//  double* data : Values of the matrix_t in a vector 1x(n*m) 
struct matrix_t {
  int  n, m;
  double * data;

  matrix_t (int n, int m) :  n(n),  m(m)
  {
    data = new double [n*m];
    for (int i=0;i<n*m;i++) data[i] = 0.0;
  }
};

// Create_matrix return a instance of matrix_t
// @IN  :
//	    n : number of rows
//		m : number of columns
// @OUT :
//		K : matrix
matrix_t * create_matrix (int n, int m){
  return new matrix_t (n,m);
}

// Delete_matrix delete a instance of matrix_t
// @IN  :
//	    K : matrix
// @OUT :
//		/
void delete_matrix (matrix_t *K) {
  delete K;
}

// Add_to_matrix add the value val at the position (n,m)
// @IN  :
//	    K : matrix
//		n : number of rows
//		m : number of columns
//		val: values to add
// @OUT :
//		/
void add_to_matrix (matrix_t *K, int n, int m, double val){
  K->data[n+K->n*m] += val;
}


// Get_form_matrix return the value val at the posotion (n,m)
// @IN  :
//	    K : matrix
//		n : number of rows
//		m : number of columns
// @OUT :
//		val: values to add
double get_from_matrix (matrix_t *K, int n, int m){
  return K->data[n+K->n*m];
}

// Solve_linear_system solve the equation Kx=f by using
//					   the LU decomposition
// @IN  :
//	    K : matrix
//		f : matrix
//		x : matrix to be fill
// @OUT :
//		result: 1 if all goes right, 0 otherwise
bool solve_linear_system(matrix_t *K, double *f, double *x)
{
  	int n = K->n;
	int m = K->m;
	int i,j,k;
	double p,r,Atemp,ftemp;

	double **A;
	A = new double*[n];
	for (i = 0; i < n; ++i)
	{
		A[i] = new double[m];
	}

	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
		{
			A[i][j] = get_from_matrix (K, i, j);
		}
	}
	
	// With pivoting
	for(k=0;k<n;k++)
	{	
		r = fabs(A[k][k]);
		m = k;    
		for(i=k+1;i<n;i++) // Finding largest element
		{
			Atemp = fabs(A[i][k]);	// Swap pivot if necessary
			if(Atemp>r)
			{
				r = Atemp;
				m = i;
			}
		}            
		if(m!=k)	// Swap rows if pivot was swapped
		{
			for(i=k;i<n;i++)
			{
				Atemp = A[k][i];
				A[k][i] = A[m][i];
				A[m][i] = Atemp; 
			}
			ftemp = f[k];
			f[k] = f[m];
			f[m] = ftemp;
		}

		p=A[k][k];

		if (p==0)
		{
			printf("A pivot is null\n");
			return 0;
		}

		//normalisation 
		for (j=k;j<n;j++) A[k][j]=A[k][j]/p;
		f[k]=f[k]/p;

		//réduction    
		for(i=0;i<n;i++)
		{
			 if (i!=k)
			 {
				  p=A[i][k];                
				  for (j=k;j<n;j++) A[i][j]=A[i][j]-p*A[k][j];	
				  f[i]=f[i]-p*f[k];
			 }
		}
	}
	for(i=0;i<n;i++)
	{
		x[i] = f[i];
	}
	return 1;
}
