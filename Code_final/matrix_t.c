#include "matrix_t.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

matrix_t *create_matrix(int n, int m)
{
	assert(n > 0);
    assert(m > 0);
	matrix_t *mat = malloc(sizeof(matrix_t));
	mat->n = n;
	mat->m = m;

	mat->mData = (double**) malloc(sizeof(double*)*n);

	for (int i = 0; i < n; ++i)
	{
		mat->mData[i] = (double*) calloc(m,sizeof(double));
	}

	return mat;
}

void delete_matrix(matrix_t *K)
{
	int n = K->n;

	for (int i = 0; i < n; ++i)
	{
		free(K->mData[i]);
	}

	free(K->mData);

	free(K);
}

void add_to_matrix(matrix_t *K, int n, int m, double val)
{
	K->mData[n][m] = val;
}

int solve_linear_system(matrix_t *K, double *f, double *x)
{
	int n = K->n;
	int m = K->m;
	int i,j,k;
	double p,r,Atemp,ftemp;

	double ** A = (double**) malloc(sizeof(double*)*n);
	for (i = 0; i < n; ++i)
	{
		A[i] = (double*) calloc(m,sizeof(double));
	}

	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
		{
			A[i][j] = K->mData[i][j];
		}
	}
	//printf("A[2][7] = %f , A[2][4] = %f , A[3][4] = %f \n", A[2][7],A[2][4],A[3][4]);

	
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
	return 1;

}
