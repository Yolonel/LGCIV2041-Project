#include "matrix_t.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

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
	//int m = K->m;

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
	//printf("K->mData[%d][%d] = %f et val = %f\n",n,m,K->mData[n][m],val );
}

int solve_linear_system(matrix_t *K, double *f, double *x)
{
	int n = K->n;
	int m = K->m;
	int i,j,k;
	double p;

	printf("Initialiser A[%u][%u]\n",n,m);
	double ** A = (double**) malloc(sizeof(double*)*n);
	for (i = 0; i < n; ++i)
	{
		A[i] = (double*) calloc(m,sizeof(double));
	}
	printf("oui oui\n");
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
		{
			A[i][j] = K->mData[i][j];
		}
	}
	//printf("A[2][7] = %f , A[2][4] = %f , A[3][4] = %f \n", A[2][7],A[2][4],A[3][4]);

	printf("Commencer solver\n");
	for(k=0;k<n;k++)
	{
		if (A[k][k]==0)
		{
			printf("\n\n * Un pivot nul ! => methode de Jordan non applicable\n\n");
			return 0;
		}                 

		p=A[k][k];
		printf("Voici p : %f\n", p);

		//normalisation 
		printf("normalisation\n");              
		for (j=k;j<n;j++) A[k][j]=A[k][j]/p;
		f[k]=f[k]/p;

		//réduction    
		printf("reduction\n");            
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
	for ( int i=0; i<n; i++) x[i] = f[i] ;
	return 1;

}

void print_m(matrix_t *K)
{
	int n = K->n;
	int m = K->m;
	
	for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				printf("%5.1f\t",K->mData[i][j] );
			}
			printf("\n");
		}	


}

void print(double *K, int n)
{
	
for (int i = 0; i < n; ++i)
{
	printf("%5.1f\n",K[i] );
}

}

// void jordan(float a[19][19],float b[19],int n)

// {
// 	float p;int i,j,k;

// 	for(k=0;k<n;k++)
// 	{
// 		if (a[k][k]==0)
// 		{
// 			printf("\n\n * Un pivot nul ! => methode de Jordan non applicable\n\n");
// 			system("PAUSE");main();
// 		}                 

// 		p=a[k][k]; 

// 		//normalisation               
// 		for (j=k;j<n;j++) a[k][j]=a[k][j]/p;
// 		b[k]=b[k]/p;

// 		//réduction                
// 		for(i=0;i<n;i++)
// 		{
// 			 if (i!=k)
// 			 {
// 				  p=a[i][k];                
// 				  for (j=k;j<n;j++) a[i][j]=a[i][j]-p*a[k][j];
// 				  b[i]=b[i]-p*b[k];
// 			 }
// 		}
// 	}
// 	zero(a,b,n);
// 	printf("\n-------------- Jordan -------------\n");
// 	printf("\n * La matrice reduite :");
// 	aff_syst(a,b,n);
// 	printf("\n * La resolution donne :\n\n");
// 	for(i=0;i<n;i++) printf(" X_%d =  %f  ;\n",i+1,b[i]);
// 	printf("\n");
// }

