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


// -------------------------------------------------------------
//       Résolution Par décomposition LU Crout
// -------------------------------------------------------------

int solve_linear_system(matrix_t *K, double *f, double *x)
{
	int n = K->n; 
	int p = K->m;
	int i,j,k,m;
	double s;

	printf("Initialiser A[%u][%u] L et U \n",n,p);
	double ** A = (double**) malloc(sizeof(double*)*n);
	double ** L = (double**) malloc(sizeof(double*)*n);
	double ** U = (double**) malloc(sizeof(double*)*n);
	double * y = (double*) malloc(sizeof(double)*n);
	for (i = 0; i < n; ++i)
	{
		A[i] = (double*) calloc(p,sizeof(double));
		L[i] = (double*) calloc(p,sizeof(double));
		U[i] = (double*) calloc(p,sizeof(double));
	}
	printf("oui oui\n");
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < p; ++j)
		{
			A[i][j] = K->mData[i][j];
		}
	}

	printf("Commencer solver\n");

	for (i=0;i<n;i++) 
	{
		for (j=0;j<p;j++) 
		{
			if(i==j) 
				U[i][j]=1; 
			else 
				U[i][j]=0;
				L[i][j]=0;
		}
	}


	for (m=0;m<n;m++)
	{    
		for (i=m;i<n;i++) 
		{
			s=0;    
			for (k=0;k<m;k++) s=s+L[i][k]*U[k][m];    
			L[i][m]=A[i][m]-s;
		}

		if (L[k][k]==0) 
		{
			printf("L[k][k] = %f\n", L[k][k]);
			printf("\n\n* Un mineur nul ! => methode de LU non applicable\n\n");
			return 0;
		}                 

		for (j=m+1;j<n;j++)   
		{
			s=0;    
			for (k=0;k<m;k++) s=s+L[m][k]*U[k][j];     
			U[m][j]=(A[m][j]-s)/L[m][m];
		}
	}

	// resolution
	for(i=0;i<n;i++)
	{
		s=0;                      
		for(j=0;j<i;j++) s=s+L[i][j]*y[j];                    
		y[i]=(f[i]-s)/L[i][i];
	}

	for(i=n-1;i>=0;i--)
	{
		s=0;                      
		for(j=i+1;j<n;j++)s=s+U[i][j]*x[j];                    
		x[i]=(y[i]-s)/U[i][i];
	}

	return 1;
}
	







// int solve_linear_system(matrix_t *K, double *f, double *x)
// {
// 	int n = K->n;
// 	int m = K->m;
// 	int i,j,k;
// 	double p;

// 	printf("Initialiser A[%u][%u]\n",n,m);
// 	double ** A = (double**) malloc(sizeof(double*)*n);
// 	for (i = 0; i < n; ++i)
// 	{
// 		A[i] = (double*) calloc(m,sizeof(double));
// 	}
// 	printf("oui oui\n");
// 	for (i = 0; i < n; ++i)
// 	{
// 		for (j = 0; j < m; ++j)
// 		{
// 			A[i][j] = K->mData[i][j];
// 		}
// 	}
// 	//printf("A[2][7] = %f , A[2][4] = %f , A[3][4] = %f \n", A[2][7],A[2][4],A[3][4]);

// 	printf("Commencer solver\n");
// 	for(k=0;k<n;k++)
// 	{
// 		if (A[k][k]==0)
// 		{
// 			printf("\n\n * Un pivot nul ! => methode de Jordan non applicable\n\n");
// 			return 0;
// 		}                 

// 		p=A[k][k];
// 		printf("Voici p : %f\n", p);

// 		//normalisation 
// 		printf("normalisation\n");              
// 		for (j=k;j<n;j++) A[k][j]=A[k][j]/p;
// 		f[k]=f[k]/p;

// 		//réduction    
// 		printf("reduction\n");            
// 		for(i=0;i<n;i++)
// 		{
// 			 if (i!=k)
// 			 {
// 				  p=A[i][k];                
// 				  for (j=k;j<n;j++) A[i][j]=A[i][j]-p*A[k][j];
// 				  f[i]=f[i]-p*f[k];
// 			 }
// 		}
// 	}
// 	for ( int i=0; i<n; i++) x[i] = f[i] ;
// 	return 1;

// }

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

// -------------------------------------------------------------
//       Résolution Par décomposition LU Crout
// -------------------------------------------------------------

// void LU_crout(float a[19][19],float b[19],int n)

// {
// float L[19][19],U[19][19],x[19],y[19],s;int i,j,k,m;

// for (i=0;i<n;i++) for (j=0;j<n;j++) 
// {
// if(i==j) U[i][j]=1; else U[i][j]=0;
// L[i][j]=0;
// }

// for (m=0;m<n;m++)
// {    
// for (i=m;i<n;i++) 
// {
// s=0;    
// for (k=0;k<m;k++) s=s+L[i][k]*U[k][m];    
// L[i][m]=a[i][m]-s;
// }

// if (L[k][k]==0) 
// {
// printf("\n\n* Un mineur nul ! => methode de LU non applicable\n\n");
// system("PAUSE");main();
// }                 

// for (j=m+1;j<n;j++)   
// {
// s=0;    
// for (k=0;k<m;k++) s=s+L[m][k]*U[k][j];     
// U[m][j]=(a[m][j]-s)/L[m][m];
// }
// }

// // resolution
// for(i=0;i<n;i++)
// {
// s=0;                      
// for(j=0;j<i;j++) s=s+L[i][j]*y[j];                    
// y[i]=(b[i]-s)/L[i][i];
// }

// for(i=n-1;i>=0;i--)
// {
// s=0;                      
// for(j=i+1;j<n;j++)s=s+U[i][j]*x[j];                    
// x[i]=(y[i]-s)/U[i][i];
// }

