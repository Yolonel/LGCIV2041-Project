#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include "matrix.h"
#include "glut.hpp"

void stiffnessMatrix (double x1, double y1, double x2, double y2, double ea , double k[4][4])
{
  double t = atan2(y2-y1,x2-x1);
  double c = cos(t);
  double s = sin(t);
  double l = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
  double f = ea/l;
  k[0][0] = k[2][2] = c*c*f;
  k[1][1] = k[3][3] = s*s*f;
  k[0][2] = k[2][0] = -c*c*f;
  k[3][1] = k[1][3] = -s*s*f;
  k[0][1] = k[1][0] = k[3][2] = k[2][3] = s*c*f;
  k[3][0] = k[0][3] = k[2][1] = k[1][2] = -s*c*f;
}

void postPro (int B, double *xy, int *ind, 
	     double *ea, double *x, double *n) {
  double stiff [4][4];
  for (int i=0;i<B;i++){
    int n1 = ind[2*i];   
    int n2 = ind[2*i+1];
    int indx[4] = {2*n1,2*n1+1,2*n2,2*n2+1};
    stiffnessMatrix ( xy[indx[0]], xy[indx[1]], 
                      xy[indx[2]], xy[indx[3]], ea[i] , stiff );
    double f [4] = {0,0,0,0};
    for (int j=0;j<4;j++) {
      for (int k=0;k<4;k++) {
	f[j] += stiff[j][k] * x[indx[k]];  
      }
    }
    double t = atan2( xy[indx[1]]- xy[indx[3]], xy[indx[0]]- xy[indx[2]]);
    n[i] = f[0] * cos(t) + f[1] * sin(t);
    double l2 = sqrt((xy[indx[1]]+x[indx[1]]- xy[indx[3]]-x[indx[3]])*
		     (xy[indx[1]]+x[indx[1]]- xy[indx[3]]-x[indx[3]])+
		     (xy[indx[0]]+x[indx[0]]- xy[indx[2]]-x[indx[2]])*
		     (xy[indx[0]]+x[indx[0]]- xy[indx[2]]-x[indx[2]]));
    double l1 = sqrt((xy[indx[1]]- xy[indx[3]])*(xy[indx[1]]- xy[indx[3]])+
		     (xy[indx[0]]- xy[indx[2]])*(xy[indx[0]]- xy[indx[2]]));
    double dl = l2 - l1;
    double eps = dl / l1;
    double N = eps*ea[i];
    
    printf("n[%d] = %g N = %g\n",i,n[i],N);
  }
}

int trussSolver (int N, int B, int M, double *xy, int *ind, 
                 double *f, int *nc, double *vc, double *ea,  
                 double *x, double *r) {
  double stiff [4][4];
  matrix_t *K = create_matrix (2*N+M, 2*N+M);
  double *F = (double*) malloc ((2*N+M)*sizeof(double));
  double *X = (double*) malloc ((2*N+M)*sizeof(double));
  for (int i=0;i<B;i++){
    int n1 = ind[2*i];   
    int n2 = ind[2*i+1];
    int indx[4] = {2*n1,2*n1+1,2*n2,2*n2+1};
    stiffnessMatrix ( xy[indx[0]], xy[indx[1]], 
                      xy[indx[2]], xy[indx[3]], ea[i] , stiff );
    for (int j=0;j<4;j++) for (int k=0;k<4;k++)
      add_to_matrix(K,indx[j],indx[k],stiff[j][k]);
  }
  for (int i=0;i<2*N;i++) F[i] = f[i];
  for (int i=0;i<M;i++) {
    int indx[2] = {2*nc[i], 2*nc[i] +1};
    double v[3] = {vc[3*i], vc[3*i+1], vc[3*i+2]};
    double norm = sqrt(v[0]*v[0] + v[1]*v[1]);
    add_to_matrix(K,2*N+i,indx[0], v[0]/norm);
    add_to_matrix(K,2*N+i,indx[1], v[1]/norm);
    add_to_matrix(K,indx[0],2*N+i, v[0]/norm);
    add_to_matrix(K,indx[1],2*N+i, v[1]/norm);
    F[2*N+i] = v[2];
  }
  bool result = solve_linear_system (K, F, X);
  if (result == false){
    printf("ERROR : the truss is not stable\n");
  }
  for (int i=0;i<2*N;i++) x[i] = X[i];
  for (int i=0;i<M;i++) {r[i] = X[2*N+i];printf("r(%d) = %g\n",i,r[i]);}

  delete_matrix(K);
  free(F);
  free(X);
  return result;
}

