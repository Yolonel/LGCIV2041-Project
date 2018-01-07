#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include "matrix_t.hpp"
#include "fonction.hpp"


// StiffnessMatrix computes of the stiffness matrix of a bar
void stiffnessMatrix (double x1, double y1, double x2, double y2, double ea, double ei, double k[6][6])
{
  double t = atan2(y2-y1,x2-x1);
  double c = cos(t);
  double s = sin(t);
  double c2 = c*c;
  double s2 = s*s;
  double cs = c*s;
  
  double l = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
  double f = ea/l;
  double f2 = ei/pow(l,3);
  
  double kf[4][4];
  
  kf[0][0] = kf[2][2] = 12*f2;
  kf[1][1] = kf[3][3] = 4*l*l*f2;
  kf[3][1] = kf[1][3] = 2*l*l*f2;
  kf[0][2] = kf[2][0] = -12*f2;
  kf[0][1] = kf[1][0] = kf[0][3] = kf[3][0] = 6*l*f2;
  kf[3][2] = kf[2][3] = kf[2][1] = kf[1][2] = -6*l*f2;
  
  k[0][0] = f*c2+kf[0][0]*s2;
  k[0][1] = k[1][0] = f*cs-kf[0][0]*cs;
  k[0][2] = -kf[0][1]*s;
  k[0][3] = -f*c2+kf[0][2]*s2;
  k[0][4] = k[1][3] = -f*cs-kf[0][2]*cs;
  k[0][5] = -kf[0][3]*s;
  
  //k[1][0] 
  k[1][1] = f*s2+kf[0][0]*c2;
  k[1][2] = kf[0][1]*c;
  //k[1][3]  
  k[1][4] = -f*s2+kf[0][2]*c2; 
  k[1][5] = kf[0][3]*c;
      
  k[2][0] = -kf[1][0]*s;
  k[2][1] = kf[1][0]*c;
  k[2][2] = kf[1][1];
  k[2][3] = -kf[1][2]*s;
  k[2][4] = kf[1][2]*c;
  k[2][5] = kf[1][3];
  
  k[3][0] = -f*c2+kf[2][0]*s2;
  k[3][1] = k[4][0] = -f*cs-kf[2][0]*cs;
  k[3][2] = -kf[2][1]*s;
  k[3][3] = f*c2+kf[2][2]*s2;
  k[3][4] = k[4][3] = f*cs-kf[2][2]*cs;
  k[3][5] = -kf[2][3]*s;
  
  //k[4][0] = 
  k[4][1] = -f*s2+kf[2][0]*c2;
  k[4][2] = kf[2][1]*c;
  //k[4][3] = 
  k[4][4] = f*s2+kf[2][2]*c2;
  k[4][5] = kf[2][3]*c;
  
  k[5][0] = -kf[3][0]*s;
  k[5][1] = kf[3][0]*c;
  k[5][2] = kf[3][1];
  k[5][3] = -kf[3][2]*s;
  k[5][4] = kf[3][2]*c;
  k[5][5] = kf[3][3];
}


int frameSolver (int N, int B, int M, double *xy, int *ind, 
                 double *f, int *nc, double *vc, double *ea,  
                 double *ei, double *x) {
  double stiff [6][6];
  matrix_t *K = create_matrix(3*N+M, 3*N+M);
  double *F = (double*) malloc ((3*N+M)*sizeof(double));
  double *X = (double*) malloc ((3*N+M)*sizeof(double));
  for (int i=0;i<B;i++){
    int n1 = ind[2*i];   // Index of starting node of bar i
    int n2 = ind[2*i+1]; // Index of ending node of bar i
    int indx[6] = {3*n1,3*n1+1,3*n1+2,3*n2,3*n2+1,3*n2+2};
    stiffnessMatrix ( xy[2*n1], xy[2*n1+1], xy[2*n2], xy[2*n2+1], ea[i] , ei[i], stiff );
    
    for (int j=0;j<6;j++) for (int k=0;k<6;k++)
      add_to_matrix(K,indx[j],indx[k],stiff[j][k]);
  }
  for (int i=0;i<3*N;i++) F[i] = f[i];
  for (int i=0;i<M;i++) {
    int indx[2] = {3*nc[i], 3*nc[i] +1};
    double v[3] = {vc[3*i], vc[3*i+1], vc[3*i+2]};
    double norm = sqrt(v[0]*v[0] + v[1]*v[1]);
    add_to_matrix(K,3*N+i,indx[0], v[0]/norm);
    add_to_matrix(K,3*N+i,indx[1], v[1]/norm);
    add_to_matrix(K,indx[0],3*N+i, v[0]/norm);
    add_to_matrix(K,indx[1],3*N+i, v[1]/norm);
    F[3*N+i] = v[2];
  }
  bool result = solve_linear_system (K, F, X);
  if (result == false){
    printf("ERROR : the truss is not stable\n");
  }
  for (int i=0;i<3*N;i++) x[i] = X[i];

  delete_matrix(K);
  free(F);
  free(X);
  return result;
}

// void postPro (int B, double *xy, int *ind, 
// 	     double *ea, double *ei, double *x, double *n) {
//   double stiff [6][6];
//   double stiffn [4][4];
//   for (int i=0;i<B;i++)
//   {
//     int n1 = ind[2*i];   
//     int n2 = ind[2*i+1];
//     int indx[4] = {3*n1,3*n1+1,3*n2,3*n2+1};
//     stiffnessMatrix ( xy[2*n1], xy[2*n1+1], 
//                       xy[2*n2], xy[2*n2+1], ea[i] , ei[i] , stiff );
//     for (int j=0;j<4;j++)
//     {
//     	for (int k=0;k<4;k++)
//     	{
//     		int jn = (j==2)*3 + (j==3)*4 + (j==0 || j==1)*j;
//     		int kn = (k==2)*3 + (k==3)*4 + (k==0 || k==1)*k;
//     		stiffn[j][k] = stiff[jn][kn];
//     	}
//     }
//     double f [4] = {0,0,0,0};
//     for (int j=0;j<4;j++) 
//     {
//       for (int k=0;k<4;k++) 
//       {
// 		f[j] += stiffn[j][k] * x[indx[k]];  
//       }
//     }
//     double t = atan2( xy[2*n1+1]- xy[2*n2+1], xy[2*n1]- xy[2*n2]);
//     n[i] = f[0] * cos(t) + f[1] * sin(t);
//   }
// }

void postPro (int B, double *xy, int *ind, 
	     double *ea, double *ei, double *x, double *n) {
  double stiff [6][6];
  for (int i=0;i<B;i++)
  {
    int n1 = ind[2*i];   
    int n2 = ind[2*i+1];
    int indx[6] = {3*n1,3*n1+1,3*n1+2,3*n2,3*n2+1,3*n2+2};
    stiffnessMatrix ( xy[2*n1], xy[2*n1+1], 
                      xy[2*n2], xy[2*n2+1], ea[i] , ei[i] , stiff );
    double f [6] = {0,0,0,0,0,0};
    for (int j=0;j<6;j++) 
    {
      for (int k=0;k<6;k++) 
      {
		f[j] += stiff[j][k] * x[indx[k]];  
      }
    }
    double t = atan2( xy[2*n1+1]- xy[2*n2+1], xy[2*n1]- xy[2*n2]);
    n[i] = f[0] * cos(t) + f[1] * sin(t);
  }
}



