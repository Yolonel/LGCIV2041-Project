#ifndef Fonction_h
#define Fonction_h

#include "matrix_t.h"

// Fonction données
void stiffnessMatrix( double x1 , double y1 , double x2 , double y2 , double ea , double k[4][4] );

int trussSolver( int N, int B , int M, double *xy , int *ind ,
 double *f , int *nc , double *vc , double *ea , double *x , double *r );

void postPro( int B , double *xy , int *ind , double *ea , double *x , double *n);




#endif 