//  main.c
//  LGCIV2125_project 2017
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "glut.hpp"
#include "matrix.h"

int main( void ) 
{
    const double L = 1.0 ;
    const double F = 1000;
    const double A = 1e-4;
    const double E = 210e9 ;

    const int N = 3 ;
    double xy[6] = {0,0,0,L,L,0} ;
    
    const int B = 2 ;
    int ind[4] = {1,2,0,2} ;
    double ea[2] = {E*A,E*A} ;
    double f[6] = {0,0,0,0,0,-F} ;
    
    const int M = 4 ;
    int nc[4] = {0,0,1,1} ;
    double vc[12] = {1,0,0,0,1,0,1,0,0,0,1,0} ;
    double x[4] ;
    double r[4] ;
    
    int result = trussSolver(N, B ,M, xy , ind , f , nc , vc , ea , x , r ) ;
    return result ;
}


