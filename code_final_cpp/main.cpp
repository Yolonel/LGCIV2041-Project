//  main.cpp
//  LGCIV2041_NUMERICAL ANALYSIS OF CIVIL ENGINEERING STRUCTURES : project 2017
//  Author : Ange ISHIMWE
//		   : Geoffroy JACQUET
//	Noma   : 39351300
//		   : 45631300
//	
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "fonction.hpp"
#include "matrix_t.hpp"


// Truss solver
int main( void ) 
{
	// Initialisation of the constante
    const double L = 1.0 ;
    const double F = 1000;
    const double A = 1e-4;
    const double E = 210e9 ;

    // 
    const int N = 3 ;
    double xy[2*N] = {0,0,0,L,L,0} ;
    
    //
    const int B = 2 ;
    int ind[2*B] = {1,2,0,2} ;
    double ea[B] = {E*A,E*A} ;
    double f[2*N] = {0,0,0,0,0,-F} ;
    
    // 
    const int M = 4 ;
    int nc[M] = {0,0,1,1} ;
    double vc[3*M] = {1,0,0,0,1,0,1,0,0,0,1,0} ;
    double x[2*N] ;
    double r[M] ;
    double n[2] = {0 ,0};
    
    // 
    int result = trussSolver(N, B ,M, xy , ind , f , nc , vc , ea , x , r ) ;
    if (result)
    {
    	postPro(B, xy, ind, ea, x, n);
    }
    
    return result;
}
