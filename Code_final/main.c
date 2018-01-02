//  main.c
//  LGCIV2125_project 2017
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "fonction.h"

int main( void ) 
{
    const double L = 1.0 ;
    const double F = 1000;
    const double A = 1e-4;
    const double E = 210e9 ;

    const int N = 3 ;
    double xy[6] = {0,0,0,L,L,0} ; // Coordonees des noeuds
    
    const int B = 2 ;
    int ind[4] = {1,2,0,2} ;
    double ea[2] = {E*A,E*A} ;
    double f[6] = {0,0,0,0,0,-F} ;
    
    const int M = 4 ; // 4 contraintees
    int nc[4] = {0,0,1,1} ; // 2 contraintes sur le noeud 0 et 2 sur le noeud 1
    double vc[12] = {1,0,0,0,1,0,1,0,0,0,1,0} ; //vc[3*j] & vc[3*j+1] : direction of displacement at node nc[j]. And vc[3*j+2] : value of displacement.
    double x[4] ; // Déplacements aux noeuds
    double r[4] ; // Réactions aux appuis
    double n[2] = {0 ,0};
    
    int result = trussSolver(N, B ,M, xy , ind , f , nc , vc , ea , x , r ) ;
    
    printf("x[1] = %f \n x[2] = %f \n x[3] = %f \n x[4] = %f \n",x[1],x[2],x[3],x[4] );
    printf("r[1] = %f \n r[2] = %f \n r[3] = %f \n r[4] = %f \n",r[1],r[2],r[3],r[4] );

    if (result)
    {
    	postPro(B, xy, ind, ea, x, n);
    }
    // ligne 19 - 22 - 23 - 24 27 28 29 30 
}


