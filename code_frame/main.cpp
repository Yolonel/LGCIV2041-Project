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
	// Initialisation of the constants
    const double L = 1.0 ;
    const double F = 1000;
    const double A = 1e-4;
    const double E = 210e9;
    const double I = 1e-8;

    const int N = 3 ; // Number of nodes 
    double xy[2*N] = {0,0,0,L,L,0} ; // Coordinates of nodes
   
    const int B = 2 ; // Number of bars
    int ind[2*B] = {1,2,0,2} ; // Index of nodes associated to the bars
    double ea[B] = {E*A,E*A} ; 
    double ei[B] = {E*I,E*I} ;
    double f[3*N] = {0,0,0,0,0,0,0,-F,0} ; // Forces/torques applied to nodes
    
    const int M = 4 ; // Number of constraints
    int nc[M] = {0,0,1,1} ; // Nodes on which constraints are applied
    double vc[3*M] = {1,0,0,0,1,0,1,0,0,0,1,0} ; // vc[3*j] & vc[3*j+1] : direction of displacement at node nc[j]. And vc[3*j+2] : value of displacement.
    
    double x[3*N] ; // Displacements and orientation at nodes
     
    // Resolution of the truss
    int result = trussSolver(N, B ,M, xy , ind , f , nc , vc , ea, ei , x) ;
    
    printf("u11 = %3.9f \tu12 = %3.9f \ttheta1 = %3.9f \n",x[0],x[1],x[2]);
    printf("u21 = %3.9f \tu22 = %3.9f \ttheta2 = %3.9f \n",x[3],x[4],x[5]);
    printf("u31 = %3.9f \tu32 = %3.9f \ttheta3 = %3.9f \n",x[6],x[7],x[8]);
    
    return result;
}
