//
// Useful Fonction
//


// StiffnessMatrix computes of the stiffness matrix of a bar
void stiffnessMatrix( double x1 , double y1 , double x2 , double y2 , double ea, double ei, double k[6][6] );


// frameSolver computes all the displacements [*x]
int frameSolver( int N, int B , int M, double *xy , int *ind ,
 double *f , int *nc , double *vc , double *ea , double *ei, double *x);


// postPro computes the normal effort of each bars in vector n
void postPro (int B, double *xy, int *ind, 
	     double *ea, double *ei, double *x, double *n);



