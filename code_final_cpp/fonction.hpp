//
// Useful Fonction
//


// StiffnessMatrix computes of the stiffness matrix of a bar
void stiffnessMatrix( double x1 , double y1 , double x2 , double y2 , double ea , double k[4][4] );


// TrussSolver computes all the displacements [*x]
// and the reaction forces [*r] corresponding to the constraints
int trussSolver( int N, int B , int M, double *xy , int *ind ,
 double *f , int *nc , double *vc , double *ea , double *x , double *r );


// PostPro computes normal efforts in all bars
void postPro( int B , double *xy , int *ind , double *ea , double *x , double *n);


