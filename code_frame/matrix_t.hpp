#ifndef _MATRIX_T_H_
#define _MATRIX_T_H_
//
// Matrix_t Structure
//


// Matrix_t contains :
//	int n 		 : number of rows
//  int m 		 : number of columns
//  double* data : Values of the matrix_t in a vector 1x(n*m) 
struct matrix_t;


// Create_matrix return a instance of matrix_t
// @IN  :
//	    n : number of rows
//		m : number of columns
// @OUT :
//		K : matrix
matrix_t * create_matrix (int n, int m);


// Delete_matrix delete a instance of matrix_t
// @IN  :
//	    K : matrix
// @OUT :
//		/
void delete_matrix (matrix_t *K);


// Add_to_matrix add the value val at the position (n,m)
// @IN  :
//	    K : matrix
//		n : number of rows
//		m : number of columns
//		val: values to add
// @OUT :
//		/
void add_to_matrix (matrix_t *K, int n, int m, double val);


// Get_form_matrix return the value val at the posotion (n,m)
// @IN  :
//	    K : matrix
//		n : number of rows
//		m : number of columns
// @OUT :
//		val: values to add
double get_from_matrix (matrix_t *K, int n, int m);


// Solve_linear_system solve the equation Kx=f by using
//					   the LU decomposition
// @IN  :
//	    K : matrix
//		f : matrix
//		x : matrix to be fill
// @OUT :
//		result: 1 if all goes right, 0 otherwise
bool solve_linear_system (matrix_t *K, double *f, double *x);

#endif
