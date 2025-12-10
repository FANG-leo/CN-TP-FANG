/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  // TODO: Fill AB with the tridiagonal Poisson operator
  int n = *la;
  int lab_v = *lab;
  int ku = *kv; if (ku <= 0) ku = 1;
  double h = 1.0/(n+1.0);
  double inv_h2 = 1.0/(h*h);
  int size = lab_v * n;
  for (int k = 0; k < size; ++k) AB[k] = 0.0;
  for (int j = 0; j < n; ++j) {
    int row_diag = ku + j - j; /* = ku */
    AB[row_diag + j*lab_v] = 2.0 * inv_h2;
    if (j+1 < n) AB[(ku + 1) + j*lab_v] = -1.0 * inv_h2;
    if (j-1 >= 0) AB[(ku - 1) + j*lab_v] = -1.0 * inv_h2;
  }
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  // TODO: Fill AB with the identity matrix
  // Only the main diagonal should have 1, all other entries are 0
  int n = *la;
  int lab_v = *lab;
  int ku = *kv;
  if (ku <= 0) ku = 1;

  int size = lab_v * n;
  for (int i = 0; i < size; ++i) AB[i] = 0.0;

  for (int j = 0; j < n; ++j) {
    AB[ku + j*lab_v] = 1.0;
  }
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  // TODO: Compute RHS vector
  int n = *la;
  double h = 1.0 / (n + 1.0);
  double inv_h2 = 1.0 / (h*h);

  for (int i = 0; i < n; ++i) RHS[i] = 0.0;
  RHS[0] += (*BC0) * inv_h2;
  RHS[n-1] += (*BC1) * inv_h2;
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  // TODO: Compute the exact analytical solution at each grid point
  // This depends on the source term f(x) used in set_dense_RHS_DBC_1D
  int n = *la;
  double T0 = *BC0, T1 = *BC1;
  for (int i = 0; i < n; ++i) EX_SOL[i] = T0 + X[i] * (T1 - T0);
}  

void set_grid_points_1D(double* x, int* la){
  // TODO: Generate uniformly spaced grid points in [0,1]
}

double relative_forward_error(double* x, double* y, int* la){
  // TODO: Compute the relative error using BLAS functions (dnrm2, daxpy or manual loop)
  return 0.0;
}

int indexABCol(int i, int j, int *lab){
  // TODO: Return the correct index formula for column-major band storage
  return 0;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  // TODO: Implement specialized LU factorization for tridiagonal matrices
  return *info;
}
