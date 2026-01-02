/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <string.h>

void eig_poisson1D(double* eigval, int *la){
  // TODO: Compute all eigenvalues for the 1D Poisson operator
  int n = *la;
  for (int k = 1; k <= n; k++){
    eigval[k-1] = 2.0 - 2.0 * cos(k * M_PI / (n + 1.0));
  }
}

double eigmax_poisson1D(int *la){
  // TODO: Compute and return the maximum eigenvalue for the 1D Poisson operator
  int n = *la;
  return 2.0 - 2.0 * cos(n * M_PI / (n + 1.0));
}

double eigmin_poisson1D(int *la){
  // TODO: Compute and return the minimum eigenvalue for the 1D Poisson operator
  int n = *la;
  return 2.0 - 2.0 * cos(M_PI / (n + 1.0));
}

double richardson_alpha_opt(int *la){
  // TODO: Compute alpha_opt
  double lmin = eigmin_poisson1D(la);
  double lmax = eigmax_poisson1D(la);
  return 2.0 / (lmin + lmax);
}

/**
 * Solve linear system Ax=b using Richardson iteration with fixed relaxation parameter alpha.
 * The iteration is: x^(k+1) = x^(k) + alpha*(b - A*x^(k))
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  // TODO: Implement Richardson iteration
  // 1. Compute residual r = b - A*x (use dgbmv for matrix-vector product)
  // 2. Update x = x + alpha*r (use daxpy)
  // 3. Check convergence: ||r||_2 < tol (use dnrm2)
  // 4. Store residual norm in resvec and repeat
  int n = *la;
  int ldab = *lab;
  double alpha = *alpha_rich;

  double *Ax = (double*) malloc(n * sizeof(double));
  double *r  = (double*) malloc(n * sizeof(double));

  double normb = cblas_dnrm2(n, RHS, 1);
  if (normb == 0.0) normb = 1.0;

  int k = 0;
  double relres = 1.0;

  while (k < *maxit && relres > *tol){

    cblas_dgbmv(CblasColMajor, CblasNoTrans,
                n, n, *kl, *ku,
                1.0, AB, ldab,
                X, 1, 0.0, Ax, 1);

    for (int i = 0; i < n; i++)
      r[i] = RHS[i] - Ax[i];

    relres = cblas_dnrm2(n, r, 1) / normb;
    resvec[k] = relres;

    cblas_daxpy(n, alpha, r, 1, X, 1);

    k++;
  }

  *nbite = k;

  free(Ax);
  free(r);
}

/**
 * Extract MB for Jacobi method from tridiagonal matrix.
 * Such as the Jacobi iterative process is: x^(k+1) = x^(k) + D^(-1)*(b - A*x^(k))
 */
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  // TODO: Extract diagonal elements from AB and store in MB
  // MB should contain only the diagonal of A
  int n = *la;
  int ldab = *lab;
  int ku_v = *ku;

  memset(MB, 0, ldab * n * sizeof(double));

  for (int j = 0; j < n; j++){
    MB[indexABCol(ku_v, j, lab)] =
      AB[indexABCol(ku_v, j, lab)];
  }
}

/**
 * Extract MB for Gauss-Seidel method from tridiagonal matrix.
 * Such as the Gauss-Seidel iterative process is: x^(k+1) = x^(k) + (D-E)^(-1)*(b - A*x^(k))
 */
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  // TODO: Extract diagonal and lower diagonal from AB
  // MB should contain the lower triangular part (including diagonal) of A
  int n = *la;
  int ldab = *lab;
  int ku_v = *ku;

  memset(MB, 0, ldab * n * sizeof(double));

  for (int j = 0; j < n; j++){
    MB[indexABCol(ku_v, j, lab)] =
      AB[indexABCol(ku_v, j, lab)];

    if (j < n-1){
      MB[indexABCol(ku_v+1, j, lab)] =
        AB[indexABCol(ku_v+1, j, lab)];
    }
  }
}

/**
 * Solve linear system Ax=b using preconditioned Richardson iteration.
 * The iteration is: x^(k+1) = x^(k) + M^(-1)*(b - A*x^(k))
 * where M is either D for Jacobi or (D-E) for Gauss-Seidel.
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  // TODO: Implement Richardson iterative method
  int n = *la;
  int ldab = *lab;
  int ku_v = *ku;

  double *Ax = (double*) malloc(n * sizeof(double));
  double *r  = (double*) malloc(n * sizeof(double));
  double *z  = (double*) malloc(n * sizeof(double));

  double normb = cblas_dnrm2(n, RHS, 1);
  if (normb == 0.0) normb = 1.0;

  int k = 0;
  double relres = 1.0;

  while (k < *maxit && relres > *tol){

    cblas_dgbmv(CblasColMajor, CblasNoTrans,
                n, n, *kl, *ku,
                1.0, AB, ldab,
                X, 1, 0.0, Ax, 1);

    for (int i = 0; i < n; i++)
      r[i] = RHS[i] - Ax[i];

    relres = cblas_dnrm2(n, r, 1) / normb;
    resvec[k] = relres;

    /* detect Jacobi or GS */
    int is_jacobi = 1;
    for (int j = 0; j < n-1; j++){
      if (MB[indexABCol(ku_v+1, j, lab)] != 0.0){
        is_jacobi = 0;
        break;
      }
    }

    if (is_jacobi){
      for (int i = 0; i < n; i++)
        z[i] = r[i] / MB[indexABCol(ku_v, i, lab)];
    }
    else{
      for (int i = 0; i < n; i++){
        double sum = r[i];
        if (i > 0)
          sum -= MB[indexABCol(ku_v+1, i-1, lab)] * z[i-1];
        z[i] = sum / MB[indexABCol(ku_v, i, lab)];
      }
    }

    cblas_daxpy(n, 1.0, z, 1, X, 1);

    k++;
  }

  *nbite = k;

  free(Ax);
  free(r);
  free(z);
}

