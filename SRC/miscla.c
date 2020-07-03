//============================================================================
// Routines for computing eigenvalues of a symmetric tridiagonal matrix.
// They are wrappers of the LAPACK routine DSTEV() or sstev_()
//============================================================================
#include "pevsl_int.h"
/**
 * @file miscla.c
 * @brief Miscellaneous la functions
 * */

/**-----------------------------------------------------------------------
 *  @brief compute all eigenvalues and eigenvectors of a symmetric tridiagonal
 *  matrix
 *  @param[in] n                The  dimension of the symmetric tridiagonal  matrix
 *  @param[in, out] pevsl pEVSL data struct (time is updated)
 *  @param[in] diag   Define the symmetric tridiagonal  matrix:  the
 *          diagonal elements are diag[0,...,n-1]
 *  @param[in] sdiag Subdiagonal elements
 *  @param[out] eigVal The output vector of length n containing all eigenvalues
 *          in ascending order
 *  @param[out] eigVec The output n-by-n matrix with columns as eigenvectors,
 *          in the order as elements in eigVal. If NULL, then no eigenvector
 *          will be computed
 *  @return The flag returned by the
 *  LAPACK routine DSTEV() (if double  is double) or stev_() (if double
 *  is float)
 * --------------------------------------------------------------------- */

int SymmTridEig(pevsl_Data *pevsl,
                double *eigVal, double *eigVec, int n,
                const double *diag, const double *sdiag) {

  double tms = pEVSL_Wtime();
  // compute eigenvalues and eigenvectors or eigvalues only
  char jobz = eigVec ? 'V' : 'N';
  int nn = n;
  int ldz = n;
  int info;  // output flag
  // copy diagonal and subdiagonal elements to alp and bet
  double *alp = eigVal;
  double *bet;
  PEVSL_MALLOC(bet, n-1, double);
  memcpy(alp, diag, n*sizeof(double));
  memcpy(bet, sdiag, (n-1)*sizeof(double));
  // allocate storage for computation
  double *sv = eigVec;
  double *work = NULL;
  if (jobz == 'V') {
    PEVSL_MALLOC(work, 2*n-2, double);
  }
  DSTEV(&jobz, &nn, alp, bet, sv, &ldz, work, &info);
  // free memory
  PEVSL_FREE(bet);
  if (work) {
    PEVSL_FREE(work);
  }
  if (info) {
    printf("DSTEV ERROR: INFO %d\n", info);
  }

  double tme = pEVSL_Wtime();
  pevsl->stats->t_eig += tme - tms;

  // return info
  return info;
}

/**-----------------------------------------------------------------------
 *  @brief compute  eigenvalues and  eigenvectors of  a symmetric  tridiagonal
 *  matrix in a slice
 *  @param[in] pevsl pEVSL data struct
 *  @param[in] n The  dimension of  the  symmetric tridiagonal  matrix
 *  @param[in] diag Diagonal elements
 *  @param[in] sdiag Sub-diagonal elements
 *  @param[out] eigVal Total number of eigenvalues found.
 *  @param[out] eigVec The first M elements contain teh selected eigenvalues in
 *  ascending oredr
 *  @param[in] vl If range='V', The lower bound of the interval to be searched
 *  for eigen values.
 *  @param[in] vu If  range='V', the upper bound of the interval to be searched
 *  for eigenvalues.
 *  @param[in] nevO If range='I', the index of the smallest eigen value to be
 *  returned.
 *
 *  This
 *  routine  computes selected  eigenvalues/vectors as  specified by  a
 *  range of values. This is a wrapper to the LAPACK routine DSTEMR().
 * ----------------------------------------------------------------------- */
int SymmTridEigS(pevsl_Data *pevsl,
                 double *eigVal, double *eigVec, int n, double vl, double vu,
                 int *nevO, const double *diag, const double *sdiag) {

  double tms = pEVSL_Wtime();
  char jobz = 'V';  // compute eigenvalues and eigenvectors
  char range = 'V'; // compute eigenvalues in an interval

  // this does not use mwlapack for mex files
  int info;
  //int idum = 0;
  //-------------------- isuppz not needed
  int *isuppz;
  PEVSL_MALLOC(isuppz, 2*n, int);
  //-------------------- real work array
  double *work;
  int lwork = 18*n;
  PEVSL_MALLOC(work, lwork, double);
  //-------------------- int work array
  int *iwork;
  int liwork = 10*n;
  PEVSL_CALLOC(iwork, liwork, int);
  //-------------------- copy diagonal + subdiagonal elements
  //                     to alp and bet
  double *alp;
  double *bet;
  PEVSL_MALLOC(bet, n, double);
  PEVSL_MALLOC(alp, n, double);
  //
  memcpy(alp, diag, n*sizeof(double));
  if (n > 1) {
    memcpy(bet, sdiag, (n-1)*sizeof(double));
  }

  //-------------------- allocate storage for computation
  /*logical tryrac = 1; */
  int tryrac = 1;
  double t0 = vl;
  double t1 = vu;

  DSTEMR(&jobz, &range, &n, alp, bet, &t0, &t1, NULL, NULL, nevO,
         eigVal, eigVec, &n, &n, isuppz, &tryrac, work, &lwork,
         iwork, &liwork, &info);

  if (info) {
    printf("dstemr_ error %d\n", info);
  }

  //-------------------- free memory
  PEVSL_FREE(bet);
  PEVSL_FREE(alp);
  PEVSL_FREE(work);
  PEVSL_FREE(iwork);
  PEVSL_FREE(isuppz);

  double tme = pEVSL_Wtime();
  pevsl->stats->t_eig += tme - tms;

  return info;
}


/**- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     @brief interface to   LAPACK SYMMETRIC EIGEN-SOLVER
 *     @param[in, out] pevsl pEVSL data struct (time is updated)
 *     @param[in] n Size of problem
 *     @param[in] A Matrix
 *     @param[in] lda Leading dimension
 *     @param[out] Q Eigenvectors
 *     @param[in] ldq Leading dimension q
 *     @param[out] lam Eigenvalues
 *- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
void SymEigenSolver(pevsl_Data *pevsl,
                    int n, double *A, int lda, double *Q, int ldq, double *lam) {

  double tms = pEVSL_Wtime();
  /* compute eigenvalues/vectors of A that n x n, symmetric
   * eigenvalues saved in lam: the eigenvalues in ascending order
   * eigenvectors saved in Q */
  char jobz='V';/* want eigenvectors */
  char uplo='U';/* use upper triangular part of the matrix */
  /*   copy A to Q */
  int i,j;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      Q[j+i*ldq] = A[j+i*lda];
    }
  }
  /*   first call: workspace query */
  double work_size;
  int lwork = -1, info;
  DSYEV(&jobz, &uplo, &n, Q, &ldq, lam, &work_size, &lwork, &info);
  if (info != 0) {
    fprintf(stdout, "DSYEV error [query call]: %d\n", info);
    exit(0);
  }
  /*   second call: do the computation */
  lwork = (int) work_size;
  double *work;
  PEVSL_MALLOC(work, lwork, double);
  DSYEV(&jobz, &uplo, &n, Q, &ldq, lam, work, &lwork, &info);
  if (info != 0) {
    fprintf(stdout, "DSYEV error [comp call]: %d\n", info);
    exit(0);
  }
  PEVSL_FREE(work);

  double tme = pEVSL_Wtime();
  pevsl->stats->t_eig += tme - tms;
}


/**
 * @brief Classical GS reortho with Daniel, Gragg, Kaufman, Stewart test
 * @param[in, out] pevsl pEVSL data struct (time is updated)
 * @param[in] k Number of vectors
 * @param[in] i_max Number iterations
 * @param[in] Q matrix
 * @param[out] nrmv norm of v
 * @param[in, out] v Input
 * @param[out] w Work vector
 **/
void CGS_DGKS(pevsl_Data *pevsl,
              int k, int i_max, pevsl_Parvecs *Q, pevsl_Parvec *v,
              double *nrmv, double *w) {


  double tms = pEVSL_Wtime();

  double eta = 1.0 / sqrt(2.0);
  double old_nrm, new_nrm;
  int i;

  pEVSL_ParvecNrm2(v, &old_nrm);

  for (i=0; i<i_max; i++) {
    pEVSL_ParvecsGemtvWithWspace(1.0, Q, k, v, 0.0, w, w+k);
    pEVSL_ParvecsGemv(-1.0, Q, k, w, 1.0, v);
    pEVSL_ParvecNrm2(v, &new_nrm);
    if (new_nrm > eta * old_nrm) {
      break;
    }
    old_nrm = new_nrm;
  }

  if (nrmv) {
    *nrmv = new_nrm;
  }

  double tme = pEVSL_Wtime();
  pevsl->stats->t_reorth += tme - tms;
}


/** JS 01/03/19 for complex cases
 * @brief Classical GS reortho with Daniel, Gragg, Kaufman, Stewart test
 * @param[in, out] pevsl pEVSL data struct (time is updated)
 * @param[in] k Number of vectors
 * @param[in] i_max Number iterations
 * @param[in] Q matrix
 * @param[out] nrmv norm of v
 * @param[in, out] v Input
 * @param[out] w Work vector
 **/
void CGS_ZDGKS(pevsl_Data *pevsl, int k, int i_max, 
              pevsl_Parvecs *Qr, pevsl_Parvecs *Qi, 
              pevsl_Parvec  *vr, pevsl_Parvec  *vi,
              double *nrmv, double *wr, double *wi) {


  double tms = pEVSL_Wtime();

  double *w0, *w1;
  double eta = 1.0 / sqrt(2.0);
  double old_nrm, new_nrm;
  int i;

  int one = 1; 
  double done = 1.0, dnone = -1.0;
  w0 = wr; 
  w1 = wi;

  
  pEVSL_ParvecZNrm2(vr, vi, &old_nrm);

  for (i=0; i<i_max; i++) {
    pEVSL_ParvecsGemtvWithWspace(1.0, Qr, k, vr, 0.0, wr, wr+k);
    pEVSL_ParvecsGemtvWithWspace(1.0, Qi, k, vi, 0.0, w0, w0+k);
    DAXPY(&k, &done, w0, &one, wr, &one);   
 
    pEVSL_ParvecsGemtvWithWspace(1.0, Qr, k, vi, 0.0, wi, wi+k);
    pEVSL_ParvecsGemtvWithWspace(1.0, Qi, k, vr, 0.0, w1, w1+k);
    DAXPY(&k, &dnone, w1, &one, wi, &one);   
    
    pEVSL_ParvecsGemv(-1.0, Qr, k, wr, 1.0, vr);
    pEVSL_ParvecsGemv( 1.0, Qi, k, wi, 1.0, vr);

    pEVSL_ParvecsGemv(-1.0, Qi, k, wr, 1.0, vi);
    pEVSL_ParvecsGemv(-1.0, Qr, k, wi, 1.0, vi);

    pEVSL_ParvecZNrm2(vr, vi, &new_nrm);
    if (new_nrm > eta * old_nrm) {
      break;
    }
    old_nrm = new_nrm;
  }

  if (nrmv) {
    *nrmv = new_nrm;
  }

  double tme = pEVSL_Wtime();
  pevsl->stats->t_reorth += tme - tms;
  
}


/**
 * @brief Modified GS reortho. No test. just do i_max times
 * used in generalized ev problems
 * vnew = v - (v, z_j)*v_j, for j=1,2,...
 **/
void CGS_DGKS2(pevsl_Data *pevsl,
               int k, int i_max, pevsl_Parvecs *V, pevsl_Parvecs *Z, pevsl_Parvec *v,
               double *w) {

  double tms = pEVSL_Wtime();

  int i;
  for (i=0; i<i_max; i++) {
    /* w = Z^T * v */
    pEVSL_ParvecsGemtvWithWspace(1.0, Z, k, v, 0.0, w, w+k);
    /* v = v - V *w */
    pEVSL_ParvecsGemv(-1.0, V, k, w, 1.0, v);
  }

  double tme = pEVSL_Wtime();
  pevsl->stats->t_reorth += tme - tms;
}

/** JS 01/03/19 
 * @brief Modified GS reortho. for complex vectors
 * No test. just do i_max times
 * used in generalized ev problems
 * vnew = v - (v, z_j)*v_j, for j=1,2,...
 **/
void CGS_ZDGKS2(pevsl_Data *pevsl, int k, int i_max, 
         pevsl_Parvecs *Vr, pevsl_Parvecs *Vi, pevsl_Parvecs *Zr, pevsl_Parvecs *Zi, 
         pevsl_Parvec  *vr, pevsl_Parvec  *vi, double *wr, double *wi) {

  double tms = pEVSL_Wtime();
  double *w0, *w1;

  int one = 1; 
  double done = 1.0, dnone = -1.0;
  w0 = wr + i_max*k; 
  w1 = wi + i_max*k;

  int i;
  for (i=0; i<i_max; i++) {
    /* w = Z^H * v */
    pEVSL_ParvecsGemtvWithWspace(1.0, Zr, k, vr, 0.0, wr, wr+k);
    pEVSL_ParvecsGemtvWithWspace(1.0, Zi, k, vi, 0.0, w0, w0+k);
    /* *wr = *w0 + *w1; */ 
    DAXPY(&k, &done, w0, &one, wr, &one); 

    pEVSL_ParvecsGemtvWithWspace(1.0, Zr, k, vi, 0.0, wi, wi+k);
    pEVSL_ParvecsGemtvWithWspace(1.0, Zi, k, vr, 0.0, w1, w1+k);
    /* *wi =  *wi - *w1; */
    DAXPY(&k, &dnone, w1, &one, wi, &one);  

    /* v = v - V *w no complex conjugate? */
    pEVSL_ParvecsGemv(-1.0, Vr, k, wr, 1.0, vr);
    pEVSL_ParvecsGemv( 1.0, Vi, k, wi, 1.0, vr);

    pEVSL_ParvecsGemv(-1.0, Vr, k, wi, 1.0, vi);
    pEVSL_ParvecsGemv(-1.0, Vi, k, wr, 1.0, vi);

  }


  double tme = pEVSL_Wtime();
  pevsl->stats->t_reorth += tme - tms;
  

}

