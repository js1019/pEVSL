#include "pevsl_int.h"
#include "pevsl_itsol.h"

/** @file pevsl_f90.c
 * FORTRAN interface of pEVSL:
 * 1. We use C type uintptr_t to save C-points and pass them to Fortran.
 *    This is an ``integer type capable of holding a value converted from
 *    a void pointer (i.e., void *) and then be converted back to that type
 *    with a value that compares equal to the original pointer.''.
 *    So, it is of 32 bits or 64 bits depending on the platform.
 *    On the Fortran side, there should exist a corresponding variable declared
 *    to communicate. Typically, it can be INTEGER*4 or INTEGER*8
 */

/** @brief Fortran interface for pEVSL_Start */
void PEVSL_FORT(pevsl_start)(MPI_Fint *Fcomm, uintptr_t *pevslf90) {

  pevsl_Data *pevsl;
  MPI_Comm comm = MPI_Comm_f2c(*Fcomm);

  pEVSL_Start(comm, &pevsl);

  /* cast the pointer */
  *pevslf90 = (uintptr_t) pevsl;
}

/** @brief Fortran interface for pEVSL_Finish */
void PEVSL_FORT(pevsl_finish)(uintptr_t *pevslf90) {

  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);

  pEVSL_Finish(pevsl);
}

/** @brief Fortran interface for pEVSL_SetAParcsr */
void PEVSL_FORT(pevsl_seta_parcsr)(uintptr_t *pevslf90, uintptr_t *Af90) {

  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);
  pevsl_Parcsr *A = (pevsl_Parcsr *) (*Af90);

  pEVSL_SetAParcsr(pevsl, A);
}

/** @brief Fortran interface for pEVSL_SetBParcsr */
void PEVSL_FORT(pevsl_setb_parcsr)(uintptr_t *pevslf90, uintptr_t *Bf90) {

  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);
  pevsl_Parcsr *B = (pevsl_Parcsr *) (*Bf90);

  pEVSL_SetBParcsr(pevsl, B);
}

/** @brief Fortran interface for pEVSL_SetProbSizes
 * @param[in] pevslf90  pevsl pointer
 * @param[in] N  global size of A
 * @param[in] n  local size of A
 * @param[in] nfirst  nfirst
 * @warning set nfirst < 0, if do not want to specify
 */
void PEVSL_FORT(pevsl_setprobsizes)(uintptr_t *pevslf90, int *N, int *n, int *nfirst) {

  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);

  pEVSL_SetProbSizes(pevsl, *N, *n, *nfirst);
}

/** @brief Fortran interface for pEVSL_SetAMatvec
 * @param[in] pevslf90  pevsl pointer
 * @param[in] func  function pointer
 * @param[in] data  associated data
 */
void PEVSL_FORT(pevsl_setamv)(uintptr_t *pevslf90, void *func, void *data) {

  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);

  pEVSL_SetAMatvec(pevsl, (MVFunc) func, data);
}


/** @brief Fortran interface for pEVSL_SetBMatvec
 * @param[in] pevslf90  pevsl pointer
 * @param[in] func  function pointer
 * @param[in] data  associated data
 */
void PEVSL_FORT(pevsl_setbmv)(uintptr_t *pevslf90, void *func, void *data) {

  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);

  pEVSL_SetBMatvec(pevsl, (MVFunc) func, data);
}

/** @brief Fortran interface for SetBsol */
void PEVSL_FORT(pevsl_setbsol)(uintptr_t *pevslf90, void *func, void *data) {

  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);

  pEVSL_SetBSol(pevsl, (SVFunc) func, data);
}


/**JS 01/20/19 for complex Hermitian Av
 * @brief Fortran interface for pEVSL_SetAMatvec
 * @param[in] pevslf90  pevsl pointer
 * @param[in] func  function pointer
 * @param[in] data  associated data
 */
void PEVSL_FORT(pevsl_setzamv)(uintptr_t *pevslf90, void *func, void *data) {

  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);

  pEVSL_SetZAMatvec(pevsl, (ZMVFunc) func, data);
}

/** @brief Fortran interface for pEVSL_SetBMatvec
 * @param[in] pevslf90  pevsl pointer
 * @param[in] func  function pointer
 * @param[in] data  associated data
 */
void PEVSL_FORT(pevsl_setzbmv)(uintptr_t *pevslf90, void *func, void *data) {

  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);

  pEVSL_SetZBMatvec(pevsl, (ZMVFunc) func, data);
}

/** @brief Fortran interface for SetBsol */
void PEVSL_FORT(pevsl_setzbsol)(uintptr_t *pevslf90, void *func, void *data) {

  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);

  pEVSL_SetZBSol(pevsl, (ZSVFunc) func, data);
}




/** @brief Fortran interface for SetStdEig */
void PEVSL_FORT(pevsl_set_stdeig)(uintptr_t *pevslf90) {

  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);

  pEVSL_SetStdEig(pevsl);
}

/** @brief Fortran interface for SetGenEig */
void PEVSL_FORT(pevsl_set_geneig)(uintptr_t *pevslf90) {

  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);

  pEVSL_SetGenEig(pevsl);
}


/** @brief Fortran interface for ParcsrCreate and ParcsrSetup
 * @param[in] nrow  number of rows    [global]
 * @param[in] ncol  number of columns [global]
 * @param[in] row_starts  row partitioning array [of size np + 1]
 * @param[in] col_starts  column partitioning array [of size np + 1]
 * @param[in] ia, ja, aa   local CSR matrix [NOTE: MUST be of C-index, starting with zero]
 * @param[in] Fcomm       MPI communicator
 * @param[out] matf90     matrix pointer
 * */
void PEVSL_FORT(pevsl_parcsrcreate)(int *nrow, int *ncol, int *row_starts, int *col_starts,
                                    int *ia, int *ja, double *aa, MPI_Fint *Fcomm,
                                    uintptr_t *matf90) {

  pevsl_Parcsr *mat;
  pevsl_Csr matloc;
  MPI_Comm comm = MPI_Comm_f2c(*Fcomm);

  PEVSL_MALLOC(mat, 1, pevsl_Parcsr);
  pEVSL_ParcsrCreate(*nrow, *ncol, row_starts, col_starts, mat, comm);

  /* matloc is the local CSR wrapper */
  matloc.nrows = mat->nrow_local;
  matloc.ncols = mat->ncol_global;
  matloc.ia = ia;
  matloc.ja = ja;
  matloc.a  = aa;

  /* setup parcsr with matloc */
  pEVSL_ParcsrSetup(&matloc, mat);

  *matf90 = (uintptr_t) mat;
}

/** @brief Fortran interface for ParcsrFree */
void PEVSL_FORT(pevsl_parcsr_free)(uintptr_t *Af90) {

  /* cast pointer */
  pevsl_Parcsr *A = (pevsl_Parcsr *) (*Af90);

  pEVSL_ParcsrFree(A);
}

/* @brief Fortran interface for ParcsrCreate and ParcsrMatvec0 */
void PEVSL_FORT(pevsl_parcsrmatvec)(double *x, double *y, uintptr_t *matf90) {
  /* cast pointer */
  pevsl_Parcsr *mat = (pevsl_Parcsr *) (*matf90);

  pEVSL_ParcsrMatvec0(x, y, (void *) mat);
}

/* @brief Fortran interface for B-Solve */
void PEVSL_FORT(pevsl_bsv)(uintptr_t *pevslf90, double *b, double *x) {

  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);

  pevsl->Bsol->func(b, x, pevsl->Bsol->data);
}




/** @brief Fortran interface for evsl_lanbounds
 * @param[in] pevslf90 pEVSL pointer
 * @param[in] mlan Krylov dimension
 * @param[in] nsteps number of steps
 * @param[in] tol stopping tol
 * @param[out] lmin lower bound
 * @param[out] lmax upper bound
 * */
void PEVSL_FORT(pevsl_lanbounds)(uintptr_t *pevslf90, int *mlan, int *nsteps, double *tol,
                                 double *lmin, double *lmax) {
  int N, n, nfirst;
  pevsl_Parvec vinit;

  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);

  N = pevsl->N;
  n = pevsl->n;
  nfirst = pevsl->nfirst;

  /*------------------- Create parallel vector: random initial guess */
  pEVSL_ParvecCreate(N, n, nfirst, pevsl->comm, &vinit);
  pEVSL_ParvecRand(&vinit);
  /*------------------- Lanczos Bounds */
  pEVSL_LanTrbounds(pevsl, *mlan, *nsteps, 1e-8, &vinit, 1, lmin, lmax, NULL);
  //pEVSL_LanTrbounds(pevsl, *mlan, *nsteps, 1e-8, &vinit, 1, lmin, lmax, stdout);

  pEVSL_ParvecFree(&vinit);
}


/** JS 01/20/19  
 * @brief Fortran interface for evsl_lanbounds
 * @param[in] pevslf90 pEVSL pointer
 * @param[in] mlan Krylov dimension
 * @param[in] nsteps number of steps
 * @param[in] tol stopping tol
 * @param[out] lmin lower bound
 * @param[out] lmax upper bound
 * */
void PEVSL_FORT(pevsl_zlanbounds)(uintptr_t *pevslf90, int *mlan, int *nsteps, double *tol,
                                 double *lmin, double *lmax) {
  int N, n, nfirst;
  pevsl_Parvec vrinit, viinit;

  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);

  N = pevsl->N;
  n = pevsl->n;
  nfirst = pevsl->nfirst;

  /*------------------- Create parallel vector: random initial guess */
  pEVSL_ParvecCreate(N, n, nfirst, pevsl->comm, &vrinit);
  pEVSL_ParvecRand_split(&vrinit);
  pEVSL_ParvecCreate(N, n, nfirst, pevsl->comm, &viinit);
  pEVSL_ParvecRand_split(&viinit);
  /*------------------- Lanczos Bounds */
  pEVSL_ZLanTrbounds(pevsl, *mlan, *nsteps, 1e-8, &vrinit, &viinit, 1, lmin, lmax, stdout);

  pEVSL_ParvecFree(&vrinit);
  pEVSL_ParvecFree(&viinit);
}


/** @brief Fortran interface for find_pol
 * @param[in] xintv  Interval of interest
 * @param[in] thresh_int Interior threshold
 * @param[in] thresh_ext Exterior threshold
 * @param[out] polf90  pointer of pol
 * @warning The pointer will be cast to uintptr_t
 *
 * uintptr_t Integer type capable of holding a value converted from
 * a void pointer and then be converted back to that type with a value
 * that compares equal to the original pointer
 */
void PEVSL_FORT(pevsl_findpol)(double *xintv, double *thresh_int,
                               double *thresh_ext, uintptr_t *polf90) {
  pevsl_Polparams *pol;
  PEVSL_MALLOC(pol, 1, pevsl_Polparams);
  pEVSL_SetPolDef(pol);
  pol->damping = 2;
  pol->thresh_int = *thresh_int;
  pol->thresh_ext = *thresh_ext;
  pol->max_deg  = 50000;
  pEVSL_FindPol(xintv, pol);

  fprintf(stdout, " polynomial deg %d, bar %e gam %e\n",
          pol->deg,pol->bar, pol->gam);

  *polf90 = (uintptr_t) pol;
}

/** @brief Fortran interface for free_pol */
void PEVSL_FORT(pevsl_freepol)(uintptr_t *polf90) {
  /* cast pointer */
  pevsl_Polparams *pol = (pevsl_Polparams *) (*polf90);
  pEVSL_FreePol(pol);
  PEVSL_FREE(pol);
}

/** @brief Fortran interface for ChebLanNr
 *  the results will be saved in the internal variables
 */
void PEVSL_FORT(pevsl_cheblannr)(uintptr_t *pevslf90, double *xintv, int *max_its, double *tol,
                                 uintptr_t *polf90) {
  int N, n, nfirst, nev2, ierr;
  double *lam, *res;
  pevsl_Parvecs *Y;
  FILE *fstats = stdout;
  pevsl_Parvec vinit;

  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);

  N = pevsl->N;
  n = pevsl->n;
  nfirst = pevsl->nfirst;
  /*-------------------- zero out stats */
  pEVSL_StatsReset(pevsl);
  /*------------------- Create parallel vector: random initial guess */
  pEVSL_ParvecCreate(N, n, nfirst, pevsl->comm, &vinit);
  pEVSL_ParvecRand(&vinit);
  /* cast pointer of pol*/
  pevsl_Polparams *pol = (pevsl_Polparams *) (*polf90);
  /* call ChebLanNr */
  ierr = pEVSL_ChebLanNr(pevsl, xintv, *max_its, *tol, &vinit, pol, &nev2, &lam, &Y, &res, fstats);

  if (ierr) {
    printf("ChebLanNr error %d\n", ierr);
  }

  pEVSL_ParvecFree(&vinit);
  /*--------------------- print stats */
  pEVSL_StatsPrint(pevsl, fstats);

  if (res) {
    PEVSL_FREE(res);
  }
  /* save pointers to the global variables */
  pevsl->nev_computed = nev2;
  pevsl->eval_computed = lam;
  pevsl->evec_computed = Y;
}


/** added by JS 020619 
 * @brief Fortran interface for ChebLanNr
 *  the results will be saved in the internal variables
 */
void PEVSL_FORT(pevsl_zcheblannr)(uintptr_t *pevslf90, double *xintv, int *max_its, double *tol,
                                 uintptr_t *polf90) {

  int N, n, nfirst, nev2, ierr;
  double *lam, *res;
  pevsl_Parvecs *Yr, *Yi;
  FILE *fstats = stdout;
  pevsl_Parvec vrinit, viinit;

  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);

  N = pevsl->N;
  n = pevsl->n;
  nfirst = pevsl->nfirst;
  /*-------------------- zero out stats */
  pEVSL_StatsReset(pevsl);
  /*------------------- Create parallel vector: random initial guess */
  pEVSL_ParvecCreate(N, n, nfirst, pevsl->comm, &vrinit);
  pEVSL_ParvecCreate(N, n, nfirst, pevsl->comm, &viinit);
  pEVSL_ParvecRand(&vrinit);
  pEVSL_ParvecRand(&viinit);
  /* cast pointer of pol*/
  pevsl_Polparams *pol = (pevsl_Polparams *) (*polf90);
  /* call ChebLanNr */
  ierr = pEVSL_ZChebLanNr(pevsl, xintv, *max_its, *tol, &vrinit, &viinit, pol, &nev2, &lam, &Yr, &Yi, &res, fstats);

  if (ierr) {
    printf("ZChebLanNr error %d\n", ierr);
  }

  pEVSL_ParvecFree(&vrinit);
  pEVSL_ParvecFree(&viinit);
  /*--------------------- print stats */
  pEVSL_StatsPrint(pevsl, fstats);

  if (res) {
    PEVSL_FREE(res);
  }
  /* save pointers to the global variables */
  pevsl->nev_computed       = nev2;
  pevsl->eval_computed      = lam;
  pevsl->evec_computed      = Yr;
  pevsl->evec_imag_computed = Yi;

}


/** JS 051919 for Lanczos vectors 
 * @brief Fortran interface for ChebLanNr
 *  the results will be saved in the internal variables
 */
void PEVSL_FORT(pevsl_lanvectors)(uintptr_t *pevslf90, double *xintv, int *max_its, double *tol,
                                  uintptr_t *polf90) {
  int N, n, nfirst, nev2, ierr;
  pevsl_Parvecs *Y;
  FILE *fstats = stdout;
  pevsl_Parvec vinit;

  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);

  N = pevsl->N;
  n = pevsl->n;
  nfirst = pevsl->nfirst;
  /*-------------------- zero out stats */
  pEVSL_StatsReset(pevsl);
  /*------------------- Create parallel vector: random initial guess */
  pEVSL_ParvecCreate(N, n, nfirst, pevsl->comm, &vinit);
  pEVSL_ParvecRand(&vinit);
  /* cast pointer of pol*/
  pevsl_Polparams *pol = (pevsl_Polparams *) (*polf90);
  /* lanvectors.c */
  ierr = pEVSL_Lanvectors(pevsl, xintv, *max_its, *tol, &vinit, pol, &nev2, &Y, fstats);

  if (ierr) {
    printf("lanvectors error %d\n", ierr);
  }
  /* save pointers to the global variables */
  pevsl->nev_computed = nev2;
  pevsl->evec_computed = Y;
  //printf("lanvectors done %d\n", ierr);
}


/** @brief Get the number of last computed eigenvalues
 */
void PEVSL_FORT(pevsl_get_nev)(uintptr_t *pevslf90, int *nev) {

  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);

  *nev = pevsl->nev_computed;
}

/** @brief copy the computed eigenvalues and vectors
 * @warning after this call the internal saved results will be freed
 * @param[in, out] pevslf90 Data to be cast to pevsl_Data
 * @param[in] ld leading dimension of vec [ld >= pevsl.n]
 * @param[out] val Eigenvalue output
 * @param[out] vec Eigenvector output
 */
void PEVSL_FORT(pevsl_copy_result)(uintptr_t *pevslf90, double *val, double *vec, int *ld) {
  int i;
  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);

  /* copy eigenvalues */
  memcpy(val, pevsl->eval_computed, pevsl->nev_computed*sizeof(double));

  /* copy eigenvectors */
  for (i=0; i<pevsl->nev_computed; i++) {
    double *dest = vec + i * (*ld);
    double *src = pevsl->evec_computed->data + i * pevsl->evec_computed->ld;
    memcpy(dest, src, pevsl->n*sizeof(double));
  }

  /* reset pointers */
  pevsl->nev_computed = 0;
  PEVSL_FREE(pevsl->eval_computed);
  pEVSL_ParvecsFree(pevsl->evec_computed);
  PEVSL_FREE(pevsl->evec_computed);
}


/** @brief copy the computed eigenvalues and vectors
 * @warning after this call the internal saved results will be freed
 * @param[in, out] pevslf90 Data to be cast to pevsl_Data
 * @param[in] ld leading dimension of vec [ld >= pevsl.n]
 * @param[out] vec Lanczos vectors output
 */
void PEVSL_FORT(pevsl_copy_vectors)(uintptr_t *pevslf90, double *vec, int *ld) {
  int i;
  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);


  /* copy eigenvectors */
  for (i=0; i<pevsl->nev_computed; i++) {
    double *dest = vec + i * (*ld);
    double *src = pevsl->evec_computed->data + i * pevsl->evec_computed->ld;
    memcpy(dest, src, pevsl->n*sizeof(double));
  }

  /* reset pointers */
  pevsl->nev_computed = 0;
  pEVSL_ParvecsFree(pevsl->evec_computed);
  PEVSL_FREE(pevsl->evec_computed);
}



/* added JS 020619 for complex Hermitian 
 * @brief copy the computed eigenvalues and vectors
 * @warning after this call the internal saved results will be freed
 * @param[in, out] pevslf90 Data to be cast to pevsl_Data
 * @param[in] ld leading dimension of vec [ld >= pevsl.n]
 * @param[out] val Eigenvalue output
 * @param[out] vecr, veci Eigenvector output
 */
void PEVSL_FORT(pevsl_copy_zresult)(uintptr_t *pevslf90, double *val, double *vecr, double *veci, int *ld) {

  int i;
  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);
  /* copy eigenvalues */
  memcpy(val, pevsl->eval_computed, pevsl->nev_computed*sizeof(double));

  /* copy eigenvectors */
  for (i=0; i<pevsl->nev_computed; i++) {
    double *destr = vecr + i * (*ld);
    double *srcr  = pevsl->evec_computed->data + i * pevsl->evec_computed->ld;
    memcpy(destr, srcr, pevsl->n*sizeof(double));
    double *desti = veci + i * (*ld);
    double *srci  = pevsl->evec_imag_computed->data + i*pevsl->evec_imag_computed->ld;
    memcpy(desti, srci, pevsl->n*sizeof(double));
  }

  /* reset pointers */
  pevsl->nev_computed = 0;
  PEVSL_FREE(pevsl->eval_computed);
  pEVSL_ParvecsFree(pevsl->evec_computed);
  PEVSL_FREE(pevsl->evec_computed);

  pEVSL_ParvecsFree(pevsl->evec_imag_computed);
  PEVSL_FREE(pevsl->evec_imag_computed);
}

void PEVSL_FORT(pevsl_setup_chebiter)(double *lmin, double *lmax, int *deg,
                                      uintptr_t *parcsrf90, uintptr_t *chebf90) {
  /* cast pointer of the matrix */
  pevsl_Parcsr *A = (pevsl_Parcsr *) (*parcsrf90);
  void *cheb;
  pEVSL_ChebIterSetup(*lmin, *lmax, *deg, A, &cheb);

  *chebf90 = (uintptr_t) cheb;
}

/** @brief Fortran interface for ChebIterSol */
void PEVSL_FORT(pevsl_chebiter)(int *type, double *b, double *x, uintptr_t *chebf90) {

  /* cast pointer */
  void *cheb = (void *) (*chebf90);

  if (*type == 1) {
    pEVSL_ChebIterSolv1(b, x, cheb);
  } else {
    pEVSL_ChebIterSolv2(b, x, cheb);
  }
}

/** @brief Fortran interface for ChebIterStatsPrint */
void PEVSL_FORT(pevsl_chebiterstatsprint)(uintptr_t *chebf90) {

  FILE *fstats = stdout;
  /* cast pointer */
  void *cheb = (void *) (*chebf90);

  pEVSL_ChebIterStatsPrint(cheb, fstats);
}

/** @brief Fortran interface for ChebIterFree */
void PEVSL_FORT(pevsl_free_chebiterb)(uintptr_t *chebf90) {

  /* cast pointer */
  void *cheb = (void *) (*chebf90);

  pEVSL_ChebIterFree(cheb);
}


/** @brief Fortran interface for SetBsol with ChebIterSol */
void PEVSL_FORT(pevsl_setbsol_chebiter)(uintptr_t *pevslf90, int *type, uintptr_t *chebf90) {
  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);
  void *cheb = (void *) (*chebf90);

  if (*type == 1) {
    pEVSL_SetBSol(pevsl, pEVSL_ChebIterSolv1, cheb);
  } else {
    pEVSL_SetBSol(pevsl, pEVSL_ChebIterSolv2, cheb);
  }
}


/** @brief Fortran interface for computing L-S polynomial for B^{-1/2} */
void PEVSL_FORT(pevsl_setup_lspolsqrt)(double *lmin, double *lmax, int *maxdeg, double *tol,
                                       uintptr_t *parcsrf90, uintptr_t *lspolf90) {
  /* cast pointer of the matrix */
  pevsl_Parcsr *A = (pevsl_Parcsr *) (*parcsrf90);
  void *lspol;
  pEVSL_SetupLSPolSqrt(*maxdeg, *tol, *lmin, *lmax, A, &lspol);

  *lspolf90 = (uintptr_t) lspol;
}

/** @brief Fortran interface for pEVSL_LSPolSol */
void PEVSL_FORT(pevsl_lspolsol)(double *b, double *x, uintptr_t *lspolf90) {

  /* cast pointer */
  void *lspol = (void *) (*lspolf90);
  pEVSL_LSPolSol(b, x, lspol);
}

/** @brief Fortran interface for pEVSL_LSPolFree */
void PEVSL_FORT(pevsl_free_lspol)(uintptr_t *lspolf90) {

  /* cast pointer */
  void *lspol = (void *) (*lspolf90);

  pEVSL_LSPolFree(lspol);
}


/** @brief Fortran interface for SetLTsol with LSPolSol */
void PEVSL_FORT(pevsl_setltsol_lspol)(uintptr_t *pevslf90, uintptr_t *lspolf90) {
  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);
  void *lspol = (void *) (*lspolf90);

  pEVSL_SetLTSol(pevsl, pEVSL_LSPolSol, lspol);
}


void PEVSL_FORT(pevsl_kpmdos_ecnt)(uintptr_t *pevslf90, int *Mdeg, int *damping, int *nvec,
                                   double *intv, double *mu, double *ecnt) {
  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);

  pEVSL_Kpmdos(pevsl, *Mdeg, *damping, *nvec, intv, 1, 0, MPI_COMM_NULL, mu, ecnt);
}

void PEVSL_FORT(pevsl_landos_ecnt)(uintptr_t *pevslf90, int *nvec, int *msteps, int *npts,
                                   double *intv, double *ecnt) {
  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);

  double *xdos, *ydos;
  PEVSL_MALLOC(xdos, *npts, double);
  PEVSL_MALLOC(ydos, *npts, double);

  pEVSL_LanDosG(pevsl, *nvec, *msteps, *npts, xdos, ydos, ecnt, intv, 1, 0, MPI_COMM_NULL);

  PEVSL_FREE(xdos);
  PEVSL_FREE(ydos);
}


void PEVSL_FORT(pevsl_set_sigmamult)(uintptr_t *pevslf90, double *mult) {

  /* cast pointer */
  pevsl_Data *pevsl = (pevsl_Data *) (*pevslf90);

  pEVSL_SetSigmaMult(pevsl, *mult);
}



























#if 0
void PEVSL_FORT(pevsl_testchebiterb)(uintptr_t *chebf90, MPI_Fint *Fcomm) {
  int i, N, n, nfirst;
  pevsl_Parvec y,b,x,d;

  ChebiterData *cheb = (ChebiterData *) (*chebf90);
  MPI_Comm comm = MPI_Comm_f2c(*Fcomm);
  N = pevsl_data.N;
  n = pevsl_data.n;
  nfirst = pevsl_data.nfirst;

  pEVSL_ParvecCreate(N, n, nfirst, comm, &b);
  pEVSL_ParvecDupl(&b, &y);
  pEVSL_ParvecDupl(&b, &x);
  pEVSL_ParvecDupl(&b, &d);

  /* y is the exact sol */
  pEVSL_ParvecRand(&y);
  /* rhs: b = B*y */
  pEVSL_MatvecB(&y, &b);

  pEVSL_SolveB(&b, &x);

  if (cheb->res) {
    printf("CHEB ITER RES\n");
    for (i=0; i<cheb->deg+1; i++) {
      printf("i %3d: %e\n", i, cheb->res[i]);
    }
  }

  pEVSL_MatvecB(&x, &d);
  pEVSL_ParvecAxpy(-1.0, &b, &d);
  pEVSL_ParvecAxpy(-1.0, &y, &x);
  double err_r, err_x;
  pEVSL_ParvecNrm2(&d, &err_r);
  pEVSL_ParvecNrm2(&x, &err_x);

  printf("||res|| %e, ||err|| %e\n", err_r, err_x);
}

void PEVSL_FORT(pevsl_parveccreate)(int *N, int *n, MPI_Fint *Fcomm, uintptr_t *xf90) {
  pevsl_Parvec *x;
  MPI_Comm comm = MPI_Comm_f2c(*Fcomm);
  PEVSL_MALLOC(x, 1, pevsl_Parvec);
  pEVSL_ParvecCreate(*N, *n, 0, comm, x);

  *xf90 = (uintptr_t) x;
}

void PEVSL_FORT(pevsl_parvecsetvals)(uintptr_t *xf90, double *vals) {
  int i, n;
  pevsl_Parvec *x = (pevsl_Parvec *) xf90;
  n = x->n_local;
  for (i=0; i<n; i++) {
    x->data[i] = vals[i];
  }
}

void PEVSL_FORT(pevsl_parvecgetvals)(uintptr_t *xf90, double *vals) {
  int i, n;
  pevsl_Parvec *x = (pevsl_Parvec *) xf90;
  n = x->n_local;
  for (i=0; i<n; i++) {
    vals[i] = x->data[i];
  }
}

void PEVSL_FORT(pevsl_parvecfree)(uintptr_t *xf90) {
  pevsl_Parvec *x = (pevsl_Parvec *) xf90;
  pEVSL_ParvecFree(x);
  PEVSL_FREE(x);
}

void PEVSL_FORT(pevsl_amv)(double *x, double *y) {
  pevsl_data.Amv->func(x, y, pevsl_data.Amv->data);
}
#endif
