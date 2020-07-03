#include "pevsl_int.h"

/**
 * @file lanvectors.c
 * @brief Polynomial Filtered no-restart Lanczos 
 * collect all the computed Lanczos vectors
 */

/**
 * if filter the initial vector
 */
#define FILTER_VINIT 1

/* -----------------------------------------------------------------------
 *  @brief Chebyshev polynomial filtering Lanczos process [NON-restarted version]
 *
 *  @param[in] pevsl  pEVSL data struct
 *  @param[in] intv     An array of length 4 \n
 *          [intv[0], intv[1]] is the interval of desired eigenvalues \n
 *          [intv[2], intv[3]] is the global interval of all eigenvalues \n
 *          Must contain all eigenvalues of A
 *
 *  @param[in] maxit    Max number of outer Lanczos steps  allowed --[max dim of Krylov
 *          subspace]
 *
 *  @param[in] tol
 *          Tolerance for convergence. The code uses a stopping criterion based
 *          on the convergence of the restricted trace. i.e., the sum of the
 *          eigenvalues of T_k that  are in the desired interval. This test  is
 *          rather simple since these eigenvalues are between `bar' and  1.0.
 *          We want the relative error on  this restricted  trace to be  less
 *          than  tol.  Note that the test  performed on filtered matrix only
 *          - *but* the actual residual norm associated with the original
 *          matrix A is returned
 *
 *  @param[in] vinit  initial  vector for Lanczos -- [optional]
 *
 *  @param[in] pol       A struct containing the parameters of the polynomial..
 *  This is set up by a call to find_deg prior to calling chenlanNr
 *
 *  @b Modifies:
 *  @param[out] nevOut    Number of eigenvalues/vectors computed
 *  @param[out] Wo        A set of eigenvectors  [n x nevOut matrix]
 *  @param[out] fstats    File stream which stats are printed to
 *
 *  @return Returns 0 on success (or if check_intv() is non-positive),  -1
 *  if |gamB| < 1
 *
 *
 * @warning memory allocation for Wo within this function
 *
 * ------------------------------------------------------------ */
int pEVSL_Lanvectors(pevsl_Data *pevsl, double *intv, int maxit, double tol,
                     pevsl_Parvec *vinit, pevsl_Polparams *pol, int *nevOut,
                     pevsl_Parvecs **Wo, FILE *fstats) {

  //-------------------- to report timings/
  double tall, tm1 = 0.0, tt;
  tall = pEVSL_Wtime();
  const int ifGenEv = pevsl->ifGenEv;
  double tr0, tr1;
  double *y, flami;
  int i, k, kdim=0, rank;
  /* handle case where fstats is NULL. Then no output. Needed for openMP */
  int do_print = 1;
  if (fstats == NULL) {
    do_print = 0;
  }
  /* MPI communicator */
  MPI_Comm comm = pevsl->comm;
  /*
  MPI_Comm_compare(comm, vinit->comm, &comp);
  if (comp != MPI_IDENT) {
    return 1;
  }
  */
  /*-------------------- rank in comm */
  MPI_Comm_rank(comm, &rank);
  /*-------------------- Ntest = when to start testing convergence */
  int Ntest = 30;
  /*-------------------- how often to test */
  int cycle = 20;
  /* size of the matrix. N: global size */
  int N;
  N = pevsl->N;
  /* max num of its */
  maxit = PEVSL_MIN(N, maxit);
  /*-------------------- polynomial filter  approximates the delta
                         function centered at 'gamB'.
                         bar: a bar value to threshold Ritz values of p(A) */
  double bar = pol->bar;
  double gamB = pol->gam;
  /*-------------------- interval [aa, bb], used for testing only at end */
  if (check_intv(intv, fstats) < 0) {
    *nevOut = 0;
    *Wo = NULL; 
    return 0;
  }
  double aa = intv[0];
  double bb = intv[1];
  int deg = pol->deg;
  if (do_print) {
    pEVSL_fprintf0(rank, fstats, " intv:[%e, %e, %e, %e]\n", intv[0], intv[1], intv[2], intv[3]);
    pEVSL_fprintf0(rank, fstats, " ** Cheb Poly of deg = %d, gam = %.15e, bar: %.15e\n",
                   deg, gamB, bar);
  }
  /*-------------------- gamB must be within [-1, 1] */
  if (gamB > 1.0 || gamB < -1.0) {
    pEVSL_fprintf0(rank, stdout, "gamB error %.15e\n", gamB);
    return -1;
  }
  /*-----------------------------------------------------------------------*
   * *Non-restarted* Lanczos iteration
   *-----------------------------------------------------------------------*/
  if (do_print) {
    pEVSL_fprintf0(rank, fstats, " ** Cheb-LanNr \n");
  }
  /*-------------------- Lanczos vectors V_m and tridiagonal matrix T_m */
  pevsl_Parvecs *V, *Z;
  double *dT, *eT;
  PEVSL_MALLOC(V, 1, pevsl_Parvecs);
  pEVSL_ParvecsDuplParvec(vinit, maxit+1, vinit->n_local, V);
  if (ifGenEv) {
    /* storage for Z = B * V */
    PEVSL_MALLOC(Z, 1, pevsl_Parvecs);
    pEVSL_ParvecsDuplParvec(vinit, maxit+1, vinit->n_local, Z);
  } else {
    /* Z and V are the same */
    Z = V;
  }
  /*-------------------- diag. subdiag of Tridiagional matrix */
  PEVSL_MALLOC(dT, maxit, double);
  PEVSL_MALLOC(eT, maxit, double);
  pevsl_Parvecs *Rvec;
  double *Lam, *res, *EvalT, *EvecT;
  /*-------------------- Lam, Rvec: the converged (locked) Ritz values vecs*/
  PEVSL_MALLOC(Lam,   maxit, double);       // holds computed Ritz values
  PEVSL_MALLOC(res,   maxit, double);       // residual norms (w.r.t. ro(A))
  PEVSL_MALLOC(EvalT, maxit, double);       // eigenvalues of tridiag matrix T
  /*-------------------- nev = current number of converged e-pairs
                         nconv = converged eigenpairs from looking at Tk alone */
  int nev, nconv = 0;
  /*-------------------- u is just a pointer. wk: work space [Parvec] */
  pevsl_Parvec *wk, *w2;
  int wk_size = ifGenEv ? 4 : 3;
  PEVSL_MALLOC(wk, wk_size, pevsl_Parvec);
  for (i=0; i<wk_size; i++) {
    pEVSL_ParvecDupl(vinit, &wk[i]);
  }
  w2 = wk + 1;
  /*-------------------- workspace [double * array] */
  double *warr;
  PEVSL_MALLOC(warr, 3*maxit, double);
  /*-------------------- lanczos vectors: Parvec referencing to Parvecs */
  pevsl_Parvec parvec[7];
  pevsl_Parvec *zold  = &parvec[0];
  pevsl_Parvec *z     = &parvec[1];
  pevsl_Parvec *znew  = &parvec[2];
  pevsl_Parvec *v     = &parvec[3];
  pevsl_Parvec *vnew  = &parvec[4];
  pevsl_Parvec *u     = &parvec[5];
#if FILTER_VINIT
  pevsl_Parvec *vrand = &parvec[6];
#endif
  /* v references the 1st columns of V */
  pEVSL_ParvecsGetParvecShell(V, 0, v);
#if FILTER_VINIT
  /*-------------------- compute w = p[(A-cc)/dd] * v */
  /*-------------------- Filter the initial vector */
  pEVSL_ChebAv(pevsl, pol, vinit, v, wk);
  pEVSL_ParvecDupl(vinit, vrand);
#else
  /*-------------------- copy initial vector to V(:,1) */
  pEVSL_ParvecCopy(vinit, v);
#endif
  /*-------------------- normalize it */
  double t, res0;
  if (ifGenEv) {
    /* z references the 1st columns of Z */
    pEVSL_ParvecsGetParvecShell(Z, 0, z);
    /* B norm */
    pEVSL_MatvecB(pevsl, v, z);
    pEVSL_ParvecDot(v, z, &t);
    t = 1.0 / sqrt(t);
    pEVSL_ParvecScal(z, t);
  } else {
    /* 2-norm */
    pEVSL_ParvecNrm2(v, &t);
    t = 1.0 / t;
  }
  /* unit B-norm or 2-norm */
  pEVSL_ParvecScal(v, t);
  /*-------------------- for ortho test */
  double wn = 0.0;
  int nwn = 0;
  /*-------------------- for stopping test [restricted trace]*/
  tr0 = 0;
  /*--------------------  Lanczos recurrence coefficients */
  double alpha, beta=0.0;
  int count = 0;
  // ---------------- main Lanczos loop
  for (k=0; k<maxit; k++) {
    /*-------------------- quick reference to Z(:,k-1) when k>0 */
    if (k > 0) {
      pEVSL_ParvecsGetParvecShell(Z, k-1, zold);
    }
    /*-------------------- a quick reference to V(:,k) */
    pEVSL_ParvecsGetParvecShell(V, k, v);
    /*-------------------- a quick reference to Z(:,k) */
    pEVSL_ParvecsGetParvecShell(Z, k, z);
    /*-------------------- next Lanczos vector V(:,k+1)*/
    pEVSL_ParvecsGetParvecShell(V, k+1, vnew);
    /*-------------------- next Lanczos vector Z(:,k+1)*/
    pEVSL_ParvecsGetParvecShell(Z, k+1, znew);
    /*-------------------- compute w = p[(A-cc)/dd] * v */
    /*-------------------- NOTE: z is used!!! [TODO: FIX ME] */
    pEVSL_ChebAv(pevsl, pol, z, znew, wk);
    /*------------------ znew = znew - beta*zold */
    if (k > 0) {
      pEVSL_ParvecAxpy(-beta, zold, znew);
    }
    /*-------------------- alpha = znew'*v */
    pEVSL_ParvecDot(v, znew, &alpha);
    /*-------------------- T(k,k) = alpha */
    dT[k] = alpha;
    wn += fabs(alpha);
    /*-------------------- znew = znew - alpha*z */
    pEVSL_ParvecAxpy(-alpha, z, znew);
    /*-------------------- FULL reortho to all previous Lan vectors */
    if (ifGenEv) {
      /* znew = znew - Z(:,1:k)*V(:,1:k)'*znew */
      CGS_DGKS2(pevsl, k+1, NGS_MAX, Z, V, znew, warr);
      /* vnew = B \ znew */
      pEVSL_SolveB(pevsl, znew, vnew);
      /*-------------------- beta = (vnew, znew)^{1/2} */
      pEVSL_ParvecDot(vnew, znew, &beta);
      beta = sqrt(beta);
    } else {
      /* vnew = vnew - V(:,1:k)*V(:,1:k)'*vnew */
      /* beta = norm(vnew) */
      CGS_DGKS(pevsl, k+1, NGS_MAX, V, vnew, &beta, warr);
    }
    wn += 2.0 * beta;
    nwn += 3;
    /*-------------------- lucky breakdown test */
    if (beta*nwn < orthTol*wn) {
      if (do_print) {
        pEVSL_fprintf0(rank, fstats, "it %4d: Lucky breakdown, beta = %.15e\n", k, beta);
      }
#if FILTER_VINIT
      /*------------------ generate a new init vector */
      pEVSL_ParvecRand(vrand);
      /*------------------  Filter the initial vector*/
      pEVSL_ChebAv(pevsl, pol, vrand, vnew, wk);
#else
      pEVSL_ParvecRand(vnew);
#endif
      if (ifGenEv) {
        /* vnew = vnew - V(:,1:k)*Z(:,1:k)'*vnew */
        CGS_DGKS2(pevsl, k+1, NGS_MAX, V, Z, vnew, warr);
        pEVSL_MatvecB(pevsl, vnew, znew);
        pEVSL_ParvecDot(vnew, znew, &beta);
        beta = sqrt(beta);
        /*-------------------- vnew = vnew / beta */
        t = 1.0 / beta;
        pEVSL_ParvecScal(vnew, t);
	/*-------------------- znew = znew / beta */
        pEVSL_ParvecScal(znew, t);
        beta = 0.0;
      } else {
        /* vnew = vnew - V(:,1:k)*V(:,1:k)'*vnew */
        /* beta = norm(vnew) */
        CGS_DGKS(pevsl, k+1, NGS_MAX, V, vnew, &beta, warr);
        /*-------------------- vnew = vnew / beta */
        t = 1.0 / beta;
        pEVSL_ParvecScal(vnew, t);
        beta = 0.0;
      }
    } else {
      /*-------------------- vnew = vnew / beta */
      t = 1.0 / beta;
      pEVSL_ParvecScal(vnew, t);
      if (ifGenEv) {
        /*-------------------- znew = znew / beta */
        pEVSL_ParvecScal(znew, t);
      }
    }
    /*-------------------- T(k+1,k) = beta */
    eT[k] = beta;
    /*---------------------- test for Ritz vectors */
    if ( (k < Ntest || (k-Ntest) % cycle != 0) && k != maxit-1 ) {
      continue;
    }
    /*---------------------- diagonalize  T(1:k,1:k)
                             vals in EvalT, vecs in EvecT  */
    kdim = k+1;
#if 1
    /*-------------------- THIS uses dsetv, do not need eig vector */
    SymmTridEig(pevsl, EvalT, NULL, kdim, dT, eT);
    count = kdim;
#else
    /*-------------------- THIS uses dstemr */
    double vl = bar - DBL_EPSILON, vu = 2.0;  /* needed by SymmTridEigS */
    SymmTridEigS(pevsl, EvalT, EvecT, kdim, vl, vu, &count, dT, eT);
#endif
    /*-------------------- restricted trace: used for convergence test */
    tr1 = 0;
    /*-------------------- get residual norms and check acceptance of Ritz
     *                     values for p(A). nconv records number of
     *                     eigenvalues whose residual for p(A) is smaller
     *                     than tol. */
    nconv = 0;
    for (i=0; i<count; i++) {
      flami = EvalT[i];
      if (flami + 5*DBL_EPSILON >= bar) {
        tr1+= flami;
        nconv++;
      }
    }

    if (do_print) {
      pEVSL_fprintf0(rank, fstats, "  k %4d:   nconv %4d  tr1 %21.15e\n",
                     k, nconv,tr1);
    }
    /* -------------------- simple test because all eigenvalues
                            are between gamB and ~1. */
    if ( (fabs(tr1-tr0) < 2e-12) || (fabs(tr1)+fabs(tr0) < 1e-10) ) {
      break;
    }
    tr0 = tr1;
  } /* end of the main loop */

  if (k >= maxit) {
     pEVSL_fprintf0(rank, fstats, "The max number of iterations [%d] has been reached. The eigenvalues computed may not have converged\n", maxit);
  }


  tm1 = pEVSL_Wtime() - tt;

  //pEVSL_fprintf0(rank, fstats, "k = %d; check 0 \n", k);
  /*-------------------- Done.  output : */
  *nevOut = kdim;
  *Wo = V;
 
  //pEVSL_fprintf0(rank, fstats, "k = %d; check 1 \n", k);
  
  /*-------------------- free arrays */
  //pEVSL_ParvecsFree(V);
  //PEVSL_FREE(V);
  PEVSL_FREE(dT);
  PEVSL_FREE(eT);
  PEVSL_FREE(EvalT);
  //PEVSL_FREE(EvecT);
  for (i=0; i<wk_size; i++) {
    pEVSL_ParvecFree(&wk[i]);
  }
  PEVSL_FREE(wk);
#if FILTER_VINIT
  pEVSL_ParvecFree(vrand);
#endif
  if (ifGenEv) {
    pEVSL_ParvecsFree(Z);
    PEVSL_FREE(Z);
  }
  PEVSL_FREE(warr);

  //pEVSL_fprintf0(rank, fstats, "k = %d; check 2 \n", k);
  /*-------------------- record stats */
  tall = pEVSL_Wtime() - tall;

  pevsl->stats->t_iter = tall;
  pevsl->stats->t_ritz = tm1;

  return 0;
}

