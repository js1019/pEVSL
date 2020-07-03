#include "pevsl_int.h"

/**
 * @file cheblanNr.c
 * @brief Polynomial Filtered no-restart Lanczos
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
 *  @param[out] reso      Associated residual norms [nev x 1 vector]
 *  @param[out] lamo      Lambda computed
 *  @param[out] fstats    File stream which stats are printed to
 *
 *  @return Returns 0 on success (or if check_intv() is non-positive),  -1
 *  if |gamB| < 1
 *
 *
 * @warning memory allocation for Wo/lamo/reso within this function
 *
 * ------------------------------------------------------------ */
int pEVSL_ChebLanNr(pevsl_Data *pevsl, double *intv, int maxit, double tol,
                    pevsl_Parvec *vinit, pevsl_Polparams *pol, int *nevOut,
                    double **lamo, pevsl_Parvecs **Wo,
                    double **reso, FILE *fstats) {

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
    *lamo = NULL; *Wo = NULL; *reso = NULL;
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

  /*-------------------- compute eig vals and vector */
  size_t kdim_l = kdim; /* just in case if kdim is > 65K */
  PEVSL_MALLOC(EvecT, kdim_l*kdim_l, double); // Eigen vectors of T
  SymmTridEig(pevsl, EvalT, EvecT, kdim, dT, eT);

  tt = pEVSL_Wtime();
  /*-------------------- done == compute Ritz vectors */
  PEVSL_MALLOC(Rvec, 1, pevsl_Parvecs);
  pEVSL_ParvecsDuplParvec(vinit, nconv, vinit->n_local, Rvec);

  nev = 0;
  for (i=0; i<count; i++) {
    flami = EvalT[i];
    //-------------------- reject eigenvalue if rho(lam)<bar
    if (flami < bar) {
      continue;
    }
    y = &EvecT[i*kdim_l];
    /*-------------------- make sure to normalize */
    /*
    t = DNRM2(&kdim, y, &one);
    t = 1.0 / t;
    DSCAL(&kdim, &t, y, &one);
    */
    /*-------------------- compute Ritz vectors */
    pEVSL_ParvecsGetParvecShell(Rvec, nev, u);
    pEVSL_ParvecsGemv(1.0, V, kdim, y, 0.0, u);
    /*-------------------- normalize u */
    if (ifGenEv) {
      /* B-norm, w2 = B*u */
      pEVSL_MatvecB(pevsl, u, w2);
      pEVSL_ParvecDot(u, w2, &t);
      t = sqrt(t); /* should be one */
    } else {
      /* 2-norm */
      pEVSL_ParvecNrm2(u, &t); /* should be one */
    }
    /*-------------------- return code 2 --> zero eigenvector found */
    if (t == 0.0) {
      return 2;
    }
    /*-------------------- scal u */
    t = 1.0 / t;
    pEVSL_ParvecScal(u, t);
    /*-------------------- scal B*u */
    if (ifGenEv) {
      /*------------------ w2 = B*u */
      pEVSL_ParvecScal(w2, t);
    }
    /*-------------------- w = A*u */
    pEVSL_MatvecA(pevsl, u, wk);
    /*-------------------- Ritz val: t = (u'*w)/(u'*u)
                                     t = (u'*w)/(u'*B*u) */
    pEVSL_ParvecDot(wk, u, &t);
    /*-------------------- if lambda (==t) is in [a,b] */
    if (t < aa - DBL_EPSILON || t > bb + DBL_EPSILON) {
      continue;
    }
    /*-------------------- compute residual wrt A for this pair */
    if (ifGenEv) {
      /*-------------------- w = w - t*B*u */
      pEVSL_ParvecAxpy(-t, w2, wk);
    } else {
      /*-------------------- w = w - t*u */
      pEVSL_ParvecAxpy(-t, u, wk);
    }
    /*-------------------- res0 = 2-norm(wk) */
    pEVSL_ParvecNrm2(wk, &res0);
    /*--------------------   accept (t, y) */
    if (res0 < tol) {
      Lam[nev] = t;
      res[nev] = res0;
      nev++;
    }
  }
  tm1 = pEVSL_Wtime() - tt;

  /*-------------------- Done.  output : */
  *nevOut = nev;
  *lamo = Lam;
  *Wo = Rvec;
  *reso = res;
  /*-------------------- free arrays */
  pEVSL_ParvecsFree(V);
  PEVSL_FREE(V);
  PEVSL_FREE(dT);
  PEVSL_FREE(eT);
  PEVSL_FREE(EvalT);
  PEVSL_FREE(EvecT);
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

  /*-------------------- record stats */
  tall = pEVSL_Wtime() - tall;

  pevsl->stats->t_iter = tall;
  pevsl->stats->t_ritz = tm1;

  return 0;
}


/* JS 02/05/2019
 * -----------------------------------------------------------------------
 *  @brief Chebyshev polynomial filtering Lanczos process [NON-restarted version]
 *  for complex Hermitian problems
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
 *  @param[in] vrinit  initial  vector for Lanczos -- [optional]
 *  @param[in] viinit  initial  vector for Lanczos -- [optional]
 *
 *  @param[in] pol       A struct containing the parameters of the polynomial..
 *  This is set up by a call to find_deg prior to calling chenlanNr
 *
 *  @b Modifies:
 *  @param[out] nevOut    Number of eigenvalues/vectors computed
 *  @param[out] Wor, Woi  A set of eigenvectors  [n x nevOut matrix]
 *  @param[out] reso      Associated residual norms [nev x 1 vector]
 *  @param[out] lamo      Lambda computed
 *  @param[out] fstats    File stream which stats are printed to
 *
 *  @return Returns 0 on success (or if check_intv() is non-positive),  -1
 *  if |gamB| < 1
 *
 *
 * @warning memory allocation for Wor/Woi/lamo/reso within this function
 *
 * ------------------------------------------------------------ */

int pEVSL_ZChebLanNr(pevsl_Data *pevsl, double *intv, int maxit, double tol,
                     pevsl_Parvec *vrinit, pevsl_Parvec *viinit, pevsl_Polparams *pol, int *nevOut,
                     double **lamo, pevsl_Parvecs **Wor, pevsl_Parvecs **Woi, 
                     double **reso, FILE *fstats) {

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
    *lamo = NULL; *Wor = NULL; *Woi = NULL; *reso = NULL;
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
  pevsl_Parvecs *Vr, *Vi, *Zr, *Zi;
  double *dT, *eT;
  PEVSL_MALLOC(Vr, 1, pevsl_Parvecs);
  PEVSL_MALLOC(Vi, 1, pevsl_Parvecs);
  pEVSL_ParvecsDuplParvec(vrinit, maxit+1, vrinit->n_local, Vr);
  pEVSL_ParvecsDuplParvec(viinit, maxit+1, viinit->n_local, Vi);

  if (ifGenEv) {
    /* storage for Z = B * V */
    PEVSL_MALLOC(Zr, 1, pevsl_Parvecs);
    PEVSL_MALLOC(Zi, 1, pevsl_Parvecs);
    pEVSL_ParvecsDuplParvec(vrinit, maxit+1, vrinit->n_local, Zr);
    pEVSL_ParvecsDuplParvec(viinit, maxit+1, viinit->n_local, Zi);
  } else {
    /* Z and V are the same */
    Zr = Vr;
    Zi = Vi;
  }

  /*-------------------- diag. subdiag of Tridiagional matrix */
  PEVSL_MALLOC(dT, maxit, double);
  PEVSL_MALLOC(eT, maxit, double);
  pevsl_Parvecs *Rvecr, *Rveci;
  double *Lam, *res, *EvalT, *EvecT;
  /*-------------------- Lam, Rvec: the converged (locked) Ritz values vecs*/
  PEVSL_MALLOC(Lam,   maxit, double);       // holds computed Ritz values
  PEVSL_MALLOC(res,   maxit, double);       // residual norms (w.r.t. ro(A))
  PEVSL_MALLOC(EvalT, maxit, double);       // eigenvalues of tridiag matrix T

  /* need to be changed JS 02/06/19 */
  int nev, nconv = 0;
  /*-------------------- u is just a pointer. wk: work space [Parvec] */
  pevsl_Parvec *wkr, *wki, *w2r, *w2i;
  int wk_size = ifGenEv ? 4 : 3;
  PEVSL_MALLOC(wkr, wk_size, pevsl_Parvec);
  PEVSL_MALLOC(wki, wk_size, pevsl_Parvec);
  for (i=0; i<wk_size; i++) {
    pEVSL_ParvecDupl(vrinit, &wkr[i]);
    pEVSL_ParvecDupl(viinit, &wki[i]);
  }
  w2r = wkr + 1;
  w2i = wki + 1;

  /*-------------------- workspace [double * array] */
  double *warr, *wari;
  PEVSL_MALLOC(warr, (2*NGS_MAX+1)*maxit, double);
  PEVSL_MALLOC(wari, (2*NGS_MAX+1)*maxit, double);
  /*-------------------- lanczos vectors: Parvec referencing to Parvecs */
  pevsl_Parvec parvec[14];
  pevsl_Parvec *zrold  = &parvec[0];
  pevsl_Parvec *ziold  = &parvec[1];
  pevsl_Parvec *zr     = &parvec[2];
  pevsl_Parvec *zi     = &parvec[3];
  pevsl_Parvec *zrnew  = &parvec[4];
  pevsl_Parvec *zinew  = &parvec[5];
  pevsl_Parvec *vr     = &parvec[6];
  pevsl_Parvec *vi     = &parvec[7];
  pevsl_Parvec *vrnew  = &parvec[8];
  pevsl_Parvec *vinew  = &parvec[9];
  pevsl_Parvec *ur     = &parvec[10];
  pevsl_Parvec *ui     = &parvec[11];
#if FILTER_VINIT
  pevsl_Parvec *vrrand = &parvec[12];
  pevsl_Parvec *virand = &parvec[13];
#endif

  /* v references the 1st columns of V */
  pEVSL_ParvecsGetParvecShell(Vr, 0, vr);
  pEVSL_ParvecsGetParvecShell(Vi, 0, vi);
#if FILTER_VINIT
  /*-------------------- compute w = p[(A-cc)/dd] * v */
  /*-------------------- Filter the initial vector */
  pEVSL_ZChebAv(pevsl, pol, vrinit, viinit, vr, vi, wkr, wki);
  pEVSL_ParvecDupl(vrinit, vrrand);
  pEVSL_ParvecDupl(viinit, virand);
#else
  /*-------------------- copy initial vector to V(:,1) */
  pEVSL_ParvecCopy(vrinit, vr);
  pEVSL_ParvecCopy(viinit, vi);
#endif

  /*-------------------- normalize it */
  double t, tr, ti, res0;
  if (ifGenEv) {
    /* z references the 1st columns of Z */
    pEVSL_ParvecsGetParvecShell(Zr, 0, zr);
    pEVSL_ParvecsGetParvecShell(Zi, 0, zi);
    /* B norm */
    pEVSL_ZMatvecB(pevsl, vr, vi, zr, zi);
    pEVSL_ParvecZDot(vr, vi, zr, zi, &tr, &ti);
    t = 1.0 / sqrt(tr);
    pEVSL_ParvecScal(zr, t);
    pEVSL_ParvecScal(zi, t);
  } else {
    /* 2-norm */
    pEVSL_ParvecZNrm2(vr, vi, &t);
    t = 1.0 / t;
  }

  /* unit B-norm or 2-norm */
  pEVSL_ParvecScal(vr, t);
  pEVSL_ParvecScal(vi, t);

  /*-------------------- for ortho test */
  double wn = 0.0;
  int nwn = 0;
  /*-------------------- for stopping test [restricted trace]*/
  tr0 = 0;
  /*--------------------  Lanczos recurrence coefficients */
  double alpha, alphar, alphai;  
  double beta=0.0, betar, betai;
  int count = 0;

  // ---------------- main Lanczos loop
  for (k=0; k<maxit; k++) {
    /*-------------------- quick reference to Z(:,k-1) when k>0 */
    if (k > 0) {
      pEVSL_ParvecsGetParvecShell(Zr, k-1, zrold);
      pEVSL_ParvecsGetParvecShell(Zi, k-1, ziold);
    }
    /*-------------------- a quick reference to V(:,k) */
    pEVSL_ParvecsGetParvecShell(Vr, k, vr);
    pEVSL_ParvecsGetParvecShell(Vi, k, vi);
    /*-------------------- a quick reference to Z(:,k) */
    pEVSL_ParvecsGetParvecShell(Zr, k, zr);
    pEVSL_ParvecsGetParvecShell(Zi, k, zi);
    /*-------------------- next Lanczos vector V(:,k+1)*/
    pEVSL_ParvecsGetParvecShell(Vr, k+1, vrnew);
    pEVSL_ParvecsGetParvecShell(Vi, k+1, vinew);
    /*-------------------- next Lanczos vector Z(:,k+1)*/
    pEVSL_ParvecsGetParvecShell(Zr, k+1, zrnew);
    pEVSL_ParvecsGetParvecShell(Zi, k+1, zinew);
    /*-------------------- compute w = p[(A-cc)/dd] * v */
    /*-------------------- NOTE: z is used!!! [TODO: FIX ME] */
    pEVSL_ZChebAv(pevsl, pol, zr, zi, zrnew, zinew, wkr, wki);
    /*------------------ znew = znew - beta*zold */
    if (k > 0) {
      pEVSL_ParvecAxpy(-beta, zrold, zrnew);
      pEVSL_ParvecAxpy(-beta, ziold, zinew);
    }
    /*-------------------- alpha = znew'*v */
    pEVSL_ParvecZDot(vr, vi, zrnew, zinew, &alphar, &alphai);
    alpha = alphar;
    /*-------------------- T(k,k) = alpha */
    dT[k] = alpha;
    wn += fabs(alpha);
    /*-------------------- znew = znew - alpha*z */
    pEVSL_ParvecAxpy(-alpha, zr, zrnew);
    pEVSL_ParvecAxpy(-alpha, zi, zinew);

    /*-------------------- FULL reortho to all previous Lan vectors */
    if (ifGenEv) {
      /* znew = znew - Z(:,1:k)*V(:,1:k)'*znew */
      CGS_ZDGKS2(pevsl, k+1, NGS_MAX, Zr, Zi, Vr, Vi, zrnew, zinew, warr, wari);
      /* vnew = B \ znew */
      pEVSL_ZSolveB(pevsl, zrnew, zinew, vrnew, vinew);
      /*-------------------- beta = (vnew, znew)^{1/2} */
      pEVSL_ParvecZDot(vrnew,vinew, zrnew, zinew, &betar, &betai);
      beta = sqrt(betar);
    } else {
      /* vnew = vnew - V(:,1:k)*V(:,1:k)'*vnew */
      /* beta = norm(vnew) */
      CGS_ZDGKS(pevsl, k+1, NGS_MAX, Vr, Vi, vrnew, vinew, &beta, warr, wari);
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
      pEVSL_ParvecRand(vrrand);
      pEVSL_ParvecRand(virand);
      /*------------------  Filter the initial vector*/
      pEVSL_ZChebAv(pevsl, pol, vrrand, virand, vrnew, vinew, wkr, wki);
#else
      pEVSL_ParvecRand(vrnew);
      pEVSL_ParvecRand(vinew);
#endif
      if (ifGenEv) {
        /* vnew = vnew - V(:,1:k)*Z(:,1:k)'*vnew */
        CGS_ZDGKS2(pevsl, k+1, NGS_MAX, Vr, Vi, Zr, Zi, vrnew, vinew, warr, wari);
        pEVSL_ZMatvecB(pevsl, vrnew, vinew, zrnew, zinew);
        pEVSL_ParvecZDot(vrnew, vinew, zrnew, zinew, &betar, &betai);
        beta = sqrt(betar);
        /*-------------------- vnew = vnew / beta */
        t = 1.0 / beta;
        pEVSL_ParvecScal(vrnew, t);
        pEVSL_ParvecScal(vinew, t);
	/*-------------------- znew = znew / beta */
        pEVSL_ParvecScal(zrnew, t);
        pEVSL_ParvecScal(zinew, t);
        beta = 0.0;
      } else {
        /* vnew = vnew - V(:,1:k)*V(:,1:k)'*vnew */
        /* beta = norm(vnew) */
        CGS_ZDGKS(pevsl, k+1, NGS_MAX, Vr, Vi, vrnew, vinew, &beta, warr, wari);
        /*-------------------- vnew = vnew / beta */
        t = 1.0 / beta;
        pEVSL_ParvecScal(vrnew, t);
        pEVSL_ParvecScal(vinew, t);
        beta = 0.0;
      }
    } else {
      /*-------------------- vnew = vnew / beta */
      t = 1.0 / beta;
      pEVSL_ParvecScal(vrnew, t);
      pEVSL_ParvecScal(vinew, t);
      if (ifGenEv) {
        /*-------------------- znew = znew / beta */
        pEVSL_ParvecScal(zrnew, t);
        pEVSL_ParvecScal(zinew, t);
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

  /*-------------------- compute eig vals and vector */
  size_t kdim_l = kdim; /* just in case if kdim is > 65K */
  PEVSL_MALLOC(EvecT, kdim_l*kdim_l, double); // Eigen vectors of T
  SymmTridEig(pevsl, EvalT, EvecT, kdim, dT, eT);

  tt = pEVSL_Wtime();
  /*-------------------- done == compute Ritz vectors */
  PEVSL_MALLOC(Rvecr, 1, pevsl_Parvecs);
  PEVSL_MALLOC(Rveci, 1, pevsl_Parvecs);
  pEVSL_ParvecsDuplParvec(vrinit, nconv, viinit->n_local, Rvecr);
  pEVSL_ParvecsDuplParvec(viinit, nconv, viinit->n_local, Rveci);

  nev = 0;
  for (i=0; i<count; i++) {
    flami = EvalT[i];
    //-------------------- reject eigenvalue if rho(lam)<bar
    if (flami < bar) {
      continue;
    }
    y = &EvecT[i*kdim_l];

    /*-------------------- compute Ritz vectors */
    pEVSL_ParvecsGetParvecShell(Rvecr, nev, ur);
    pEVSL_ParvecsGetParvecShell(Rveci, nev, ui);
    pEVSL_ParvecsGemv(1.0, Vr, kdim, y, 0.0, ur);
    pEVSL_ParvecsGemv(1.0, Vi, kdim, y, 0.0, ui);

    if (ifGenEv) {
      /* B-norm, w2 = B*u */
      pEVSL_ZMatvecB(pevsl, ur, ui, w2r, w2i);
      pEVSL_ParvecZDot(ur, ui, w2r, w2i, &tr, &ti);
      t = sqrt(tr); /* should be one */
    } else {
      /* 2-norm */
      pEVSL_ParvecZNrm2(ur, ui, &t); /* should be one */
    }

    /*-------------------- return code 2 --> zero eigenvector found */
    if (t == 0.0) {
      return 2;
    }
    /*-------------------- scal u */
    t = 1.0 / t;
    pEVSL_ParvecScal(ur, t);
    pEVSL_ParvecScal(ui, t);
    /*-------------------- scal B*u */
    if (ifGenEv) {
      /*------------------ w2 = B*u */
      pEVSL_ParvecScal(w2r, t);
      pEVSL_ParvecScal(w2i, t);
    }
    /*-------------------- w = A*u */
    pEVSL_ZMatvecA(pevsl, ur, ui, wkr, wki);
    /*-------------------- Ritz val: t = (u'*w)/(u'*u)
                                     t = (u'*w)/(u'*B*u) */
    pEVSL_ParvecZDot(wkr, wki, ur, ui, &tr, &ti);
    t = tr;
    /*-------------------- if lambda (==t) is in [a,b] */
    if (t < aa - DBL_EPSILON || t > bb + DBL_EPSILON) {
      continue;
    }
    /*-------------------- compute residual wrt A for this pair */
    if (ifGenEv) {
      /*-------------------- w = w - t*B*u */
      pEVSL_ParvecAxpy(-t, w2r, wkr);
      pEVSL_ParvecAxpy(-t, w2i, wki);
    } else {
      /*-------------------- w = w - t*u */
      pEVSL_ParvecAxpy(-t, ur, wkr);
      pEVSL_ParvecAxpy(-t, ui, wki);
    }

    /*-------------------- res0 = 2-norm(wk) */
    pEVSL_ParvecZNrm2(wkr, wki, &res0);
    /*--------------------   accept (t, y) */
    if (res0 < tol) {
      Lam[nev] = t;
      res[nev] = res0;
      nev++;
    }
  }
  tm1 = pEVSL_Wtime() - tt;

  /*-------------------- Done.  output : */
  *nevOut = nev;
  *lamo = Lam;
  *Wor = Rvecr;
  *Woi = Rveci;
  *reso = res;

  /*-------------------- free arrays */
  pEVSL_ParvecsFree(Vr);
  pEVSL_ParvecsFree(Vi);
  PEVSL_FREE(Vr);
  PEVSL_FREE(Vi);
  PEVSL_FREE(dT);
  PEVSL_FREE(eT);
  PEVSL_FREE(EvalT);
  PEVSL_FREE(EvecT);
  for (i=0; i<wk_size; i++) {
    pEVSL_ParvecFree(&wkr[i]);
    pEVSL_ParvecFree(&wki[i]);
  }
  PEVSL_FREE(wkr);
  PEVSL_FREE(wki);
#if FILTER_VINIT
  pEVSL_ParvecFree(vrrand);
  pEVSL_ParvecFree(virand);
#endif
  if (ifGenEv) {
    pEVSL_ParvecsFree(Zr);
    pEVSL_ParvecsFree(Zi);
    PEVSL_FREE(Zr);
    PEVSL_FREE(Zi);
  }
  PEVSL_FREE(warr);
  PEVSL_FREE(wari);

  /*-------------------- record stats */
  tall = pEVSL_Wtime() - tall;

  pevsl->stats->t_iter = tall;
  pevsl->stats->t_ritz = tm1;

  return 0; 
}

