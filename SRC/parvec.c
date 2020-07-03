#include "pevsl_int.h"
/**
 * @file parvec.c
 * @brief Parallel vector related functions
 */


/*!
 * @brief Create a parallel vector struct without allocating memory for data
 * @param[in] nglobal Number of global elements
 * @param[in] nlocal Number local elements
 * @param[in] nfirst Index of first local element
 * @param[in] comm Communicator
 * @param[out] x Vector
 * @param[in] data Data
 */
void pEVSL_ParvecCreateShell(int nglobal, int nlocal, int nfirst, MPI_Comm comm,
                             pevsl_Parvec *x, double *data) {

  x->comm = comm;
  x->n_global = nglobal;
  x->n_local = nlocal;
  x->n_first = nfirst < 0 ? PEVSL_NOT_DEFINED : nfirst;
  x->data = data;
}

/*!
 * @brief Create a parallel vector struct
 * @param[in] nglobal Number of global elements
 * @param[in] nlocal Number local elements
 * @param[in] nfirst Index of first local element
 * @param[in] comm Communicator
 * @param[out] x Vector
 */
void pEVSL_ParvecCreate(int nglobal, int nlocal, int nfirst, MPI_Comm comm,
                        pevsl_Parvec *x) {

  pEVSL_ParvecCreateShell(nglobal, nlocal, nfirst, comm, x, NULL);
  PEVSL_MALLOC(x->data, nlocal, double);
}

/*!
 * @brief Create a parallel vector struct
 * @param[in] x Input vector
 * @param[out] y Output vector
 */
void pEVSL_ParvecDupl(pevsl_Parvec *x, pevsl_Parvec *y) {

  pEVSL_ParvecCreate(x->n_global, x->n_local, x->n_first, x->comm, y);
}

/*! JS 01/02/19 create x*H change the imag part xi to -xi
 * @param[xi] xi Input vector 
 * @param[xiH] xi^H Output vector 
 */
void pEVSL_ParvecCreateH(pevsl_Parvec *xi, pevsl_Parvec *xiH) {

  pEVSL_ParvecDupl(xi, xiH); 

  int i;
  for (i=0; i<xi->n_local; i++) {
    xiH->data[i] = - xi->data[i];
  }
}

/*!
 * @brief Destroy a Parvec struct
 */
void pEVSL_ParvecFree(pevsl_Parvec *x) {

  PEVSL_FREE(x->data);
}


/*!
 * @brief generate a random parallel vector, each element in [-1,1]
 * @param[in,out] x Vector
 */
void pEVSL_ParvecRand(pevsl_Parvec *x) {
  int i;
  double t = ((double) RAND_MAX) / 2.0;
  for (i=0; i<x->n_local; i++) {
    double z = (rand() - t) / t;
    x->data[i] = z;
  }
}

/* JS to split vector initialization
 * @brief generate a random parallel vector, each element in [-1,1]
 * @param[in,out] x Vector
 */
void pEVSL_ParvecRand_split(pevsl_Parvec *x) {
  int i;
  double t = ((double) RAND_MAX) / 2.0;
  for (i=0; i<x->n_local/2; i++) {
    double z = (rand() - t) / t;
    x->data[i] = z;
    x->data[x->n_local/2+i] = z;
  }
}

/*!
 * @brief generate a random parallel vector
 * @warning This is generated in such a way that it is the same as the
 * sequential vector with the same starting seed. NOT scalable, ONLY FOR DEBUG
 *
 * @param[in,out] x Vector
 */
void pEVSL_ParvecRand2(pevsl_Parvec *x) {
  int i;
  double t = ((double) RAND_MAX) / 2.0;
  /* parallel random vector is generated in a way such that
   * it is the sames as sequential vector with the same starting seed
   * In this case, it is also independent of number of MPI ranks.
   * [NOT scalable] ONLY for DEBUG
   */
  int np, rank, *nlocal_all, nfirst;
  nfirst = x->n_first;
  if (nfirst == PEVSL_NOT_DEFINED) {
    /* if nfirst has not been defined, compute it */
    MPI_Comm comm = x->comm;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    PEVSL_MALLOC(nlocal_all, np, int);
    MPI_Allgather(&x->n_local, 1, MPI_INT, nlocal_all, 1, MPI_INT, comm);
    nfirst=0;
    /* compute nfirst */
    for (i=0; i<rank; i++) {
      nfirst += nlocal_all[i];
    }
    PEVSL_FREE(nlocal_all);
  }
  /* populate the local portion of the vector */
  for (i=0; i<x->n_global; i++) {
    double z = (rand() - t) / t;
    if (i >= nfirst && i < nfirst + x->n_local) {
      x->data[i-nfirst] = z;
    }
  }
}



/**
 * @brief Parallel dot product
 * @param[in] x Input vector
 * @param[in] y Input vector
 * @param[out] t dot product
 *
 */
void pEVSL_ParvecDot(pevsl_Parvec *x, pevsl_Parvec *y, double *t) {
  PEVSL_CHKERR(x->n_global != y->n_global);
  PEVSL_CHKERR(x->n_local != y->n_local);
  PEVSL_CHKERR(x->n_first != y->n_first);
  double tlocal;
  int one = 1;
  tlocal = DDOT(&(x->n_local), x->data, &one, y->data, &one);
  MPI_Allreduce(&tlocal, t, 1, MPI_DOUBLE, MPI_SUM, x->comm);
}

/** JS 12/26/18 
 * computes x^H * y = tr + ti * i
 * @brief Parallel complex dot product
 * @param[in] xr, xi Input vector
 * @param[in] yr, yi Input vector
 * @param[out] tr, ti dot product
 *
 */
void pEVSL_ParvecZDot(pevsl_Parvec *xr, pevsl_Parvec *xi, pevsl_Parvec *yr, 
                      pevsl_Parvec *yi, double *tr, double *ti) {
  PEVSL_CHKERR(xr->n_global != yr->n_global);
  PEVSL_CHKERR(xr->n_local  != yr->n_local);
  PEVSL_CHKERR(xr->n_first  != yr->n_first);
  PEVSL_CHKERR(xi->n_global != yi->n_global);
  PEVSL_CHKERR(xi->n_local  != yi->n_local);
  PEVSL_CHKERR(xi->n_first  != yi->n_first);
  double tlocal0,tlocal1,t0,t1;
  /* make x complex conjugate*/
  /*pevsl_Parvec *xiH;
  pEVSL_ParvecCreateH(xi,xiH); 
  */

  int one = 1;
  tlocal0 = DDOT(&(xr->n_local), xr->data, &one, yr->data, &one);
  tlocal1 = DDOT(&(xi->n_local), xi->data, &one, yi->data, &one);
  t0 = tlocal0 + tlocal1;
  MPI_Allreduce(&t0, tr, 1, MPI_DOUBLE, MPI_SUM, xr->comm);
  tlocal0 = DDOT(&(xr->n_local), xr->data, &one, yi->data, &one);
  tlocal1 = DDOT(&(xi->n_local), xi->data, &one, yr->data, &one);
  t1 = tlocal0 - tlocal1;
  MPI_Allreduce(&t1, ti, 1, MPI_DOUBLE, MPI_SUM, xr->comm);
  /* pEVSL_ParvecFree(xiH); */
}


/**
 * @brief Euclidean norm
 * @param[in] x Input vector
 * @param[out] t Euclidean norm
 * */
void pEVSL_ParvecNrm2(pevsl_Parvec *x, double *t) {
  double t2;
  pEVSL_ParvecDot(x, x, &t2);
  *t = sqrt(t2);
}

/** JS 12/26/18 
 * @brief Euclidean norm
 * @param[in] xr, xi Input vector
 * @param[out] t Euclidean norm
 * */
void pEVSL_ParvecZNrm2(pevsl_Parvec *xr, pevsl_Parvec *xi, double *t) {
  double t1, t2;

  /* pevsl_Parvec *xiH;   
  pEVSL_ParvecCreateH(xi,xiH); */ 

  pEVSL_ParvecZDot(xr, xi, xr, xi, &t1, &t2);
  *t = sqrt(t1);

  /* pEVSL_ParvecFree(xiH); */
}


// y := x
/**
 * @brief Copies x to y
 * @param[in] x Input vector
 * @param[out] y Output vector
 */
void pEVSL_ParvecCopy(pevsl_Parvec *x, pevsl_Parvec *y) {
  PEVSL_CHKERR(x->n_global != y->n_global);
  PEVSL_CHKERR(x->n_local  != y->n_local);
  PEVSL_CHKERR(x->n_first  != y->n_first);
  int one = 1;
  DCOPY(&(x->n_local), x->data, &one, y->data, &one);
}

/**
 * @brief Sums up each element of the vector
 * @param[in] x Input vector
 * @param[out] t Sum of elements
 * */
void pEVSL_ParvecSum(pevsl_Parvec *x, double *t) {
  double localsum = 0.0;
  int i;
  for (i=0; i<x->n_local; i++) {
    localsum += x->data[i];
  }
  MPI_Allreduce(&localsum, t, 1, MPI_DOUBLE, MPI_SUM, x->comm);
}

/**
 * @brief Scales vector
 * @param[in,out] x Vector to scale
 * @param[in] t scale
 * */
void pEVSL_ParvecScal(pevsl_Parvec *x, double t) {
  int one = 1;
  DSCAL(&(x->n_local), &t, x->data, &one);
}

/*
   void pEVSL_ParvecAddScalar(pevsl_Parvec *x, double t) {
   int i;
   for (i=0; i<x->n_local; i++) {
   x->data[i] += t;
   }
   }
   */

/**
 * @brief Sets all elements of a vector to a scalar
 * @param[in, out] x Vector to set to constant
 * @param[in] t scalar to set vector to
 *
 */
void pEVSL_ParvecSetScalar(pevsl_Parvec *x, double t) {
  int one = 1, zero = 0;
  DCOPY(&(x->n_local), &t, &zero, x->data, &one);
}

/**
 * @brief Sets a vector to 0
 * @param[in,out] x Vector to set to 0
 * */
void pEVSL_ParvecSetZero(pevsl_Parvec *x) {
  memset(x->data, 0, x->n_local*sizeof(double));
}

/**
 * @brief Parallel y=ax+y
 * @param[in] a scalar
 * @param[in] x input vector
 * @param[in, out] y input and output vector
 *
 */
void pEVSL_ParvecAxpy(double a, pevsl_Parvec *x, pevsl_Parvec *y) {
  PEVSL_CHKERR(x->n_global != y->n_global);
  PEVSL_CHKERR(x->n_local != y->n_local);
  PEVSL_CHKERR(x->n_first != y->n_first);
  int one = 1;
  DAXPY(&(x->n_local), &a, x->data, &one, y->data, &one);
}


/**
 * @brief Checks if two vectors are the same size
 * @param[in] x Input vector
 * @param[in] y Input vector
 * @return 1 if same size, else 0
 *
 */
int pEVSL_ParvecSameSize(pevsl_Parvec *x, pevsl_Parvec *y) {
  if (x->n_global != y->n_global || x->n_local != y->n_local ||
      x->n_first != y->n_first) {
    return 0;
  }

  return 1;
}

/** @brief print Parvec on comm (squential in comm)
 * @param[in] x Input vector
 *  @param[in] fn filename (if NULL, print it to stdout) */
int pEVSL_ParvecWrite(pevsl_Parvec *x, const char *fn) {
  MPI_Comm comm = x->comm;
  /* if fn == NULL, print to stdout */
  FILE *fp = fn ? NULL : stdout;
  /* sequential loop with all ranks in comm */
  PEVSL_SEQ_BEGIN(comm)
  {
    /* if print to file, rank 0 open it */
    if (fn) {
      if (rank) {
        fp = fopen(fn, "a");
      } else {
        fp = fopen(fn, "w");
      }
    }
    /* if fp != NULL */
    if (fp) {
      /* rank 0 writes the global size */
      if (rank == 0) {
        fprintf (fp, "%% %d\n", x->n_global);
      }
      /* write local data */
      int j;
      for (j=0; j<x->n_local; j++) {
        fprintf(fp, "%.15e\n", x->data[j]);
      }
      if (fn) {
        fclose(fp);
      }
    }
  }
  PEVSL_SEQ_END(comm)

    return 0;
}

/*
   int pEVSL_ParvecRead(pevsl_Parvec *x, const char* fn) {
   char str[PEVSL_MAX_LINE];
   int k;
   FILE *fp = fopen(fn, "r");
   do {
   fgets(str, PEVSL_MAX_LINE, fp);
   } while (str[0] == '%');

   for (k=0; k<x->n_first+x->n_local; k++) {
   if (k >= x->n_first) {
   x->data[k-x->n_first] = atof(str);
   }
   fscanf(fp, "%s\n", str);
   }
   fclose(fp);

   return 0;
   }
   */
