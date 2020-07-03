#include "pevsl_int.h"

/**
 * @file pevsl.c
 * @brief pEVSL interface functions
 */
/**
 * @brief Initialize pEVSL
 * @param[in] comm MPI comm
 * @param[out] data pevsl data struct
 *
 * */
int pEVSL_Start(MPI_Comm comm, pevsl_Data **data) {

  pevsl_Data *pevsl_data;
  PEVSL_MALLOC(pevsl_data, 1, pevsl_Data);

  /* Initialize pevsl_data */
  pevsl_data->comm    = comm;
  pevsl_data->N       = PEVSL_NOT_DEFINED;
  pevsl_data->n       = PEVSL_NOT_DEFINED;
  pevsl_data->nfirst  = PEVSL_NOT_DEFINED;
  pevsl_data->ifGenEv = 0;
  pevsl_data->Amv     = NULL;
  pevsl_data->Bmv     = NULL;
  pevsl_data->Bsol    = NULL;
  pevsl_data->LTsol   = NULL;

  pevsl_data->sigma_mult = 1.0;

  /* Use MPI_COMM_WORLD rank as rand seed,
   * so each proc will have a different seed,
   * cf. parvec.c: pEVSL_ParvecRand */
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  srand(rank);

  /* rest stats */
  PEVSL_MALLOC(pevsl_data->stats, 1, pevsl_Stat);
  pEVSL_StatsReset(pevsl_data);

  *data = pevsl_data;

  return 0;
}

/**
 * @brief Finish pEVSL
 * @param[in,out] pevsl_data pEVSL data struct
 *
 * */
int pEVSL_Finish(pevsl_Data *pevsl_data) {

  if (pevsl_data->Amv) {
    PEVSL_FREE(pevsl_data->Amv);
  }
  if (pevsl_data->Bmv) {
    PEVSL_FREE(pevsl_data->Bmv);
  }
  if (pevsl_data->Bsol) {
    PEVSL_FREE(pevsl_data->Bsol);
  }
  if (pevsl_data->LTsol) {
    PEVSL_FREE(pevsl_data->LTsol);
  }
  if (pevsl_data->stats) {
    PEVSL_FREE(pevsl_data->stats);
  }

  PEVSL_FREE(pevsl_data);

  return 0;
}

/**
 * @brief Set the matrix A as Parcsr
 * @param[in,out] pevsl_data pevsl data struct
 * @param[in] A Matrix
 *
 * */
int pEVSL_SetAParcsr(pevsl_Data *pevsl_data, pevsl_Parcsr *A) {

  /* set N */
  pevsl_data->N = A->ncol_global;
  /* set n */
  pevsl_data->n = A->ncol_local;
  /* set nfirst */
  pevsl_data->nfirst = A->first_col;
  /* set A matvec */
  if (!pevsl_data->Amv) {
    PEVSL_CALLOC(pevsl_data->Amv, 1, pevsl_Matvec);
  }
  pevsl_data->Amv->func = pEVSL_ParcsrMatvec0;
  pevsl_data->Amv->data = (void *) A;

  return 0;
}

/**
 * @brief Set the B matrix as Parcsr
 * @param[in,out] pevsl_data pevsl data struct
 * @param[in] B Matrix
 *
 * */
int pEVSL_SetBParcsr(pevsl_Data *pevsl_data, pevsl_Parcsr *B) {

  /* set N */
  pevsl_data->N = B->ncol_global;
  /* set n */
  pevsl_data->n = B->ncol_local;
  /* set nfirst */
  pevsl_data->nfirst = B->first_col;
  /* set B matvec */
  if (!pevsl_data->Bmv) {
    PEVSL_CALLOC(pevsl_data->Bmv, 1, pevsl_Matvec);
  }
  pevsl_data->Bmv->func = pEVSL_ParcsrMatvec0;
  pevsl_data->Bmv->data = (void *) B;

  return 0;
}

/**
 * @brief Set problem sizes
 * @param[in,out] pevsl_data pevsl data struct
 * @param[in] N global size
 * @param[in] n local size
 * @param[in] nfirst first row/col
 * @warning if nfirst < 0, nfirst = PEVSL_NOT_DEFINED
 * */
int pEVSL_SetProbSizes(pevsl_Data *pevsl_data, int N, int n, int nfirst) {

  pevsl_data->N = N;
  pevsl_data->n = n;
  nfirst = nfirst < 0 ? PEVSL_NOT_DEFINED : nfirst;
  pevsl_data->nfirst = nfirst;

  return 0;
}



/**
 * @brief Set the user-input matvec routine and the associated data for A.
 * Save them in pevsl_data
 * @param[in, out] pevsl_data pEVSL data struct
 * @param[in] func Matrix vector function
 * @param[in] data data required for matvec
 * */
int pEVSL_SetAMatvec(pevsl_Data *pevsl_data, MVFunc func, void *data) {

  if (!pevsl_data->Amv) {
    PEVSL_CALLOC(pevsl_data->Amv, 1, pevsl_Matvec);
  }
  pevsl_data->Amv->func = func;
  pevsl_data->Amv->data = data;

  return 0;
}


/**
 * @brief Set the user-input matvec routine and the associated data for B.
 * Save them in pevsl_data
 * @param[in, out] pevsl_data pEVSL data struct
 * @param[in] func Matrix vector function
 * @param[in] data data required for matvec
 * */
int pEVSL_SetBMatvec(pevsl_Data *pevsl_data, MVFunc func, void *data) {

  if (!pevsl_data->Bmv) {
    PEVSL_CALLOC(pevsl_data->Bmv, 1, pevsl_Matvec);
  }
  pevsl_data->Bmv->func = func;
  pevsl_data->Bmv->data = data;

  return 0;
}

/**
 * @brief Set the solve routine and the associated data for B
 * @param[in, out] pevsl_data pEVSL data structure
 * @param[in] func Function for BSol
 * @param[in] data Data required for func
 * */
int pEVSL_SetBSol(pevsl_Data *pevsl_data, SVFunc func, void *data) {

  if (!pevsl_data->Bsol) {
    PEVSL_CALLOC(pevsl_data->Bsol, 1, pevsl_Bsol);
  }
  pevsl_data->Bsol->func = func;
  pevsl_data->Bsol->data = data;

  return 0;
}


/** JS 01/20/19 for complex Hermitian cases
 * @brief Set the user-input matvec routine and the associated data for A.
 * Save them in pevsl_data
 * @param[in, out] pevsl_data pEVSL data struct
 * @param[in] func Matrix vector function
 * @param[in] data data required for matvec
 * */
int pEVSL_SetZAMatvec(pevsl_Data *pevsl_data, ZMVFunc func, void *data) {

  /*if (!pevsl_data->ZAmv) {
    PEVSL_CALLOC(pevsl_data->ZAmv, 1, pevsl_ZMatvec);
  }*/
  PEVSL_CALLOC(pevsl_data->ZAmv, 1, pevsl_ZMatvec);
  pevsl_data->ZAmv->func = func;
  pevsl_data->ZAmv->data = data;

  return 0;
}


/**
 * @brief Set the user-input matvec routine and the associated data for B.
 * Save them in pevsl_data
 * @param[in, out] pevsl_data pEVSL data struct
 * @param[in] func Matrix vector function
 * @param[in] data data required for matvec
 * */
int pEVSL_SetZBMatvec(pevsl_Data *pevsl_data, ZMVFunc func, void *data) {

  /*if (!pevsl_data->ZBmv) {
    PEVSL_CALLOC(pevsl_data->ZBmv, 1, pevsl_ZMatvec);
  }*/
  PEVSL_CALLOC(pevsl_data->ZBmv, 1, pevsl_ZMatvec);
  pevsl_data->ZBmv->func = func;
  pevsl_data->ZBmv->data = data;

  return 0;
}


/**
 * @brief Set the solve routine and the associated data for B
 * @param[in, out] pevsl_data pEVSL data structure
 * @param[in] func Function for BSol
 * @param[in] data Data required for func
 * */
int pEVSL_SetZBSol(pevsl_Data *pevsl_data, ZSVFunc func, void *data) {

  /*if (!pevsl_data->ZBsol) {
    PEVSL_CALLOC(pevsl_data->ZBsol, 1, pevsl_ZBsol);
  }*/
  PEVSL_CALLOC(pevsl_data->ZBsol, 1, pevsl_ZBsol);
  pevsl_data->ZBsol->func = func;
  pevsl_data->ZBsol->data = data;

  return 0;
}


/**
 * @brief Set the solve routine and the associated data for LT
 * @param[in, out] pevsl_data pEVSL data structure
 * @param[in] func Function for LTSol
 * @param[in] data Data required for func
 * */
int pEVSL_SetLTSol(pevsl_Data *pevsl_data, SVFunc func, void *data) {

  if (!pevsl_data->LTsol) {
    PEVSL_CALLOC(pevsl_data->LTsol, 1, pevsl_Ltsol);
  }
  pevsl_data->LTsol->func = func;
  pevsl_data->LTsol->data = data;

  return 0;
}

/**
 * @brief Set the problem to standard eigenvalue problem
 * @param[in, out] pevsl_data pEVSL data struct
 * */
int pEVSL_SetStdEig(pevsl_Data *pevsl_data) {

  pevsl_data->ifGenEv = 0;

  return 0;
}

/**
 * @brief Set the problem to generalized eigenvalue problem
 * @param[in, out] pevsl_data pEVSL data struct
 *
 * */
int pEVSL_SetGenEig(pevsl_Data *pevsl_data) {

  pevsl_data->ifGenEv = 1;

  return 0;
}


int pEVSL_SetSigmaMult(pevsl_Data *pevsl_data, double mult) {

  pevsl_data->sigma_mult = mult;

  return 0;
}

