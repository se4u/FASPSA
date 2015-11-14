/*=========================================================
 * rank_two_update_v2_fast.c -
 *
 * C = matrixMultiply(Bbar_by_a, b, delta_tilda_k, delta_k)
 * Compile using
 * mex('-largeArrayDims', '-lmwblas', 'rank_two_update_v2_fast.c')
 * USAGE
 * rank_two_update_v2_fast(Bbar
 *=======================================================*/

#include "mex.h"
#include "blas.h"
// Constants
double one = 1.0;
double zero = 0.0;
mwSignedIndex inc = 1;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Allocate space.
  double *Bbar, *delta_tilda_k, *delta_k, *tmp_stage, *Binv, *tmp_stage_b;
  double a, b;
  mwSignedIndex p;

  // Read Input
  Bbar = mxGetPr(prhs[0]);
  a = *(mxGetPr(prhs[1]));
  b = *(mxGetPr(prhs[2]));
  delta_tilda_k = mxGetPr(prhs[2]);
  delta_k = mxGetPr(prhs[3]);
  tmp_stage = mxGetPr(prhs[4]);
  Binv = mxGetPr(prhs[5]);
  tmp_stage_b = mxGetPr(prhs[6]);
  p = (mwSignedIndex)mxGetM(prhs[0]);

  /* Stage 1: Bbar is symmetric.
   */
  // tmp_stage = Bbar * delta_tilda_k;
  dsymv("U", &p, &one, Bbar, &p, delta_tilda_k, &inc, &zero, tmp_stage, &inc);

  // tmp_deno = a + b * (delta_k' * tmp_stage);
  double tmp_deno = a + b * ddot(&p, delta_k, &inc, tmp_stage, &inc);

  // tmp_stage_b = Bbar' * delta_k; == Bbar * delta_k (since Bbar is symmetric)
  dsymv("U", &p, &one, Bbar, &p, delta_k, &inc, &zero, tmp_stage_b, &inc);

  // Bbar = Bbar + (-b/tmp_deno) * tmp_stage * tmp_stage_b';
  double b_by_tmp = -b / tmp_deno;
  dger(&p, &p, &b_by_tmp, tmp_stage, &inc, tmp_stage_b, &inc, Bbar, &p);

  /* Stage 2: Bbar is no longer symmetric.
     Need to use dgemv instead of dsymv.
   */
  // tmp_stage = Bbar * delta_k;
  dgemv("N", &p, &p, &one, Bbar, &p, delta_k, &inc, &zero, tmp_stage, &inc);

  // tmp_deno = a + b * (delta_tilda_k' * tmp_stage);
  tmp_deno = a + b * ddot(&p, delta_tilda_k, &inc, tmp_stage, &inc);

  // tmp_stage_b = Bbar' * delta_tilda_k;
  dgemv("T", &p, &p, &one, Bbar, &p, delta_tilda_k, &inc, &zero, tmp_stage_b, &inc);

  // Bbar = Bbar + (-b/tmp_deno) * tmp_stage * tmp_stage_b';
  b_by_tmp = -b / tmp_deno;
  dger(&p, &p, &b_by_tmp, tmp_stage, &inc, tmp_stage_b, &inc, Bbar, &p);

  // In place edit Bbar and return the value.
  plhs[0] = prhs[0];
}
