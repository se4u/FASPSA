/*=========================================================
 * rank_two_update_v2_fast.c -
 *
 * C = matrixMultiply(Bbar_by_a, b, delta_tilda_k, delta_k, tmp_v1, tmp_v2)
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
void init_to_zero(double* ptr, int n){
  mwSignedIndex i;
  for (i=0; i < n; i++){
    ptr[i] = 0;
  }
}
void print_arr(double* ptr, int n, char* name){
  mwSignedIndex i;
  for (i=0; i < n; i++){
    mexPrintf("\n i %d", i);
    mexPrintf(name);
    mexPrintf("[i] %f", ptr[i]);
  }
  fflush(stdout);
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Allocate space.
  double *Bbar, *delta_tilda_k, *delta_k;
  double *tmp_stage, *tmp_stage_b;
  double a, b;
  mwSignedIndex p;
  mwSignedIndex i;
  #ifdef SUPERSAFE
  mxArray* Bbar_copy = mxDuplicateArray(prhs[0]); // Safe
  #else
  mxArray* Bbar_copy = (prhs[0]); // Fast and Dangerous
  #endif
  // Read Input
  p = mxGetM(prhs[0]);
  Bbar = mxGetPr(Bbar_copy);
  a = *(mxGetPr(prhs[1]));
  b = *(mxGetPr(prhs[2]));
  delta_tilda_k = mxGetPr(prhs[3]);
  delta_k = mxGetPr(prhs[4]);
  /* mexPrintf("a %f b %f Bbar(1,1) %f Bbar(1,3) %f", a, b, Bbar[0], Bbar[2]); */
  /* print_arr(delta_tilda_k, p, "delta_tilda_k"); */
  /* print_arr(delta_k, p, "delta_k"); */

  #ifdef SAFE
  tmp_stage = mxCalloc(p , sizeof(double));
  tmp_stage_b = mxCalloc(p , sizeof(double));
  #else
  tmp_stage = mxGetPr(prhs[5]);
  tmp_stage_b = mxGetPr(prhs[6]);
  #endif
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
  plhs[0] = Bbar_copy;

  #ifdef SAFE
  mxFree(tmp_stage);
  mxFree(tmp_stage_b);
  #endif
}
