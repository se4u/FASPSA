/*=========================================================
 * quartic_loss_fast.c -
 *
 * C = matrixMultiply(Bbar_by_a, b, delta_tilda_k, delta_k, tmp_v1, tmp_v2)
 * Compile using
 * mex('quartic_loss_fast.c')
 * USAGE
 * quartic_loss_fast(Bbar
 *=======================================================*/

#include "mex.h"
// #include "blas.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Read Input
  double* t = mxGetPr(prhs[0]);
  mwSignedIndex p = mxGetM(prhs[0]);
  if (mxGetN(prhs[0]) != 1){
    mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
                      "Input is not a vector.");
  }


  // Allocate the accumulators.
  // Since the quartic loss function uses a special triangular matrix that
  // contains all ones therefore B*t is just a cumulative sum of the entries
  // of t. acc contains that cumulative sum of entries.
  double acc=0;
  // Also note, we don't need to actually allocate storage for a vector since
  // The loss function is just
  // 1 * sum(Bt.^2) + 0.1 * sum(Bt.^3) + 0.01 * sum(Bt.^4)
  // val contains the final value returned.
  double val=0;
  // Bt_i, Bt_i_sq contain Bt[i] and Bt[i]^2 respectively.
  double Bt_i, Bt_i_sq;

  // The main loop.
  for (int i=p-1; i >= 0; i--){
    Bt_i = (t[i] + acc)/p;
    Bt_i_sq = Bt_i*Bt_i;
    val +=  Bt_i_sq * 1 \
            + (Bt_i * Bt_i_sq * 0.1) \
            + (Bt_i_sq * Bt_i_sq * 0.01);
    acc += t[i];
  }
  plhs[0] = mxCreateDoubleScalar(val);
}
