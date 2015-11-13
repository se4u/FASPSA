/*=========================================================
 * rank_two_update.c - 
 *
 * C = matrixMultiply(Bbar_by_a, b, delta_tilda_k, delta_k) 
 *
 
 *=======================================================*/

#if !defined(_WIN32)
#define dgemm dgemm_
#endif
#define ONE 1
#include "mex.h"
#include "blas.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Read Input */
    double *Bbar_by_a, *delta_tilda_k, *delta_k; 
    mwSignedIndex p;      
    double *pb = mxGetPr(prhs[1]);
    double b = *pb;
    Bbar_by_a = mxGetPr(prhs[0]); 
    delta_tilda_k = mxGetPr(prhs[2]); 
    delta_k = mxGetPr(prhs[3]);
    p = (mwSignedIndex)mxGetM(prhs[0]);
    
    /* form of op(A) & op(B) to use in matrix multiplication */
    char *chn = "N";
    double one = 1.0, zero = 0.0;    
    
    /* create intermediate vector and matrix */
    double *tmp_state = mxCreateDoubleMatrix(p, ONE, mxREAL);
    double *Binv = mxCreateDoubleMatrix(p, p, mxREAL);
    
    /* tmp_stage_1 = Bbar_by_a * delta_tilda_k; */
    dsymv('U', &p, &one, Bbar_by_a, &p, delta_tilda_k, &one, &zero, tmp_state, &one);

    /* stage1_deno = 1 + b * (delta_k' * tmp_stage_1); */
    double stage1_deno = 1 + b * ddot(&p, delta_k, &one, tmp_state, &one);
    
    /* Binv = Bbar_by_a - tmp_stage_1 * ((delta_k' *(b/stage1_deno)) * Bbar_by_a);*/
    double *bbar_delta_k = mxCreateDoubleMatrix(p, ONE, mxREAL);
    
    
    
    C = mxGetPr(plhs[0]);

    /* Pass arguments to Fortran by reference */
    dgemm(chn, chn, &m, &n, &p, &one, A, &m, B, &p, &zero, C, &m);
}
