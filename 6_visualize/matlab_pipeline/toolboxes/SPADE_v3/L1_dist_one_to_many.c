/* L1_dist_one_to_many(data(:,1),data) */


#include "math.h"
#include "matrix.h"
#include "mex.h"   /* This one is required */
#define	MIN(A, B)	((A) < (B) ? (A) : (B))


double L1_dist(double* A, double*B, int dim)
{
    double dist;
    int i;
    dist = 0;
    for(i=0;i<dim;i++)
        dist += fabs(*(A++)-*(B++));
    return(dist);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
/* calculate the L1 distance between each item in prhs[0] and its nearest neighbor in prhs[1]     */
    
    int n, N, dim, i, j;
    double *ref_data, *ref_data_i;  
    double *new_data, *new_data_j;
    double *all_dist; 
    
    
    n = mxGetN(prhs[0]);    /* n : number of new points */
    N = mxGetN(prhs[1]);    /* N : number of ref points */
    dim = mxGetM(prhs[0]);  /* d : dimension */
    new_data = mxGetPr(prhs[0]);  /* pointer that points to matrix that has the new poitns, each col is a point */
    ref_data = mxGetPr(prhs[1]);  /* pointer to matrix that has all the ref points, each col is a point */

    plhs[0] = mxCreateDoubleMatrix(N, n, mxREAL);
    all_dist = mxGetPr(plhs[0]);
    
   
    for(j=0;j<n;j++){
        new_data_j = &new_data[j*dim];
        for(i=0;i<N;i++){
            *(all_dist++) = L1_dist(new_data_j, &ref_data[i*dim], dim);
        }
    }
}