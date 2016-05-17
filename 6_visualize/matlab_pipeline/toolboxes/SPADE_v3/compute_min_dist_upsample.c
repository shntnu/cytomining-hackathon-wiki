#include "math.h"
#include "matrix.h"
#include "mex.h"   /* This one is required */
#define	MIN(A, B)	((A) < (B) ? (A) : (B))


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
/* calculate the L1 distance between each item in prhs[0] and its nearest neighbor in prhs[1]    */
    
    int n, N, dim, i, j, d;
    double *ref_data, *ref_data_i;  
    double *new_data, *new_data_j, *new_data_j_tmp, *min_dist;
    double *all_dist; 
    double dist_tmp;
    double *NN_index;
    
    
    n = mxGetN(prhs[0]);    /* n : number of new points  */
    N = mxGetN(prhs[1]);    /* N : number of ref points  */
    dim = mxGetM(prhs[0]);  /* d : dimension  */
    new_data = mxGetPr(prhs[0]);  /* pointer that points to matrix that has the new poitns, each col is a point  */
    ref_data = mxGetPr(prhs[1]);  /* pointer to matrix that has all the ref points, each col is a point  */

    plhs[0] = mxCreateDoubleMatrix(1, n, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, n, mxREAL);
    min_dist = mxGetPr(plhs[0]);
    NN_index = mxGetPr(plhs[1]);
    
    new_data_j = (double*)mxMalloc(sizeof(double)*dim);
    
    mexPrintf ("%4d%%", 0);
	mexEvalString("drawnow;");

    new_data_j = new_data;
    for(j=0;j<n;j++){
        *NN_index=0;                    /* for the j'th new point, the NN and min_dist has never been assigned yet */
        ref_data_i = ref_data;
        {   
            i=0;
            dist_tmp=0;
            for(d=0;d<dim;d++)
                dist_tmp += fabs(new_data_j[d] - *(ref_data_i++));
            *min_dist = dist_tmp;
            *NN_index = i+1;
            i++;
        }
        for(i=1;i<N;i++){                     /* visit the remaining points in the ref set, and see whether there is anyone closer  */
            dist_tmp=0;
            for(d=0;d<dim;d++)
                dist_tmp += fabs(new_data_j[d] - *(ref_data_i++));
            if(dist_tmp<*min_dist)
            {
                *min_dist = dist_tmp;
                *NN_index=i+1;
            }
        }
        min_dist++;
        NN_index++;
        new_data_j += dim;
   		if ((j % 100) == 1)
		{
			mexPrintf ("\b\b\b\b\b%4d%%", (j * 100) / n);
			mexEvalString("drawnow;");
		}
    }
   	mexPrintf ("\b\b\b\b\b%4d%%", 100);
	mexEvalString("drawnow;");
}


        