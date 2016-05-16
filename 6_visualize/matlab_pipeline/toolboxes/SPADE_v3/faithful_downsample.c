
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
/* calculate the L1 distance between each item in prhs[0] and its nearest neighbor in prhs[1]    */
    
    int N, dim, i, j;
    double *data, *second_input;  
    double kernel_width;
    double *data_j;
    double *weights, *group_assign;
    double dist_tmp; 
    
    
    N = mxGetN(prhs[0]);    /* n : number of new points */
    dim = mxGetM(prhs[0]);  /* d : dimension */
    data = mxGetPr(prhs[0]);  /* pointer that points to data matrix, each col is a point */
    second_input = mxGetPr(prhs[1]);  
    kernel_width = *second_input; 
    

    plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);
    weights = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1, N, mxREAL);
    group_assign = mxGetPr(plhs[1]);
    
    for(i=0;i<N;i++)
    {
        weights[i]=1;
        group_assign[i]=0;
    }
    
    mexPrintf ("%4d%%", 0);
    mexEvalString("drawnow;");
    
    for(j=0;j<N;j++){
        if (weights[j]!=0)
        {
            data_j = &data[j*dim];
            group_assign[j]=j+1;
            for(i=0;i<j;i++){
                if(weights[i]!=0 && L1_dist(data_j, &data[i*dim], dim) < kernel_width)
                {
                    weights[j] += weights[i]; weights[i]=0;
                    group_assign[i] = j+1;
                }
            }
            for(i=j+1;i<N;i++){
                if(weights[i]!=0 && L1_dist(data_j, &data[i*dim], dim) < kernel_width)
                {
                    weights[j] += weights[i]; weights[i]=0;
                    group_assign[i] = j+1;
                }
            }
        }
        if ((j % 500) == 1)
		{
			mexPrintf ("\b\b\b\b\b%4d%%", (j * 100) / N);
			mexEvalString("drawnow;");
		}
    }
    mexPrintf ("\b\b\b\b\b%4d%%\n", 100);
    mexEvalString("drawnow;");
}