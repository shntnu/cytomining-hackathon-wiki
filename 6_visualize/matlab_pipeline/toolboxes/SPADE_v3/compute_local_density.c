#include "stdio.h"
#include "math.h"
#include "matrix.h"
#include "mex.h"   /* This one is required */


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
/* calculate the local density of each point in prhs[0], neighborhood size is defined in prhs[1], optimization parameter is in prhs[2] */

    int N, dim, i, j, d;
    double *data, *data_point_i, *data_point_j, *data_point_j_tmp;
    double *kernel_width, *optimization_para;
    double *local_density_j, *local_density_i, *local_density;
    double *dist_tmp;
    double **optimized_LD_index;
    int optimized_LD_index_counter;
    
    /* input variables */
    N = mxGetN(prhs[0]);                /* N : number of points */
    dim = mxGetM(prhs[0]);              /* d : dimension */
    data = mxGetPr(prhs[0]);            /* pointer that points to matrix that has all the data points, each column */
    kernel_width = mxGetPr(prhs[1]);
    optimization_para = mxGetPr(prhs[2]);
    /* output variables */
    plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);
    local_density = mxGetPr(plhs[0]);
    local_density_j = mxGetPr(plhs[0]);
    for(j=0;j<N;j++) {*(local_density_j++) = 0;}  /* initialize local density of all cells to be 0 */
    
    /* an intermediate variable */
    optimized_LD_index = (double**)mxMalloc(sizeof(double*)*N);
    dist_tmp = (double*)mxMalloc(sizeof(double));
    
    data = mxGetPr(prhs[0]);
    local_density_j = local_density;
    data_point_j = data;
    mexPrintf ("%4d%%", 0);
	mexEvalString("drawnow;");
    for(j=0;j<N;j++){
        if (*local_density_j == 0) {
            /* calculate the distance between point j and all points */
            data_point_i = data;
            local_density_i = local_density;
            optimized_LD_index_counter=0;
            for(i=0;i<N;i++){
                *dist_tmp = 0;
                data_point_j_tmp = data_point_j;
                for(d=0;d<dim;d++){
                   *dist_tmp += fabs(*(data_point_j_tmp++)-*(data_point_i++));
                }
                if (*dist_tmp<*kernel_width)
                    *local_density_j += 1;
                if (*dist_tmp<*optimization_para){
                    *(optimized_LD_index++)=local_density_i;
                    optimized_LD_index_counter++;
                }
                local_density_i++;
            }
            /* use the optimization parameter */
            for(i=0;i<optimized_LD_index_counter;i++){
                **(--optimized_LD_index) = *local_density_j;
            }
            local_density_j++;
            data_point_j = data_point_j_tmp;
        }
        else
        {
            local_density_j++;
            data_point_j += dim; /* move to the next data point */
        }
		if ((j % 500) == 1)
		{
			mexPrintf ("\b\b\b\b\b%4d%%", (j * 100) / N);
			mexEvalString("drawnow;");
		}
    }
	mexPrintf ("\b\b\b\b\b%4d%%", 100);
	mexEvalString("drawnow;");
    mxFree(optimized_LD_index);
    mxFree(dist_tmp);
}
        