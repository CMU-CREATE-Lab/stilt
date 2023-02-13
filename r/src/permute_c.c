#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>

double *footprint;
int footprint_nrows, footprint_ncols;

void create_footprint(int *nrows_in, int *ncols_in) {
    footprint_nrows = *nrows_in;
    footprint_ncols = *ncols_in;
    if (footprint) free(footprint);
    footprint = calloc(footprint_nrows*footprint_ncols, sizeof(double));
    printf("Created footprint size %d x %d at %p\n", footprint_nrows, footprint_ncols, footprint);
}

void permute(double *k, int *nkx_in, int *nky_in, // kernel and its dimensions
             int *nl_in, int *lai, int *loi, double *l_weights) // number of locations to write and the locations
{
    int nkx = *nkx_in;
    int nky = *nky_in;
    int nl = *nl_in;

    // TODO: try clipping this to the correct area to make sure we never overwrite the buffer

    for (int i = 0; i < nky; i++) { // loop over kernel Y
        int ys = i - (nky + 1) / 2; // subtract half the width of the kernel
        for (int j = 0; j < nkx; j++) { // loop over kernel X
            int xs = j - (nkx + 1) / 2;
            double ks = k[i+j*nky]; // col-major [r,c] = r+c*nrows
            for (int n = 0; n < nl; n++) {
                int row = loi[n]+ys;
                int col = lai[n]+xs;
                footprint[row + col * footprint_nrows] += l_weights[n] * ks; // col-major [r,c] = r+c*nrows
            }
        }
    }
}

SEXP read_footprint()
{
    double sum = 0;
    int len = footprint_nrows * footprint_ncols;
    SEXP out = PROTECT(allocVector(REALSXP, len));
    for(int i=0; i<len; ++i){
        REAL(out)[i]=footprint[i];
        sum += footprint[i];
    }
    UNPROTECT(1);
    printf("read_footprint returning sum %lf\n", sum);
    return out;
}