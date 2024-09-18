#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <assert.h>

double *footprint;
int footprint_nrows, footprint_ncols;

void create_footprint(int *nrows_in, int *ncols_in) {
    footprint_nrows = *nrows_in;
    footprint_ncols = *ncols_in;
    if (footprint) free(footprint);
    footprint = calloc(footprint_nrows*footprint_ncols, sizeof(double));
    //fprintf(stderr, "Created footprint size %d x %d at %p\n", footprint_nrows, footprint_ncols, footprint);
}

void permute(double *sigma_in, int *nkx_in, int *nky_in, // kernel and its dimensions
             int *nl_in, int *lai, int *loi, double *l_weights) // number of locations to write and the locations
{
    int nkx = *nkx_in;
    int nky = *nky_in;
    double sigma = *sigma_in;
    //fprintf(stderr, "Computing 1-d kernel with sigma %lg into array of size %d\n", sigma, nkx);
    assert(nkx == nky);
    assert(nkx % 2 == 1);  // always odd dimension so the gaussian is centered on a sample
    int nl = *nl_in;

    double la_ctr = 0, lo_ctr = 0;
    for (int i = 0; i < nl; i++) {
        la_ctr += lai[i];
        lo_ctr += loi[i];
    }
    la_ctr /= nl;
    lo_ctr /= nl;

    //fprintf(stderr, "permute called with sigma %lg, center of mass %lg,%lg, nkx %d, nky %d, nl %d\n", sigma, la_ctr, lo_ctr, nkx, nky, nl);
    double gaussian[nkx];
    int kernel_midpoint = nkx/2; // midpoint is 2,2 for 5x5;  3,3 for 7x7, etc
    for (int x = 0; x <= nkx/2; x++) { // x is dist from center of gaussian
        gaussian[kernel_midpoint + x] = exp(-(x * x) / (2 * sigma * sigma));
        gaussian[kernel_midpoint - x] = gaussian[kernel_midpoint + x];
    }

    // Normalize sum to 1
    double sum = 0;
    for (int i = 0; i < nkx; i++) {
        sum += gaussian[i];
    }
    for (int i = 0; i < nkx; i++) {
        gaussian[i] /= sum;
    }

    // Compute 1-d kernel of sigma *sigma_in

    for (int i = 0; i < nky; i++) { // loop over kernel Y
        int ys = i - (nky + 1) / 2; // subtract half the width of the kernel
        for (int j = 0; j < nkx; j++) { // loop over kernel X
            int xs = j - (nkx + 1) / 2;
            double ks = gaussian[i] * gaussian[j]; // 2d kernel comes from multplying the 1d kernels
            for (int n = 0; n < nl; n++) {
                int row = loi[n]+ys;
                int col = lai[n]+xs;

                // Boundary checks
                if (row < 0 || row >= footprint_nrows || col < 0 || col >= footprint_ncols) {
                    continue; // Skip out-of-bounds indices
                }

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
    //printf("read_footprint returning sum %lf\n", sum);
    return out;
}