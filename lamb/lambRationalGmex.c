/* lambRationalGmex.c  Inner loop of lambDepthCorr>rationalG (numerical law), in C.
 *   g = lambRationalGmex(Pt, P2t, At, Et, id, K, y)
 *   Pt, P2t, At: L x nd complex (pole l, dk entry j), P = sqrt(dk^2 - s_l), P2 = P^2, A = r_l/(2P).
 *   Et: L x nd complex e^{-P y} for a common y, or [] (then computed per query).
 *   id: nq x 1 double, dk entry (1-based) of each query;  K: nq x 1 complex;  y: nq x 1 real.
 *   g(q) = sum_l A [P + iK + (1 - e^{-Py} e^{-iKy}) (P - iK)]/(P2 + K^2), series where |(P + iK) y| < 1e-4.
 * Build: mex -R2018a -O lambRationalGmex.c   (interleaved complex API); without an accepted Xcode
 *   licence, clang directly (ML = matlabroot):
 *   clang -O3 -fcx-limited-range -bundle -fPIC -DMATLAB_DEFAULT_RELEASE=R2018a -DMATLAB_MEX_FILE -I$ML/extern/include \
 *     lambRationalGmex.c $ML/extern/version/c_mexapi_version.c -L$ML/bin/maca64 -lmx -lmex \
 *     -Wl,-exported_symbols_list,$ML/extern/lib/maca64/c_exportsmexfileversion.map -o lambRationalGmex.mexmaca64
 */
#include "mex.h"
#include <complex.h>

static inline double complex cdiv(double complex a, double complex b)  /* a/b, finite b: no __divdc3 */
{
    double br = creal(b), bi = cimag(b), d = br*br + bi*bi;
    return (a*(br - I*bi))/d;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 7) mexErrMsgIdAndTxt("lambRationalGmex:nrhs", "7 inputs: Pt, P2t, At, Et, id, K, y");
    for (int i = 0; i < 3; i++)
        if (!mxIsComplex(prhs[i])) mexErrMsgIdAndTxt("lambRationalGmex:type", "Pt, P2t, At must be complex");
    if (!mxIsComplex(prhs[5])) mexErrMsgIdAndTxt("lambRationalGmex:type", "K must be complex");
    mwSize L = mxGetM(prhs[0]), nd = mxGetN(prhs[0]), nq = mxGetNumberOfElements(prhs[5]);
    const double complex *P = (const double complex *)mxGetComplexDoubles(prhs[0]);
    const double complex *P2 = (const double complex *)mxGetComplexDoubles(prhs[1]);
    const double complex *A = (const double complex *)mxGetComplexDoubles(prhs[2]);
    const double complex *E = NULL;
    if (!mxIsEmpty(prhs[3])) {
        if (!mxIsComplex(prhs[3])) mexErrMsgIdAndTxt("lambRationalGmex:type", "Et must be complex or []");
        E = (const double complex *)mxGetComplexDoubles(prhs[3]);
    }
    const double *id = mxGetDoubles(prhs[4]);
    const double complex *K = (const double complex *)mxGetComplexDoubles(prhs[5]);
    const double *y = mxGetDoubles(prhs[6]);
    if (mxGetNumberOfElements(prhs[4]) != nq || mxGetNumberOfElements(prhs[6]) != nq)
        mexErrMsgIdAndTxt("lambRationalGmex:size", "id, K, y must have the same length");
    plhs[0] = mxCreateDoubleMatrix(nq, 1, mxCOMPLEX);
    double complex *g = (double complex *)mxGetComplexDoubles(plhs[0]);
    for (mwSize q = 0; q < nq; q++) {
        mwSize j = (mwSize)id[q] - 1;
        if (j >= nd) mexErrMsgIdAndTxt("lambRationalGmex:id", "id out of range");
        double complex k = K[q], iK = I*k, K2 = k*k, eK = cexp(-iK*y[q]);
        double yq = y[q], lim = 1e-8/(yq*yq);                       /* |t|^2 < lim <=> |t y| < 1e-4 */
        const double complex *Pj = P + j*L, *P2j = P2 + j*L, *Aj = A + j*L, *Ej = E ? E + j*L : NULL;
        double complex s = 0;
        for (mwSize l = 0; l < L; l++) {
            double complex p = Pj[l], t = p + iK, m = p - iK;
            double ta = creal(t)*creal(t) + cimag(t)*cimag(t);
            if (yq == 0 || ta < lim) {                              /* small (P + iK) y: series */
                double complex x = t*yq;
                s += Aj[l]*(cdiv(1.0, m) + yq*(1.0 - x/2.0 + x*x/6.0));
            } else {
                double complex e = Ej ? Ej[l] : cexp(-p*yq);
                s += cdiv(Aj[l]*(t + (1.0 - e*eK)*m), P2j[l] + K2);
            }
        }
        g[q] = s;
    }
}
