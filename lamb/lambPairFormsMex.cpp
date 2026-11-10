// lambPairFormsMex.cpp  lambRandomGram>pairForms with the rational g of lambDepthCorr (numerical law), in C++.
//   Q = lambPairFormsMex(A, e, Bs, Kh, dh, zz, kap, k, phi, s_fit, r_fit)
//   A: N x 9 x npc complex;  e: N x npc complex (IW: 1 x nt);  Bs: 9 x nph x 3 x nu complex;  Kh: nu complex;
//   dh, zz, k: real scalars;  kap: N real;  phi: nph real (2 pi (0:nph-1)/nph);  s_fit, r_fit: L complex (R.s, R.r).
//   Q(a, b, n) = 1/2 sum_{phi, j} T_a conj(T_b) [g_a + conj g_b]  (na x na x N, na = nu npc, as the MATLAB code),
//   T_(u,c)(phi, j, n) = A(n, :, c) . Bs(:, phi, j, u),  g_(u,c)(phi, n) = g(Kh(u) + Kc(n, c), |zz|, dk(n, phi)),
//   Kc = -i e dh (sign flipped for zz < 0). Kernel of cpp_rtmclib/bench_tabela/bench_pairforms.cpp
//   Threads: env LAMB_MEX_THREADS (default 1); run with VECLIB_MAXIMUM_THREADS=1.
// Build (ML = matlabroot; clang directly, mex needs an accepted Xcode licence). OpenMP headers from Homebrew,
// runtime = MATLAB's own $ML/bin/maca64/libomp.dylib (already in the process: no second libomp):
//   clang -O2 -c -DMATLAB_DEFAULT_RELEASE=R2018a -DMATLAB_MEX_FILE -I$ML/extern/include \
//     $ML/extern/version/c_mexapi_version.c -o c_mexapi_version.o
//   clang++ -std=c++17 -O3 -DNDEBUG -fno-fast-math -fcx-limited-range -mcpu=native -DACCELERATE_NEW_LAPACK \
//     -Xclang -fopenmp -I/opt/homebrew/opt/libomp/include -bundle -fPIC -DMATLAB_DEFAULT_RELEASE=R2018a \
//     -DMATLAB_MEX_FILE -I$ML/extern/include lambPairFormsMex.cpp c_mexapi_version.o \
//     -L$ML/bin/maca64 -lmx -lmex -lomp -Wl,-rpath,$ML/bin/maca64 -framework Accelerate \
//     -Wl,-exported_symbols_list,$ML/extern/lib/maca64/c_exportsmexfileversion.map -o lambPairFormsMex.mexmaca64
#include "mex.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <vector>
#include <omp.h>
#include <Accelerate/Accelerate.h>

using cplx = std::complex<double>;
static inline cplx cdiv(cplx a, cplx b) { double d = std::norm(b); return a*std::conj(b)/d; }

static const cplx* cin(const mxArray* a, const char* name) {
    if (!mxIsDouble(a) || !mxIsComplex(a)) mexErrMsgIdAndTxt("lambPairFormsMex:type", "%s must be complex double", name);
    return reinterpret_cast<const cplx*>(mxGetComplexDoubles(a));
}
static const double* rin(const mxArray* a, const char* name) {
    if (!mxIsDouble(a) || mxIsComplex(a)) mexErrMsgIdAndTxt("lambPairFormsMex:type", "%s must be real double", name);
    return mxGetDoubles(a);
}

void mexFunction(int, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 11) mexErrMsgIdAndTxt("lambPairFormsMex:nrhs", "11 inputs: A, e, Bs, Kh, dh, zz, kap, k, phi, s_fit, r_fit");
    const cplx *A = cin(prhs[0], "A"), *E = cin(prhs[1], "e"), *Bs = cin(prhs[2], "Bs"), *Kh = cin(prhs[3], "Kh"),
               *s = cin(prhs[9], "s_fit"), *r = cin(prhs[10], "r_fit");
    const double dh = mxGetScalar(prhs[4]), zz = mxGetScalar(prhs[5]), k = mxGetScalar(prhs[7]),
                 *kapv = rin(prhs[6], "kap"), *phi = rin(prhs[8], "phi");
    const mwSize* da = mxGetDimensions(prhs[0]);  const mwSize* db = mxGetDimensions(prhs[2]);
    const size_t N = da[0], npc = mxGetNumberOfDimensions(prhs[0]) > 2 ? da[2] : 1, nph = mxGetNumberOfElements(prhs[8]),
                 nu = mxGetNumberOfDimensions(prhs[2]) > 3 ? db[3] : 1, na = nu*npc, ncol = nph*3,
                 L = mxGetNumberOfElements(prhs[9]), nmir = nph/2 + 1;
    if (da[1] != 9 || db[0] != 9 || db[1] != nph || mxGetNumberOfElements(prhs[2]) != 9*nph*3*nu
        || mxGetNumberOfElements(prhs[1]) != N*npc || mxGetNumberOfElements(prhs[3]) != nu
        || mxGetNumberOfElements(prhs[6]) != N || mxGetNumberOfElements(prhs[10]) != L)
        mexErrMsgIdAndTxt("lambPairFormsMex:size", "inconsistent sizes");
    const mwSize dq[3] = {na, na, N};
    plhs[0] = mxCreateNumericArray(3, dq, mxDOUBLE_CLASS, mxCOMPLEX);
    cplx* Q = reinterpret_cast<cplx*>(mxGetComplexDoubles(plhs[0]));
    const char* env = std::getenv("LAMB_MEX_THREADS");
    const int nth = std::max(1, env ? std::atoi(env) : 1);
    const double y = std::abs(zz);
    auto Aat = [&](size_t n, size_t i, size_t c) { return A[n + N*(i + 9*c)]; };
    auto Bat = [&](size_t i, size_t ph, size_t j, size_t u) { return Bs[i + 9*(ph + nph*(j + 3*u))]; };
#pragma omp parallel num_threads(nth)
    {
        std::vector<cplx> T(na*ncol), GT(na*ncol), G(na*nmir), P(L), Al(L);
        std::vector<double> Pr(L), Pi(L), Qr(L), Qi(L), Ar(L), Ai(L), Er(L), Ei(L);
#pragma omp for schedule(static, 1)   // not dynamic: MATLAB's libomp lacks __kmpc_dispatch_deinit
        for (size_t n = 0; n < N; ++n) {
            const double kap = kapv[n];
            for (size_t c = 0; c < npc; ++c)                          // T(al, phi + nph j), al = u + nu c
                for (size_t u = 0; u < nu; ++u)
                    for (size_t col = 0; col < ncol; ++col) {
                        size_t ph = col % nph, j = col/nph;  cplx t = 0;
                        for (size_t i = 0; i < 9; ++i) t += Aat(n, i, c)*Bat(i, ph, j, u);
                        T[(u + nu*c) + na*col] = t;
                    }
            for (size_t m = 0; m < nmir; ++m) {                       // g on phi in [0, pi]
                double dk = std::sqrt(std::max(k*k + kap*kap - 2*k*kap*std::cos(phi[m]), 0.0));
                for (size_t l = 0; l < L; ++l) {
                    cplx p2 = dk*dk - s[l], p = std::sqrt(p2), al = r[l]/(2.0*p), e = std::exp(-p*y);
                    Pr[l] = p.real();  Pi[l] = p.imag();  Qr[l] = p2.real();  Qi[l] = p2.imag();
                    Ar[l] = al.real();  Ai[l] = al.imag();  Er[l] = e.real();  Ei[l] = e.imag();
                    P[l] = p;  Al[l] = al;
                }
                for (size_t c = 0; c < npc; ++c) {
                    cplx Kc = cplx(0, -1)*E[n + N*c]*dh;  if (zz < 0) Kc = -Kc;
                    for (size_t u = 0; u < nu; ++u) {
                        const cplx K = Kh[u] + Kc, eK = std::exp(cplx(0, -1)*K*y);
                        const double kr = -K.imag(), ki = K.real(), ekr = eK.real(), eki = eK.imag(), lim = 1e-8/(y*y);  // iK = kr + i ki
                        const cplx K2 = K*K;  const double k2r = K2.real(), k2i = K2.imag();
                        double sr = 0, si = 0;  bool ser = false;
                        // term = A (t + (1 - E eK) m)/(P^2 + K^2),  t = P + iK, m = P - iK: one real division
                        for (size_t l = 0; l < L; ++l) {
                            double tr = Pr[l] + kr, ti = Pi[l] + ki, mr = Pr[l] - kr, mi = Pi[l] - ki;
                            ser |= tr*tr + ti*ti < lim;
                            double wr = 1 - (Er[l]*ekr - Ei[l]*eki), wi = -(Er[l]*eki + Ei[l]*ekr);
                            double nr = tr + wr*mr - wi*mi, ni = ti + wr*mi + wi*mr;
                            double dr = Qr[l] + k2r, di = Qi[l] + k2i, dn = 1/(dr*dr + di*di);
                            double xr = (nr*dr + ni*di)*dn, xi = (ni*dr - nr*di)*dn;
                            sr += Ar[l]*xr - Ai[l]*xi;  si += Ar[l]*xi + Ai[l]*xr;
                        }
                        cplx acc(sr, si);
                        if (ser || y == 0) {                                    // rare: redo with the series where |t y| < 1e-4
                            acc = 0;  const cplx iK = cplx(0, 1)*K;
                            for (size_t l = 0; l < L; ++l) {
                                cplx t = P[l] + iK, mm = P[l] - iK;
                                if (y == 0 || std::norm(t) < lim) { cplx x = t*y;  acc += Al[l]*(cdiv(1.0, mm) + y*(1.0 - x/2.0 + x*x/6.0)); }
                                else acc += Al[l]*(cdiv(1.0, mm) + cdiv(1.0 - cplx(Er[l], Ei[l])*eK, t));
                            }
                        }
                        G[(u + nu*c) + na*m] = acc;
                    }
                }
            }
            cplx* Qn = Q + na*na*n;                                    // Q = (GT T' + T GT')/2: Hermitian rank-2k update
            for (size_t col = 0; col < ncol; ++col) {
                size_t ph = col % nph, m = std::min(ph, nph - ph);
                for (size_t a = 0; a < na; ++a) GT[a + na*col] = T[a + na*col]*G[a + na*m];
            }
            const cplx half(0.5, 0);
            cblas_zher2k(CblasColMajor, CblasUpper, CblasNoTrans, int(na), int(ncol), &half, GT.data(), int(na), T.data(), int(na), 0.0, Qn, int(na));
            for (size_t b = 0; b < na; ++b) for (size_t a = b + 1; a < na; ++a) Qn[a + na*b] = std::conj(Qn[b + na*a]);
        }
    }
}
