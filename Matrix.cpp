#include "stdafx.h"
#include "Matrix.h"
#pragma intrinsic (fabs, sqrt)

#define D_MACHINE_EPS   2.22044605e-016
#define F_MACHINE_EPS   1.1920929e-007f

void
Matrix::Print(int verbose)
{
    int i, j;

    if (pf == 0) {
        printf("[null]\n");
    } else {
        if (!verbose) {
            printf("[ ");
            for (i = 0; i < N && i < 8; i++)
                printf("%6.3f ", (*this)[i][i]);
            printf(" ... ]\n");
        } else {
            for (i = 0; i < M; i++) {
                printf("[ ");
                for (j = 0; j < N; j++)
                    printf("%6.3f ", (*this)[i][j]);
                printf("]\n");
            }
        }
        if (p) {
            printf("( ");
            for (i = 0; i < M; i++)
                printf("%6d ", p[i]);
            printf(" )\n");
        }
        if (q) {
            printf("( ");
            for (i = 0; i < N; i++)
                printf("%6d ", q[i]);
            printf(" )\n");
        }
    }
}
//
//  Naive      - 704.5 seconds
//  Transpose   - 73.1 seconds
//  LoopTile    - 63.3          TILE = 32 (L1 cache size)
//  LoopTile    - 65.5          TILE = 128 (L2 cache size)
//  LoopTile    - 69.4          TILE = 512 (L3 cache size)
//  Unwind 2    - 58.8
//  Unwind 4    - 55.5, 51.6
//  Unwind 8    - 57.6
//  No Trans    - 93.0
//  SSE         - 29.5          4x1 multiplcation using dp
//  SSE         - 28.0          4x1 using set1 and mul, more accurate
//  2xLoopTile  - 29.6          L1 and L2 cache tiling, not a win
//  AVX         - 29.0          no faster than SSE
//  SSE 4x4     - 19.0          full 4x4 matrix multiple in registers in inner loop
//  Skylake     - 11.8          new computer
//
//  Efficent matrix multiply for N divisible by 4, using cache blocking and SSE5 vector instructions
//
Matrix
Matrix::operator*(Matrix B)
{
    Matrix &A(*this);
    Matrix C(M, B.N);
    int i, j, k, iB, jB, kB, NB, TILE;
    __m128 rC0, rC1, rC2, rC3, rB0, rB1, rB2, rB3;

    if (N != B.M)
        return C.Delete();
    if (C.M != C.N || N & 0x3 || M & 0x3) {
        for (i = 0; i < C.M; i++) {
            for (j = 0; j < C.N; j++) {
                C[i][j] = 0.0;
                for (k = 0; k < B.M; k++)
                    C[i][j] += A[i][k] * B[k][j];
            }
        }
        return C;
    }
    NB = A.N;
    TILE = 32;
    for (i = 0; i < NB; i++) {
        for (j = 0; j < NB; j++) {
            C[i][j] = 0.0;
        }
    }
    for (iB = 0; iB < NB; iB += TILE) {
        for (jB = 0; jB < NB; jB += TILE) {
            for (kB = 0; kB < NB; kB += TILE) {
                for (i = iB; i < iB+TILE && i < NB; i+=4) {
                    for (j = jB; j < jB+TILE && j < NB; j+=4) {
                        rC0 = _mm_load_ps(C[i  ]+j);
                        rC1 = _mm_load_ps(C[i+1]+j);
                        rC2 = _mm_load_ps(C[i+2]+j);
                        rC3 = _mm_load_ps(C[i+3]+j);
                        for (k = kB; k < kB+TILE && k < NB; k+=4) {
                            rB0 = _mm_load_ps(B[k  ]+j);
                            rB1 = _mm_load_ps(B[k+1]+j);
                            rB2 = _mm_load_ps(B[k+2]+j);
                            rB3 = _mm_load_ps(B[k+3]+j);
                            rC0 = _mm_add_ps(rC0, _mm_mul_ps(rB0, _mm_set1_ps(A[i  ][k  ])));
                            rC0 = _mm_add_ps(rC0, _mm_mul_ps(rB1, _mm_set1_ps(A[i  ][k+1])));
                            rC0 = _mm_add_ps(rC0, _mm_mul_ps(rB2, _mm_set1_ps(A[i  ][k+2])));
                            rC0 = _mm_add_ps(rC0, _mm_mul_ps(rB3, _mm_set1_ps(A[i  ][k+3])));

                            rC1 = _mm_add_ps(rC1, _mm_mul_ps(rB0, _mm_set1_ps(A[i+1][k  ])));
                            rC1 = _mm_add_ps(rC1, _mm_mul_ps(rB1, _mm_set1_ps(A[i+1][k+1])));
                            rC1 = _mm_add_ps(rC1, _mm_mul_ps(rB2, _mm_set1_ps(A[i+1][k+2])));
                            rC1 = _mm_add_ps(rC1, _mm_mul_ps(rB3, _mm_set1_ps(A[i+1][k+3])));

                            rC2 = _mm_add_ps(rC2, _mm_mul_ps(rB0, _mm_set1_ps(A[i+2][k  ])));
                            rC2 = _mm_add_ps(rC2, _mm_mul_ps(rB1, _mm_set1_ps(A[i+2][k+1])));
                            rC2 = _mm_add_ps(rC2, _mm_mul_ps(rB2, _mm_set1_ps(A[i+2][k+2])));
                            rC2 = _mm_add_ps(rC2, _mm_mul_ps(rB3, _mm_set1_ps(A[i+2][k+3])));

                            rC3 = _mm_add_ps(rC3, _mm_mul_ps(rB0, _mm_set1_ps(A[i+3][k  ])));
                            rC3 = _mm_add_ps(rC3, _mm_mul_ps(rB1, _mm_set1_ps(A[i+3][k+1])));
                            rC3 = _mm_add_ps(rC3, _mm_mul_ps(rB2, _mm_set1_ps(A[i+3][k+2])));
                            rC3 = _mm_add_ps(rC3, _mm_mul_ps(rB3, _mm_set1_ps(A[i+3][k+3])));
                        }
                        _mm_store_ps(C[i  ]+j, rC0);
                        _mm_store_ps(C[i+1]+j, rC1);
                        _mm_store_ps(C[i+2]+j, rC2);
                        _mm_store_ps(C[i+3]+j, rC3);
                    }
                }
            }
        }
    }
    return C;
}
//
//  Gaussian Elimiation to put matrix into LU form (in-place)
//
inline static void
ExchangeRows(int iRow1, int iRow2, Matrix mat)
{
    float f;
    int j;

    for (j = 0; j < mat.N; j++) {
        f = mat[iRow1][j];  mat[iRow1][j] = mat[iRow2][j];  mat[iRow2][j] = f;
    }
}

inline static void
Row(int iRow, float f, Matrix mat)
{
    int j;

    for (j = 0; j < mat.N; j++)
        mat[iRow][j] *= f;
}

inline static void
Addd(int iRow1, float f, int iRow2, Matrix mat)
{
    int j;

    for (j = 0; j < mat.N; j++)
        mat[iRow2][j] += f*mat[iRow1][j];
}

double ML_fMatrixCondition;  // internal report from GECP and QRCP

Matrix
Matrix::GaussianLU()
{
    int i, j, k, n, iPivot, jPivot, *p, *q;
    float fPivot, fTmp, f, fLastPivot, fEps; // *AI, *AJ;
    Matrix mat(M, N);

    mat.p = p = new int[M];
    mat.q = q = new int[N];
    n = N;
    for (i = 0; i < M * N; i++)
        mat.pf[i] = pf[i];
    for (k = 0; k < n; k++)             // Start with identity permutations of rows and columns
        p[k] = q[k] = k;
    fEps = 8.0f*n*F_MACHINE_EPS;
    for (k = 0; k < n - 1; k++) {       //REVIEW: if non-square, loop to k < n
        //
        //  Find maximum element of submatrix
        //
        fPivot = 0.0;
        for (i = k; i < n; i++) {
            for (j = k; j < n; j++) {
                if ((f = fabs(mat[i][j])) > fPivot) {
                    fPivot = f; 
                    iPivot = i;
                    jPivot = j;
                }
            }
        }
        //
        //  Rank determination
        //
        if (fPivot == 0.0) {
            ML_fMatrixCondition = fPivot;
            mat.kRank = k;
            return mat;
        }
        if (k == 0) {
            fLastPivot = fPivot;
        } else if ( fPivot < fLastPivot*fEps) {
            ML_fMatrixCondition = fPivot/fLastPivot;
            mat.kRank = k;
            return mat;
        }
        //
        //  Exchange rows and columns to move pivot to k,k
        //
        p[k] = iPivot;
        q[k] = jPivot;
        if (jPivot != k) {
            for (i = 0; i < n; i++) {
                fTmp = mat[i][k];
                mat[i][k] = mat[i][jPivot];
                mat[i][jPivot] = fTmp;
            }
        }
        if (iPivot != k) {
            for (j = 0; j < n; j++) {
                fTmp = mat[k][j];
                mat[k][j] = mat[iPivot][j];
                mat[iPivot][j] = fTmp;
            }
        }
        //
        //  Elimination
        //
        f = 1.0f/mat[k][k];
        for (i = k + 1; i < n; i++) {
            mat[i][k] *= f;
            for (j = k + 1; j < n; j++)
                mat[i][j] -= mat[i][k] * mat[k][j];
        }
    }
    fPivot = fabs(mat[k][k]);
    if (k > 0)
        ML_fMatrixCondition = fPivot / fLastPivot;
    else
        ML_fMatrixCondition = fPivot;
    if (fPivot == 0.0 || (k > 0 && fPivot < fLastPivot*fEps))
        mat.kRank = k;   // rank = n - 1
    else
        mat.kRank = k+1; // full rank = n
    return mat;
}

Matrix
Matrix::SolveLU(Matrix &v)
{
    int i, j, n;
    float f;
    Matrix LU(0), x(v.M, v.N);

    LU = *this;
    n = LU.N;
    for (i = 0; i < x.M*x.N; i++)
        x.pf[i] = v.pf[i];
    for (i = 0; i < n - 1; i++) {   // Permuation P
        f = x[p[i]][0];
        x[p[i]][0] = x[i][0];
        x[i][0] = f;
    }
    for (i = 1; i < n; i++) {       // Forward Substitution
        f = x[i][0];
        for (j = 0; j < i; j++)
            f -= LU[i][j] * x[j][0];      // LU[n*i + j]
        x[i][0] = f;
    }
    x[n-1][0] /= LU[n - 1][n - 1];
    for (i = n - 2; i >= 0; --i) {    // Back Substitution
        f = x[i][0];
        for (j = i + 1; j < n; j++)
            f -= LU[i][j] * x[j][0];
        x[i][0] = f / LU[i][i];
    }
    for (i = n - 2; i >= 0; --i) {    // Inverse Q
        f = x[q[i]][0];
        x[q[i]][0] = x[i][0];
        x[i][0] = f;
    }
    return x;
}

//
//  QR Decomposition With Complete Pivoting: AP = QR
//  where: R is m x n.  V is n x m
//  j-th Householder transform: I - pfBeta[j]*Outer(V[m*j], V[m*j])
//
//  Input:  matrix A is input in the array R.
//  Output: Q is encoded in pfBeta and V.
//          R is written over storage used to input A.  //REVIEW: bad API design
//
//
//  Vector L2 Norm
//
float
ML_VectorNorm(float *pv, int n, int nSpan)
{
    float f, fMax, fNorm;
    int i, nExp;

    fMax = 0;
    n = nSpan*(n - 1);
    for (i = n; i >= 0; i -= nSpan) {
        f = fabs(pv[i]);
        if (f > fMax)
            fMax = f;
    }
    //
    //   to avoid overflow.  Use close power of two
    //  to avoid any roundoff error.
    //
    frexp(fMax, &nExp);
    f = ldexp(1.0f, -nExp);
    fNorm = 0.0;
    for (i = n; i >= 0; i -= nSpan) {
        f = pv[i] * f;
        fNorm += f*f;
    }
    return sqrt(fNorm)/f;
}

int
ML_HouseholderQR(float *V, float pfBeta[], float *R, int p[], int q[], int m, int n)
{
    int i, j, k, iMax, jMax, nExp;
    float fScale, fNorm, fAvj, fLastNorm, f, fEps, *pfNorms;

    if (m < n)
        return 0;
    pfNorms = V + m*(n - 1);
    for (k = 0; k < n; k++)
        q[k] = p[k] = k;
    fEps = 4.0f*sqrt(float(n))*F_MACHINE_EPS;
    for (k = 0; k < n; k++) {
        //
        //  Downdate column norms, recalculate periodically or if bad cancellation
        //
        for (j = k; j < n; j++) {
            if ((k & 7) && fabs(f = R[n*(k - 1) + j]/pfNorms[j]) < 0.999)
                pfNorms[j] *= sqrt(1.0f - f*f);
            else
                pfNorms[j] = ML_VectorNorm(R + n*k + j, m - k, n);
        }
        //
        //  Pivot best column
        //
        jMax = k;
        for (j = k + 1; j < n; j++)
            if (pfNorms[j] > pfNorms[jMax])
                jMax = j;
        p[k] = jMax;
        if (jMax != k) {
            f = pfNorms[k];
            pfNorms[k] = pfNorms[jMax];
            pfNorms[jMax] = f;
            for (i = 0; i < m; i++) {
                f = R[n*i + k];
                R[n*i + k] = R[n*i + jMax];
                R[n*i + jMax] = f;
            }
        }
        //
        //  Pivot best row
        //
        iMax = k;
        for (i = k + 1; i < m; i++)
            if (fabs(R[n*i + k]) > fabs(R[n*iMax + k]))
                iMax = i;
        q[k] = iMax;
        if (iMax != k) {
            for (j = 0; j < n; j++) {
                f = R[n*k + j];
                R[n*k + j] = R[n*iMax + j];
                R[n*iMax + j] = f;
            }
        }
        fNorm = pfNorms[k] = ML_VectorNorm(R + n*k + k, m - k, n);
        //
        //  Rank determination
        //
        if (fNorm == 0.0) {
            ML_fMatrixCondition = fNorm;
            return k;
        }
        if (k == 0) {
            fLastNorm = fNorm;
        } else if (fNorm < fLastNorm*fEps) {
            ML_fMatrixCondition = fNorm/fLastNorm;
            return k;
        }
        //
        //  Householder: I - fBeta*Outer(V, V)
        //
        frexp(fNorm, &nExp);
        fScale = ldexp(1.0f, -nExp);
        for (i = 0; i < k; i++)
            V[m*k + i] = 0.0;
        for (i = k; i < m; i++)
            V[m*k + i] = R[n*i + k] * fScale;   // overwrites pfNorms on last round
        if (R[n*k + k] < 0.0) {
            V[m*k + k] -= fNorm*fScale;
            R[n*k + k] = +fNorm;
        } else {
            V[m*k + k] += fNorm*fScale;
            R[n*k + k] = -fNorm;
        }
        pfBeta[k] = 1.0f/fabs(fNorm * fScale * V[m*k + k]);
        for (i = k + 1; i < m; i++)
            R[n*i + k] = 0.0;
        for (j = k + 1; j < n; j++) {
            fAvj = 0.0;
            for (i = k; i < m; i++)
                fAvj += V[m*k + i] * R[n*i + j];
            fAvj *= pfBeta[k];
            for (i = k; i < m; i++)
                R[n*i + j] -= fAvj * V[m*k + i];
        }
    }
    ML_fMatrixCondition = fNorm / fLastNorm;
    return k;
}
