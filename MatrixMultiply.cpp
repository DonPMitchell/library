#include "stdafx.h"
#include "Timing.h"
#include "Montecarlo.h"
#include "smmintrin.h"
#include "immintrin.h"

//
//  Results: 1013.354126 1026.571045 1025.290894
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
#define TILE_L1      32
#define TILE_L2     128
#define TILE_L3     512
#define NB          4096

static void
MatrixMultiply_Naive(float C[NB][NB], float A[NB][NB], float B[NB][NB])
{
    int i, j, k;

    for (i = 0; i < NB; i++) {
        for (j = 0; j < NB; j++) {
            C[i][j] = 0.0;
            for (k = 0; k < NB; k++)
                C[i][j] += A[i][k]*B[k][j];
        }
    }
}

static void
MatrixMultiply_Transpose(float C[NB][NB], float A[NB][NB], float B[NB][NB])
{
    int i, j, k;
    static float D[NB][NB];

    for (i = 0; i < NB; i++) {
        for (j = 0; j < NB; j++) {
            D[i][j] = B[j][i];
        }
    }
    for (i = 0; i < NB; i++) {
        for (j = 0; j < NB; j++) {
            C[i][j] = 0.0;
            for (k = 0; k < NB; k++)
                C[i][j] += A[i][k]*D[j][k];
        }
    }
}

#define TILE    32
#define UNWIND  8

static void
MatrixMultiply_LoopTilingNoTrans(float C[NB][NB], float A[NB][NB], float B[NB][NB])
{
    int i, j, k, iB, jB, kB;
    float f;
    static float D[NB][NB];

    for (i = 0; i < NB; i++) {
        for (j = 0; j < NB; j++) {
            C[i][j] = 0.0;
        }
    }
    for (iB = 0; iB < NB; iB += TILE) {
        for (jB = 0; jB < NB; jB += TILE) {
            for (kB = 0; kB < NB; kB += TILE) {
                for (i = iB; i < iB+TILE; i++) {
                    for (j = jB; j < jB+TILE; j++) {
                        f = 0.0;
                        for (k = kB; k < kB+TILE; k+=UNWIND)
                            f       += A[i][k  ]*B[k  ][j]
                                     + A[i][k+1]*B[k+1][j]
                                     + A[i][k+2]*B[k+2][j]
                                     + A[i][k+3]*B[k+3][j]
                                     ;
                            C[i][j] += f;
                    }
                }
            }
        }
    }
}

static void
MatrixMultiply_LoopTiling(float C[NB][NB], float A[NB][NB], float B[NB][NB])
{
    int i, j, k, iB, jB, kB;
    __declspec( align( 16 ) )static float D[NB][NB];

    for (i = 0; i < NB; i++) {
        for (j = 0; j < NB; j++) {
            D[i][j] = B[j][i];
            C[i][j] = 0.0;
        }
    }
    for (iB = 0; iB < NB; iB += TILE) {
        for (jB = 0; jB < NB; jB += TILE) {
            for (kB = 0; kB < NB; kB += TILE) {
                for (i = iB; i < iB+TILE; i++) {
                    for (j = jB; j < jB+TILE; j++) {
                        for (k = kB; k < kB+TILE; k+=UNWIND)
                            C[i][j] += A[i][k  ]*D[j][k  ]
                                     + A[i][k+1]*D[j][k+1]
                                     + A[i][k+2]*D[j][k+2]
                                     + A[i][k+3]*D[j][k+3]
                                     ;
                    }
                }
            }
        }
    }
}

static void
MatrixMultiply_LoopTilingSSE(float C[NB][NB], float A[NB][NB], float B[NB][NB])
{
    int i, j, k, iB, jB, kB;
    __declspec( align( 16 ) )static float D[NB][NB];
    __m128 rC, rA, rD;
    float fC;

    for (i = 0; i < NB; i++) {
        for (j = 0; j < NB; j++) {
            D[i][j] = B[j][i];
            C[i][j] = 0.0;
        }
    }
    for (iB = 0; iB < NB; iB += TILE) {
        for (jB = 0; jB < NB; jB += TILE) {
            for (kB = 0; kB < NB; kB += TILE) {
                for (i = iB; i < iB+TILE; i++) {
                    for (j = jB; j < jB+TILE; j++) {
                        fC = 0.0;
                        for (k = kB; k < kB+TILE; k+=4) {
                            rA = _mm_load_ps(A[i] + k    );
                            rD = _mm_load_ps(D[j] + k    );
                            rC = _mm_dp_ps(rA, rD, 0xF1);
                            fC += rC.m128_f32[0];
                            //rA = _mm_load_ps(A[i] + k + 4);
                            //rD = _mm_load_ps(D[j] + k + 4);
                            //rC = _mm_dp_ps(rA, rD, 0xF1);
                            //fC += rC.m128_f32[0];
                        }
                        C[i][j] += fC;
                    }
                }
            }
        }
    }
}

static void
MatrixMultiply_LoopTilingSSE2(float C[NB][NB], float A[NB][NB], float B[NB][NB])
{
    int i, j, k, iB, jB, kB;
    __declspec( align( 16 ) )static float D[NB][NB];
    __m128 rC, rA, rDot, rD0, rD1, rD2, rD3;

    for (i = 0; i < NB; i++) {
        for (j = 0; j < NB; j++) {
            D[i][j] = B[j][i];
            C[i][j] = 0.0;
        }
    }
    for (iB = 0; iB < NB; iB += TILE) {
        for (jB = 0; jB < NB; jB += TILE) {
            for (kB = 0; kB < NB; kB += TILE) {
                for (i = iB; i < iB+TILE; i++) {
                    for (j = jB; j < jB+TILE; j+=4) {
                        rC = _mm_load_ps(C[i]+j);
                        for (k = kB; k < kB+TILE; k+=4) {
                            rA = _mm_load_ps(A[i]+k);
                            rD0 = _mm_load_ps(D[j  ]+k);
                            rD1 = _mm_load_ps(D[j+1]+k);
                            rD2 = _mm_load_ps(D[j+2]+k);
                            rD3 = _mm_load_ps(D[j+3]+k);
                            rC = _mm_add_ps(rC, _mm_dp_ps(rD0, rA, 0xF1));
                            rC = _mm_add_ps(rC, _mm_dp_ps(rD1, rA, 0xF2));
                            rC = _mm_add_ps(rC, _mm_dp_ps(rD2, rA, 0xF4));
                            rC = _mm_add_ps(rC, _mm_dp_ps(rD3, rA, 0xF8));
                        }
                        _mm_store_ps(C[i]+j, rC);
                    }
                }
            }
        }
    }
}

static void
MatrixMultiply_LoopTilingSSE3(float C[NB][NB], float A[NB][NB], float B[NB][NB])
{
    int i, j, k, iB, jB, kB;
    __m128 rC;

    for (i = 0; i < NB; i++) {
        for (j = 0; j < NB; j++) {
            C[i][j] = 0.0;
        }
    }
    for (iB = 0; iB < NB; iB += TILE) {
        for (jB = 0; jB < NB; jB += TILE) {
            for (kB = 0; kB < NB; kB += TILE) {
                for (i = iB; i < iB+TILE; i++) {
                    for (j = jB; j < jB+TILE; j+=4) {
                        rC = _mm_load_ps(C[i]+j);
                        for (k = kB; k < kB+TILE; k+=4) {
                            rC = _mm_add_ps(rC, _mm_mul_ps(_mm_load_ps(B[k  ]+j), _mm_set1_ps(A[i][k  ])));
                            rC = _mm_add_ps(rC, _mm_mul_ps(_mm_load_ps(B[k+1]+j), _mm_set1_ps(A[i][k+1])));
                            rC = _mm_add_ps(rC, _mm_mul_ps(_mm_load_ps(B[k+2]+j), _mm_set1_ps(A[i][k+2])));
                            rC = _mm_add_ps(rC, _mm_mul_ps(_mm_load_ps(B[k+3]+j), _mm_set1_ps(A[i][k+3])));
                        }
                        _mm_store_ps(C[i]+j, rC);
                    }
                }
            }
        }
    }
}

static void
MatrixMultiply_LoopTilingSSE4(float C[NB][NB], float A[NB][NB], float B[NB][NB])
{
    int i, j, k, iB, jB, kB, iB2, jB2, kB2;
    __m128 rC;

    for (i = 0; i < NB; i++) {
        for (j = 0; j < NB; j++) {
            C[i][j] = 0.0;
        }
    }
    for (iB2 = 0; iB2 < NB; iB2 += TILE_L2) {
        for (jB2 = 0; jB2 < NB; jB2 += TILE_L2) {
            for (kB2 = 0; kB2 < NB; kB2 += TILE_L2) {
                for (iB = iB2; iB < iB2+TILE_L2; iB += TILE_L1) {
                    for (jB = jB2; jB < jB2+TILE_L2; jB += TILE_L1) {
                        for (kB = kB2; kB < kB2+TILE_L2; kB += TILE_L1) {
                            for (i = iB; i < iB+TILE_L1; i++) {
                                for (j = jB; j < jB+TILE_L1; j+=4) {
                                    rC = _mm_load_ps(C[i]+j);
                                    for (k = kB; k < kB+TILE_L1; k+=4) {
                                        rC = _mm_add_ps(rC, _mm_mul_ps(_mm_load_ps(B[k  ]+j), _mm_set1_ps(A[i][k  ])));
                                        rC = _mm_add_ps(rC, _mm_mul_ps(_mm_load_ps(B[k+1]+j), _mm_set1_ps(A[i][k+1])));
                                        rC = _mm_add_ps(rC, _mm_mul_ps(_mm_load_ps(B[k+2]+j), _mm_set1_ps(A[i][k+2])));
                                        rC = _mm_add_ps(rC, _mm_mul_ps(_mm_load_ps(B[k+3]+j), _mm_set1_ps(A[i][k+3])));
                                    }
                                    _mm_store_ps(C[i]+j, rC);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

static void
MatrixMultiply_LoopTilingAVX(float C[NB][NB], float A[NB][NB], float B[NB][NB])
{
    int i, j, k, iB, jB, kB;
    __m256 rC;

    for (i = 0; i < NB; i++) {
        for (j = 0; j < NB; j++) {
            C[i][j] = 0.0;
        }
    }
    for (iB = 0; iB < NB; iB += TILE) {
        for (jB = 0; jB < NB; jB += TILE) {
            for (kB = 0; kB < NB; kB += TILE) {
                for (i = iB; i < iB+TILE; i++) {
                    for (j = jB; j < jB+TILE; j+=8) {
                        rC = _mm256_load_ps(C[i]+j);
                        for (k = kB; k < kB+TILE; k+=8) {
                            rC = _mm256_add_ps(rC, _mm256_mul_ps(_mm256_load_ps(B[k  ]+j), _mm256_set1_ps(A[i][k  ])));
                            rC = _mm256_add_ps(rC, _mm256_mul_ps(_mm256_load_ps(B[k+1]+j), _mm256_set1_ps(A[i][k+1])));
                            rC = _mm256_add_ps(rC, _mm256_mul_ps(_mm256_load_ps(B[k+2]+j), _mm256_set1_ps(A[i][k+2])));
                            rC = _mm256_add_ps(rC, _mm256_mul_ps(_mm256_load_ps(B[k+3]+j), _mm256_set1_ps(A[i][k+3])));
                            rC = _mm256_add_ps(rC, _mm256_mul_ps(_mm256_load_ps(B[k+4]+j), _mm256_set1_ps(A[i][k+4])));
                            rC = _mm256_add_ps(rC, _mm256_mul_ps(_mm256_load_ps(B[k+5]+j), _mm256_set1_ps(A[i][k+5])));
                            rC = _mm256_add_ps(rC, _mm256_mul_ps(_mm256_load_ps(B[k+6]+j), _mm256_set1_ps(A[i][k+6])));
                            rC = _mm256_add_ps(rC, _mm256_mul_ps(_mm256_load_ps(B[k+7]+j), _mm256_set1_ps(A[i][k+7])));
                        }
                        _mm256_store_ps(C[i]+j, rC);
                    }
                }
            }
        }
    }
}

static void
MatrixMultiply_LoopTilingSSE5(float C[NB][NB], float A[NB][NB], float B[NB][NB])
{
    int i, j, k, iB, jB, kB;
    __m128 rC0, rC1, rC2, rC3, rB0, rB1, rB2, rB3;

    for (i = 0; i < NB; i++) {
        for (j = 0; j < NB; j++) {
            C[i][j] = 0.0;
        }
    }
    for (iB = 0; iB < NB; iB += TILE) {
        for (jB = 0; jB < NB; jB += TILE) {
            for (kB = 0; kB < NB; kB += TILE) {
                for (i = iB; i < iB+TILE; i+=4) {
                    for (j = jB; j < jB+TILE; j+=4) {
                        rC0 = _mm_load_ps(C[i  ]+j);
                        rC1 = _mm_load_ps(C[i+1]+j);
                        rC2 = _mm_load_ps(C[i+2]+j);
                        rC3 = _mm_load_ps(C[i+3]+j);
                        for (k = kB; k < kB+TILE; k+=4) {
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
}

void
TestSSE()
{
    __m128 a, b;
    const int mask = 0xF1;

    printf("Test SSE\n");
    a.m128_f32[0] = 1.5;
    a.m128_f32[1] = 10.25;
    a.m128_f32[2] = -11.0625;
    a.m128_f32[3] = 81.0;
    b.m128_f32[0] = -1.5;
    b.m128_f32[1] = 3.125;
    b.m128_f32[2] = -50.5;
    b.m128_f32[3] = 100.0;

    __m128 res = _mm_dp_ps(a, b, mask);

    printf_s("Original a: %f\t%f\t%f\t%f\nOriginal b: %f\t%f\t%f\t%f\n",
                a.m128_f32[0], a.m128_f32[1], a.m128_f32[2], a.m128_f32[3],
                b.m128_f32[0], b.m128_f32[1], b.m128_f32[2], b.m128_f32[3]);
    printf_s("Result res: %f\t%f\t%f\t%f\n",
                res.m128_f32[0], res.m128_f32[1], res.m128_f32[2], res.m128_f32[3]);
}
void
TestMatrixMultiply()
{
    TimingInfo ti;
    static __declspec( align( 32 ) ) float C[NB][NB], A[NB][NB], B[NB][NB];
    int i, j;

    //TestSSE();
    printf("Testing Matrix Multiplication\n");
    for (j = 0; j < NB; j++) {
        for (i = 0; i < NB; i++) {
            A[i][j] = RandomDouble();
            B[i][j] = RandomDouble();
        }
    }
    StartTiming(ti);
    MatrixMultiply_LoopTilingSSE5(C, A, B);
    StopTiming(ti);
    printf(" Correct: 1013.354126 1026.571045 1025.290894\n");
    printf(" Results: %f %f %f\n", C[117][200], C[0][0], C[3999][2999]);
    ReportTiming(ti);
}








#define NBLOCK 4096

static void
MatrixMultiply(float *pmA, const float *pmB, const float *pmC)
{
    int i, j, k;
    const float *pmSaveC;
    double f;

    if (pmA == pmB || pmA == pmC)
        return;
    pmSaveC = pmC;
    for (i = 0; i < NBLOCK; i++) {              // 7 times slower than loop-tiled multiply
        for (j = 0; j < NBLOCK; j++) {
            f = 0.0;
            for (k = 0; k < NBLOCK; k++) {
                f += pmB[k] * pmC[j];
                pmC += NBLOCK;
            }
            *pmA++ = f;
            pmC = pmSaveC;
        }
        pmB += NBLOCK;
    }
}

static
void matrix_mult_wiki_block(const float*A , const float* B, float* C, const int N, const int M, const int K) {
   const int block_size = 8;  //I have tried several different block sizes
   for(int i=0; i<N; i++) {
       for(int j=0; j<K; j++) {
           C[K*i + j] = 0;
       }
    }
    for(int l2=0; l2<M; l2+=block_size) {
        for(int j2=0; j2<K; j2+=block_size) {
        #pragma omp parallel for
            for(int i=0; i<N; i++) {
                for(int l=l2; l<min(M, l2+block_size); l++) {
                    for(int j=j2; j<min(K, j2+block_size); j++) {
                        C[K*i + j] += A[M*i+l]*B[K*l+j];
                    }
                }
            }
        }
    }
}

#define LOOPTILE 32

void 
matrix_mult_looptiled(float* C, const float*A , const float* B) 
{
   for(int i=0; i<NBLOCK; i++) {
       for(int j=0; j<NBLOCK; j++) {
           C[NBLOCK*i + j] = 0;
       }
    }
    for(int l2=0; l2<NBLOCK; l2+=LOOPTILE) {
        for(int j2=0; j2<NBLOCK; j2+=LOOPTILE) {
            for(int i=0; i<NBLOCK; i++) {
                for(int l=l2; l<l2+LOOPTILE; l++) {
                    for(int j=j2; j<j2+LOOPTILE; j++) {
                        C[NBLOCK*i + j  ] += A[NBLOCK*i+l]*B[NBLOCK*l+j  ];
                    }
                }
            }
        }
    }
}

//
//  Gaussian Elimination with Complete Pivoting: PAQ = LU
//
#define D_MACHINE_EPS   2.22044605e-016
#define F_MACHINE_EPS   1.1920929e-007f

float ML_fMatrixCondition;  // internal report from GECP and QRCP

int
ML_GaussianElimination(float *A, int p[], int q[], int n)
{
    int i, j, k, iPivot, jPivot;
    float fPivot, fTmp, fLastPivot, fEps, *AI, *AJ;

    for (k = 0; k < n; k++)
        p[k] = q[k] = k;
    fEps = 8.0*n*D_MACHINE_EPS;
    for (k = 0; k < n - 1; k++) {     //REVIEW: if non-square, loop to k < n
        //
        //  Find maximum element of submatrix
        //
        fPivot = -1.0;
        AI = A + n*k;
        for (i = k; i < n; i++) {
            for (j = k; j < n; j++) {
                if (fabs(AI[j]) > fPivot) {
                    fPivot = fabs(AI[j]);   // A[n*i + j]
                    iPivot = i;
                    jPivot = j;
                }
            }
            AI += n;
        }
        //
        //  Rank determination
        //
        if (fPivot == 0.0) {
            ML_fMatrixCondition = fPivot;
            return k;
        }
        if (k == 0) {
            fLastPivot = fPivot;
        } else if ( fPivot < fLastPivot*fEps) {
            ML_fMatrixCondition = fPivot/fLastPivot;
            return k;
        }
        //
        //  Exchange rows and columns to move pivot to k,k
        //
        p[k] = iPivot;
        q[k] = jPivot;
        if (jPivot != k) {
            AI = A;
            for (i = 0; i < n; i++) {
                fTmp = AI[k];
                AI[k] = AI[jPivot];
                AI[jPivot] = fTmp;
                AI += n;
            }
        }
        AI = A + n*k;
        if (iPivot != k) {
            AJ = A + n*iPivot;
            for (j = 0; j < n; j++) {
                fTmp = AI[j];
                AI[j] = AJ[j];
                AJ[j] = fTmp;
            }
        }
        //
        //  Elimination
        //
        fPivot = 1.0/AI[k];
        AJ = AI;
        for (i = k + 1; i < n; i++) {
            AJ += n;
            AJ[k] *= fPivot;
            for (j = k + 1; j < n; j++)
                AJ[j] -= AJ[k] * AI[j];     // A[n*i+j] -= A[n*i+k] * A[n*k + j]
        }
    }
    fPivot = fabs(A[n*k + k]);
    if (k > 0)
        ML_fMatrixCondition = fPivot / fLastPivot;
    else
        ML_fMatrixCondition = fPivot;
    if (fPivot == 0.0 || (k > 0 && fPivot < fLastPivot*fEps))
        return k;   // rank = n - 1
    else
        return k+1; // full rank = n
}

void
ML_SolveLinearSystem(float *x, float *LU, int p[], int q[], int n)
{
    int i, j;
    float f;

    for (i = 0; i < n - 1; i++) {   // Permuation P
        f = x[p[i]];
        x[p[i]] = x[i];
        x[i] = f;
    }
    for (i = 1; i < n; i++) {       // Forward Substitution
        f = x[i];
        LU += n;
        for (j = 0; j < i; j++)
            f -= LU[j] * x[j];      // LU[n*i + j]
        x[i] = f;
    }
    x[n-1] /= LU[n - 1];
    for (i = n - 2; i >= 0; --i) {    // Back Substitution
        f = x[i];
        LU -= n;
        for (j = i + 1; j < n; j++)
            f -= LU[j] * x[j];
        x[i] = f / LU[i];
    }
    for (i = n - 2; i >= 0; --i) {    // Inverse Q
        f = x[q[i]];
        x[q[i]] = x[i];
        x[i] = f;
    }
}

int
ML_InvertMatrix(float *A, float *Y, int n)
{
    int *p, *q, kRank, i, j;
    float *x;

    p = new int[n];
    q = new int[n];
    x = new float[n];
    kRank = ML_GaussianElimination(A, p, q, n);
    if (kRank != n)
        return kRank;
    for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++)
            x[i] = 0.0;
        x[j] = 1.0;
        ML_SolveLinearSystem(x, A, p, q, n);
        for (i = 0; i < n; i++)
            Y[n*i + j] = x[i];
    }
    delete [] p;
    delete [] q;
    delete [] x;
    return n;
}