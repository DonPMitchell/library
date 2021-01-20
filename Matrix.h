//
//  General N x N matrix class
//  D.P. Mitchell 2018/08/30.
//
#pragma once
#pragma intrinsic (fabs)

struct Matrix {
    float  *pf;
    int     *p, *q;
    int     M, N, kRank;

            Matrix(int n) : M(n), N(n), p(0), q(0) { if (n) pf = new float[n*n]; };
            Matrix(int m, int n) : M(m), N(n), p(0), q(0) { if (n && m) pf = new float[m*n]; };
            //~Matrix() { if(--nCount > 0) return; printf("deleting %f\n", pf[0]);  delete pf; pf = 0; };

    float* operator [](int j) { return pf + j*N; }
    Matrix  operator +(Matrix &m) { Matrix m2(M, N); int i; 
                                    for (i = 0; i < M*N; i++) m2.pf[i] = pf[i] + m.pf[i]; 
                                    return m2;
                                }
    Matrix  operator -(Matrix &m) { Matrix m2(M, N); int i; 
                                    for (i = 0; i < M*N; i++) m2.pf[i] = pf[i] - m.pf[i]; 
                                    return m2;
                                }
    Matrix  operator *(float f) { Matrix m2(M, N); int i; 
                                    for (i = 0; i < M*N; i++) m2.pf[i] = pf[i] * f; 
                                    return m2;
                                }
    double  MaxNorm() { double max; int i; max = fabs(*pf); 
                        for (i = 1; i < M*N; i++)
                            if (fabs(pf[i]) > max) max = fabs(pf[i]);  return max; }
    Matrix  operator *(Matrix B);

    Matrix  GaussianLU();       // return an LU decomposition of matrix
    Matrix  SolveLU(Matrix &x);

    void    Print(int verbose = 1);
    Matrix  Delete() { M = N = 0; delete pf; pf = 0; return *this; }
};
