//
//  Pseudo-random number generators
//  D.P. Mitchell  2018/09/11.
//
#pragma once

extern "C"
{
   unsigned __int64 __emulu(unsigned int a, unsigned int b);
}

#pragma intrinsic(__emulu)

extern unsigned     ML_nMarsagliaX, ML_nMarsagliaC;
extern unsigned     ML_TrueRandomUnsigned();
extern unsigned     ML_RandomLessThanN(unsigned N);

extern double       ML_RandomNormal();
extern double       ML_RandomExponential(double mu);
extern double       ML_RandomGamma(double a);
extern double       ML_RandomBeta(double a, double b);
extern double       ML_RandomChiSquare(double nu);
extern double       ML_RandomFDistribution(double nu1, double nu2);
extern double       ML_RandomTDistribution(double nu);
extern int          ML_RandomGeometric(double p);
extern int          ML_RandomBinomial(double p, int trials);
extern int          ML_RandomPoisson(double mu);

extern double       ML_RandomDisk(double x[2]);             // random point on unit disk
extern double       ML_RandomCircle(double x[2]);           // random point on unit circle
extern void         ML_RandomSphere(double x[3]);           // random point on unit sphere
extern void         ML_RandomHyperSphere(double x[4]);      // random point on unit 4D sphere
extern void         ML_RandomRotation(double m[3][3]);      // random rotation matrix

inline unsigned
ML_RandomUnsigned()
{
    unsigned __int64 n;

    n = __emulu(1965537969, ML_nMarsagliaX);
    n = n + ML_nMarsagliaC;
    ML_nMarsagliaX = unsigned(n & 0xFFFFFFFF);
    ML_nMarsagliaC = unsigned(n >> 32);
    return ML_nMarsagliaX;
}

inline float
ML_RandomFloat()
{
    return float(ML_RandomUnsigned()) * (1.0f/4294967296.0f);
}

inline double
ML_RandomDouble()
{
    return (double(ML_RandomUnsigned()) + double(ML_RandomUnsigned())*(1.0/4294967296.0))*(1.0/4294967296.0);
}

inline void
ML_InitializeRandom(unsigned nX = 886459, unsigned nC = 361290869)
{
    if (nX == 0 && nC == 0) {
        nX = ML_TrueRandomUnsigned();
        nC = ML_TrueRandomUnsigned();
    }
    ML_nMarsagliaX = nX;
    ML_nMarsagliaC = nC;
}
//
//  Randomly shuffle an array
//
template <class T>
void
ML_RandomShuffle(T rgtItem[], int nItems)
{
    int i, j;
    T tTemp;

    for (j = nItems - 1; j > 0; --j) {
        i = ML_RandomLessThanN(j + 1);
        tTemp = rgtItem[i];
        rgtItem[i] = rgtItem[j];
        rgtItem[j] = tTemp;
    }
}
//
//  Simple but not very good mixed-congruential RNG, as described by Knuth
//
/*
double
ML_SimpleRand()
{
    static unsigned x = 1;

    x = 1099087573*x + 2654435761;
    return double(x) / 4294967296.0;
}
*/