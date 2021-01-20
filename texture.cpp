#include "stdafx.h"
#include "texture.h"
#define F_PI    3.1415926535897932384626433832795028842f
#define F_2PI   (2.0f*F_PI)
#pragma intrinsic(_lrotl, _lrotr, sin, cos)
//
//  The lattice value is based on a fast, good hashing function that combines
//  multiplication and data-dependant rotation (inspired by current block ciphers)
//REVIEW: spectrum should be more isotropic
//
static inline unsigned
LatticeValue(unsigned n)
{
    n = 1099087573 * n;
    return _lrotr(n, n);
}
//
//  Recursive N-dimensional interpolation of the lattice values
//
struct Vector4 {
    float x, y, z, w;
};

static inline double
DotProduct(Vector4 v1, Vector4 v2)
{
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z + v1.w*v2.w;
}

static Vector4  rgvGrad[MAX_TEXTURE_DIMENSION];
static Vector4  rgvCoef[MAX_TEXTURE_DIMENSION];
static unsigned rgnFloor[MAX_TEXTURE_DIMENSION];

static float
Spline(int nDim, unsigned nLattice)
{
    unsigned nFloor, nH0, nH1, nH2, nH3;
    Vector4 vKnots;

    nFloor = rgnFloor[nDim];
    nH0 = LatticeValue(nFloor - 1 ^ nLattice);
    nH1 = LatticeValue(nFloor     ^ nLattice);
    nH2 = LatticeValue(nFloor + 1 ^ nLattice);
    nH3 = LatticeValue(nFloor + 2 ^ nLattice);
    if (nDim == 1) {
        vKnots.x = float(nH0) - (0.5f*4294967296.0f);
        vKnots.y = float(nH1) - (0.5f*4294967296.0f);
        vKnots.z = float(nH2) - (0.5f*4294967296.0f);
        vKnots.w = float(nH3) - (0.5f*4294967296.0f);
    } else {
        --nDim;
        vKnots.x = Spline(nDim, nH0);
        vKnots.y = Spline(nDim, nH1);
        vKnots.z = Spline(nDim, nH2);
        vKnots.w = Spline(nDim, nH3);
        nDim++;
    }
    return float(DotProduct(rgvCoef[nDim], vKnots));
}
//
//  Precompute cubic b-spline filter coefficients in O(N) instead of doing it
//  in the exponential recursion above.  Note that cubic b-spline is much more
//  isotropic than the catmul-rom spline often suggested.
//

float
BandNoise(const float rgf[], int nDim)
{
    int i, iDim, nFloor;
    float x, fGain;

    fGain = 1.0f;
    for (i = 0; i < nDim; i++) {
        //
        //  fGain prevents the noise from flattening out in higher dimension
        //
        fGain *= 1.27f;
        iDim = nDim - i;
        nFloor = int(rgf[i]) - (rgf[i] < 0.0f && rgf[i] != int(rgf[i]));
        x = rgf[i] - nFloor;
        rgnFloor[iDim] = nFloor;
        rgvCoef[iDim].x = 1.0f/6.0f + x*(x*(0.5f - 1.0f/6.0f*x) - 0.5f);
        rgvCoef[iDim].y = 4.0f/6.0f + x*x*(0.5f*x - 1.0f);
        rgvCoef[iDim].z = 1.0f/6.0f + 0.5f*x*(1.0f + x*(1.0f - x));
        rgvCoef[iDim].w = 1.0f/6.0f*x*x*x;
    }
    return fGain * Spline(nDim, 123456789) * (0.9f/4294967296.0f) + 0.5f;
}
//
//  Compute BandNoise and its N-dimensional gradient
//
static float
GradientSpline(float rgfGradient[], int nDim, unsigned nLattice, int nStride)
{
    unsigned nFloor, nH0, nH1, nH2, nH3;
    Vector4 vKnots, rgvGradientKnots[MAX_TEXTURE_DIMENSION];
    int i, iStride;


    nFloor = rgnFloor[nDim];
    nH0 = LatticeValue(nFloor - 1 ^ nLattice);
    nH1 = LatticeValue(nFloor     ^ nLattice);
    nH2 = LatticeValue(nFloor + 1 ^ nLattice);
    nH3 = LatticeValue(nFloor + 2 ^ nLattice);
    if (nDim == 1) {
        vKnots.x = float(nH0) - (0.5f*4294967296.0f);
        vKnots.y = float(nH1) - (0.5f*4294967296.0f);
        vKnots.z = float(nH2) - (0.5f*4294967296.0f);
        vKnots.w = float(nH3) - (0.5f*4294967296.0f);
    } else {
        --nDim;
        vKnots.x = GradientSpline(&rgvGradientKnots[0].x, nDim, nH0, 4);
        vKnots.y = GradientSpline(&rgvGradientKnots[0].y, nDim, nH1, 4);
        vKnots.z = GradientSpline(&rgvGradientKnots[0].z, nDim, nH2, 4);
        vKnots.w = GradientSpline(&rgvGradientKnots[0].w, nDim, nH3, 4);
        nDim++;
    }
    for (iStride = nStride, i = 1; i < nDim; i++, iStride += nStride)
        rgfGradient[iStride] = float(DotProduct(rgvCoef[nDim], rgvGradientKnots[i - 1]));
    rgfGradient[0] = float(DotProduct(rgvGrad[nDim], vKnots));
    return float(DotProduct(rgvCoef[nDim], vKnots));
}

float
GradientBandNoise(float rgfGradient[], const float rgf[], int nDim)
{
    int i, iDim, nFloor;
    float x, fGain, fValue;

    fGain = 1.0f;
    for (i = 0; i < nDim; i++) {
        //
        //  fGain prevents the noise from flattening out in higher dimension
        //
        fGain *= 1.27f;
        iDim = nDim - i;
        nFloor = int(rgf[i]) - (rgf[i] < 0.0f && rgf[i] != int(rgf[i]));
        x = rgf[i] - nFloor;
        rgnFloor[iDim] = nFloor;
        rgvCoef[iDim].x = 1.0f/6.0f + x*(x*(0.5f - 1.0f/6.0f*x) - 0.5f);
        rgvCoef[iDim].y = 4.0f/6.0f + x*x*(0.5f*x - 1.0f);
        rgvCoef[iDim].z = 1.0f/6.0f + 0.5f*x*(1.0f + x*(1.0f - x));
        rgvCoef[iDim].w = 1.0f/6.0f*x*x*x;

        rgvGrad[iDim].x = -0.5f + x*(1.0f - 0.5f*x);
        rgvGrad[iDim].y = x*(1.5f*x - 2.0f);
        rgvGrad[iDim].z = 0.5f + x*(1.0f - 1.5f*x);
        rgvGrad[iDim].w = 0.5f*x*x;
    }
    fValue = fGain * GradientSpline(rgfGradient, nDim, 123456789, 1) * (0.9f/4294967296.0f) + 0.5f;
    for (i = 0; i < nDim; i++)
        rgfGradient[i] *= fGain * (0.9f/4294967296.0f);
    return fValue;
}
//
//  Musgrave's refinement of fractal noise.  Worley's rotation trick would be nice, but
//  it makes calculation of the gradient too expensive.
//

#define F_FIBSIN    (0.93203242381323f)
#define F_FIBCOS   (-0.36237489008048f)

float
FractalNoise(const float rgf[], int nDim, float fLevelDetail,
       float fFractalDimension, float fLacunarity, int bAbsoluteValue)
{
    float fWaveLength, fOmega, fNoise, fB, fTotalGain, fRemainder;
    float rgfScaled[MAX_TEXTURE_DIMENSION];
    int i, j, nLevelDetail;

    fOmega = powf(fLacunarity, -fFractalDimension);
    fWaveLength = 1.0f;
    fNoise = 0.0f;
    fTotalGain = 0.0f;
    nLevelDetail = int(fLevelDetail);
    for (j = 0; j < nDim; j++)
        rgfScaled[j] = rgf[j];
    for (i = 0; i < nLevelDetail; i++) {
        fB = fWaveLength * (BandNoise(rgfScaled, nDim) - 0.5f);
        if (bAbsoluteValue)
            fB = fabsf(fB);
        fNoise += fB;
        fTotalGain += fWaveLength;
        fWaveLength *= fOmega;
        for (j = 0; j < nDim; j++)
            rgfScaled[j] *= fLacunarity;
    }
    //
    //  If level of detail is animated, blending in the remainder will prevent snapping
    //
    //REVIEW: rescale to last wavelength
    fRemainder = fLevelDetail - float(nLevelDetail);
    if (fRemainder) {
        fB = fRemainder * fWaveLength * (BandNoise(rgfScaled, nDim) - 0.5f);
        if (bAbsoluteValue)
            fB = fabsf(fB);
        fNoise += fB;
        fTotalGain += fWaveLength*fRemainder;
    }
    return fNoise / fTotalGain + 0.5f;
}

float
GradientFractalNoise(float rgfGradient[], const float rgf[], int nDim, float fLevelDetail,
       float fFractalDimension, float fLacunarity, int bAbsoluteValue)
{
    float fWaveLength, fScale, fSign, fOmega, fNoise, fB, fTotalGain, fRemainder;
    float rgfScaled[MAX_TEXTURE_DIMENSION], rgfBGrad[MAX_TEXTURE_DIMENSION];
    int i, j, nLevelDetail;

    fOmega = powf(fLacunarity, -fFractalDimension);
    fWaveLength = 1.0f;
    fScale = 1.0f;
    fNoise = 0.0f;
    fTotalGain = 0.0f;
    nLevelDetail = int(fLevelDetail);
    for (j = 0; j < nDim; j++){
        rgfScaled[j] = rgf[j];
        rgfGradient[j] = 0.0f;
    }
    for (i = 0; i < nLevelDetail; i++) {
        fB = fWaveLength * (GradientBandNoise(rgfBGrad, rgfScaled, nDim) - 0.5f);
        fSign = 1.0f;
        if (bAbsoluteValue) {
            if (fB < 0.0)
                fSign = -1.0f;
            fB = fabsf(fB);
        }
        fNoise += fB;
        for (j = 0; j < nDim; j++) {
            rgfScaled[j] *= fLacunarity;
            rgfGradient[j] += rgfBGrad[j]*fScale*fWaveLength*fSign;
        }
        fTotalGain += fWaveLength;
        fWaveLength *= fOmega;
        fScale *= fLacunarity;
    }
    //
    //  If level of detail is animated, blending in the remainder will prevent snapping
    //
    //REVIEW: rescale to last wavelength
    fRemainder = fLevelDetail - float(nLevelDetail);
    if (fRemainder) {
        fB = fRemainder * fWaveLength * (BandNoise(rgfScaled, nDim) - 0.5f);    //BUGBUG GradientBandnoise?
        fSign = 1.0f;
        if (bAbsoluteValue) {
            if (fB < 0.0)
                fSign = -1.0f;
            fB = fabsf(fB);
        }
        fNoise += fB;
        for (j = 0; j < nDim; j++)
            rgfGradient[j] += rgfBGrad[j]*fScale*fRemainder*fWaveLength*fSign;
        fTotalGain += fWaveLength*fRemainder;
    }
    fTotalGain = 1.0f/fTotalGain;
    for (j = 0; j < nDim; j++)
        rgfGradient[j] *= fTotalGain;
    return fNoise * fTotalGain + 0.5f;
}

//
//  Generate periodic texture by mapping an intrinsically flat torus into
//  2N dimentions.
//
float
ML_PeriodicTexture(const float rgf[], int nDim, float fBandlimit, ML_TextureType *pt)
{
    float fOver2Pi, rgf2[MAX_TEXTURE_DIMENSION];
    int i;

    fOver2Pi = 0.5f/F_PI;
    for (i = 0; i < nDim; i++) {
        rgf2[2*i] = sinf(F_2PI*rgf[i]) * fOver2Pi;
        rgf2[2*i+1] = cosf(F_2PI*rgf[i]) * fOver2Pi;
    }
    return (*pt)(rgf2, 2*nDim, fBandlimit);
}
