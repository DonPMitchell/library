#include "stdafx.h"
#include "Polyhedra.h"

#define D_SQRT3         1.73205080756887729352744634150587236694280525381038
#define D_PHI           1.61803398874989484820458683436563811772030917980576
#define D_CBRT4         1.5874010519681994747517056392723
#define SIN2PI5         0.95105651629515357211643933337938
#define D_GOLD D_PHI
#define T   (1.0/D_SQRT3)
#define P1  (D_PHI/D_SQRT3)
#define P2  ((D_PHI - 1.0)/D_SQRT3)
#define I1  (0.5/SIN2PI5)               // 0.52573111211913     
#define I2  (D_PHI/(2.0*SIN2PI5))       //0.85065080835204	

const Vector3 ML_rgvSimpleCubic[3] = {
    Vector3(1.0, 0.0, 0.0),
    Vector3(0.0, 1.0, 0.0),
    Vector3(0.0, 0.0, 1.0)
};

const Vector3 ML_rgvFaceCenteredCubic[3] = {
    Vector3(0.0, 0.5*D_CBRT4, 0.5*D_CBRT4),
    Vector3(0.5*D_CBRT4, 0.0, 0.5*D_CBRT4),
    Vector3(0.5*D_CBRT4, 0.5*D_CBRT4, 0.0),
};

const Vector3 ML_rgvTetrahedron[4] = {
    Vector3( T,  T,  T),
    Vector3( T, -T, -T),
    Vector3(-T,  T, -T),
    Vector3(-T, -T,  T)
};

const Vector3 ML_rgvCube[8] = {
    Vector3( T,  T,  T),
    Vector3( T,  T, -T),
    Vector3( T, -T,  T),
    Vector3( T, -T, -T),
    Vector3(-T,  T,  T),
    Vector3(-T,  T, -T),
    Vector3(-T, -T,  T),
    Vector3(-T, -T, -T)
};

const Vector3 ML_rgvOctahedron[6] = {
    Vector3( 1.0,  0.0,  0.0),
    Vector3(-1.0,  0.0,  0.0),
    Vector3( 0.0,  1.0,  0.0),
    Vector3( 0.0, -1.0,  0.0),
    Vector3( 0.0,  0.0,  1.0),
    Vector3( 0.0,  0.0, -1.0)
};

const Vector3 ML_rgvDodecahedron[20] = {
    Vector3( T,  T,  T),
    Vector3( T,  T, -T),
    Vector3( T, -T,  T),
    Vector3( T, -T, -T),
    Vector3(-T,  T,  T),
    Vector3(-T,  T, -T),
    Vector3(-T, -T,  T),
    Vector3(-T, -T, -T),

    Vector3(0.0,  P2,  P1),
    Vector3(0.0,  P2, -P1),
    Vector3(0.0, -P2,  P1),
    Vector3(0.0, -P2, -P1),

    Vector3( P1, 0.0,  P2),
    Vector3( P1, 0.0, -P2),
    Vector3(-P1, 0.0,  P2),
    Vector3(-P1, 0.0, -P2),

    Vector3( P2,  P1, 0.0),
    Vector3( P2, -P1, 0.0),
    Vector3(-P2,  P1, 0.0),
    Vector3(-P2, -P1, 0.0)
};

const Vector3 ML_rgvIcosahedron[12] = {
    Vector3( I1, 0.0,  I2),
    Vector3( I1, 0.0, -I2),
    Vector3(-I1, 0.0,  I2),
    Vector3(-I1, 0.0, -I2),

    Vector3(0.0,  I2,  I1),
    Vector3(0.0,  I2, -I1),
    Vector3(0.0, -I2,  I1),
    Vector3(0.0, -I2, -I1),

    Vector3( I2,  I1, 0.0),
    Vector3(-I2,  I1, 0.0),
    Vector3( I2, -I1, 0.0),
    Vector3(-I2, -I1, 0.0)
};

#undef T
#undef P1
#undef P2
#undef I1
#undef I2

const Matrix3 ML_rgmTetrahedralGroup[12] = {
    Matrix3( 1,     0,     0,
                  0,     1,     0,
                  0,     0,     1)
    ,
    Matrix3( 0,     0,     1,
                  1,     0,     0,
                  0,     1,     0)
    ,
    Matrix3( 0,     1,     0,
                  0,     0,     1,
                  1,     0,     0)
    ,
    Matrix3( 0,     0,    -1,
                 -1,     0,     0,
                  0,     1,     0)
    ,
    Matrix3( 0,     1,     0,
                  0,     0,    -1,
                 -1,     0,     0)
    ,
    Matrix3(-1,     0,     0,
                  0,     1,     0,
                  0,     0,    -1)
    ,
    Matrix3( 0,     0,     1,
                 -1,     0,     0,
                  0,    -1,     0)
    ,
    Matrix3( 0,    -1,     0,
                  0,     0,     1,
                 -1,     0,     0)
    ,
    Matrix3(-1,     0,     0,
                  0,    -1,     0,
                  0,     0,     1)
    ,
    Matrix3( 0,     0,    -1,
                  1,     0,     0,
                  0,    -1,     0)
    ,
    Matrix3( 0,    -1,     0,
                  0,     0,    -1,
                  1,     0,     0)
    ,
    Matrix3( 1,     0,     0,
                  0,    -1,     0,
                  0,     0,    -1)
};

const Matrix3 ML_rgmOctahedralGroup[24] = {
    Matrix3( 1,     0,     0,
                  0,     1,     0,
                  0,     0,     1)
    ,
    Matrix3( 0,     0,     1,
                  1,     0,     0,
                  0,     1,     0)
    ,
    Matrix3( 0,     1,     0,
                  0,     0,     1,
                  1,     0,     0)
    ,
    Matrix3(-1,     0,     0,
                  0,     1,     0,
                  0,     0,    -1)
    ,
    Matrix3( 0,     0,    -1,
                 -1,     0,     0,
                  0,     1,     0)
    ,
    Matrix3( 0,     1,     0,
                  0,     0,    -1,
                 -1,     0,     0)
    ,
    Matrix3( 1,     0,     0,
                  0,     0,    -1,
                  0,     1,     0)
    ,
    Matrix3( 0,     1,     0,
                  1,     0,     0,
                  0,     0,    -1)
    ,
    Matrix3( 0,     0,    -1,
                  0,     1,     0,
                  1,     0,     0)
    ,
    Matrix3(-1,     0,     0,
                  0,     0,    -1,
                  0,    -1,     0)
    ,
    Matrix3( 0,    -1,     0,
                 -1,     0,     0,
                  0,     0,    -1)
    ,
    Matrix3( 0,     0,    -1,
                  0,    -1,     0,
                 -1,     0,     0)
    ,
    Matrix3( 1,     0,     0,
                  0,    -1,     0,
                  0,     0,    -1)
    ,
    Matrix3( 0,     0,    -1,
                  1,     0,     0,
                  0,    -1,     0)
    ,
    Matrix3( 0,    -1,     0,
                  0,     0,    -1,
                  1,     0,     0)
    ,
    Matrix3(-1,     0,     0,
                  0,    -1,     0,
                  0,     0,     1)
    ,
    Matrix3( 0,     0,     1,
                 -1,     0,     0,
                  0,    -1,     0)
    ,
    Matrix3( 0,    -1,     0,
                  0,     0,     1,
                 -1,     0,     0)
    ,
    Matrix3( 1,     0,     0,
                  0,     0,     1,
                  0,    -1,     0)
    ,
    Matrix3( 0,    -1,     0,
                  1,     0,     0,
                  0,     0,     1)
    ,
    Matrix3( 0,     0,     1,
                  0,    -1,     0,
                  1,     0,     0)
    ,
    Matrix3(-1,     0,     0,
                  0,     0,     1,
                  0,     1,     0)
    ,
    Matrix3( 0,     1,     0,
                 -1,     0,     0,
                  0,     0,     1)
    ,
    Matrix3( 0,     0,     1,
                  0,     1,     0,
                 -1,     0,     0)
};

const Matrix3 ML_rgmIcosahedralGroup[60] = {
Matrix3( 0,     1,     0,
              0,     0,     1,
              1,     0,     0)
,
Matrix3( 1,     0,     0,
              0,     1,     0,
              0,     0,     1)
,
Matrix3( 0,     0,     1,
              1,     0,     0,
              0,     1,     0)
,
Matrix3( 0,     0,    -1,
             -1,     0,     0,
              0,     1,     0)
,
Matrix3( 0,     1,     0,
              0,     0,    -1,
             -1,     0,     0)
,
Matrix3(-1,     0,     0,
              0,     1,     0,
              0,     0,    -1)
,
Matrix3( 1,     0,     0,
              0,    -1,     0,
              0,     0,    -1)
,
Matrix3( 0,     0,    -1,
              1,     0,     0,
              0,    -1,     0)
,
Matrix3( 0,    -1,     0,
              0,     0,    -1,
              1,     0,     0)
,
Matrix3(-1,     0,     0,
              0,    -1,     0,
              0,     0,     1)
,
Matrix3( 0,     0,     1,
             -1,     0,     0,
              0,    -1,     0)
,
Matrix3( 0,    -1,     0,
              0,     0,     1,
             -1,     0,     0)
,
Matrix3(0.5*D_GOLD, -0.5*(D_GOLD-1.0f),   0.5,
          -0.5*(D_GOLD-1.0f),   0.5, 0.5*D_GOLD,
           -0.5, -0.5*D_GOLD, 0.5*(D_GOLD-1.0f))
,
Matrix3( -0.5, -0.5*D_GOLD, 0.5*(D_GOLD-1.0f),
          0.5*D_GOLD, -0.5*(D_GOLD-1.0f),   0.5,
          -0.5*(D_GOLD-1.0f),   0.5, 0.5*D_GOLD)
,
Matrix3(-0.5*(D_GOLD-1.0f),   0.5, 0.5*D_GOLD,
           -0.5, -0.5*D_GOLD, 0.5*(D_GOLD-1.0f),
          0.5*D_GOLD, -0.5*(D_GOLD-1.0f),   0.5)
,
Matrix3(-0.5*D_GOLD, -0.5*(D_GOLD-1.0f),   0.5,
          0.5*(D_GOLD-1.0f),   0.5, 0.5*D_GOLD,
           -0.5, 0.5*D_GOLD, -0.5*(D_GOLD-1.0f))
,
Matrix3( -0.5, 0.5*D_GOLD, -0.5*(D_GOLD-1.0f),
          -0.5*D_GOLD, -0.5*(D_GOLD-1.0f),   0.5,
          0.5*(D_GOLD-1.0f),   0.5, 0.5*D_GOLD)
,
Matrix3(0.5*(D_GOLD-1.0f),   0.5, 0.5*D_GOLD,
           -0.5, 0.5*D_GOLD, -0.5*(D_GOLD-1.0f),
          -0.5*D_GOLD, -0.5*(D_GOLD-1.0f),   0.5)
,
Matrix3(-0.5*D_GOLD, -0.5*(D_GOLD-1.0f),  -0.5,
          -0.5*(D_GOLD-1.0f),  -0.5, 0.5*D_GOLD,
           -0.5, 0.5*D_GOLD, 0.5*(D_GOLD-1.0f))
,
Matrix3( -0.5, 0.5*D_GOLD, 0.5*(D_GOLD-1.0f),
          -0.5*D_GOLD, -0.5*(D_GOLD-1.0f),  -0.5,
          -0.5*(D_GOLD-1.0f),  -0.5, 0.5*D_GOLD)
,
Matrix3(-0.5*(D_GOLD-1.0f),  -0.5, 0.5*D_GOLD,
           -0.5, 0.5*D_GOLD, 0.5*(D_GOLD-1.0f),
          -0.5*D_GOLD, -0.5*(D_GOLD-1.0f),  -0.5)
,
Matrix3(0.5*D_GOLD, -0.5*(D_GOLD-1.0f),  -0.5,
          0.5*(D_GOLD-1.0f),  -0.5, 0.5*D_GOLD,
           -0.5, -0.5*D_GOLD, -0.5*(D_GOLD-1.0f))
,
Matrix3( -0.5, -0.5*D_GOLD, -0.5*(D_GOLD-1.0f),
          0.5*D_GOLD, -0.5*(D_GOLD-1.0f),  -0.5,
          0.5*(D_GOLD-1.0f),  -0.5, 0.5*D_GOLD)
,
Matrix3(0.5*(D_GOLD-1.0f),  -0.5, 0.5*D_GOLD,
           -0.5, -0.5*D_GOLD, -0.5*(D_GOLD-1.0f),
          0.5*D_GOLD, -0.5*(D_GOLD-1.0f),  -0.5)
,
Matrix3(0.5*(D_GOLD-1.0f),  -0.5, -0.5*D_GOLD,
            0.5, 0.5*D_GOLD, -0.5*(D_GOLD-1.0f),
          0.5*D_GOLD, -0.5*(D_GOLD-1.0f),   0.5)
,
Matrix3(0.5*D_GOLD, -0.5*(D_GOLD-1.0f),   0.5,
          0.5*(D_GOLD-1.0f),  -0.5, -0.5*D_GOLD,
            0.5, 0.5*D_GOLD, -0.5*(D_GOLD-1.0f))
,
Matrix3(  0.5, 0.5*D_GOLD, -0.5*(D_GOLD-1.0f),
          0.5*D_GOLD, -0.5*(D_GOLD-1.0f),   0.5,
          0.5*(D_GOLD-1.0f),  -0.5, -0.5*D_GOLD)
,
Matrix3(-0.5*(D_GOLD-1.0f),  -0.5, 0.5*D_GOLD,
            0.5, -0.5*D_GOLD, -0.5*(D_GOLD-1.0f),
          0.5*D_GOLD, 0.5*(D_GOLD-1.0f),   0.5)
,
Matrix3(0.5*D_GOLD, 0.5*(D_GOLD-1.0f),   0.5,
          -0.5*(D_GOLD-1.0f),  -0.5, 0.5*D_GOLD,
            0.5, -0.5*D_GOLD, -0.5*(D_GOLD-1.0f))
,
Matrix3(  0.5, -0.5*D_GOLD, -0.5*(D_GOLD-1.0f),
          0.5*D_GOLD, 0.5*(D_GOLD-1.0f),   0.5,
          -0.5*(D_GOLD-1.0f),  -0.5, 0.5*D_GOLD)
,
Matrix3(-0.5*(D_GOLD-1.0f),  -0.5, -0.5*D_GOLD,
           -0.5, 0.5*D_GOLD, -0.5*(D_GOLD-1.0f),
          0.5*D_GOLD, 0.5*(D_GOLD-1.0f),  -0.5)
,
Matrix3(0.5*D_GOLD, 0.5*(D_GOLD-1.0f),  -0.5,
          -0.5*(D_GOLD-1.0f),  -0.5, -0.5*D_GOLD,
           -0.5, 0.5*D_GOLD, -0.5*(D_GOLD-1.0f))
,
Matrix3( -0.5, 0.5*D_GOLD, -0.5*(D_GOLD-1.0f),
          0.5*D_GOLD, 0.5*(D_GOLD-1.0f),  -0.5,
          -0.5*(D_GOLD-1.0f),  -0.5, -0.5*D_GOLD)
,
Matrix3(  0.5, 0.5*D_GOLD, -0.5*(D_GOLD-1.0f),
          -0.5*D_GOLD, 0.5*(D_GOLD-1.0f),  -0.5,
          -0.5*(D_GOLD-1.0f),   0.5, 0.5*D_GOLD)
,
Matrix3(-0.5*(D_GOLD-1.0f),   0.5, 0.5*D_GOLD,
            0.5, 0.5*D_GOLD, -0.5*(D_GOLD-1.0f),
          -0.5*D_GOLD, 0.5*(D_GOLD-1.0f),  -0.5)
,
Matrix3(-0.5*D_GOLD, 0.5*(D_GOLD-1.0f),  -0.5,
          -0.5*(D_GOLD-1.0f),   0.5, 0.5*D_GOLD,
            0.5, 0.5*D_GOLD, -0.5*(D_GOLD-1.0f))
,
Matrix3(  0.5, 0.5*D_GOLD, 0.5*(D_GOLD-1.0f),
          0.5*D_GOLD, -0.5*(D_GOLD-1.0f),  -0.5,
          -0.5*(D_GOLD-1.0f),   0.5, -0.5*D_GOLD)
,
Matrix3(-0.5*(D_GOLD-1.0f),   0.5, -0.5*D_GOLD,
            0.5, 0.5*D_GOLD, 0.5*(D_GOLD-1.0f),
          0.5*D_GOLD, -0.5*(D_GOLD-1.0f),  -0.5)
,
Matrix3(0.5*D_GOLD, -0.5*(D_GOLD-1.0f),  -0.5,
          -0.5*(D_GOLD-1.0f),   0.5, -0.5*D_GOLD,
            0.5, 0.5*D_GOLD, 0.5*(D_GOLD-1.0f))
,
Matrix3(0.5*(D_GOLD-1.0f),  -0.5, -0.5*D_GOLD,
           -0.5, -0.5*D_GOLD, 0.5*(D_GOLD-1.0f),
          -0.5*D_GOLD, 0.5*(D_GOLD-1.0f),  -0.5)
,
Matrix3(-0.5*D_GOLD, 0.5*(D_GOLD-1.0f),  -0.5,
          0.5*(D_GOLD-1.0f),  -0.5, -0.5*D_GOLD,
           -0.5, -0.5*D_GOLD, 0.5*(D_GOLD-1.0f))
,
Matrix3( -0.5, -0.5*D_GOLD, 0.5*(D_GOLD-1.0f),
          -0.5*D_GOLD, 0.5*(D_GOLD-1.0f),  -0.5,
          0.5*(D_GOLD-1.0f),  -0.5, -0.5*D_GOLD)
,
Matrix3(0.5*(D_GOLD-1.0f),   0.5, -0.5*D_GOLD,
           -0.5, 0.5*D_GOLD, 0.5*(D_GOLD-1.0f),
          0.5*D_GOLD, 0.5*(D_GOLD-1.0f),   0.5)
,
Matrix3(0.5*D_GOLD, 0.5*(D_GOLD-1.0f),   0.5,
          0.5*(D_GOLD-1.0f),   0.5, -0.5*D_GOLD,
           -0.5, 0.5*D_GOLD, 0.5*(D_GOLD-1.0f))
,
Matrix3( -0.5, 0.5*D_GOLD, 0.5*(D_GOLD-1.0f),
          0.5*D_GOLD, 0.5*(D_GOLD-1.0f),   0.5,
          0.5*(D_GOLD-1.0f),   0.5, -0.5*D_GOLD)
,
Matrix3(-0.5*(D_GOLD-1.0f),   0.5, -0.5*D_GOLD,
           -0.5, -0.5*D_GOLD, -0.5*(D_GOLD-1.0f),
          -0.5*D_GOLD, 0.5*(D_GOLD-1.0f),   0.5)
,
Matrix3(-0.5*D_GOLD, 0.5*(D_GOLD-1.0f),   0.5,
          -0.5*(D_GOLD-1.0f),   0.5, -0.5*D_GOLD,
           -0.5, -0.5*D_GOLD, -0.5*(D_GOLD-1.0f))
,
Matrix3( -0.5, -0.5*D_GOLD, -0.5*(D_GOLD-1.0f),
          -0.5*D_GOLD, 0.5*(D_GOLD-1.0f),   0.5,
          -0.5*(D_GOLD-1.0f),   0.5, -0.5*D_GOLD)
,
Matrix3(  0.5, -0.5*D_GOLD, 0.5*(D_GOLD-1.0f),
          0.5*D_GOLD, 0.5*(D_GOLD-1.0f),  -0.5,
          0.5*(D_GOLD-1.0f),   0.5, 0.5*D_GOLD)
,
Matrix3(0.5*(D_GOLD-1.0f),   0.5, 0.5*D_GOLD,
            0.5, -0.5*D_GOLD, 0.5*(D_GOLD-1.0f),
          0.5*D_GOLD, 0.5*(D_GOLD-1.0f),  -0.5)
,
Matrix3(0.5*D_GOLD, 0.5*(D_GOLD-1.0f),  -0.5,
          0.5*(D_GOLD-1.0f),   0.5, 0.5*D_GOLD,
            0.5, -0.5*D_GOLD, 0.5*(D_GOLD-1.0f))
,
Matrix3(-0.5*D_GOLD, -0.5*(D_GOLD-1.0f),   0.5,
          -0.5*(D_GOLD-1.0f),  -0.5, -0.5*D_GOLD,
            0.5, -0.5*D_GOLD, 0.5*(D_GOLD-1.0f))
,
Matrix3(  0.5, -0.5*D_GOLD, 0.5*(D_GOLD-1.0f),
          -0.5*D_GOLD, -0.5*(D_GOLD-1.0f),   0.5,
          -0.5*(D_GOLD-1.0f),  -0.5, -0.5*D_GOLD)
,
Matrix3(-0.5*(D_GOLD-1.0f),  -0.5, -0.5*D_GOLD,
            0.5, -0.5*D_GOLD, 0.5*(D_GOLD-1.0f),
          -0.5*D_GOLD, -0.5*(D_GOLD-1.0f),   0.5)
,
Matrix3(0.5*(D_GOLD-1.0f),  -0.5, 0.5*D_GOLD,
            0.5, 0.5*D_GOLD, 0.5*(D_GOLD-1.0f),
          -0.5*D_GOLD, 0.5*(D_GOLD-1.0f),   0.5)
,
Matrix3(-0.5*D_GOLD, 0.5*(D_GOLD-1.0f),   0.5,
          0.5*(D_GOLD-1.0f),  -0.5, 0.5*D_GOLD,
            0.5, 0.5*D_GOLD, 0.5*(D_GOLD-1.0f))
,
Matrix3(  0.5, 0.5*D_GOLD, 0.5*(D_GOLD-1.0f),
          -0.5*D_GOLD, 0.5*(D_GOLD-1.0f),   0.5,
          0.5*(D_GOLD-1.0f),  -0.5, 0.5*D_GOLD)
,
Matrix3(  0.5, -0.5*D_GOLD, -0.5*(D_GOLD-1.0f),
          -0.5*D_GOLD, -0.5*(D_GOLD-1.0f),  -0.5,
          0.5*(D_GOLD-1.0f),   0.5, -0.5*D_GOLD)
,
Matrix3(0.5*(D_GOLD-1.0f),   0.5, -0.5*D_GOLD,
            0.5, -0.5*D_GOLD, -0.5*(D_GOLD-1.0f),
          -0.5*D_GOLD, -0.5*(D_GOLD-1.0f),  -0.5)
,
Matrix3(-0.5*D_GOLD, -0.5*(D_GOLD-1.0f),  -0.5,
          0.5*(D_GOLD-1.0f),   0.5, -0.5*D_GOLD,
            0.5, -0.5*D_GOLD, -0.5*(D_GOLD-1.0f))
};
//
//  Test/Verify the polyhedra and groups
//
#define ML_Assert(__expr)                                       \
      do {                                                      \
        if (!(__expr)) {                                        \
            printf("ML_Assert failed in %s\n",__FILE__);        \
            printf(" line %d: (%s)\n", __LINE__, #__expr);      \
             __debugbreak();                                    \
        }                                                       \
      } while (0) 

static void
TestPolyhedron(int c, const Vector3 rgpnt[])
{
    int i;
    Vector3 vec;
    Vector3 pntCent(0.0, 0.0, 0.0);

    for (i = 0; i < c; i++) {
        vec = rgpnt[i];
        ML_Assert(fabs(Dot(vec, vec) - 1.0) < 0.00001); // unit distant verts
        pntCent = pntCent + vec;
    }
    vec = pntCent / double(c);
    ML_Assert(Dot(vec, vec) < 0.00001);     // centroid is origin
}

static float fEps = 0.00001f;

static int
Approx(const Vector3 &vec1, const Vector3 &vec2)
{
    return Dot(vec1 - vec2, vec1 - vec2) < 0.00001;
}

static int
Approx(const Matrix3 &mat1, const Matrix3 &mat2)
{
    int i, j;
    double fErr = 0.0;

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            fErr += fabs(mat1.m[i][j] - mat2.m[i][j]);
    return fErr <= fEps;
}
//
//  The regular isometry groups have to be tested very carefully
//  since constructing them is tricky.
//
static void
IsaGroup(int c, const Matrix3 rgmat[])
{
    Matrix3 mat;
    int i, j, k;
    int cIdents, cInverts;

    cIdents = 0;
    for (i = 0; i < c; i++)
        if (Approx(rgmat[i], Matrix3(1.0)))      // Identity exists
            cIdents++;
    ML_Assert(cIdents == 1);
    for (i = 0; i < c; i++) {
        cInverts = 0;
        for (j = 0; j < c; j++)
            if (Approx(rgmat[i]*rgmat[j], Matrix3(1.0)))     // Inverses
                cInverts++;
        ML_Assert(cInverts == 1);
    }
    for (i = 0; i < c; i++) {
        for (j = 0; j < c; j++) {
            mat = rgmat[i]*rgmat[j];
            for (k = 0; k < c; k++)
                if (Approx(mat, rgmat[k]))      // Closure
                    goto FoundIt;
            ML_Assert("not a group" == 0);
FoundIt:
            continue;
        }
    }
}

static void
IsaSymmetry(int cMats, const Matrix3 rgmat[], int cPnts, const Vector3 rgpnt[])
{
    Vector3 pnt;
    int iPnt, jPnt, iMat;

    for (iMat = 0; iMat < cMats; iMat++) {
        for (iPnt = 0; iPnt < cPnts; iPnt++) {
            pnt = rgmat[iMat]*rgpnt[iPnt];
            for (jPnt = 0; jPnt < cPnts; jPnt++)
                if (Approx(pnt, rgpnt[jPnt]))       // Permutation of vertices
                    goto FoundIt;
            ML_Assert("not a symmetry" == 0);
FoundIt:
            continue;
        }
    }
}

void
TestPolyhedra()
{
    printf("Testing polyhedra and symmetry groups\n");
    printf("  Tetrahedral Group\n");
    TestPolyhedron(4, ML_rgvTetrahedron);
    fEps = 0;   // should be exact, just 0, 1, -1 entries in matrices
    IsaGroup(12, ML_rgmTetrahedralGroup);
    fEps = 0.000001f;
    IsaSymmetry(12, ML_rgmTetrahedralGroup, 4, ML_rgvTetrahedron);
    printf("  Octahedral Group\n");
    TestPolyhedron(8, ML_rgvCube);
    fEps = 0;
    IsaGroup(24, ML_rgmOctahedralGroup);
    fEps = 0.000001f;
    IsaSymmetry(24, ML_rgmOctahedralGroup, 8, ML_rgvCube);
    TestPolyhedron(6, ML_rgvOctahedron);
    IsaSymmetry(24, ML_rgmOctahedralGroup, 6, ML_rgvOctahedron);
    printf("  Icosahedral Group\n");
    TestPolyhedron(20, ML_rgvDodecahedron);
    IsaGroup(60, ML_rgmIcosahedralGroup);
    IsaSymmetry(60, ML_rgmIcosahedralGroup, 20, ML_rgvDodecahedron);
    TestPolyhedron(12, ML_rgvIcosahedron);
    IsaSymmetry(60, ML_rgmIcosahedralGroup, 12, ML_rgvIcosahedron);
}
