//
//  3D affine space: points, vectors, co-vectors, transformations
//  D.P. Mitchell  2018/09/29.
//
#pragma once
#pragma intrinsic (sqrt, sin, cos)
//
//  Vector space
//
struct Vector3 {
    union {
        struct { double x, y, z; };
        double rgf[3];
    };

            Vector3() {};
            Vector3(double f) : x(f), y(f), z(f) {};
            Vector3(double f1, double f2, double f3) : x(f1), y(f2), z(f3) {};

    Vector3  operator +(const Vector3 &v) const { return Vector3(x+v.x, y+v.y, z+v.z); }
    Vector3  operator -(const Vector3 &v) const { return Vector3(x-v.x, y-v.y, z-v.z); }
    Vector3  operator -()         const { return Vector3(-x, -y, -z); }
    Vector3  operator *(double f) const { return Vector3(x*f, y*f, z*f); }
    Vector3  operator /(double f) const { return Vector3(x/f, y/f, z/f); }      // the compiler only computes 1/f once

    void    Print() const { printf("(%f %f %f)\n", x, y, z); }
};


inline double
Dot(const Vector3 &u, const Vector3 &v)
{
    return u.x*v.x + u.y*v.y + u.z*v.z;
}

inline Vector3
Cross(const Vector3 &u, const Vector3 &v)
{
    return Vector3(u.y*v.z-u.z*v.y, u.z*v.x-u.x*v.z, u.x*v.y-u.y*v.x);
}

inline double
Norm(const Vector3 &v)
{
    return sqrt(Dot(v, v));
}

inline Vector3
Normalize(const Vector3 &v)
{
    return v / Norm(v);
}
//
//  Matrix, linear transformations
//
struct Matrix3 {
    union {
        double   m[3][3];            // m[iRow][iColumn]
        struct {
            double f00, f01, f02;
            double f10, f11, f12;
            double f20, f21, f22;
        };
    };
                Matrix3() {};
                Matrix3(double x) :
                                f00(x),   f01(0.0), f02(0.0),
                                f10(0.0), f11(x),   f12(0.0),
                                f20(0.0), f21(0.0), f22(x)      {};
                Matrix3(double x00, double x01, double x02,
                       double x10, double x11, double x12,
                       double x20, double x21, double x22) :
                                f00(x00), f01(x01), f02(x02),
                                f10(x10), f11(x11), f12(x12),
                                f20(x20), f21(x21), f22(x22)    {};

    Matrix3  operator +(const Matrix3 &m2) const;
    Matrix3  operator -(const Matrix3 &m2) const;
    Matrix3  operator *(double f) const;
    Matrix3  operator /(double f) const;
    Matrix3  operator *(const Matrix3 &m2) const;
    Vector3  operator *(const Vector3 &v) const;        // v' = M*v

    double  Determinant() const;
    void    Print() const { printf("[%9f %9f %9f]\n[%9f %9f %9f]\n[%9f %9f %9f]\n", 
                                f00, f01, f02, f10, f11, f12, f20, f21, f22); }
};

inline Matrix3
Matrix3::operator +(const Matrix3 &m2) const {
    return Matrix3 ( f00 + m2.f00, f01 + m2.f01, f02 + m2.f02,
                    f10 + m2.f10, f11 + m2.f11, f12 + m2.f12,
                    f20 + m2.f20, f21 + m2.f21, f22 + m2.f22);
}

inline Matrix3
Matrix3::operator -(const Matrix3 &m2) const {
    return Matrix3( f00 - m2.f00, f01 - m2.f01, f02 - m2.f02,
                    f10 - m2.f10, f11 - m2.f11, f12 - m2.f12,
                    f20 - m2.f20, f21 - m2.f21, f22 - m2.f22);
}

inline Matrix3
Matrix3::operator *(double f) const {
    return Matrix3( f00 * f, f01 * f, f02 * f,
                    f10 * f, f11 * f, f12 * f,
                    f20 * f, f21 * f, f22 * f);
}

inline Matrix3
Matrix3::operator /(double f) const {
    return *this * (1.0/f);
}

inline Matrix3
Matrix3::operator *(const Matrix3 &m2) const
{
    return Matrix3(
        m[0][0]*m2.m[0][0] + m[0][1]*m2.m[1][0] + m[0][2]*m2.m[2][0],
        m[0][0]*m2.m[0][1] + m[0][1]*m2.m[1][1] + m[0][2]*m2.m[2][1],
        m[0][0]*m2.m[0][2] + m[0][1]*m2.m[1][2] + m[0][2]*m2.m[2][2],

        m[1][0]*m2.m[0][0] + m[1][1]*m2.m[1][0] + m[1][2]*m2.m[2][0],
        m[1][0]*m2.m[0][1] + m[1][1]*m2.m[1][1] + m[1][2]*m2.m[2][1],
        m[1][0]*m2.m[0][2] + m[1][1]*m2.m[1][2] + m[1][2]*m2.m[2][2],

        m[2][0]*m2.m[0][0] + m[2][1]*m2.m[1][0] + m[2][2]*m2.m[2][0],
        m[2][0]*m2.m[0][1] + m[2][1]*m2.m[1][1] + m[2][2]*m2.m[2][1],
        m[2][0]*m2.m[0][2] + m[2][1]*m2.m[1][2] + m[2][2]*m2.m[2][2]
    );
}

inline Vector3
Matrix3::operator *(const Vector3 &vec) const
{
return Vector3(vec.x*m[0][0] + vec.y*m[0][1] + vec.z*m[0][2],
               vec.x*m[1][0] + vec.y*m[1][1] + vec.z*m[1][2],
               vec.x*m[2][0] + vec.y*m[2][1] + vec.z*m[2][2]);
}

inline double
Matrix3::Determinant() const
{
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
         + m[0][1] * (m[1][2] * m[2][0] - m[1][0] * m[2][2])
         + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

inline Matrix3
InverseMatrix(const Matrix3 &mat)
{
    double fDet, fDetInv, f1, f2, f3;

    fDet = mat.m[0][0] * ( f1 = mat.m[1][1] * mat.m[2][2] - mat.m[1][2] * mat.m[2][1] )
         + mat.m[0][1] * ( f2 = mat.m[1][2] * mat.m[2][0] - mat.m[1][0] * mat.m[2][2] )
         + mat.m[0][2] * ( f3 = mat.m[1][0] * mat.m[2][1] - mat.m[1][1] * mat.m[2][0] );
    if (fDet)
        fDetInv = (1.0/fDet);
    else
        fDetInv = 0.0;
    return Matrix3(
          f1 * fDetInv,
        - ( mat.m[0][1] * mat.m[2][2] - mat.m[0][2] * mat.m[2][1] ) * fDetInv,
          ( mat.m[0][1] * mat.m[1][2] - mat.m[0][2] * mat.m[1][1] ) * fDetInv,

          f2 * fDetInv,
          ( mat.m[0][0] * mat.m[2][2] - mat.m[0][2] * mat.m[2][0] ) * fDetInv,
        - ( mat.m[0][0] * mat.m[1][2] - mat.m[0][2] * mat.m[1][0] ) * fDetInv,

          f3 * fDetInv,
        - ( mat.m[0][0] * mat.m[2][1] - mat.m[0][1] * mat.m[2][0] ) * fDetInv,
          ( mat.m[0][0] * mat.m[1][1] - mat.m[0][1] * mat.m[1][0] ) * fDetInv
    );
}


inline Matrix3
ScaleMatrix(double x, double y, double z)
{
    return Matrix3(  x, 0.0, 0.0,
                   0.0,   y, 0.0,
                   0.0, 0.0,   z);
}
//
//  Rotation by Euler angles
//
inline void
RotateBasisAxis(int i1, int i2, double fRadians, Matrix3 &mat)
{
    double s, c, x;
    int j;

    s = sin(fRadians);
    c = cos(fRadians);
    //
    //  efficient multiplication by an axis-rotation matrix
    //
    for (j = 0; j < 3; j++) {
        x =            c * mat.m[i1][j] - s * mat.m[i2][j];
        mat.m[i2][j] = s * mat.m[i1][j] + c * mat.m[i2][j];
        mat.m[i1][j] = x;
    }    
}

inline void
InitialBasisAxis(int i1, int i2, double fRadians, Matrix3 &mat)
{
    double s, c;

    s = sin(fRadians);
    c = cos(fRadians);
    mat.m[i1][i1] = c;
    mat.m[i1][i2] = -s;
    mat.m[i2][i1] = s;
    mat.m[i2][i2] = c;
}

inline Matrix3
EulerMatrix(double fRadians0, double fRadians1, double fRadians2, char *szConvention = "xyz")
{
    Matrix3 mat(1.0);

    switch(szConvention[0]) {
        case 'x':   InitialBasisAxis(1, 2, fRadians0, mat);
                    break;
        case 'y':   InitialBasisAxis(2, 0, fRadians0, mat);
                    break;
        case 'z':   InitialBasisAxis(0, 1, fRadians0, mat);
                    break;
    }
    switch(szConvention[1]) {
        case 'x':   RotateBasisAxis(1, 2, fRadians1, mat);
                    break;
        case 'y':   RotateBasisAxis(2, 0, fRadians1, mat);
                    break;
        case 'z':   RotateBasisAxis(0, 1, fRadians1, mat);
                    break;
    }
    switch(szConvention[2]) {
        case 'x':   RotateBasisAxis(1, 2, fRadians2, mat);
                    break;
        case 'y':   RotateBasisAxis(2, 0, fRadians2, mat);
                    break;
        case 'z':   RotateBasisAxis(0, 1, fRadians2, mat);
                    break;
    }
    return mat;
}
//
//  Dual vector space
//
struct CoVector3 {
    union {
        struct { double x, y, z; };
        double rgf[3];
    };

            CoVector3() {};
            CoVector3(double f) : x(f), y(f), z(f) {};
            CoVector3(double f1, double f2, double f3) : x(f1), y(f2), z(f3) {};

    CoVector3  operator +(const CoVector3 &v) const { return CoVector3(x+v.x, y+v.y, z+v.z); }
    CoVector3  operator -(const CoVector3 &v) const { return CoVector3(x-v.x, y-v.y, z-v.z); }
    CoVector3  operator -()          const { return CoVector3(-x, -y, -z); }
    CoVector3  operator *(double f)  const { return CoVector3(x*f, y*f, z*f); }
    CoVector3  operator /(double f)  const { return CoVector3(x/f, y/f, z/f); }
    CoVector3  operator *(const Matrix3 &m2) const;       // v' = v*M
    double     operator *(const Vector3 &v) const { return x*v.x + y*v.y + z*v.z; }

    void    Print() const { printf("(%f %f %f)\n", x, y, z); }
};


inline double
Dot(const CoVector3 &u, const CoVector3 &v)
{
    return u.x*v.x + u.y*v.y + u.z*v.z;
}

inline CoVector3
Cross(const CoVector3 &u, const CoVector3 &v)
{
    return CoVector3(u.y*v.z-u.z*v.y, u.z*v.x-u.x*v.z, u.x*v.y-u.y*v.x);
}

inline double
Norm(const CoVector3 &v)
{
    return sqrt(Dot(v, v));
}

inline CoVector3
Normalize(const CoVector3 &v)
{
    return v / Norm(v);
}

inline CoVector3
CoVector3::operator *(const Matrix3 &m) const
{
return CoVector3( x*m.m[0][0] + y*m.m[1][0] + z*m.m[2][0],
                  x*m.m[0][1] + y*m.m[1][1] + z*m.m[2][1],
                  x*m.m[0][2] + y*m.m[1][2] + z*m.m[2][2] );
}
//
//  Affine point space
//
struct Point3 {
    union {
        struct { double x, y, z; };
        double rgf[3];
    };
            Point3() : x(0.0), y(0.0), z(0.0) {};
            Point3(double a, double b, double c) : x(a), y(b), z(c) {};

    Point3   operator +(const Vector3 &v) const { return Point3(x+v.x, y+v.y, z+v.z); }
    Vector3  operator -(const Point3 &p)  const { return Vector3 (x-p.x, y-p.y, z-p.z); }
    void    Print() const { printf("(%f %f %f)\n", x, y, z); }
};

inline Point3
Blend(double u, const Point3 &p1, const Point3 &p2)
{
    return p1 + (p2 - p1)*u;
}
