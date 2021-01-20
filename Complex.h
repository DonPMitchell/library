#pragma once
//
//  Complex number template (for float or double precision)
//  D.P.Mitchell 04/18/2009.
//
#pragma intrinsic(fabs, sqrt, log, exp, sin, cos, tan, atan2, sinh, cosh, asin, acos)
#define D_PI_COMPLEX    3.14159265358979323846264338327950288419716939937510

template <class T>
class ML_Imaginary {
public:
    ML_Imaginary(const T imaginary) : im(imaginary) {}
    ML_Imaginary() {}

    T im;
};

template <class T>
class ML_Complex {
public:
        ML_Complex(const T real, const T imaginary);
        ML_Complex(const T real);
        ML_Complex() {}

    ML_Complex& operator +=(const ML_Complex<T>&);
    ML_Complex& operator -=(const ML_Complex<T>&);
    ML_Complex& operator *=(const ML_Complex<T>&);
    ML_Complex& operator /=(const ML_Complex<T>&);
    ML_Complex& operator +=(T);
    ML_Complex& operator -=(T);
    ML_Complex& operator *=(T);
    ML_Complex& operator /=(T);

    T   re;
    T   im;
};

typedef ML_Complex<float>  ML_FComplex;
typedef ML_Complex<double> ML_DComplex;

template <class T>
inline
ML_Complex<T>::ML_Complex(T real, T imaginary)
: re(real), im(imaginary)
{
}

template <class T>
inline
ML_Complex<T>::ML_Complex(T real)
: re(real), im(T(0.0))
{
}

template <class T>
inline void
ML_Print(ML_Complex<T> c)
{
    printf("( ");
    ML_Print(c.re);
    printf(", ");
    ML_Print(c.im);
    printf(" )");
}

void inline
ML_Print(float f)
{
    printf("%g", f);
}

void inline
ML_Print(double f)
{
    printf("%lg", f);
}

template<class T>
inline ML_Complex<T>
ML_Conjugate(ML_Complex<T> c)
{
    return ML_Complex<T>(c.re, -c.im);
}

template<class T>
inline T
ML_Real(ML_Complex<T> c)
{
    return c.re;
}

template<class T>
inline ML_Complex<T>&
ML_Complex<T>::operator +=(const ML_Complex<T> &a)
{
    re += a.re;
    im += a.im;
    return *this;
}

template<class T>
inline ML_Complex<T>&
ML_Complex<T>::operator -=(const ML_Complex<T> &a)
{
    re -= a.re;
    im -= a.im;
    return *this;
}

template<class T>
inline ML_Complex<T>&
ML_Complex<T>::operator *=(const ML_Complex<T> &a)
{
    T tReal;

    tReal = a.re*re - a.im*im;
    im = a.re*im + a.im*re;
    re = tReal;
    return *this;
}

template<class T>
inline ML_Complex<T>&
ML_Complex<T>::operator /=(const ML_Complex<T> &b)
{
    T fDenom, Q, tReal;

    if (fabs(b.re) < fabs(b.im)) {              // overflow-avoiding method
        Q = b.re/b.im;
        fDenom = b.re*Q + b.im;
        tReal = (re*Q + im)/fDenom;
        im = (im*Q - re)/fDenom;
        re = tReal;
    } else {
        Q = b.im/b.re;
        fDenom = b.im*Q + b.re;
        tReal = (im*Q + re)/fDenom;
        im = (im - re*Q)/fDenom;
        re = tReal;
    }
    return *this;
}

template<class T>
inline ML_Complex<T>&
ML_Complex<T>::operator +=(T t)
{
    re += t;
    return *this;
}

template<class T>
inline ML_Complex<T>&
ML_Complex<T>::operator -=(T t)
{
    re -= t;
    return *this;
}

template<class T>
inline ML_Complex<T>&
ML_Complex<T>::operator *=(T t)
{
    re *= t;
    im *= t;
    return *this;
}

template<class T>
inline ML_Complex<T>&
ML_Complex<T>::operator /=(T t)
{
    re /= t;
    im /= t;
    return *this;
}

template<class T>
inline ML_Complex<T>
operator +(const ML_Complex<T> &a, const ML_Complex<T> &b)
{
    return ML_Complex<T>(a.re + b.re, a.im + b.im);
}

template<class T>
inline ML_Complex<T>
operator +(const ML_Complex<T> &a, const T &t)
{
    return ML_Complex<T>(a.re + t, a.im);
}

template<class T>
inline ML_Complex<T>
operator +(const T &t, const ML_Complex<T> &b)
{
    return ML_Complex<T>(t + b.re, b.im);
}

template<class T>
inline ML_Complex<T>
operator -(const ML_Complex<T> &a, const ML_Complex<T> &b)
{
    return ML_Complex<T>(a.re - b.re, a.im - b.im);
}

template<class T>
inline ML_Complex<T>
operator -(const ML_Complex<T> &a, const T &t)
{
    return ML_Complex<T>(a.re - t, a.im);
}

template<class T>
inline ML_Complex<T>
operator -(const T &t, const ML_Complex<T> &b)
{
    return ML_Complex<T>(t - b.re, -b.im);
}

template<class T>
inline ML_Complex<T>
operator *(const ML_Complex<T> &a, const ML_Complex<T> &b)
{
    return ML_Complex<T>(a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re);
}

template<class T>
inline ML_Complex<T>
operator *(const ML_Complex<T> &a, const T &t)
{
    return ML_Complex<T>(a.re*t, a.im*t);
}

template<class T>
inline ML_Complex<T>
operator *(const T &t, const ML_Complex<T> &b)
{
    return ML_Complex<T>(t*b.re, t*b.im);
}

template<class T>
inline ML_Complex<T>
operator /(const ML_Complex<T> &a, const ML_Complex<T> &b)
{
    T fDenom, Q;

    if (fabs(b.re) < fabs(b.im)) {              // overflow-avoiding method
        Q = b.re/b.im;
        fDenom = b.re*Q + b.im;
        return ML_Complex<T>((a.re*Q + a.im)/fDenom, (a.im*Q - a.re)/fDenom);
    } else {
        Q = b.im/b.re;
        fDenom = b.im*Q + b.re;
        return ML_Complex<T>((a.im*Q + a.re)/fDenom, (a.im - a.re*Q)/fDenom);
    }
}

template<class T>
inline ML_Complex<T>
operator /(const ML_Complex<T> &a, const T &t)
{
    return ML_Complex<T>(a.re/t, a.im/t);
}

template<class T>
inline ML_Complex<T>
operator /(const T &t, const ML_Complex<T> &b)
{
    T fDenom;

    fDenom = b.re*b.re + b.im*b.im;
    return ML_Complex<T>(t*b.re/fDenom, -t*b.im/fDenom);
}

template<class T>
inline ML_Complex<T>
operator -(const ML_Complex<T> &a)
{
    return ML_Complex<T>(-a.re, -a.im);
}

template<class T>
inline ML_Complex<T>
operator +(const ML_Complex<T> &a)
{
    return a;
}
//
//  Comparison Operators
//
template<class T>
inline int
operator ==(const ML_Complex<T> &a, const ML_Complex<T> &b)
{
    return a.re == b.re && a.im == b.im;
}

template<class T>
inline int
operator ==(const ML_Complex<T> &a, T t)
{
    return a.re == t && a.im == 0.0;
}

template<class T>
inline int
operator ==(T t, const ML_Complex<T> &b)
{
    return b.re == t && b.im == 0.0;
}

template<class T>
inline int
operator !=(const ML_Complex<T> &a, const ML_Complex<T> &b)
{
    return a.re != b.re || a.im != b.im;
}

template<class T>
inline int
operator !=(const ML_Complex<T> &a, T t)
{
    return a.re != t || a.im != 0.0;
}

template<class T>
inline int
operator !=(T t, const ML_Complex<T> &b)
{
    return b.re != t || b.im != 0.0;
}

template<class T>
inline int
subset(const ML_Complex<T> &a, const ML_Complex<T> &b)
{
    return ML_Subset(a.re, b.re) && ML_Subset(a.im, b.im);
}
//
//  Standard math functions
//
template<class T>
inline T
fabs(const ML_Complex<T> &z)
{
    return _hypot(z.re, z.im);
}

template<class T>
ML_Complex<T>
sqrt(const ML_Complex<T> &z)
{
    T t;

    t = sqrt(T(0.5)*(fabs(z) + fabs(z.re)));        // avoid subtraction cancellation error
    if (z.re >= 0.0) {
        return ML_Complex<T>(t, T(0.5)*z.im/t);
    } else {
        if (z.im < 0.0)
            return ML_Complex<T>(T(-0.5)*z.im/t, -t);
        else
            return ML_Complex<T>(T(0.5)*z.im/t, t);
    }
}

template<class T>
ML_Complex<T>
log(const ML_Complex<T> &z)
{
    T fPhase;

    if (z.re < 0.0)
        fPhase = T(-D_PI_COMPLEX);
    else
        fPhase = atan2(z.im, z.re);
    return ML_Complex<T>(log(fabs(z)), fPhase);
}

template<class T>
ML_Complex<T>
exp(const ML_Complex<T> &z)
{
    T fRadius;

    fRadius = exp(z.re);
    return ML_Complex<T>(fRadius*cos(z.im), fRadius*sin(z.im));
}

template<class T>
ML_Complex<T>
cos(const ML_Complex<T> &z)
{
    return ML_Complex<T>(cos(z.re)*cosh(z.im), -sin(z.re)*sinh(z.im));
}

template<class T>
ML_Complex<T>
sin(const ML_Complex<T> &z)
{
    return ML_Complex<T>(sin(z.re)*cosh(z.im), cos(z.re)*sinh(z.im));
}
//REVIEW: why is tan less accurate than sin and cos?
template<class T>
ML_Complex<T>
tan(const ML_Complex<T> &z)
{
    T fDenom;

    fDenom = cos(z.re + z.re) + cosh(z.im + z.im);
    return ML_Complex<T>(sin(z.re + z.re)/fDenom, sinh(z.im + z.im)/fDenom);
}

template<class T>
ML_Complex<T>
sinh(const ML_Complex<T> &z)
{
    return T(0.5)*(exp(z) - exp(-z));
}

template<class T>
ML_Complex<T>
cosh(const ML_Complex<T> &z)
{
    return T(0.5)*(exp(z) + exp(-z));
}

template<class T>
ML_Complex<T>
tanh(const ML_Complex<T> &z)
{
    ML_Complex<T> z1, z2;

    z1 = exp(z);
    z2 = exp(-z);
    return (z1 - z2) / (z1 + z2);
}

template<class T>
ML_Complex<T>
atanh(const ML_Complex<T> &z)
{
    return T(0.5)*(log(T(1.0) + z) - log(T(1.0) - z));   // Kahan's proposal
}

template<class T>
ML_Complex<T>
asinh(const ML_Complex<T> &z)
{
    return log(z + sqrt(T(1.0) + z*z));
}

template<class T>
ML_Complex<T>
acosh(const ML_Complex<T> &z)
{
    return T(2.0)*log(sqrt(T(0.5)*(z + T(1.0))) + sqrt(T(0.5)*(z - T(1.0))));  // Kahan
}

template<class T>
ML_Complex<T>
asin(const ML_Complex<T> &z)
{
    ML_Complex<T> iz, y;

    iz = ML_Complex<T>(-z.im, z.re);
    y = asinh(iz);
    return ML_Complex<T>(y.im, -y.re);    // Kahan's recommended asinh(i*z) / i
}

template<class T>
ML_Complex<T>
acos(const ML_Complex<T> &z)
{
    return T(0.5*D_PI_COMPLEX) - asin(z);
}

template<class T>
ML_Complex<T>
atan(const ML_Complex<T> &z)
{
    ML_Complex<T> iz, y;

    iz = ML_Complex<T>(-z.im, z.re);
    y = atanh(iz);
    return ML_Complex<T>(y.im, -y.re);    // Kahan's recommended atanh(i*z) / i
}

template<class T>
ML_Complex<T>
atan2(const ML_Complex<T> &z, const ML_Complex<T> &x)
{
    ML_Complex<T> iz, y, r;

    r = z/x;
    iz = ML_Complex<T>(-r.im, r.re);
    y = atanh(iz);
    return ML_Complex<T>(y.im, -y.re);    // Kahan's recommended atanh(i*z) / i
}

template<class T>
ML_Complex<T>
pow(const ML_Complex<T> &b, T x)
{
    if (x == 0.0)
        return ML_Complex<T>(1.0);        // IEEE
    return exp(b*log(x));
}

template<class T>
ML_Complex<T>
cis(T x)
{
    return ML_DComplex(cos(x), sin(x));
}
//
//  Operations with pure real and imaginary argument
//
template<class T>
T
operator *(ML_Imaginary<T> im1, ML_Imaginary<T> im2)
{
    return -im1.im*im2.im;
}

template<class T>
inline ML_Complex<T>
operator *(const ML_Complex<T> &a, const ML_Imaginary<T> &b)
{
    return ML_Complex<T>(- a.im*b.im, a.re*b.im);
}

template<class T>
inline ML_Complex<T>
operator *(const ML_Imaginary<T> &a, const ML_Complex<T> &b)
{
    return ML_Complex<T>(- a.im*b.im, a.im*b.re);
}
