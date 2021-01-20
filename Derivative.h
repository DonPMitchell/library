#pragma once
////
//  derivative - automatic-differentiation arithmetic
//  D. P. Mitchell
//
#pragma intrinsic (fabs, sqrt, log, exp, sin, cos, tan, asin, acos, atan, atan2, sinh, cosh, tanh)

template <class S, class V>
class Derivative {
public:
                        Derivative() {}
                        //Derivative(double f) : x((S)f) { Dx = ZeroValue(Dx); }
                        Derivative(S s) : x(s), Dx(0.0f) {  }
                        Derivative(S s, V v) : x((S)s), Dx(v) {}

    Derivative<S,V>&    operator +=(Derivative<S,V>&);
    Derivative<S,V>&    operator +=(S);
    Derivative<S,V>&    operator -=(Derivative<S,V>&);
    Derivative<S,V>&    operator -=(S);
    Derivative<S,V>&    operator *=(Derivative<S,V>&);
    Derivative<S,V>&    operator *=(S);
    Derivative<S,V>&    operator /=(Derivative<S,V>&);
    Derivative<S,V>&    operator /=(S);

    S   x;      // scalar value
    V   Dx;     // vector derivative value
};
//
//  Some common instances of the template class
//
#include "MatrixVector.h"

typedef Derivative<float, float>            DERIV_FF;
typedef Derivative<float, ML_FVector2>      DERIV_FV2;
typedef Derivative<float, ML_FVector3>      DERIV_FV3;

typedef Derivative<double, double>          DERIV_DD;
typedef Derivative<double, ML_DVector2>     DERIV_DV2;
typedef Derivative<double, ML_DVector3>     DERIV_DV3;

//extern float   ZeroValue(float);
//extern ML_FVector2 ZeroValue(ML_FVector2);
//extern ML_FVector3 ZeroValue(ML_FVector3);

//inline double subset(double a, double b) { return a == b; } // compatability with intervals

template <class S, class V>
inline void
ML_Print(Derivative<S, V> c)
{
    printf("< ");
    ML_Print(c.x);
    printf(", ");
    ML_Print(c.Dx);
    printf(" >");
}

template <class S, class V>
int
subset(Derivative<S,V> a, Derivative<S,V> b)
{
    return subset(a.x, b.x);
}

template <class S, class V>
int
subset(S a, Derivative<S,V> b)
{
    return subset(a, b.x);
}

template<class S, class V>
inline Derivative<S,V>
operator +(Derivative<S,V> &a, Derivative<S,V> &b)
{
	return Derivative<S,V>(a.x + b.x, a.Dx + b.Dx);
}

template<class S, class V>
inline Derivative<S,V>
operator +(Derivative<S,V> &a, S b)
{
	return Derivative<S,V>(a.x + b, a.Dx);
}

template<class S, class V>
inline Derivative<S,V>
operator +(S a, Derivative<S,V> &b)
{
	return Derivative<S,V>(a + b.x, b.Dx);
}

template<class S, class V>
inline Derivative<S,V>
operator -(Derivative<S,V> &a, Derivative<S,V> &b)
{
	return Derivative<S,V>(a.x - b.x, a.Dx - b.Dx);
}

template<class S, class V>
inline Derivative<S,V>
operator -(Derivative<S,V> &a, S b)
{
	return Derivative<S,V>(a.x - b, a.Dx);
}

template<class S, class V>
inline Derivative<S,V>
operator -(S a, Derivative<S,V> &b)
{
	return Derivative<S,V>(a - b.x, b.Dx);
}

template<class S, class V>
inline Derivative<S,V>
operator -(Derivative<S,V> &a)
{
	return Derivative<S,V>(-a.x, -a.Dx);
}

template<class S, class V>
inline Derivative<S,V>
operator +(Derivative<S,V> &a)
{
	return a;
}

template<class S, class V>
inline Derivative<S,V>
operator *(Derivative<S,V> &a, Derivative<S,V> &b)
{
	return Derivative<S,V>(a.x*b.x, a.x*b.Dx + b.x*a.Dx);
}

template<class S, class V>
inline Derivative<S,V>
operator *(Derivative<S,V> &a, S b)
{
	return Derivative<S,V>(a.x*b, a.Dx*b);
}

template<class S, class V>
inline Derivative<S,V>
operator *(S a, Derivative<S,V> &b)
{
	return Derivative<S,V>(a*b.x, a*b.Dx);
}

template<class S, class V>
inline Derivative<S,V>
operator /(Derivative<S,V> &a, Derivative<S,V> &b)
{
	return Derivative<S,V>(a.x/b.x, (a.Dx*b.x - a.x*b.Dx)/(b.x*b.x));
}

template<class S, class V>
inline Derivative<S,V>
operator /(Derivative<S,V> &a, S b)
{
	return Derivative<S,V>(a.x/b, a.Dx/b);
}

template<class S, class V>
inline Derivative<S,V>
operator /(S a, Derivative<S,V> &b)
{
	return Derivative<S,V>(a/b.x, -a*b.Dx/(b.x*b.x));
}

template<class S, class V>
inline Derivative<S,V>&
Derivative<S,V>::operator +=(Derivative<S,V> &a)
{
	x += a.x;
	Dx += a.Dx;
    return *this;
}

template<class S, class V>
inline Derivative<S,V>&
Derivative<S,V>::operator +=(S s)
{
	x += s;
    return *this;
}

template<class S, class V>
inline Derivative<S,V>&
Derivative<S,V>::operator -=(Derivative<S,V> &a)
{
	x -= a.x;
	Dx -= a.Dx;
    return *this;
}

template<class S, class V>
inline Derivative<S,V>&
Derivative<S,V>::operator -=(S s)
{
	x -= s;
    return *this;
}

template<class S, class V>
inline Derivative<S,V>&
Derivative<S,V>::operator *=(Derivative<S,V> &a)
{
	Dx = a.x*Dx + x*a.Dx;
    x *= a.x;
    return *this;
}

template<class S, class V>
inline Derivative<S,V>&
Derivative<S,V>::operator *=(S s)
{
	x *= s;
    Dx *= s;
    return *this;
}

template<class S, class V>
inline Derivative<S,V>&
Derivative<S,V>::operator /=(Derivative<S,V> &a)
{
	Dx = (Dx*a.x - x*a.Dx)/(a.x*a.x);
    x /= a.x;
    return *this;
}

template<class S, class V>
inline Derivative<S,V>&
Derivative<S,V>::operator /=(S s)
{
	x /= s;
    Dx /= s;
    return *this;
}

template<class S, class V>
inline int
operator ==(Derivative<S,V> &a, Derivative<S,V> &b)
{
    return a.x == b.x;
}

template<class S, class V>
inline int
operator ==(S s, Derivative<S,V> &b)
{
    return s == b.x;
}

template<class S, class V>
inline int
operator ==(Derivative<S,V> &a, S s)
{
    return a.x == s;
}

template<class S, class V>
inline int
operator !=(Derivative<S,V> &a, Derivative<S,V> &b)
{
    return a.x != b.x;
}

template<class S, class V>
inline int
operator !=(S s, Derivative<S,V> &b)
{
    return s != b.x;
}

template<class S, class V>
inline int
operator !=(Derivative<S,V> &a, S s)
{
    return a.x != s;
}
template<class S, class V>
inline int
operator <(Derivative<S,V> &a, Derivative<S,V> &b)
{
    return a.x < b.x;
}

template<class S, class V>
inline int
operator <(S s, Derivative<S,V> &b)
{
    return s < b.x;
}

template<class S, class V>
inline int
operator <(Derivative<S,V> &a, S s)
{
    return a.x < s;
}
template<class S, class V>
inline int
operator <=(Derivative<S,V> &a, Derivative<S,V> &b)
{
    return a.x <= b.x;
}

template<class S, class V>
inline int
operator <=(S s, Derivative<S,V> &b)
{
    return s <= b.x;
}

template<class S, class V>
inline int
operator <=(Derivative<S,V> &a, S s)
{
    return a.x <= s;
}
template<class S, class V>
inline int
operator >(Derivative<S,V> &a, Derivative<S,V> &b)
{
    return a.x > b.x;
}

template<class S, class V>
inline int
operator >(S s, Derivative<S,V> &b)
{
    return s > b.x;
}

template<class S, class V>
inline int
operator >(Derivative<S,V> &a, S s)
{
    return a.x > s;
}
template<class S, class V>
inline int
operator >=(Derivative<S,V> &a, Derivative<S,V> &b)
{
    return a.x >= b.x;
}

template<class S, class V>
inline int
operator >=(S s, Derivative<S,V> &b)
{
    return s >= b.x;
}

template<class S, class V>
inline int
operator >=(Derivative<S,V> &a, S s)
{
    return a.x >= s;
}

template <class S, class V>
Derivative<S,V>
acos(Derivative<S,V> &a)
{
	return Derivative<S,V>((S)acos(a.x), -a.Dx/(S)sqrt((S)1.0 - a.x*a.x));
}

template <class S, class V>
Derivative<S,V>
asin(Derivative<S,V> &a)
{
	return Derivative<S,V>((S)asin(a.x), a.Dx/(S)sqrt((S)1.0 - a.x*a.x));
}

template <class S, class V>
Derivative<S,V>
atan(Derivative<S,V> &a)
{
	return Derivative<S,V>((S)atan(a.x), a.Dx/((S)1.0 + a.x*a.x));
}
//
//	Note: imaginary part of complex log(x(t) + iy(t))
//
template <class S, class V>
Derivative<S,V>
atan2(Derivative<S,V> &a, Derivative<S,V> &b)
{
	return Derivative<S,V>((S)atan2(a.x, b.x), (b.x*a.Dx - a.x*b.Dx)/(a.x*a.x + b.x*b.x));
}

template <class S, class V>
Derivative<S,V>
ceil(Derivative<S,V> &a)
{
	S c;

	c = (S)ceil(a.x);
	if (c == a.x)
		return Derivative<S,V>(c, a.Dx * (S)F_MACHINE_MAX);
	else
		return Derivative<S,V>(c, a.Dx * (S)0.0);
}

template <class S, class V>
Derivative<S,V>
cos(Derivative<S,V> &a)
{
	return Derivative<S,V>((S)cos(a.x), -(S)sin(a.x)*a.Dx);
}

template <class S, class V>
Derivative<S,V>
cosh(Derivative<S,V> &a)
{
	return Derivative<S,V>((S)cosh(a.x), (S)sinh(a.x)*a.Dx);
}

template <class S, class V>
Derivative<S,V>
exp(Derivative<S,V> &a)
{
	S e;

	e = (S)exp(a.x);
	return Derivative<S,V>(e, e*a.Dx);
}
//REVIEW: this doesn't work well with interval derivative types
template <class S, class V>
Derivative<S,V>
fabs(Derivative<S,V> &a)
{
	if (subset(0.0, a.x))
		return Derivative<S,V>(fabs(a.x), a.Dx * (S)F_MACHINE_MAX);
	else if (a.x < 0.0)
		return -a;
	else
		return a;
}

template <class S, class V>
Derivative<S,V>
floor(Derivative<S,V> &a)
{
        S f;

        f = (S)floor(a.x);
        if (f == a.x)
                return Derivative<S,V>(f, a.Dx * (S)F_MACHINE_MAX);
        else
                return Derivative<S,V>(f, a.Dx * (S)0.0);
}

template <class S, class V>
Derivative<S,V>
hypot(Derivative<S,V> &a, Derivative<S,V> &b)
{
	S r;

	r = (S)hypot(a.x, b.x);
	return Derivative<S,V>(r, (a.x*a.Dx + b.x*b.Dx)/r);
}

template <class S, class V>
Derivative<S,V>
log(Derivative<S,V> &a)
{
	return Derivative<S,V>((S)log(a.x), a.Dx/a.x);
}

template <class S, class V>
Derivative<S,V>
pow(Derivative<S,V> &a, float b)
{
	return Derivative<S,V>((S)pow(a.x, b), b*(S)pow(a.x, b - (S)1.0)*a.Dx);
}

template <class S, class V>
Derivative<S,V>
sin(Derivative<S,V> &a)
{
	return Derivative<S,V>((S)sin(a.x), (S)cos(a.x)*a.Dx);
}

template <class S, class V>
Derivative<S,V>
sinh(Derivative<S,V> &a)
{
	return Derivative<S,V>((S)sinh(a.x), (S)cosh(a.x)*a.Dx);
}

template <class S, class V>
Derivative<S,V>
sqr(Derivative<S,V> &a)
{
	return Derivative<S,V>(a.x*a.x, (S)2.0*a.x*a.Dx);
}

template <class S, class V>
Derivative<S,V>
sqrt(Derivative<S,V> &a)
{
	S s;

	s = (S)sqrt(a.x);
	return Derivative<S,V>(s, (S)0.5*a.Dx/s);
}

template <class S, class V>
Derivative<S,V>
tan(Derivative<S,V> &a)
{
	S t;

	t = (S)tan(a.x);
	return Derivative<S,V>(t, a.Dx*((S)1.0 + t*t));
}

template <class S, class V>
Derivative<S,V>
tanh(Derivative<S,V> &a)
{
	S t;

	t = (S)tanh(a.x);
	return Derivative<S,V>(t, a.Dx*((S)1.0 - t*t));
}

template <class S, class V>
Derivative<S,V>
asinh(Derivative<S,V> &a)
{
    return Derivative<S,V>(asinh(a.x), a.Dx*(S)1.0/sqrt((S)1.0 + a.x*a.x));
}

template <class S, class V>
Derivative<S,V>
acosh(Derivative<S,V> &a)
{
    return Derivative<S,V>(acosh(a.x), a.Dx*(S)1.0/sqrt((a.x + (S)1.0)*(a.x - (S)1.0)));
}

template <class S, class V>
Derivative<S,V>
atanh(Derivative<S,V> &a)
{
    return Derivative<S,V>(atanh(a.x), a.Dx*(S)1.0/((S)1.0 - a.x*a.x));
}
