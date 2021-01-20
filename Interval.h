#pragma once
//
//  Interval arithmetic template
//  D.P.Mitchell 04/18/2009.
//
#pragma intrinsic(fabs, sqrt, log, exp, sin, cos, tan, atan2, sinh, cosh, asin, acos)
#define D_INTERVAL_PI            3.14159265358979323846264338327950288419716939937510
#define D_INTERVAL_2PI           6.28318530717958647692528676655900576839433879875021
#define D_INTERVAL_PIINV         0.31830988618379067153776752674502872406891929148091
#define D_INTERVAL_MACHINE_EPS   2.22044605e-016
#define D_INTERVAL_MACHINE_MAX   1.79e+308
#define D_INTERVAL_MACHINE_MIN   2.23e-308
#define F_INTERVAL_MACHINE_EPS   1.1920929e-007f
#define F_INTERVAL_MACHINE_MAX   3.40e+38f
#define F_INTERVAL_MACHINE_MIN   1.18e-38f
#define F_INTERVAL_PI            3.1415926535897932384626433832795028842f
#define F_INTERVAL_PIINV         0.3183098861837906715377675267450287241f

template <class T>
class ML_Interval {
public:
        ML_Interval(const T real, const T imaginary);
        ML_Interval(const T real);
        ML_Interval() {}

    ML_Interval& operator +=(const ML_Interval<T>&);
    ML_Interval& operator -=(const ML_Interval<T>&);
    ML_Interval& operator *=(const ML_Interval<T>&);
    ML_Interval& operator /=(const ML_Interval<T>&);
    ML_Interval& operator +=(T);
    ML_Interval& operator -=(T);
    ML_Interval& operator *=(T);
    ML_Interval& operator /=(T);

    T   lo;
    T   hi;
};

typedef ML_Interval<float> ML_FInterval;
typedef ML_Interval<double> ML_DInterval;

template <class T>
inline
ML_Interval<T>::ML_Interval(T low, T high)
: lo(low), hi(high)
{
}

template <class T>
inline
ML_Interval<T>::ML_Interval(T real)
: lo(real), hi(real)
{
}

template <class T>
inline void
ML_Print(ML_Interval<T> c)
{
    printf("[ ");
    ML_Print(c.lo);
    printf(", ");
    ML_Print(c.hi);
    printf(" ]");
}

template <class T>
inline ML_Interval<T>
Reciprocol(const ML_Interval<T> &b)
{
	if (b.lo*b.hi <= 0.0) {
		if (b.lo == 0.0 && b.hi > 0.0)
			return ML_Interval<T>(T(1.0)/b.hi, F_INTERVAL_MACHINE_MAX);
		else if (b.lo < 0.0 && b.hi == 0.0)
			return ML_Interval<T>(T(-F_INTERVAL_MACHINE_MAX), T(1.0)/b.lo);
		else
			return ML_Interval<T>(T(-F_INTERVAL_MACHINE_MAX), T(F_INTERVAL_MACHINE_MAX));
	} else
		return ML_Interval<T>(T(1.0)/b.hi, T(1.0)/b.lo);
}

template<class T>
inline ML_Interval<T>&
ML_Interval<T>::operator +=(const ML_Interval<T> &a)
{
    lo += a.lo;
    hi += a.hi;
    return *this;
}

template<class T>
inline ML_Interval<T>&
ML_Interval<T>::operator -=(const ML_Interval<T> &a)
{
    lo -= a.lo;
    hi -= a.hi;
    return *this;
}

template <class T>
inline ML_Interval<T>&
ML_Interval<T>::operator *=(const ML_Interval<T> &a)
{
    T f0, f1, f2, f3, fTemp;
    //
    //  Faster to just do four multiplications, then optimal min/max selection.
    //
    f0 = a.lo * lo;
    f1 = a.lo * hi;
    f2 = a.hi * lo;
    f3 = a.hi * hi;
    if (f0 > f1) {
        fTemp = f0;
        f0 = f1;
        f1 = fTemp;
    }
    if (f2 > f3) {
        fTemp = f2;
        f2 = f3;
        f3 = fTemp;
    }
    if (f1 > f3)
        f3 = f1;
    if (f2 < f0)
        f0 = f2;
    lo = f0;
    hi = f3;
    return *this;
}

template <class T>
inline ML_Interval<T>&
ML_Interval<T>::operator /=(const ML_Interval<T> &a)
{
    return *this *= Reciprocol(a);
}

template<class T>
inline ML_Interval<T>&
ML_Interval<T>::operator +=(T t)
{
    lo += t;
    hi += t;
    return *this;
}

template<class T>
inline ML_Interval<T>&
ML_Interval<T>::operator -=(T t)
{
    lo -= t;
    hi -= t;
    return *this;
}

template<class T>
inline ML_Interval<T>&
ML_Interval<T>::operator *=(T t)
{
    T tmp;

    lo *= t;
    hi *= t;
    if (t < 0.0) {
        tmp = lo;
        lo = hi;
        hi = tmp;
    }
    return *this;
}

template<class T>
inline ML_Interval<T>&
ML_Interval<T>::operator /=(T t)
{
    T tmp;

    lo /= t;
    hi /= t;
    if (t < 0.0) {
        tmp = lo;
        lo = hi;
        hi = tmp;
    }
    return *this;
}

template <class T>
inline int
operator ==(const ML_Interval<T> &a, const ML_Interval<T> &b)
{
	return a.lo == b.lo && a.hi == b.hi;
}

template <class T>
inline int
operator ==(const ML_Interval<T> &a, T f)
{
	return a.lo == f && a.hi == f;
}

template <class T>
inline int
operator ==(T f, const ML_Interval<T> &b)
{
	return f == b.lo && f == b.hi;
}

template <class T>
inline int
operator !=(const ML_Interval<T> &a, const ML_Interval<T> &b)
{
	return a.lo != b.lo || a.hi != b.hi;
}

template <class T>
inline int
operator !=(const ML_Interval<T> &a, T f)
{
	return a.lo != f || a.hi != f;
}

template <class T>
inline int
operator !=(T f, const ML_Interval<T> &b)
{
	return f != b.lo || f != b.hi;
}

template <class T>
inline int
operator <(const ML_Interval<T> &a, const ML_Interval<T> &b)
{
	return a.hi < b.lo;
}

template <class T>
inline int
operator <(const ML_Interval<T> &a, T f)
{
	return a.hi < f;
}

template <class T>
inline int
operator <(T f, const ML_Interval<T> &b)
{
	return f < b.lo;
}

template <class T>
inline int
operator <=(const ML_Interval<T> &a, const ML_Interval<T> &b)
{
	return a.hi <= b.lo;
}

template <class T>
inline int
operator <=(const ML_Interval<T> &a, T f)
{
	return a.hi <= f;
}

template <class T>
inline int
operator <=(T f, const ML_Interval<T> &b)
{
	return f <= b.lo;
}

template <class T>
inline int
operator >(const ML_Interval<T> &a, const ML_Interval<T> &b)
{
	return a.lo > b.hi;
}

template <class T>
inline int
operator >(const ML_Interval<T> &a, T b)
{
	return a.lo > b;
}

template <class T>
inline int
operator >(T a, const ML_Interval<T> &b)
{
	return a > b.hi;
}

template <class T>
inline int
operator >=(const ML_Interval<T> &a, const ML_Interval<T> &b)
{
	return a.lo >= b.hi;
}

template <class T>
inline int
operator >=(const ML_Interval<T> &a, T b)
{
	return a.lo >= b;
}

template <class T>
inline int
operator >=(T a, const ML_Interval<T> &b)
{
	return a >= b.hi;
}

template <class T>
inline int
subset(const ML_Interval<T> &a, const ML_Interval<T> &b)
{
	return a.lo >= b.lo && a.hi <= b.hi;
}

template <class T>
inline int
subset(T f, const ML_Interval<T> &b)
{
	return f >= b.lo && f <= b.hi;
}

template <class T>
inline ML_Interval<T>
ML_Conjugate(const ML_Interval<T> &a)
{
    return a;
}

template <class T>
inline ML_Interval<T>
ML_Real(const ML_Interval<T> &a)
{
    return a;       // T cannot be complex, it has to have a total ordering
}

template <class T>
inline ML_Interval<T>
operator +(const ML_Interval<T> &a, const ML_Interval<T> &b)
{
	return ML_Interval<T>(a.lo + b.lo , a.hi + b.hi);
}

template <class T>
inline ML_Interval<T>
operator +(const ML_Interval<T> &a, T f)
{
	return ML_Interval<T>(a.lo + f, a.hi + f);
}

template <class T>
inline ML_Interval<T>
operator +(T f, const ML_Interval<T> &a)
{
	return ML_Interval<T>(a.lo + f, a.hi + f);
}

template <class T>
inline ML_Interval<T>
operator -(const ML_Interval<T> &a, const ML_Interval<T> &b)
{
	return ML_Interval<T>(a.lo - b.hi, a.hi - b.lo);
}

template <class T>
inline ML_Interval<T>
operator -(const ML_Interval<T> &a, T f)
{
	return ML_Interval<T>(a.lo - f, a.hi - f);
}

template <class T>
inline ML_Interval<T>
operator -(T f, const ML_Interval<T> &b)
{
	return ML_Interval<T>(f - b.hi, f - b.lo);
}

template <class T>
inline ML_Interval<T>
operator -(const ML_Interval<T> &a)
{
	return ML_Interval<T>(-a.hi, -a.lo);
}

template <class T>
inline ML_Interval<T>
operator +(const ML_Interval<T> &a)
{
	return a;
}

template <class T>
inline ML_Interval<T>
operator *(const ML_Interval<T> &a, const T f)
{
	return (f < 0) ? ML_Interval<T>(a.hi*f, a.lo*f) : ML_Interval<T>(a.lo*f, a.hi*f);
}

template <class T>
inline ML_Interval<T>
operator *(const T f, const ML_Interval<T> &a)
{
	return (f < 0) ? ML_Interval<T>(a.hi*f, a.lo*f) : ML_Interval<T>(a.lo*f, a.hi*f);
}

template <class T>
inline ML_Interval<T>
operator /(const ML_Interval<T> &a, T f)
{
	return (f < 0) ? ML_Interval<T>(a.hi/f, a.lo/f) : ML_Interval<T>(a.lo/f, a.hi/f);
}

template <class T>
inline ML_Interval<T>
operator *(const ML_Interval<T> &a, const ML_Interval<T> &b)
{
    T f0, f1, f2, f3, fTemp;

    //
    //  Faster to just do four multiplications, then optimal min/max selection.
    //
    f0 = a.lo * b.lo;
    f1 = a.lo * b.hi;
    f2 = a.hi * b.lo;
    f3 = a.hi * b.hi;
    if (f0 > f1) {
        fTemp = f0;
        f0 = f1;
        f1 = fTemp;
    }
    if (f2 > f3) {
        fTemp = f2;
        f2 = f3;
        f3 = fTemp;
    }
    if (f1 > f3)
        f3 = f1;
    if (f2 < f0)
        f0 = f2;
    return ML_Interval<T>(f0, f3);
}

template <class T>
inline ML_Interval<T>
operator /(const ML_Interval<T> &a, const ML_Interval<T> &b)
{
	return a * Reciprocol(b);
}

template <class T>
inline ML_Interval<T>
operator /(T a, const ML_Interval<T> &b)
{
	return a * Reciprocol(b);
}

template <class T>
inline ML_Interval<T>
fabs(const ML_Interval<T>& a)
{

	if (a.lo >= 0.0)
		return a;
	else if (a.hi <= 0.0)
		return ML_Interval<T>(-a.hi, -a.lo);
	else
		return ML_Interval<T>(0.0, (-a.lo > a.hi) ? -a.lo : a.hi);
}

template <class T>
inline ML_Interval<T>
pow(const ML_Interval<T> &a, T y)
{
	ML_Interval<T> aa;

	aa = fabs(a);
	return ML_Interval<T>(pow(aa.lo, y), pow(aa.hi, y));
}

template <class T>
inline ML_Interval<T>
sqrt(const ML_Interval<T> &a)
{
	ML_Interval<T> aa;

	aa = fabs(a);
	return ML_Interval<T>(sqrt(aa.lo), sqrt(aa.hi));
}

template <class T>
inline ML_Interval<T>
exp(const ML_Interval<T> &a)
{

	return ML_Interval<T>(exp(a.lo), exp(a.hi));
}

template <class T>
inline ML_Interval<T>
log(const ML_Interval<T> &a)
{
	return ML_Interval<T>(log(a.lo), log(a.hi));
}

template <class T>
inline ML_Interval<T>
cos(const ML_Interval<T> &a)
{
    T fLow, fHigh, fTmp;
    ML_Interval<T> b;
    
    if (a.hi - a.lo >= T(2.0*F_INTERVAL_PI))
        return ML_Interval<T>(-1.0, 1.0);
    b = a - T(2.0*D_INTERVAL_PI * floor(0.5*D_INTERVAL_PIINV * double(a.lo)));
    fLow = cos(a.lo);
    fHigh = cos(a.hi);
    if (fLow > fHigh) {
        fTmp = fLow;
        fLow = fHigh;
        fHigh = fTmp;
    }
    if (subset(T(0.0), b) || subset(T(2.0*D_INTERVAL_PI), b))
        fHigh = 1.0;
    if (subset(T(D_INTERVAL_PI), b) || subset(T(3.0*D_INTERVAL_PI), b))
        fLow = -1.0;
    return ML_Interval<T>(fLow, fHigh);
}

template <class T>
inline ML_Interval<T>
sin(const ML_Interval<T> &a)
{
    T fLow, fHigh, fTmp;
    ML_Interval<T> b;
    
    if (a.hi - a.lo >= T(2.0*D_INTERVAL_PI))
        return ML_Interval<T>(-1.0, 1.0);
    b = a - T(2.0*D_INTERVAL_PI * floor(0.5*D_INTERVAL_PIINV * double(a.lo)));
    fLow = sin(a.lo);
    fHigh = sin(a.hi);
    if (fLow > fHigh) {
        fTmp = fLow;
        fLow = fHigh;
        fHigh = fTmp;
    }
    if (subset(T(1.5*D_INTERVAL_PI), b) || subset(T(3.5*D_INTERVAL_PI), b))
        fLow = -1.0;
    if (subset(T(0.5*D_INTERVAL_PI), b) || subset(T(2.5*D_INTERVAL_PI), b))
        fHigh = 1.0;
    return ML_Interval<T>(fLow, fHigh);
}

template <class T>
inline ML_Interval<T>
tan(const ML_Interval<T> &a)
{

	if (a.hi*T(D_INTERVAL_PIINV) + 0.5 < ceil(a.lo*T(D_INTERVAL_PIINV) + 0.5))
		return ML_Interval<T>(tan(a.lo), tan(a.hi));
	else
		return ML_Interval<T>(T(-F_INTERVAL_MACHINE_MAX), T(F_INTERVAL_MACHINE_MAX));
}

template <class T>
inline ML_Interval<T>
_hypot(const ML_Interval<T> &x, const ML_Interval<T> &y)
{
	T mindx, mindy, maxdx, maxdy;

	if (x.lo > 0.0) {
		mindx = x.lo;
		maxdx = x.hi;
	} else if (x.hi < 0.0) {
		mindx = x.hi;
		maxdx = x.lo;
	} else {
		mindx = 0.0;
		maxdx = (x.hi > -x.lo) ? x.hi : x.lo;
	}
	if (y.lo > 0.0) {
		mindy = y.lo;
		maxdy = y.hi;
	} else if (y.hi < 0.0) {
		mindy = y.hi;
		maxdy = y.lo;
	} else {
		mindy = 0.0;
		maxdy = (y.hi > -y.lo) ? y.hi : y.lo;
	}
	return ML_Interval<T>(_hypot(mindx, mindy), _hypot(maxdx, maxdy));
}

template <class T>
inline ML_Interval<T>
atan2(const ML_Interval<T> &y, const ML_Interval<T> &x)
{

	if (subset(T(0.0), x)) {
		if (subset(T(0.0), y))
			return ML_Interval<T>(T(-D_INTERVAL_PI), T(D_INTERVAL_PI));
		else if (y > T(0.0))
			return ML_Interval<T>(atan2(y.lo,x.hi), atan2(y.lo,x.lo));
		else
			return ML_Interval<T>(atan2(y.hi,x.lo), atan2(y.hi,x.hi));
	} else if (subset(T(0.0), y)) {
		if (x > T(0.0))
			return ML_Interval<T>(atan2(y.lo,x.lo), atan2(y.hi,x.lo));
		else
		    return ML_Interval<T>(T(-D_INTERVAL_PI), T(D_INTERVAL_PI));
	} else if (x > T(0.0)) {
		if (y > T(0.0))
			return ML_Interval<T>(atan2(y.lo,x.hi), atan2(y.hi,x.lo));
		else
			return ML_Interval<T>(atan2(y.lo,x.lo), atan2(y.hi,x.hi));
	} else {
		if (y > T(0.0))
			return ML_Interval<T>(atan2(y.hi,x.hi), atan2(y.lo,x.lo));
		else
			return ML_Interval<T>(atan2(y.hi,x.lo), atan2(y.lo,x.hi));
	}
}

template <class T>
inline ML_Interval<T>
atan(const ML_Interval<T> &x)
{

	return ML_Interval<T>(atan(x.lo), atan(x.hi));
}

template <class T>
inline ML_Interval<T>
asin(const ML_Interval<T> &x)
{

	return ML_Interval<T>(asin(x.lo), asin(x.hi));
}

template <class T>
inline ML_Interval<T>
acos(const ML_Interval<T> &x)
{

	return ML_Interval<T>(acos(x.hi), acos(x.lo));
}

template <class T>
inline ML_Interval<T>
sinh(const ML_Interval<T> &x)
{

	return ML_Interval<T>(sinh(x.lo), sinh(x.hi));
}

template <class T>
inline ML_Interval<T>
tanh(const ML_Interval<T> &x)
{

	return ML_Interval<T>(tanh(x.lo), tanh(x.hi));
}

template <class T>
inline ML_Interval<T>
cosh(const ML_Interval<T> &x)
{

	if (x.lo >= 0.0)
		return ML_Interval<T>(cosh(x.lo), cosh(x.hi));
	else if (x.hi <= 0.0)
		return ML_Interval<T>(cosh(x.hi), cosh(x.lo));
	else
		return ML_Interval<T>(1.0, (-x.lo > x.hi) ? cosh(x.lo) : cosh(x.hi));
}

template <class T>
inline ML_Interval<T>
asinh(const ML_Interval<T> &x)
{

	return ML_Interval<T>(asinh(x.lo), asinh(x.hi));
}

template <class T>
inline ML_Interval<T>
acosh(const ML_Interval<T> &x)
{

	return ML_Interval<T>(acosh(x.lo), acosh(x.hi));        // monotonic on x > 1
}

template <class T>
inline ML_Interval<T>
atanh(const ML_Interval<T> &x)
{

	return ML_Interval<T>(atanh(x.lo), atanh(x.hi));
}

template <class T>
inline ML_Interval<T>
ceil(const ML_Interval<T> &x)
{

	return ML_Interval<T>(ceil(x.lo), ceil(x.hi));
}

template <class T>
inline ML_Interval<T>
floor(const ML_Interval<T> &x)
{

	return ML_Interval<T>(floor(x.lo), floor(x.hi));
}
