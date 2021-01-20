#include "stdafx.h"
#include "MonteCarlo.h"

#pragma intrinsic (sqrt, log, exp, sin, cos, pow)

#define D_PI            3.14159265358979323846264338327950288419716939937510
#define D_E             2.71828182845904523536028747135266249775724709369995

unsigned ML_nMarsagliaX  = 886459;
unsigned ML_nMarsagliaC  = 361290869;
//
//  True random numbers generated from jitter in clock interrupt
//
unsigned
ML_TrueRandomUnsigned()
{
    static LARGE_INTEGER nPerformanceCount;
    static unsigned n;
    BOOL bResult;
    int i;


    for (i = 0; i < 32; i += 8) {
        bResult = QueryPerformanceCounter(&nPerformanceCount);
        Sleep(1L);
        n = (n << 8) | (n >> 24);
        n ^= nPerformanceCount.LowPart;
    }
    return n;
}
//
//  Uniform distribution from 0 to N-1
//
unsigned
ML_RandomLessThanN(unsigned nModulus)
{
    unsigned nRemainder, n;

    nRemainder = 0xFFFFFFFF % nModulus;
    do {
        n = ML_RandomUnsigned();
    } while (n <= nRemainder);
    return n % nModulus;
}
//
//  Various random distributions
//
double
ML_RandomDisk(double x[2])
{
    double r2;

    do {
        x[0] = 2.0*ML_RandomDouble() - 1.0;
        x[1] = 2.0*ML_RandomDouble() - 1.0;
    } while ((r2 = x[0]*x[0] + x[1]*x[1]) > 1.0);
    return r2;
}

double
ML_RandomCircle(double x[2])            // twice as fast as sin/cos of random angle
{
    double d[2], r2;

    r2 = ML_RandomDisk(d);
    x[0] = 2.0*d[0]*d[1]/r2;            // von Neuman's trick to avoid sqrt
    x[1] = (d[0]*d[0] - d[1]*d[1])/r2;
    return r2;                          
}

double
ML_RandomNormal()
{
    double s, x[2];
    static double fNormal;
    static int nSecondVariate = 1;

    nSecondVariate = !nSecondVariate;
    if (nSecondVariate) {
        return fNormal;
    } else {
        s = ML_RandomDisk(x);
        s = sqrt(-2.0*log(s)/s);
        fNormal = x[0]*s;
        return x[1]*s;
    }
}

void
ML_RandomSphere(double x[3])
{
    double s, fCos, fSin, v[2];

    s = ML_RandomCircle(v);
    fCos = 2.0*s - 1.0;
    fSin = sqrt(1.0 - fCos*fCos);
    x[0] = v[0]*fSin;
    x[1] = v[1]*fSin;
    x[2] = fCos;
}

void
ML_RandomHyperSphere(double x[4])
{
    double r, r1, r2;

    ML_RandomCircle(x);
    ML_RandomCircle(x+2);
    r = ML_RandomDouble();
    r1 = sqrt(1.0 - r);
    r2 = sqrt(r);
    x[0] *= r1;
    x[1] *= r1;
    x[2] *= r2;
    x[3] *= r2;
}

void
ML_RandomRotation(double m[3][3])
{
    union {
        double a[4];
        struct { double x, y, z, t; } v;
    };
    double xx, xy, xz, xw, yy, yz, yw, zz, zw;

    ML_RandomHyperSphere(a);
    xx = 2.0 * v.x * v.x;                   // quanterion to matrix transform
    yy = 2.0 * v.y * v.y;
    zz = 2.0 * v.z * v.z;
    xy = 2.0 * v.x * v.y;
    xz = 2.0 * v.x * v.z;
    yz = 2.0 * v.y * v.z;
    xw = 2.0 * v.x * v.t;
    yw = 2.0 * v.y * v.t;
    zw = 2.0 * v.z * v.t;
    m[0][0] = 1.0 - yy - zz;
    m[0][1] = xy - zw;
    m[0][2] = xz + yw;
    m[1][0] = xy + zw;
    m[1][1] = 1.0 - xx - zz;
    m[1][2] = yz - xw;
    m[2][0] = xz - yw;
    m[2][1] = yz + xw;
    m[2][2] = 1.0 - xx - yy;
}

double
ML_RandomExponential(double mu)
{
	return -mu * log(ML_RandomDouble());
}

double
ML_RandomGamma(double a)
{
	double u, x, y, v, v1, v2, s, p, q, t;
	int k;

	k = int(a);
	if (a < 20.0 && a == (double)k) {
		for (u = 1.0; k > 0; --k)
			u *= ML_RandomDouble();
		return - log(u);
	}
	if (a > 1.0) {
		do {
			do {
				v1 = 2.0 * ML_RandomDouble() - 1.0;
				v2 = 2.0 * ML_RandomDouble() - 1.0;
				s = v1*v1 + v2*v2;
			} while (s >= 1.0);
			y = v2/v1;
			x = sqrt(2.0 * a - 1.0) * y + a - 1.0;
			if (x <= 0.0)
				continue;
			t = (1.0 + y*y) * exp((a - 1.0)*log(x/(a - 1.0))
						- x + a - 1.0);
		} while (ML_RandomDouble() > t);
	} else {
		p = D_E/(a + D_E);
		do {
			u = ML_RandomDouble();
			v = ML_RandomDouble();
			if (u < p) {
				x = exp(log(v) / a);
				q = exp(-x);
			} else {
				x = 1 - log(v);
				q = exp(log(x) * (a - 1.0));
			}
		} while (ML_RandomDouble() >= q);
	}
	return x;
}

double
ML_RandomBeta(double a, double b)
{
	double x1, x2;

	x1 = ML_RandomGamma(a);
	x2 = ML_RandomGamma(b);
	return x1 / (x1 + x2);
}

double
ML_RandomChiSquare(double nu)
{

	return 2.0*ML_RandomGamma(0.5 * nu);
}

double
ML_RandomFDistribution(double nu1, double nu2)
{
	double y1, y2;

	y1 = ML_RandomChiSquare(nu1);
	y2 = ML_RandomChiSquare(nu2);
	return y1*nu2 / (y2*nu1);
}

double
ML_RandomTDistribution(double nu)
{
	double y1, y2;

	y1 = ML_RandomNormal();
	y2 = ML_RandomChiSquare(nu);
	return y1/sqrt(y2/nu);
}

int
ML_RandomGeometric(double p)
{
	return (int)ceil(log(ML_RandomDouble()) / log(1.0 - p));
}

int
ML_RandomBinomial(double p, int trials)
{
	int i, n;
	int a, b;
	double x;

	if (trials < 10) {
		n = 0;
		for (i = 0; i < trials; i++)
			n += ML_RandomDouble() < p;
		return n;
	} else {
		a = (int)(1.0 + floor(0.5*(double)trials));
		b = trials - a + 1;
		x = ML_RandomBeta(a, b);
		if (x >= p)
			return ML_RandomBinomial(p/x, a - 1);
		else
			return a + ML_RandomBinomial((p - x)/(1.0 - x), b - 1);
	}
}

int
ML_RandomPoisson(double mu)
{
	int n;
	double m;
	double x, u;

	if (mu < 10) {
		u = exp(-mu);
		x = 1.0;
		n = 0;
		do {
			n++;
			x *= ML_RandomDouble();
		} while (x > u);
		return n - 1;
	}
	m = floor(7.0*mu/8.0);
	x = ML_RandomGamma(m);
	if (x < mu)
		return (int)m + ML_RandomPoisson(mu - x);
	else
		return ML_RandomBinomial(mu/x, (int)(m - 1.0));
}


//
//  A collection of Knuth's and Marsaglia's statistical randomness tests
//
static int
Report(double x2, double x97, const char *szTest)
{
	if (x2 > x97) {
		printf("       FAIL: %8d.%1d chi-square above 97th percentile in %s\n",
			int(x2), int(10.0*x2)%10, szTest);
		return 0;
	} else {
		printf("       pass: %8d.%1d chi-square in %s\n",
			int(x2), int(10.0*x2)%10, szTest);
		return 1;
	}
}
//
//	A.  Frequency Test
//
static double freqhistogram[256];

static void
Frequency(unsigned long c)
{
	freqhistogram[c&0xFF] += 1.0;
}

static int
FrequencyReport(double count)
{
	double x2, n, p;
	int i, nResult;

	n = count;
	p = 1.0/256.0;
	x2 = 0;
	for (i = 0; i < 256; i++)
		x2 += (freqhistogram[i] - n*p)*(freqhistogram[i] - n*p)/(n*p);
	nResult = Report(x2, 298.953724, "frequency test");
	for (i = 0; i < 256; i++)
		freqhistogram[i] = 0.0;	// Report functions clean-up for the next run
	return nResult;
}
//
//	B.  Serial Test
//
static double paircount[65536];

static void
Serial(unsigned long c)
{
	static unsigned toggle=0, word;

	c &= 0xFF;
	toggle = !toggle;
	if (toggle)
		word = c;
	else
		paircount[(word << 8) + c] += 1.0;
}

static int
SerialReport(double count)
{
	double x2, n, p;
	int i;

	n = floor(count/2.0);
	p = 1.0/65536.0;
	x2 = 0;
	for (i = 0; i < 65536; i++)
		x2 += (paircount[i] - n*p)*(paircount[i] - n*p)/(n*p);
	if (n*p < 5.0)
		return 1; // not really enough statistics (return optimistic guess)
	return Report(x2, 66181.302788, "serial test");
	//
	//  Poker Report cleans up Serial state !!!
	//
}
//
//	C.  Gap Test
//
#define GAP_T 16
#define ALPHA 0.25
#define BETA  0.666

static double gaphistogram[GAP_T+1];
static double ngaps = 0.0, grun = 0.0;

static void
Gap(unsigned long c)
{
	double f = double(c) * (1.0/4294967296.0);

	if (ALPHA <= f && f < BETA) {
		if (grun >= GAP_T)
			gaphistogram[GAP_T] += 1.0;
		else
			gaphistogram[(int)grun] += 1.0;
		grun = 0.0;
		ngaps += 1.0;
	} else
		grun += 1.0;
}

static int
GapReport(double count)
{
	double x2, n, p, q;
	int i, nResult;

	n = ngaps;
	p = (BETA - ALPHA);
	q = 1.0;
	x2 = 0.0;
	for (i = 0; i < GAP_T; i++) {
		x2 += (gaphistogram[i] - n*p*q)*(gaphistogram[i] - n*p*q)/(n*p*q);
		q *= (1.0 - p);
	}
	x2 += (gaphistogram[GAP_T] - n*q)*(gaphistogram[GAP_T] - n*q)/(n*q);
	nResult = Report(x2, 28.186043, "gap test");
	for (i = 0; i < GAP_T+1; i++)
		gaphistogram[i] = 0.0;
	ngaps = 0.0;
	grun = 0.0;
	return nResult;
}
//
//	D. Poker Test
//
static int
PokerReport(double count)
{
	int a, b, c, d;
	int i, nResult;
	static double alldiff, onepair, twopair, threekind, fourkind;
	double p, x2, npairs;

	alldiff = onepair = twopair = threekind = fourkind = 0;
	npairs = floor(count/2.0);
	for (i = 0; i < 65536; i++) {
		a = i >> 12 & 0xf;
		b = i >>  8 & 0xf;
		c = i >>  4 & 0xf;
		d = i >>  0 & 0xf;
		if (a == b && a == c && a == d)
			fourkind += paircount[i];
		else if ((a == b && a == c)||(a == b && a == d)||
			 (a == c && a == d)||(b == c && b == d))
			threekind += paircount[i];
		else if ((a == b && c == d)||(a == d && b == c)||
			 (a == c && b == d))
			twopair += paircount[i];
		else if (a == b || a == c || a == d || b == c || b == d || c == d)
			onepair += paircount[i];
		else
			alldiff += paircount[i];
	}
	p = 16.0 * 15.0 * 14.0 * 13.0 / 65536.0;
	x2 = (alldiff - npairs*p)*(alldiff - npairs*p)/(npairs*p);
	p = 16.0 * 15.0 * 14.0 * 6.0 / 65536.0;
	x2 += (onepair - npairs*p)*(onepair - npairs*p)/(npairs*p);
	p = 16.0 * 15.0 * 3.0 / 65536.0;
	x2 += (twopair - npairs*p)*(twopair - npairs*p)/(npairs*p);
	p = 16.0 * 15.0 * 4.0 / 65536.0;
	x2 += (threekind - npairs*p)*(threekind - npairs*p)/(npairs*p);
	p = 16.0 / 65536.0;
	x2 += (fourkind - npairs*p)*(fourkind - npairs*p)/(npairs*p);
	nResult = Report(x2, 10.704621, "poker test");
	for (i = 0; i < 65536; i++)
		paircount[i] = 0.0;		// clean up serial-test tables now
	return nResult;
}
//
//	E.  Coupon Collector Test
//
#define COUPON_T 40
#define D 8

static void CouponInit();

static double couponhistogram[COUPON_T + 1];
static int occurs[D];
static double ncoupons = 0.0;
static double rCoupon = 0.0, qCoupon = 0.0;
static double p[COUPON_T + 1];

static void
Coupon(unsigned long c)
{
	int u, i;

	u = c % D;
	rCoupon += 1.0;
	if (occurs[u])
		return;
	else {
		occurs[u]++;
		qCoupon += 1.0;
		if (qCoupon < D)
			return;
		if (rCoupon >= COUPON_T)
			couponhistogram[COUPON_T] += 1.0;
		else
			couponhistogram[(int)rCoupon] += 1.0;
		ncoupons += 1.0;
		qCoupon = 0.0;
		rCoupon = 0.0;
		for (i = 0; i < D; i++)
			occurs[i] = 0;
	}
}

static int
CouponReport(double count)
{
	double x2, n;
	int i, nResult;

	CouponInit();
	x2 = 0.0;
	n = ncoupons;
	for (i = D; i < COUPON_T; i++) {
		x2 += (couponhistogram[i] - n*p[i])*(couponhistogram[i] - n*p[i])/(n*p[i]);
	}
	x2 += (couponhistogram[COUPON_T] - n*p[i])*(couponhistogram[COUPON_T] - n*p[i])/(n*p[i]);
	nResult = Report(x2, 48.611826, "coupon-collection test");
	for (i = 0; i < COUPON_T+1; i++)
		couponhistogram[i] = 0.0;
	for (i = 0; i < D; i++)
		occurs[i] = 0;
	ncoupons = 0.0;
	rCoupon = qCoupon = 0.0;
	return nResult;
}

static double q[COUPON_T+1][D+1];

static void
CouponProb(int d, int t)
{
	int r;
	double ratio, power;

	if (d == 1) {
		q[0][d] = 0.0;
		for (r = 1; r < t; r++)
			q[r][d] = 1.0;
		return;
	}
	q[0][d] = 0.0;
	power = 1.0;
	ratio = (double)(d-1)/(double)d;
	for (r = 1; r < t; r++) {
		q[r][d] = power*q[r-1][d-1] + q[r-1][d];
		power *= ratio;
	}
}

static void
CouponInit()
{
	int d, r;

	for (d = 1; d <= D; d++)
		CouponProb(d, COUPON_T);
	for (r = D-1; r < COUPON_T; r++)
		q[r][D] = 1.0 - q[r][D];
	for (r = D; r < COUPON_T; r++)
		p[r] = q[r-1][D] - q[r][D];
	p[COUPON_T] = q[COUPON_T-1][D];
}
//
//	F.  Permutation Test
//
static double permhistogram[120+1];
static double nperms;
static int iPermute = 0;

static int perm(unsigned long, unsigned long, unsigned long,
				unsigned long, unsigned long);

static void
Permute(unsigned long c)
{
	static unsigned long cp[5];
	int i;

	cp[iPermute++] = c;
	if (iPermute > 4) {
		i = perm(cp[0], cp[1], cp[2], cp[3], cp[4]);
		permhistogram[i] += 1.0;
		if (i < 120)
			nperms += 1.0;
		iPermute = 0;
	}
}

static int
PermuteReport(double count)
{
	double x2, n, p;
	int i, nResult;

	x2 = 0;
	n = nperms;
	p = 1.0/(5.0 * 4.0 * 3.0 * 2.0 * 1.0);
	for (i = 0; i < 120; i++)
		x2 += (permhistogram[i] - n*p)*(permhistogram[i] - n*p)/(n*p);
	nResult = Report(x2, 149.652965, "permutations test");
	for (i = 0; i < 120+1; i++)
		permhistogram[i] = 0.0;
	nperms = 0.0;
	iPermute = 0;
	return nResult;
}

static int
perm(unsigned long a, unsigned long b, unsigned long c,
	 unsigned long d, unsigned long e)
{
	unsigned long s, i;
	unsigned long *u;
	unsigned long r, f;
	struct { unsigned long a, b, c, d, e; } arg;

	arg.a = a;
	arg.b = b;
	arg.c = c;
	arg.d = d;
	arg.e = e;
	u = (unsigned long *)&arg;
	f = 0;
	s = 0;
	r = 5;
	for (i = 1; i < r; i++)
		if (u[i] > u[s])
			s = i;
	f = r * f + s;
	--r;
	i = u[s];
	u[s] = u[r];
	u[r] = i;
	s = 0;
	for (i = 1; i < r; i++)
		if (u[i] > u[s])
			s = i;
	f = r * f + s;
	--r;
	i = u[s];
	u[s] = u[r];
	u[r] = i;
	s = 0;
	for (i = 1; i < r; i++)
		if (u[i] > u[s])
			s = i;
	f = r * f + s;
	--r;
	i = u[s];
	u[s] = u[r];
	u[r] = i;
	s = 0;
	for (i = 1; i < r; i++)
		if (u[i] > u[s])
			s = i;
	f = r * f + s;
	--r;
	i = u[s];
	u[s] = u[r];
	u[r] = i;
	if (arg.a == arg.b||arg.b == arg.c||arg.c == arg.d||arg.d == arg.e)
		f = 120;
	return f;
}
//
//	G.  Runs Up Test
//
#define RUN_T 6

static double runshistogram[RUN_T + 1];
static double nruns;
static int qRun, rRun, toggleRun, rejectRun;

static void
Runs(unsigned long c)
{
	static  int q, r, toggle, reject;
	static unsigned long u, lastu;

	toggleRun = !toggleRun;
	if (toggleRun) {
		u = c;
		return;
	} else
		u = (u << 8) + c;
	if (qRun++ == 0) {
		lastu = u;
		return;
	}
	rRun++;
	if (u == lastu)
		rejectRun++;
	if (u > lastu) {
		lastu = u;
		return;
	}
	if (rejectRun == 0) {
		if (rRun >= RUN_T)
			runshistogram[RUN_T] += 1.0;
		else
			runshistogram[rRun] += 1.0;
		nruns += 1.0;
	}
	rRun = 0;
	qRun = 0;		//skip next item
	rejectRun = 0;
}
static double
Factorial(unsigned n)
{
	double f;

	f = 1.0;
	while (n >= 2) {
		f *= double(n);
		--n;
	}
	return f;
}

static int
RunsReport(double count)
{
	double x2, n, p;
	int i, nResult;

	n = nruns;
	x2 = 0.0;
	for (i = 1; i < RUN_T; i++) {
		p = 1.0/Factorial(i) - 1.0/Factorial(i + 1);
		x2 += (runshistogram[i] - n*p)*(runshistogram[i] - n*p)/(n*p);
	}
	p = 1.0/Factorial(RUN_T);
	x2 += (runshistogram[RUN_T] - n*p)*(runshistogram[RUN_T] - n*p)/(n*p);
	nResult = Report(x2, 12.370127, "runs-up test");
	for (i = 0; i < RUN_T+1; i++)
		runshistogram[i] = 0.0;
	nruns = 0.0;
	rRun = qRun = toggleRun = rejectRun = 0;
	return nResult;
}
//
//	H.  Maximum of t Test
//
#include <stdio.h>
#define MAX_T 8

static double maxhistogram[256];
static double nmax;
static unsigned long maxMaximum, countMaximum;

static void
Maximum(unsigned long c)
{

	c &= 0xFF;
	if (c > maxMaximum)
		maxMaximum = c;
	countMaximum++;
	if (countMaximum >= MAX_T) {
		maxhistogram[maxMaximum] += 1.0;
		nmax += 1.0;
		maxMaximum = 0;
		countMaximum = 0;
	}
}

static int
MaximumReport(double count)
{
	double x2, n, p;
	int i, nResult;

	n = nmax;
	x2 = 0;
	for (i = 200; i < 256; i++) {
		p = pow((double)(i + 1)/256.0, (double)MAX_T)
		  - pow((double)i/256.0, (double)MAX_T);
		x2 += (maxhistogram[i] - n*p)*(maxhistogram[i] - n*p)/(n*p);
	}
	nResult = Report(x2, 76.338817, "max-of-8 test");
	for (i = 0; i < 256; i++)
		maxhistogram[i] = 0.0;
	nmax = 0.0;
	maxMaximum = countMaximum = 0;
	return nResult;
}
//
//	I. mtuple test
//
static double triples[512];
static double pairs[64];
static double ntuples;
static int aTuple, bTuple, cTuple, nbytesTuple;

static void
MTuple(unsigned long c)
{
	int i;

	aTuple = bTuple;
	bTuple = cTuple;
	cTuple = c & 7;
	if (nbytesTuple++ >= 3) {
		nbytesTuple = 4;
		i = (aTuple) | (bTuple << 3);
		pairs[i] += 1.0;
		i = i | (cTuple << 6);
		triples[i] += 1.0;
		ntuples += 1.0;
	}
}

static int
MTupleReport(double count)
{
	double q3, q2, mean, n;
	int i, nResult;

	n = ntuples;
	q3 = q2 = 0;
	mean = n/512.0;
	for (i = 0; i < 512; i++)
		q3 += (triples[i] - mean)*(triples[i] - mean)/mean;
	mean = n/64.0;
	for (i = 0; i < 64; i++)
		q2 += (pairs[i] - mean)*(pairs[i] - mean)/mean;
	nResult = Report(q3 - q2, 505.513052, "lapped m-tuple test");
	for (i = 0; i < 512; i++)
		triples[i] = 0.0;
	for (i = 0; i < 64; i++)
		pairs[i] = 0.0;
	ntuples = 0.0;
	aTuple = bTuple = cTuple = nbytesTuple = 0;
	return nResult;
}
//
//	A2.  Frequency Test of ML_RandomLessThanN
//
static double freq360histogram[360];

static void
Frequency360(unsigned long c)
{
	freq360histogram[c] += 1.0;
}

static int
Frequency360Report(double count)
{
	double x2, n, p;
	int i, nResult;

	n = count;
	p = 1.0/double(360);
	x2 = 0;
	for (i = 0; i < 360; i++)
		x2 += (freq360histogram[i] - n*p)*(freq360histogram[i] - n*p)/(n*p);
	nResult = Report(x2, 411.055634, "RandomeLessThan(360) frequency test");
	for (i = 0; i < 360; i++)
		freq360histogram[i] = 0;
	return nResult;
}
//
//  Test nonuniform variates with KS
//
static double gs_rgfX[100000];

static void
QuickSort(double rgtItem[], int nItems)
{
    int i, j;
    double tTemp, tPivot;

    nItems--;   // iLast
    while (nItems >= 30) {
        i = nItems/2;
        //
        //  Sort first, middle and last elements.  This provides median-of-3
        //  partitioning, limits partitioning to only N - 3 remaining items,
        //  and creates sentinals to simplify the inner loop.
        //
        if (rgtItem[i] < rgtItem[0]) {      // 2.48 compares on average
            tTemp = rgtItem[0];
            rgtItem[0] = rgtItem[i];
            rgtItem[i] = tTemp;
        }
        if (rgtItem[nItems] < rgtItem[i]) {
            tTemp = rgtItem[nItems];
            rgtItem[nItems] = rgtItem[i];
            if (tTemp < rgtItem[0]) {
                rgtItem[i] = rgtItem[0];
                rgtItem[0] = tTemp;
            } else
                rgtItem[i] = tTemp;
        }
        j = nItems - 1;
        tPivot = rgtItem[i];
        rgtItem[i] = rgtItem[j];
        rgtItem[j] = tPivot;
        i = 0;
        //
        //  Partition, using Sedgewick's "j < i" suggestion.  Oddly, it is
        //  faster to loop on i before looping on j (on the Pentium 4).
        //
        for(;;) {
            while(rgtItem[++i] < tPivot)
                ;
            while(tPivot < rgtItem[--j])
                ;
            if (j < i)
                break;
            tTemp = rgtItem[i];
            rgtItem[i] = rgtItem[j];
            rgtItem[j] = tTemp;
        }
        tTemp = rgtItem[nItems - 1];
        rgtItem[nItems - 1] = rgtItem[i];
        rgtItem[i] = tTemp;
        //
        //  Recursing on smaller partition yields O(log N) stack growth.
        //
        if (j < nItems - i - 1) {
            QuickSort(rgtItem ,j + 1);
            rgtItem += i + 1;
            nItems -= i + 1;
        } else {
            QuickSort(rgtItem + i + 1, nItems - i);
            nItems = j;
        }
    }
    //
    //  Small partitions are insertion sorted.  Distribution is wedge
    //  shaped, with only about 3.8 comparisons done on average, and
    //  benefit gained from structuring the loop for quick out.
    //
    for (i = 1; i <= nItems; i++) {
        j = i;
        tTemp = rgtItem[j];
        if (tTemp <rgtItem[j - 1]) {
            do {
                rgtItem[j] = rgtItem[j - 1];
                j = j - 1;
            } while (j > 0 && tTemp < rgtItem[j - 1]);    //REVIEW: forloop  faster with MSVC
            rgtItem[j] = tTemp;
        }
    }
}

static int
KolmogorovSmirnovTest(double (*Random)(), double (*CDF)(double), char *szName, char *szLabel)
{
	int i, j, N;
	double fKPlus, fKMinus, fK, F, yp, K, fRootN;

	printf("%s Kolmogorov-Smirnov Test of %s\n", szLabel, szName);
	fKPlus = fKMinus = 0.0;
	N = sizeof(gs_rgfX)/sizeof(double);
	fRootN = sqrtf(float(N));
	for (i = 0; i < N; i++)
		gs_rgfX[i] = Random();
	QuickSort(gs_rgfX, N);
	for (j = 1; j <= N; j++) {
		F = CDF(gs_rgfX[j-1]);
		fK = float(j)/float(N) - F;
		if (fK > fKPlus)
			fKPlus = fK;
		fK = F - float(j - 1)/float(N);
		if (fK > fKMinus)
			fKMinus = fK;
	}
	fKPlus *= fRootN;
	fKMinus *= fRootN;
	yp = sqrt(0.5*log(1.0/(1.0 - 0.99)));   // 95th percentile
	K = yp - 1.0/(6.0*fRootN);
	if (fKPlus < K && fKMinus < K) {
		printf("       pass: %d %d %d\n", int(1000.0*K), int(1000.0*fKPlus), int(1000.0*fKMinus));
		return 1;
	} else {
		printf("       FAIL: %d %d %d\n", int(1000.0*K), int(1000.0*fKPlus), int(1000.0*fKMinus));
		return 0;
	}
}

void
TestRandomCircle()
{
    static int rgn[36];
    int i;
    double x[2], theta;

    for (i = 0; i < 10000000; i++) {
        ML_RandomCircle(x);
        theta = atan2(x[1], x[0])*18.0/D_PI + 18.0;
        if (i % 1000000 == 0)
            printf("RC: %f\n", theta);
        rgn[int(theta)]++;
    }
    for (i = 0; i < 18; i++)
        printf("%2d %7d       %2d %7d\n", i, rgn[i], i+18, rgn[i+18]);
}

void
TestRandomSphere()
{
    double freqOctants[8];
    double x[3], y[3], m[3][3], p, x2;
    int i, n;

    printf("Test random sphere:\n");
    //
    //  Test random points on sphere by historgram of randomly oriented octants
    //
    for (i = 0; i < 8; i++)
        freqOctants[i] = 0.0;
    ML_RandomRotation(m);
    for (n = 0; n < 1000000000; n++) {
        ML_RandomSphere(x);
        y[0] = x[0]*m[0][0] + x[1]*m[0][1] + x[2]*m[0][2];
        y[1] = x[0]*m[1][0] + x[1]*m[1][1] + x[2]*m[1][2];
        y[2] = x[0]*m[2][0] + x[1]*m[2][1] + x[2]*m[2][2];
        i = (y[0] < 0.0) | ((y[1] < 0.0) << 1) | ((y[2] < 0.0) << 2);
        freqOctants[i] += 1.0;
    }
	p = 1.0/8.0;
	x2 = 0;
	for (i = 0; i < 8; i++)
		x2 += (freqOctants[i] - n*p)*(freqOctants[i] - n*p)/(n*p);      // Chi-Square
	Report(x2, 16.0128, "sphere points");
    //
    //  Test random rotations
    //
    for (i = 0; i < 8; i++)
        freqOctants[i] = 0.0;
    ML_RandomSphere(x);
    for (n = 0; n < 1000000000; n++) {
        ML_RandomRotation(m);
        y[0] = x[0]*m[0][0] + x[1]*m[0][1] + x[2]*m[0][2];
        y[1] = x[0]*m[1][0] + x[1]*m[1][1] + x[2]*m[1][2];
        y[2] = x[0]*m[2][0] + x[1]*m[2][1] + x[2]*m[2][2];
        i = (y[0] < 0.0) | ((y[1] < 0.0) << 1) | ((y[2] < 0.0) << 2);
        freqOctants[i] += 1.0;
    }
	p = 1.0/8.0;
	x2 = 0;
	for (i = 0; i < 8; i++)
		x2 += (freqOctants[i] - n*p)*(freqOctants[i] - n*p)/(n*p);
	Report(x2, 16.0128, "rotations");
}
//
//  Returns number of tests failed.
//
int
RandomnessTests(unsigned (*RandomSource)(), double fTries, char *szName, char *szLabel)
{
	unsigned nR1, nR2, nLastRand, nRand;
	double fc;
	int nFail;

	printf("%s Report on %u calls to %s\n",szLabel, unsigned(fTries), szName);
	fflush(stdout);
	nR1 = RandomSource();
	nR2 = RandomSource();
	nLastRand = RandomSource();
	for (fc = 0.0; fc < fTries; fc += 1.0) {
		nRand = RandomSource();
		if (nLastRand == nR1 && nRand == nR2) {
			printf("Failed: too short of a period, %d\n", int(fc));
		}
		nLastRand = nRand;
		Frequency(nRand);
		Serial(nRand);
		Gap(nRand);
		Coupon(nRand);
		Permute(nRand);
		Runs(nRand);
		Maximum(nRand);
		MTuple(nRand);
		Frequency360(nRand % 360);
	}
	//
	//  Note, poker report uses serial-test data and then cleans up
	//
	nFail = 0;
	nFail += !FrequencyReport(fc);
	nFail += !SerialReport(fc);
	nFail += !GapReport(fc);
	nFail += !PokerReport(fc);
	nFail += !CouponReport(fc);
	nFail += !PermuteReport(fc);
	nFail += !RunsReport(fc);
	nFail += !MaximumReport(fc);
	nFail += !MTupleReport(fc);
	nFail += !Frequency360Report(fc);
	return nFail;
}
