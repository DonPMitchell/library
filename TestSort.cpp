#include"stdafx.h"
#include "Sorting.h"
#include "Montecarlo.h"
#include "Utilities.h"
#include "IsLess.h"
#include "Heap.h"
#include <algorithm>

#define ML_TESTING      1
#define D_MACHINE_EPS   2.22044605e-016
#define D_MACHINE_MAX   1.79e+308
#define NITEMS          400
#define NBIG            10000

static unsigned gs_nCompares = 0;

template <class T>
static void
VerifySort(T rgtItem[], int nItems)
{
    int i;
    unsigned nSave;

    nSave = gs_nCompares;
    for (i = 0; i < nItems - 1; i++)
        ML_Assert(!ML_IsLess(rgtItem[i + 1], rgtItem[i]));
    gs_nCompares = nSave;
}
//
//  Foo counts comparisons while being sorted and has a fairly expensive
//  compare function (not inlined).
//
class Foo {
public:
    Foo() {}
    Foo(int i) : f(double(i)) {}
    Foo(double d) : f(d) {}

    double f;
};

inline int operator  <(const Foo &a, const Foo &b) { gs_nCompares++; return (a.f)  < (b.f); }
int operator  >(Foo &a, Foo &b) { gs_nCompares++; return (a.f)  > (b.f); }
inline int operator <=(Foo &a, Foo &b) { gs_nCompares++; return (a.f) <= (b.f); }
int operator >=(Foo &a, Foo &b) { gs_nCompares++; return (a.f) >= (b.f); }
int operator ==(Foo &a, Foo &b) { return (a.f) == (b.f); }
int operator !=(Foo &a, Foo &b) { return (a.f) != (b.f); }
int operator ==(Foo &a, int i)  { return (a.f) == double(i); }
inline int ML_IsLess(Foo &a, Foo &b) { gs_nCompares++; return (a.f)  < (b.f); }
inline int ML_IsLessOrEqual(Foo &a, Foo &b) { gs_nCompares++; return (a.f)  <= (b.f); }

static int _cdecl
CompareFoo(const void *a, const void *b)
{ 
    gs_nCompares++; 
    return (((Foo *)a)->f > ((Foo *)b)->f) - (((Foo *)a)->f < ((Foo *)b)->f); 
}
//
//  Do a more intensive test of sorting, identity, stability, etc.
//
#define NMAX 20000

struct Bar {
    unsigned    n;
    unsigned    nID;
};

int inline operator  <(const Bar &a, const Bar &b) { gs_nCompares++; return a.n < b.n; }
int inline operator ==(const Bar &a, const Bar &b) { return a.n == b.n; }
int inline operator <=(const Bar &a, const Bar &b) { return a.n <= b.n; }
int inline ML_IsLess(const Bar &a, const Bar &b) { gs_nCompares++; return a.n < b.n; }
int inline ML_IsLessOrEqual(const Bar &a, const Bar &b) { gs_nCompares++; return a.n <= b.n; }

static unsigned char rgBitTable[NMAX/8 + 1];
static Bar gs_rgb[NMAX];
//
//  Make sure an array is sorted and that all the elements are present.
//  Don't change gs_nCompares.
//
static void
AssertIsSorted(Bar rgb[], int nItems)
{
    int i;

    for (i = 0; i <= nItems/8; i++)
        rgBitTable[i] = 0;
    for (i = 0; i < nItems - 1; i++)
        ML_Assert(rgb[i] <= rgb[i + 1]);
    for (i = 0; i < nItems; i++)
        rgBitTable[rgb[i].nID >> 3] |= 1 << (rgb[i].nID & 7);
    for (i = 0; i < nItems/8; i++)
        ML_Assert(rgBitTable[i] == 0xFF);
    for (i = nItems & 0xFFFFFFF8; i < nItems; i++)
        ML_Assert(rgBitTable[i >> 3] & (1 << (i & 7)));
}

static int
IsStable(Bar rgb[], int nItems)
{
    int i;

    for (i = 0; i < nItems - 1; i++)
        if (rgb[i] == rgb[i + 1] && rgb[i].nID > rgb[i + 1].nID)
            return 0;
    return 1;
}

static int gs_nPrepared;

static void
PrepareData(Bar rgb[], int nItems)
{
    int i, n;

    switch (ML_RandomUnsigned() & 31) {

    case 0:
        for (i = 0; i < nItems; i++)
            rgb[i].n = 0;                   // constant
        gs_nPrepared = 0;
        break;
    case 1:
        for (i = 0; i < nItems; i++)
            rgb[i].n = i + 17;              // sorted
        gs_nPrepared = 1;
        break;
    case 2:
        for (i = 0; i < nItems; i++)
            rgb[i].n = nItems + 17 - i;     // reverse sorted
        gs_nPrepared = 2;
        break;
    default:
        if (ML_RandomUnsigned() & 4096) {
            n = 3 + int(sqrt(double(nItems)));
            for (i = 0; i < nItems; i++)
                rgb[i].n = ML_RandomLessThanN(n);   // many duplicates
            gs_nPrepared = 3;
        } else {
            for (i = 0; i < nItems; i++)
                rgb[i].n = i;
            ML_RandomShuffle(rgb, nItems);          // a permutation
            gs_nPrepared = 4;
        }
    }
    for (i = 0; i < nItems; i++)
        rgb[i].nID = i;
}

static void
PrepareCombination(Bar rgb[], int nItems, int iComb)
{
    int i;

    for (i = 0; i < nItems; i++) {
        rgb[i].n = iComb % nItems;
        iComb = iComb / nItems;
        rgb[i].nID = i;
    }
}

static void
PreparePermutation(Bar rgb[], int nItems, int iPerm)
{
    int i, j;
    Bar bTemp;

    for (i = 0; i < nItems; i++)
        rgb[i].n = i;
    for (j = nItems - 1; j > 0; --j) {
        i = iPerm % (j + 1);
        iPerm /= j + 1;
        bTemp = rgb[i];
        rgb[i] = rgb[j];
        rgb[j] = bTemp;
    }
    for (i = 0; i < nItems; i++)
        rgb[i].nID = i;
}

extern char *ML_rgszTestNouns[10636];

//
//  Test ML_Quicksort
//
static void
TestQuickSort()
{
    int i, N, logN, nTrials;
    unsigned nX, nC, nSum, nMin, nMax;
    ML_TimingInfo ti;
    double rgfBig[NBIG];
    Foo rgFoo[NBIG];
    char  *rgsz[NBIG];

    printf("  Intensive Testing of ML_QuickSort\n");
    printf("      a. sort various arrangements of monitored Bar type\n");
    nX = ML_nMarsagliaX;
    nC = ML_nMarsagliaC;
    ML_InitializeRandom();
    for (N = 5, logN = 1; N <= NMAX; N *= 10, logN += 1) {
        nTrials = 1 + 20000000/(N*logN);
        nSum = 0;
        nMax = 0;
        nMin = 4000000000;
        ML_StartTiming(ti);
        for (i = 0; i < nTrials; i++) {
            gs_nCompares = 0;
            PrepareData(gs_rgb, N);
            ML_QuickSort(gs_rgb, N);
            nSum += gs_nCompares;
            if (gs_nCompares > nMax)
                nMax = gs_nCompares;
            if (gs_nCompares < nMin)
                nMin = gs_nCompares;
            AssertIsSorted(gs_rgb, N);
            //ML_TestAssert(IsStable(gs_rgb, N)); it is not stable, of course
        }
        ML_StopTiming(ti);
        printf("       %8d N %8d ave %8d min %8d max ", N,
                nSum/nTrials, nMin, nMax);
        ML_ReportEventsPerSecond(ti, nTrials, "sorts");
    }
    //
    //  Test and time some other data types
    //
    for (i = 0; i < NBIG; i++)
        rgfBig[i] = i;
    ML_InitializeRandom();
    ML_StartTiming(ti);
    for (i = 0; i < 500; i++) {
        ML_RandomShuffle(rgfBig, NBIG);
        ML_QuickSort(rgfBig, NBIG);
    }
    ML_StopTiming(ti);
    printf("      b. sort 10000 doubles: ");
    ML_ReportEventsPerSecond(ti, i, "sorts");
    VerifySort(rgfBig, NBIG);

    for (i = 0; i < NBIG; i++)
        rgFoo[i].f = i;
    ML_InitializeRandom();
    gs_nCompares = 0;
    ML_StartTiming(ti);
    for (i = 0; i < 500; i++) {
        ML_RandomShuffle(rgFoo, NBIG);
        ML_QuickSort(rgFoo, NBIG);
    }
    ML_StopTiming(ti);
    printf("      c. sort 10000 Foo:     ");
    ML_ReportEventsPerSecond(ti, i, "sorts");
    VerifySort(rgFoo, NBIG);

    ML_Assert(NBIG <= 10636);   // size of string table
    for (i = 0; i < NBIG; i++)
        rgsz[i] = ML_rgszTestNouns[i];
    ML_InitializeRandom();
    ML_StartTiming(ti);
    for (i = 0; i < 500; i++) {
        ML_RandomShuffle(rgsz, NBIG);
        ML_QuickSort(rgsz, NBIG);
    }
    ML_StopTiming(ti);
    printf("      d. sort 10000 strings: ");
    ML_ReportEventsPerSecond(ti, i, "sorts");
    VerifySort(rgsz, NBIG);

    ML_InitializeRandom(nX, nC);
}
//
//  Insertion sort template
//
static void
TestInsertionSort()
{
    Bar rgb[NBIG];
    int i, j, n, N;
    ML_TimingInfo ti;
    unsigned nSum, nMin, nMax;

    printf("  InsertionSort\n");
    nSum = 0;
    nMax = 0;
    nMin = 10000;
    for (i = 0; i < 2*2; i++) {
        gs_nCompares = 0;
        PrepareCombination(rgb, 2, i);
        ML_InsertionSort(rgb, 2);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        if (gs_nCompares < nMin)
            nMin = gs_nCompares;
        AssertIsSorted(rgb, 2);
        ML_Assert(IsStable(rgb, 2));
    }
    n = (10000*nSum)/i;
    printf("      a. Two:       %d.%0.4d average compares %d min %d max\n",
        n/10000, n%10000, nMin, nMax);
    nSum = 0;
    nMax = 0;
    nMin = 10000;
    for (i = 0; i < 3*3*3; i++) {
        gs_nCompares = 0;
        PrepareCombination(rgb, 3, i);
        ML_InsertionSort(rgb, 3);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        if (gs_nCompares < nMin)
            nMin = gs_nCompares;
        AssertIsSorted(rgb, 3);
        ML_Assert(IsStable(rgb, 3));
    }
    n = (10000*nSum)/i;
    printf("      b. Three:     %d.%0.4d average compares %d min %d max\n",
        n/10000, n%10000, nMin, nMax);
    nSum = 0;
    nMax = 0;
    nMin = 10000;
    for (i = 0; i < 4*4*4*4; i++) {
        gs_nCompares = 0;
        PrepareCombination(rgb, 4, i);
        ML_InsertionSort(rgb, 4);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        if (gs_nCompares < nMin)
            nMin = gs_nCompares;
        AssertIsSorted(rgb, 4);
        ML_Assert(IsStable(rgb, 4));
    }
    n = (10000*nSum)/i;
    printf("      c. Four:      %d.%0.4d average compares %d min %d max\n",
        n/10000, n%10000, nMin, nMax);
    N = 10;
    for (i = 0; i < 10000; i++) {
        PrepareData(rgb, N);
        ML_InsertionSort(rgb, N);
        AssertIsSorted(rgb, N);
        ML_Assert(IsStable(rgb, N));
    }
    N = 10000;
    ML_StartTiming(ti);
    for (i = 0; i < 5; i++) {
        for (j = 0; j < N; j++)
            rgb[j].n = ML_RandomUnsigned();
        ML_InsertionSort(rgb, N);
    }
    ML_StopTiming(ti);
    printf("      d. 10000: ");
    ML_ReportEventsPerSecond(ti, i, "sorts");
}
//
//  Test the little sorts
//
static void
TestSortTwo()
{
    int i, n;
    Bar rgb[2];

    gs_nCompares = 0;
    for (i = 0; i < 2*2; i++) {
        PrepareCombination(rgb, 2, i);
        ML_SortTwo(rgb);
        AssertIsSorted(rgb, 2);
        ML_Assert(IsStable(rgb, 2));
    }
    n = (10000*gs_nCompares)/i;
    printf("      a. SortTwo:   %d.%0.4d average compares\n", n/10000, n%10000);
}

static void
TestSortThree()
{
    int i, n;
    Bar rgb[3];
    unsigned nSum, nMin, nMax;

    nSum = 0;
    nMax = 0;
    nMin = 10000;
    for (i = 0; i < 3*3*3; i++) {
        gs_nCompares = 0;
        PrepareCombination(rgb, 3, i);
        ML_SortThree(rgb);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        if (gs_nCompares < nMin)
            nMin = gs_nCompares;
        AssertIsSorted(rgb, 3);
        ML_Assert(IsStable(rgb, 3));
    }
    n = (10000*nSum)/i;
    printf("      b. SortThree: %d.%0.4d average compares %d min %d max\n",
                    n/10000, n%10000, nMin, nMax);
}
//
//  This has the same average compare count as ML_SortFour, but a min of 3.
//  It is more expensive overall.
//
template <class T>
inline static void
TEST_SortFour(T rgtItem[])
{
    T tTemp;
    int i;

    if (rgtItem[1] < rgtItem[0]) {  // sort first three
        tTemp = rgtItem[0];
        rgtItem[0] = rgtItem[1];
        rgtItem[1] = tTemp;
    }
    if (rgtItem[2] < rgtItem[1]) {
        tTemp = rgtItem[2];
        rgtItem[2] = rgtItem[1];
        if (tTemp < rgtItem[0]) {
            rgtItem[1] = rgtItem[0];
            rgtItem[0] = tTemp;
        } else
            rgtItem[1] = tTemp;
    } else {
        //
        //  We can get a best-case of 3 compares by doing insertion sort here,
        //  but this ends up being less time efficient on the average.
        //
        tTemp = rgtItem[3];
        i = 3;
        if (tTemp < rgtItem[2]) {
            rgtItem[3] = rgtItem[2];
            i = 2;
            if (tTemp < rgtItem[1]) {
                rgtItem[2] = rgtItem[1];
                i = 1;
                if (tTemp < rgtItem[0]) {
                    rgtItem[1] = rgtItem[0];
                    i = 0;
                }
            }
        }
        rgtItem[i] = tTemp;
        return;
    }
    //
    //  Could "else insertion sort" here, best case 3, worst case 5.
    //  But it turns out to be an average-case lose.  Instead, do
    //  binary insertion sort.
    //
    if (rgtItem[3] < rgtItem[1]) {
        tTemp = rgtItem[3];
        rgtItem[3] = rgtItem[2];
        rgtItem[2] = rgtItem[1];
        if (tTemp < rgtItem[0]) {
            rgtItem[1] = rgtItem[0];
            rgtItem[0] = tTemp;
        } else {
            rgtItem[1] = tTemp;
        }
    } else {
        if (rgtItem[3] < rgtItem[2]) {
            tTemp = rgtItem[3];
            rgtItem[3] = rgtItem[2];
            rgtItem[2] = tTemp;
        }
    }
}

static int
TestSortFour()
{
    int i, j, n;
    Bar rgb[4];
    int rgn[4];
    ML_TimingInfo ti;
    unsigned nSum, nMin, nMax;

    nSum = 0;
    nMax = 0;
    nMin = 10000;
    for (i = 0; i < 4*4*4*4; i++) {
        gs_nCompares = 0;
        PrepareCombination(rgb, 4, i);
        ML_SortFour(rgb);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        if (gs_nCompares < nMin)
            nMin = gs_nCompares;
        AssertIsSorted(rgb, 4);
        ML_Assert(IsStable(rgb, 4));
    }
    n = (10000*nSum)/i;
    printf("      b. SortFour:  %d.%0.4d average compares %d min %d max\n",
                    n/10000, n%10000, nMin, nMax);
    n = 0;
    ML_StartTiming(ti);
    for (j = 0; j < 1000000; j++) {
        for (i = 0; i < 4; i++)
            n += ML_RandomUnsigned();
    }
    ML_NullLoopTiming(ti);
    for (j = 0; j < 1000000; j++) {
        for (i = 0; i < 4; i++)
            rgn[i] = ML_RandomUnsigned();
        ML_InsertionSort(rgn, 4);
    }

    ML_StopTiming(ti);
    printf("        1. ML_InsertionSort: ");
    ML_ReportEventsPerSecond(ti, j, "sorts");

    ML_StartTiming(ti);
    for (j = 0; j < 1000000; j++) {
        for (i = 0; i < 4; i++)
            n += ML_RandomUnsigned();
    }
    ML_NullLoopTiming(ti);
    for (j = 0; j < 1000000; j++) {
        for (i = 0; i < 4; i++)
            rgn[i] = ML_RandomUnsigned();
        ML_SortFour(rgn);
    }
    ML_StopTiming(ti);
    printf("        2. ML_SortFour:      ");
    ML_ReportEventsPerSecond(ti, j, "sorts");
    return n;
}

static void
TestBinarySearch()
{
    int j, k, n;
    Foo rgf[16], x;

    printf("  BinarySearch\n");
    x.f = 0.0;
    for (n = 0; n < 10000; n++) {
        for (k = 0; k < 16; k++) {
            for (j = 0; j < k; j++)
                rgf[j].f = ML_RandomDouble();
            ML_QuickSort(rgf, k);
            x.f = 0.0;
            gs_nCompares = 0;
            ML_Assert(ML_BinarySearch(rgf, k, x) == 0);
            ML_Assert((1 << (gs_nCompares - 1)) <= k);
            for (j = 0; j < k; j++) {
                x.f = rgf[j].f + 2.0*D_MACHINE_EPS;
                gs_nCompares = 0;
                ML_Assert(ML_BinarySearch(rgf, k, x) == j + 1);
                ML_Assert((1 << (gs_nCompares - 1)) <= k);
            }
        }
    }
}
//
//  Evaluate a number of other sorting routines for curiosity's sake
//

//
//  Vornberger's small quicksort routine
//
static void QuicksortV(Foo a[], int unten, int oben)
{   
    Foo tmp ;
    int i = unten;
    int j = oben;
    Foo x = a[(unten+oben) / 2];               // Pivotelement, willkuerlich
  
    do {
        while (a[i] < x) i++;                     // x fungiert als Bremse
        while (a[j] > x) j--;                     // x fungiert als Bremse
        if ( i<=j )  {
            tmp  = a[i];                          // Hilfsspeicher
            a[i] = a[j];                          // a[i] und 
            a[j] = tmp;                           // a[j] werden getauscht
            i++;                  
            j--;  
        }                        
    } while (i <= j); 
    //
    // alle Elemente der linken Haelfte sind kleiner
    // als alle Elemente der rechten Haelfte 
    //
    if (unten < j)  QuicksortV(a, unten, j);      // sortiere linke Haelfte
    if (i < oben )  QuicksortV(a, i, oben );      // sortiere rechte Haelfte
}
//
//  Sedgewick's sorting programs
//
static inline void exch(Foo &a, Foo &b)
{
    Foo t;

    t = a;
    a = b;
    b = t;
}

static inline void compexch(Foo &a, Foo &b)
{
    Foo t;

    if (a > b) {
        t = a;
        a = b;
        b = t;
    }
}

static void selection(Foo a[], int l, int r)
  { for (int i = l; i < r; i++)
      { int min = i;
        for (int j = i+1; j <= r; j++) 
            if (a[j] < a[min]) min = j;
        exch(a[i], a[min]);
      } 
  }

static void insertion(Foo a[], int l, int r)
  { int i;
    for (i = r; i > l; i--) compexch(a[i-1], a[i]);
    for (i = l+2; i <= r; i++)
      { int j = i; Foo v = a[i]; 
        while (v < a[j-1])
          { a[j] = a[j-1]; j--; }
        a[j] = v; 
      } 
  }

static void bubble(Foo a[], int l, int r)
  { for (int i = l; i < r; i++)
      for (int j = r; j > i; j--)
        compexch(a[j-1], a[j]);
  }

static void shellsort(Foo a[], int l, int r)
  { int h, i, j;
    for (h = 1; h <= (r-l)/9; h = 3*h+1) ;
    for ( ; h > 0; h /= 3)
      for (i = l+h; i <= r; i++)
        { j = i; Foo v = a[i]; 
          while (j >= l+h && v < a[j-h])
            { a[j] = a[j-h]; j -= h; }
          a[j] = v; 
        } 
  }

void shellsortP(Foo a[], int l, int r)
  { int i, j, k; 
    static int incs[16] = { 1391376, 463792, 198768, 
                            86961, 33936, 13776, 4592, 
                            1968, 861, 336, 112, 48, 
                            21, 7, 3, 1 };
    for ( k = 0; k < 16; k++)
      { int h = incs[k];
        for (i = l+h; i <= r; i++)
          { Foo v = a[i]; 
            j = i;
            while (j >= h && (v < a[j-h]))
              { a[j] = a[j-h]; j -= h; }
            a[j] = v; 
          } 
      }
  }


static int partition(Foo a[], int l, int r)
  { int i = l-1, j = r; Foo v = a[r];
    for (;;)
      { 
        while (a[++i] < v) ;
        while (v < a[--j]) if (j == l) break;
        if (i >= j) break;
        exch(a[i], a[j]);
      }
    exch(a[i], a[r]);
    return i;
  }
static void quicksort1(Foo a[], int l, int r)
  {
    if (r <= l) return;
    int i = partition(a, l, r);
    quicksort1(a, l, i-1);
    quicksort1(a, i+1, r);
  }
//
//  Gosling and Ahren's quicksort
//
static inline void 
swap(Foo a[], int i, int j)
{
        Foo T;
        T = a[i]; 
        a[i] = a[j];
        a[j] = T;
}

 static inline void InsertionSort(Foo a[], int lo0, int hi0)
    {
            int i;
            int j;
            Foo v;
    
            for (i=lo0+1;i<=hi0;i++)
            {
                    v = a[i];
                    j=i;
                    while ((j>lo0) && (a[j-1]>v))
                    {
                            a[j] = a[j-1];
                            j--;
                    }
                    a[j] = v;
            }
    }

void QuickSort6(Foo a[], int l, int r)
   {
        int M = 4;
        int i;
        int j;
        Foo v;

        if ((r-l)>=M)
        {
                i = (r+l)/2;
                if (a[l]>a[i]) swap(a,l,i);     // Tri-Median Methode!
                if (a[l]>a[r]) swap(a,l,r);
                if (a[i]>a[r]) swap(a,i,r);

                j = r-1;
                swap(a,i,j);
                i = l;
                v = a[j];
                for(;;)
                {
                        while(a[++i]<v);
                        while(a[--j]>v);
                        if (j<i) break;
                        swap (a,i,j);
                 }
                swap(a,i,r-1);
                QuickSort6(a,l,j);
                QuickSort6(a,i+1,r);
        } else
            InsertionSort(a, l, r);
    }
    
    static void AhrensQuicksort2(Foo a[], int L, int R)
    {
            QuickSort6(a, L, R);
            //InsertionSort(a,L,R);
    }
/*
int partition(arr, i, j){
  int pivot = arr[(i+j)/2]; // use some way to pick a pivot
  arr[(i+j)/2]= arr[j];
  arr[j] = pivot;           // swap pivot with the rightmost element
  
  // now do the Dutch flag on the rest of the array
  red = i;  // set at the left border of the range
  blue = j; //set at the right border where the pivot sits
  while(red < blue) {
     if (arr[red] < pivot) red++;
     else {
        blue--;
        swap(arr, red, blue);
    }
  }
  swap(arr, blue, j); // put the pivot on the border
  return blue;
}


void stackQuickSort(int[] arr, int left, int right){
   Stack stack = new Stack(); //
   int border;
   stack.push(left); // assume the stack can keep ints
   stack.push(right);
   
   while(! stack.isempty()){
      int j = stack.pop();
      int i = stack.pop();
      border = partition(arr, i, j);
      if((border-1) - i > 0) {
         stack.push(i); 
         stack.push(border-1);
      }
      if(j - (border+1) > 0) {
         stack.push(border+1);
         stack.push(j);
      }
   }
}
*/

//
//  Merge Sort
//
static void mergesort(Foo a[], int lo, int hi, Foo scratch[])  {
    int k, mid, t_lo, t_hi;

        if (lo >= hi) {
            return;                 /* a[lo] is sorted already   */
        }

        mid = (lo+hi) / 2;
        mergesort(a, lo, mid, scratch);      /* Sort sublist a[lo..mid]   */
        mergesort(a, mid+1, hi, scratch);    /* Sort sublist a[mid+1..hi] */

        t_lo = lo; t_hi = mid+1;  
        for (k = lo; k <= hi; k++)            /* Merge sorted sublists    */
            if ((t_lo <= mid) && ((t_hi > hi) || (a[t_lo] < a[t_hi]))) {
                scratch[k] = a[t_lo++];
            }
            else {
                scratch[k] = a[t_hi++];
            }
        for (k = lo; k <= hi; k++)
            a[k] = scratch[k];
}
//
//  Bottom-up merge sort (doesn't work)
//
static void
merge4(Foo arr[], int iLo, int iUp, Foo temp[], int iSpan)
{
    int iEnd, i1, i2, j, k;

    i1 = iLo;
    i2 = iLo + iSpan;
    if (i2 > iUp)
        i2 = iUp;
    for (j = 1; j <= ((iUp - iLo + 1)/(iSpan*2)) + 1; j++) {
        iEnd = i2 + iSpan - 1;
        if (iEnd > iUp)
            iEnd = iUp;
        for (k = i1; k <= iEnd; k++) {
            if (arr[i1] <= arr[i2]) {
                temp[k] = arr[i1];
                arr[i1] = D_MACHINE_MAX;
                i1++;
            } else {
                temp[k] = arr[i2];
                arr[i2] = D_MACHINE_MAX;
                i2++;
                if (i2 > iEnd)
                    i2 = iEnd;
            }
        }
        i1 = iEnd + 1;
        i2 = i1 + 1;
        if (i2 > iUp)
            i2 = iUp;
    }
}

static void
mergesort4(Foo arr[], int iLo, int iUp)
{
    Foo temp[NBIG];
    int iSpan;

    for (iSpan = 1; iSpan <= iUp - iLo + 1; iSpan = 4*iSpan) {
        merge4(arr, iLo, iUp, temp, iSpan);
        merge4(temp, iLo, iUp, arr, iSpan*2);
    }
}
//
//  Floyd's Heapsort code (Algorithm 245, Treesort 3)
//  This does his O(N) buildheap invention, but does not do his
//  faster siftup trick.
//
static void
siftup(Foo M[], int i, int n)
{
    Foo copy;
    int j;

    copy = M[i];
loop: j = 2 * i;
    if (j <= n) {
        if (j < n) {
            if (M[j] < M[j+1])
                j = j+1;
        }
        if (copy < M[j]) {
            M[i] = M[j];
            i = j;
            goto loop;
        }
    }
    M[i] = copy;
}

static void
treesort3(Foo a[], int n)
{
    Foo t, M[NBIG + 1];
    int i;

    for (i = 0; i < n; i++)
        M[i+1] = a[i];      // Floyd's algol base-1 indexing!
    
    for (i = n/2; i >= 2; --i)
        siftup(M, i, n);
    for (i = n; i >= 2; --i) {
        siftup(M, 1, i);
        t = M[1];
        M[1] = M[i];
        M[i] = t;
    }

    for (i = 0; i < n; i++)
        a[i] = M[i+1];
}
//
//  Chvatal's accelerated heapsort
//
static int floor_of_lg(int n, int* power)
{
  /***********************************************************
  *  computes the largest power of 2 that is at most n and   *
  *  returns  the exponent in this power:                    *
  *            n =   1   2   3   4   5   6   7   8   9  ...  *
  *        power =   1   2   2   4   4   4   4   8   8  ...  *
  *   exponent k =   0   1   1   2   2   2   2   3   3  ...  *
  ***********************************************************/

  int k=0;
  for(*power=1; 2*(*power)<=n; (*power)*=2)
    k+=1;
  return k;
}

static void siftdown(Foo a[], int n, int vacant, Foo missing, int drop)
{
  int memo=vacant;
  int child, parent;
  int count, next_peek;

  count=0, next_peek=(drop+1)/2;

  child=2*(vacant+1);
  while(child<n)
    {
      if(ML_IsLess(a[child],a[child-1]))
        child--;
      a[vacant]=a[child], vacant=child, child=2*(vacant+1);

      count++;
      if (count==next_peek)
	{
	  if(ML_IsLessOrEqual(a[(vacant-1)/2],missing))
	    break;
	  else
	    next_peek=(count+drop+1)/2;	      
	}
    }

  if(child==n)
    {
      a[vacant]=a[n-1], vacant=n-1;
    }
  
  parent=(vacant-1)/2;
  while(vacant>memo)    
    {
      if(ML_IsLess(a[parent],missing))        
	{
          a[vacant]=a[parent], vacant=parent, parent=(vacant-1)/2;
	}
      else        
	break;
    }
  a[vacant]=missing;
}

static void makeheap(Foo a[], int n)
{
  int k, drop, first;
  
  drop=1, first=n/2-1;
  for(k=first; k>=0; k--)
    { 
      if(k==(first-1)/2)
	{
	  drop++, first=k;
	}
      siftdown(a, n, k, a[k], drop);
    }
}

static void hsort_c(Foo a[], int n)
{
  int k, drop, last;
  Foo temp;

  makeheap(a,n);

  drop = floor_of_lg(n-1, &last);
  for(k=n-1; k>0; k--) 
    {
      temp=a[k], a[k]=a[0];
      siftdown(a, k, 0, temp, drop);
      if (k==last)
	    {
	    drop--, last/=2;
	    }
    }
}

template <class T>
static void
BinaryInsertion(T r[], int lo, int up )
{
    int i, j, h, l;
    T tempr;

    for ( i=lo+1; i<=up; i++ ) {
        tempr = r[i];
        for ( l=lo-1, h=i; h-l > 1; ) {
            j = (h+l)/2;
            if ( tempr < r[j] )
                h = j;
            else
                l = j;
        }
        for (j = i; j > h + 4; j = j - 4) { //DPM a 10% optimization
            r[j] = r[j - 1];
            r[j - 1] = r[j - 2];
            r[j - 2] = r[j - 3];
            r[j - 3] = r[j - 4];
        }
        for (; j > h; j = j - 1)
            r[j] = r[j - 1];
        r[h] = tempr;
    }
}

static void
TestBinaryInsertionSort()
{
    Bar rgb[4];
    int i, n;
    unsigned nSum, nMin, nMax;

    printf("  BinaryInsertionSort\n");
    nSum = 0;
    nMax = 0;
    nMin = 10000;
    for (i = 0; i < 2*2; i++) {
        gs_nCompares = 0;
        PrepareCombination(rgb, 2, i);
        BinaryInsertion(rgb, 0, 1);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        if (gs_nCompares < nMin)
            nMin = gs_nCompares;
        AssertIsSorted(rgb, 2);
        ML_Assert(IsStable(rgb, 2));
    }
    n = (10000*nSum)/i;
    printf("      a. Two:       %d.%0.4d average compares %d min %d max\n",
        n/10000, n%10000, nMin, nMax);
    nSum = 0;
    nMax = 0;
    nMin = 10000;
    for (i = 0; i < 3*3*3; i++) {
        gs_nCompares = 0;
        PrepareCombination(rgb, 3, i);
        BinaryInsertion(rgb, 0, 2);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        if (gs_nCompares < nMin)
            nMin = gs_nCompares;
        AssertIsSorted(rgb, 3);
        ML_Assert(IsStable(rgb, 3));
    }
    n = (10000*nSum)/i;
    printf("      b. Three:     %d.%0.4d average compares %d min %d max\n",
        n/10000, n%10000, nMin, nMax);
    nSum = 0;
    nMax = 0;
    nMin = 10000;
    for (i = 0; i < 4*4*4*4; i++) {
        gs_nCompares = 0;
        PrepareCombination(rgb, 4, i);
        BinaryInsertion(rgb, 0, 3);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        if (gs_nCompares < nMin)
            nMin = gs_nCompares;
        AssertIsSorted(rgb, 4);
        ML_Assert(IsStable(rgb, 4));
    }
    n = (10000*nSum)/i;
    printf("      c. Four:      %d.%0.4d average compares %d min %d max\n",
        n/10000, n%10000, nMin, nMax);
}
//
//  Ford-Johnson Algorithm
//
struct FooBar {
    Foo foo;
    int iLabel;
};

int operator < (FooBar &a, FooBar &b) { return a.foo < b.foo; }

static void
MergeInsert(FooBar a[], FooBar b[], int N, int rgn[])
{
    int i, j, iLo, iHi, k, tk, iChain, i1, i2, iM, NB, nDelta;
    FooBar t;
    unsigned short *rgi;

    if (N <= 1)
        return;
    //
    //  Pairwise round, move winners (larger) to the front half.
    //
    for (i = 0; i < N/2; i++) {
        j = i + N/2;
        if (a[i].foo < a[j].foo) {
               t = a[j];
            a[j] = a[i];
            a[i] = t;
        }
    }
    //
    //  Sort the winners of the first round, a[0]...a[N/2 - 1].  Apply
    //  the same permutation to the losers, a[N/2]...a{N - 1], as they
    //  are moved into b[0]...b[N/2 - 1].
    //
    rgi = (unsigned short *) _alloca(sizeof(unsigned short) * N/2);
    for (i = 0; i < N/2; i++)
        rgi[i] = a[i].iLabel;
    MergeInsert(a, b, N/2, rgn);
    for (i = 0; i < N/2; i++)
        rgn[a[i].iLabel] = i;
    for (i = 0; i < N/2; i++)
        b[rgn[rgi[i]]] = a[i + N/2];
    if (N & 1) {
        b[N/2] = a[N - 1];  // the odd element
        NB = N/2 + 1;
    } else
        NB = N/2;
    for (i = 0; i < N/2; i++)
        ML_Assert(b[i].foo <= a[i].foo);
    gs_nCompares -= N/2;
    for (i = 0; i < N/2-1; i++)
        ML_Assert(a[i].foo < a[i + 1].foo);
    gs_nCompares -= N/2 - 1;
    //
    //  Binary insertion of b[] into a[] according to a schedule
    //
    for (i = 0; i < N/2; i++)
        a[i + N - N/2] = a[i];
    iChain = N - N/2 - 1;
    a[iChain] = b[0];      // the main chain is in place
    for (i = iChain; i < N - 1; i++) {
        ML_Assert(a[i] < a[i + 1]);
        gs_nCompares--;
    }
    tk = 1;
    k = 2;
    for (;;) {
        //
        //  insert b[iHi] to b[iLo] into a[iChain] to a[iChain + k - 2]
        //  tk is computed such that the binary insertions are into sets
        //  of size 2**n - 1, making an optimal use of n comparisons.
        //
        k = 2*k;
        iLo = tk;
        tk = k - tk;
        iHi = tk - 1;
        if (iHi >= NB)
            iHi = NB - 1;
        if (iLo > iHi)
            break;
        i1 = iChain - 1;
        i2 = iChain + k - 1;
        if (iChain + k - 1 > N)
            nDelta = N - iChain + 1;
        else
            nDelta = k; // try to do binary insertion into 2**n - 1 elements
        for (i = iHi; i >= iLo; --i) {
            t = b[i];
            i1 = iChain - 1;
            ML_Assert(iChain >= 0);
            i2 = i1 + nDelta;
            ML_Assert(i2 <= N);
            while(i2 - i1 > 1) {
                iM = (i1 + i2)/2;
                if (t.foo < a[iM].foo)
                    i2 = iM;
                else
                    i1 = iM;
            }
            --iChain;
            for (j = iChain; j < i1; j++)
                a[j] = a[j + 1];
            a[j] = t;
        }
    }
    ML_Assert(iChain == 0);
    //
    //  Insert the odd element if there was one
    //
    /*
    if (N & 1) {
        t = b[N/2];
        i1 = 0;
        i2 = N;
        while(i2 - i1 > 1) {
            iM = (i1 + i2)/2;
            if (t.foo < a[iM].foo)
                i2 = iM;
            else
                i1 = iM;
        }
        for (j = 0; j < i1; j++)
            a[j] = a[j + 1];
        a[i1] = t;
    }
    */
}

static void
FordJohnson(Foo a[], int N)
{
    FooBar b[NBIG], c[NBIG];
    int i, rgn[NBIG];

    for (i = 0; i < N; i++) {
        b[i].foo = a[i];
        b[i].iLabel = i;
    }
    MergeInsert(b, c, N, rgn);
    for (i = 0; i < N; i++)
        a[i] = b[i].foo;
}

#include <stdlib.h>

/*
    6. Selection:    (49995000 49995000 cmps)   1.708 sorts/sec.
    7. Insertion:    (12507283 17445075 cmps)   6.617 sorts/sec.
    8. Bubble:       (49995000 49995000 cmps)   0.695 sorts/sec.
    9. Shellsort:    (  235155   251170 cmps) 133.086 sorts/sec.
    10. Treesort3:   (  235374   235685 cmps) 230.540 sorts/sec.
    11. libc qsort:  (  154538   170226 cmps) 145.758 sorts/sec.
    12. Quicksort S: (  167212   201192 cmps) 301.860 sorts/sec.
    13. Quicksort M: (  160429   179562 cmps) 303.601 sorts/sec.
    14. Quicksort A: (  141961   160666 cmps) 234.870 sorts/sec.
    15. Heapsort:    (  138965   139183 cmps) 269.903 sorts/sec.
    16. Mergesort:   (  120450   120667 cmps) 250.205 sorts/sec.
    17. BinaryInsert: ( 118996   119101 cmps)  11.417 sorts/sec.
    18. FordJohnson: (  118817   118863 cmps)   8.187 sorts/sec.

VC++ .NET 2003:
    6. Selection:    (49995000 49995000 cmps)   3.529 sorts/sec
    7. Insertion:    (24988846 25241798 cmps)   8.148 sorts/sec
    8. Bubble:       (49995000 49995000 cmps)   0.683 sorts/sec
    9. Shellsort:    (  209759   215161 cmps) 126.330 sorts/sec
    10. Treesort3:   (  235379   235685 cmps) 226.210 sorts/sec
    11. libc qsort:  (  154353   165092 cmps) 143.976 sorts/sec
    12. Quicksort S: (  166870   192518 cmps) 291.647 sorts/sec
    13. Quicksort M: (  160044   170990 cmps) 319.419 sorts/sec
    14. Quicksort A: (  141613   152587 cmps) 292.894 sorts/sec
    15. Heapsort:    (  138965   139173 cmps) 255.580 sorts/sec
    16. Mergesort:   (  120447   120591 cmps) 247.144 sorts/sec
    17. BinaryInsert: ( 119003   119061 cmps)  11.371 sorts/sec
    18. FordJohnson: (  118816   118846 cmps)   8.389 sorts/sec

Time/P4/Intrinsics:
    6. Selection:    (49995000 49995000 cmps)   3.966 sorts/sec
    7. Insertion:    (24988846 25241798 cmps)   8.340 sorts/sec
    8. Bubble:       (49995000 49995000 cmps)   0.691 sorts/sec
    9. Shellsort:    (  209759   215161 cmps) 123.454 sorts/sec
    10. Treesort3:   (  235379   235685 cmps) 240.578 sorts/sec
    11. libc qsort:  (  154353   165092 cmps) 145.651 sorts/sec
    12. Quicksort S: (  166870   192518 cmps) 299.015 sorts/sec
    13. Quicksort M: (  160044   170990 cmps) 328.407 sorts/sec
    14. Quicksort A: (  141613   152587 cmps) 297.932 sorts/sec
    15. Heapsort:    (  138965   139173 cmps) 259.841 sorts/sec
    16. Mergesort:   (  120447   120591 cmps) 252.135 sorts/sec
    17. BinaryInsert: ( 119003   119061 cmps)  11.428 sorts/sec
    18. FordJohnson: (  118816   118846 cmps)   8.404 sorts/sec

Full Optimization: Ox
Pentium 4, 1.8 GHz:
    6. Selection:    (49995000 49995000 cmps)   4.001 sorts/sec
    7. Insertion:    (24988846 25241798 cmps)   8.395 sorts/sec
    8. Bubble:       (49995000 49995000 cmps)   0.690 sorts/sec
    9. Shellsort:    (  209759   215161 cmps) 123.051 sorts/sec
    10. Treesort3:   (  235379   235685 cmps) 240.869 sorts/sec
    11. libc qsort:  (  154353   165092 cmps) 146.959 sorts/sec
    12. Quicksort S: (  166870   192518 cmps) 297.076 sorts/sec
    13. Quicksort M: (  160044   170990 cmps) 329.523 sorts/sec
    14. Quicksort A: (  141613   152587 cmps) 297.366 sorts/sec
    15. Heapsort:    (  138965   139173 cmps) 260.039 sorts/sec
    16. Mergesort:   (  120447   120591 cmps) 254.872 sorts/sec
    17. BinaryInsert: ( 119003   119061 cmps)  11.386 sorts/sec
    18. FordJohnson: (  118816   118846 cmps)   8.402 sorts/sec

Skylake Xeon E3-1275v5 3.6 GHz:
  Selection:   (49995000 49995000 cmps)       47.658 sorts/sec
  Insertion:   (24988846 25241798 cmps)       60.346 sorts/sec
  Bubble:      (49995000 49995000 cmps)        7.105 sorts/sec
  Shellsort:   (  209759   215161 cmps)      956. 59 sorts/sec
  Treesort3:   (  235379   235685 cmps)     1314.652 sorts/sec
  libc qsort:  (  154353   165092 cmps)      919.708 sorts/sec
  Quicksort S: (  166870   192518 cmps)     1567.819 sorts/sec
  Quicksort M: (  160044   170990 cmps)     1804.696 sorts/sec
  Quicksort A: (  141613   152587 cmps)     1601.398 sorts/sec
  Heapsort:    (  138965   139173 cmps)     1325.391 sorts/sec
  Mergesort:   (  120447   120591 cmps)     1329.687 sorts/sec
  BinaryInsert: ( 119003   119061 cmps)      116.919 sorts/sec
  FordJohnson: (   88844    88874 cmps)      201.641 sorts/sec
    std::sort: (  208913   213901 cmps)     1577.393 sorts/sec
  std::stable: (  172593   173621 cmps)     1555.961 sorts/sec

*/
static void
TestOtherSorts()
{
    Foo rgn[NBIG], scratch[NBIG];
    int i;
    ML_TimingInfo ti;
    unsigned nX, nC, nSum, nMax;

    nX = ML_nMarsagliaX;
    nC = ML_nMarsagliaC;
    ML_InitializeRandom();
    for (i = 0; i < NBIG; i++)
        rgn[i].f = double(i);
    ML_RandomShuffle(rgn, NBIG);

    printf("  Selection:   ");
    gs_nCompares = 0;
    nMax = 0;
    nSum = 0;
    ML_StartTiming(ti);
    for (i = 0; i < 2; i++) {
        ML_RandomShuffle(rgn, NBIG);
        selection(rgn, 0, NBIG-1);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        gs_nCompares = 0;
    }
    ML_StopTiming(ti);
    printf("(%8d %8d cmps) ", nSum/i, nMax);
    VerifySort(rgn, NBIG);
    ML_ReportEventsPerSecond(ti, i, "sorts");

    printf("  Insertion:   ");
    gs_nCompares = 0;
    nMax = 0;
    nSum = 0;
    ML_StartTiming(ti);
    for (i = 0; i < 7; i++) {
        ML_RandomShuffle(rgn, NBIG);
        ML_InsertionSort(rgn, NBIG);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        gs_nCompares = 0;
    }
    ML_StopTiming(ti);
    printf("(%8d %8d cmps) ", nSum/i, nMax);
    VerifySort(rgn, NBIG);
    ML_ReportEventsPerSecond(ti, i, "sorts");

    printf("  Bubble:      ");
    gs_nCompares = 0;
    nMax = 0;
    nSum = 0;
    ML_StartTiming(ti);
    for (i = 0; i < 1; i++) {
        ML_RandomShuffle(rgn, NBIG);
        bubble(rgn, 0, NBIG-1);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        gs_nCompares = 0;
    }
    ML_StopTiming(ti);
    printf("(%8d %8d cmps) ", nSum/i, nMax);
    VerifySort(rgn, NBIG);
    ML_ReportEventsPerSecond(ti, i, "sorts");

    printf("  Shellsort:   ");
    gs_nCompares = 0;
    nMax = 0;
    nSum = 0;
    ML_StartTiming(ti);
    for (i = 0; i < 40; i++) {
        ML_RandomShuffle(rgn, NBIG);
        shellsortP(rgn, 0, NBIG-1);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        gs_nCompares = 0;
    }
    ML_StopTiming(ti);
    printf("(%8d %8d cmps) ", nSum/i, nMax);
    VerifySort(rgn, NBIG);
    ML_ReportEventsPerSecond(ti, i, "sorts");

    printf("  Treesort3:   ");
    gs_nCompares = 0;
    nMax = 0;
    nSum = 0;
    ML_StartTiming(ti);
    for (i = 0; i < 80; i++) {
        ML_RandomShuffle(rgn, NBIG);
        treesort3(rgn, NBIG);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        gs_nCompares = 0;
    }
    ML_StopTiming(ti);
    printf("(%8d %8d cmps) ", nSum/i, nMax);
    VerifySort(rgn, NBIG);
    ML_ReportEventsPerSecond(ti, i, "sorts");

    printf("  libc qsort:  ");
    ML_InitializeRandom();
    gs_nCompares = 0;
    nMax = 0;
    nSum = 0;
    for (i = 0; i < NBIG; i++)
        rgn[i].f = double(i);
    ML_StartTiming(ti);
    for (i = 0; i < 50; i++) {
        ML_RandomShuffle(rgn, NBIG);
        qsort(rgn, NBIG, sizeof(double), CompareFoo);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        gs_nCompares = 0;
    }
    ML_StopTiming(ti);
    printf("(%8d %8d cmps) ", nSum/i, nMax);
    VerifySort(rgn, NBIG);
    ML_ReportEventsPerSecond(ti, i, "sorts");

    printf("  Quicksort S: ");
    ML_InitializeRandom();
    gs_nCompares = 0;
    nMax = 0;
    nSum = 0;
    for (i = 0; i < NBIG; i++)
        rgn[i].f = double(i);
    ML_StartTiming(ti);
    for (i = 0; i < 100; i++) {
        ML_RandomShuffle(rgn, NBIG);
        quicksort1(rgn, 0, NBIG-1);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        gs_nCompares = 0;
    }
    ML_StopTiming(ti);
    printf("(%8d %8d cmps) ", nSum/i, nMax);
    VerifySort(rgn, NBIG);
    ML_ReportEventsPerSecond(ti, i, "sorts");

    printf("  Quicksort M: ");
    ML_InitializeRandom();
    gs_nCompares = 0;
    nMax = 0;
    nSum = 0;
    for (i = 0; i < NBIG; i++)
        rgn[i].f = double(i);
    ML_StartTiming(ti);
    for (i = 0; i < 100; i++) {
        ML_RandomShuffle(rgn, NBIG);
        ML_QuickSort(rgn, NBIG);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        gs_nCompares = 0;
    }
    ML_StopTiming(ti);
    printf("(%8d %8d cmps) ", nSum/i, nMax);
    VerifySort(rgn, NBIG);
    ML_ReportEventsPerSecond(ti, i, "sorts");

    printf("  Quicksort A: ");
    ML_InitializeRandom();
    gs_nCompares = 0;
    nMax = 0;
    nSum = 0;
    for (i = 0; i < NBIG; i++)
        rgn[i].f = double(i);
    ML_StartTiming(ti);
    for (i = 0; i < 100; i++) {
        ML_RandomShuffle(rgn, NBIG);
        AhrensQuicksort2(rgn, 0, NBIG-1);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        gs_nCompares = 0;
    }
    ML_StopTiming(ti);
    printf("(%8d %8d cmps) ", nSum/i, nMax);
    VerifySort(rgn, NBIG);
    ML_ReportEventsPerSecond(ti, i, "sorts");

    printf("  Heapsort:    ");
    ML_InitializeRandom();
    gs_nCompares = 0;
    nMax = 0;
    nSum = 0;
    for (i = 0; i < NBIG; i++)
        rgn[i].f = double(i);
    ML_StartTiming(ti);
    for (i = 0; i < 100; i++) {
        ML_RandomShuffle(rgn, NBIG);
        ML_HeapSort(rgn, NBIG);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        gs_nCompares = 0;
    }
    ML_StopTiming(ti);
    printf("(%8d %8d cmps) ", nSum/i, nMax);
    VerifySort(rgn, NBIG);
    ML_ReportEventsPerSecond(ti, i, "sorts");

    printf("  Mergesort:   ");
    ML_InitializeRandom();
    gs_nCompares = 0;
    nMax = 0;
    nSum = 0;
    for (i = 0; i < NBIG; i++)
        rgn[i].f = double(i);
    ML_StartTiming(ti);
    for (i = 0; i < 100; i++) {
        ML_RandomShuffle(rgn, NBIG);
        mergesort(rgn, 0, NBIG-1, scratch);
        //mergesort4(rgn, 0, NBIG-1);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        gs_nCompares = 0;
    }
    ML_StopTiming(ti);
    printf("(%8d %8d cmps) ", nSum/i, nMax);
    VerifySort(rgn, NBIG);
    ML_ReportEventsPerSecond(ti, i, "sorts");

    printf("  BinaryInsert: ");
    ML_InitializeRandom();
    gs_nCompares = 0;
    nMax = 0;
    nSum = 0;
    for (i = 0; i < NBIG; i++)
        rgn[i].f = double(i);
    ML_StartTiming(ti);
    for (i = 0; i < 10; i++) {
        ML_RandomShuffle(rgn, NBIG);
        BinaryInsertion(rgn, 0, NBIG-1);
        VerifySort(rgn, NBIG);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        gs_nCompares = 0;
    }
    ML_StopTiming(ti);
    printf("(%7d %8d cmps) ", nSum/i, nMax);
    VerifySort(rgn, NBIG);
    ML_ReportEventsPerSecond(ti, i, "sorts");

    printf("  FordJohnson: ");
    ML_InitializeRandom();
    gs_nCompares = 0;
    nMax = 0;
    nSum = 0;
    for (i = 0; i < NBIG; i++)
        rgn[i].f = double(i);
    ML_StartTiming(ti);
    for (i = 0; i < 10; i++) {
        ML_RandomShuffle(rgn, NBIG);
        FordJohnson(rgn, NBIG);
        VerifySort(rgn, NBIG);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        gs_nCompares = 0;
    }
    ML_StopTiming(ti);
    printf("(%8d %8d cmps) ", nSum/i, nMax);
    //for (i = 0; i < NBIG; i++)
    //    printf("  %d\n", int(rgn[i].f));
    VerifySort(rgn, NBIG);
    ML_ReportEventsPerSecond(ti, i, "sorts");

    printf("    std::sort: ");
    ML_InitializeRandom();
    gs_nCompares = 0;
    nMax = 0;
    nSum = 0;
    for (i = 0; i < NBIG; i++)
        rgn[i].f = double(i);
    ML_StartTiming(ti);
    for (i = 0; i < 10; i++) {
        ML_RandomShuffle(rgn, NBIG);
        std::sort(rgn, rgn+NBIG-1);
        VerifySort(rgn, NBIG);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        gs_nCompares = 0;
    }
    ML_StopTiming(ti);
    printf("(%8d %8d cmps) ", nSum/i, nMax);
    //for (i = 0; i < NBIG; i++)
    //    printf("  %d\n", int(rgn[i].f));
    VerifySort(rgn, NBIG);
    ML_ReportEventsPerSecond(ti, i, "sorts");

    printf("  std::stable: ");
    ML_InitializeRandom();
    gs_nCompares = 0;
    nMax = 0;
    nSum = 0;
    for (i = 0; i < NBIG; i++)
        rgn[i].f = double(i);
    ML_StartTiming(ti);
    for (i = 0; i < 10; i++) {
        ML_RandomShuffle(rgn, NBIG);
        std::stable_sort(rgn, rgn+NBIG-1);
        VerifySort(rgn, NBIG);
        nSum += gs_nCompares;
        if (gs_nCompares > nMax)
            nMax = gs_nCompares;
        gs_nCompares = 0;
    }
    ML_StopTiming(ti);
    printf("(%8d %8d cmps) ", nSum/i, nMax);
    //for (i = 0; i < NBIG; i++)
    //    printf("  %d\n", int(rgn[i].f));
    VerifySort(rgn, NBIG);
    ML_ReportEventsPerSecond(ti, i, "sorts");

    ML_InitializeRandom(nX, nC);
}
/*
Examine Quicksort minimum-item parameter
  QS M (   2): (  159981   170898 cmps)     1831.319 sorts/sec
  QS M (   4): (  159979   170896 cmps)     1879.417 sorts/sec
  QS M (   8): (  159979   170896 cmps)     1864. 31 sorts/sec
  QS M (  16): (  159988   170914 cmps)     1854. 73 sorts/sec
  QS M (  32): (  160053   170990 cmps)     1869.324 sorts/sec
  QS M (  64): (  160356   171393 cmps)     1856.158 sorts/sec
  QS M ( 128): (  161980   172107 cmps)     1814. 18 sorts/sec
  QS M ( 256): (  169054   181178 cmps)     1858.876 sorts/sec
  QS M ( 512): (  198988   223654 cmps)     1825.291 sorts/sec
  QS M (1024): (  323745   412545 cmps)     1596.511 sorts/sec
  
void
TestQuickSort2()
{
    Foo rgn[NBIG], scratch[NBIG];
    int i, k;
    ML_TimingInfo ti;
    unsigned nSum, nMax;

    printf("Examine Quicksort parameters:\n");
    for (k = 2; k <= 1024; k = k+k) {
        printf("  QS M (%4d): ", k);
        ML_InitializeRandom();
        gs_nCompares = 0;
        nMax = 0;
        nSum = 0;
        for (i = 0; i < NBIG; i++)
            rgn[i].f = double(i);
        ML_StartTiming(ti);
        for (i = 0; i < 100; i++) {
            ML_RandomShuffle(rgn, NBIG);
            ML_QuickSort(rgn, NBIG, k);         // temporarily add nMinQuickItems parameter
            nSum += gs_nCompares;
            if (gs_nCompares > nMax)
                nMax = gs_nCompares;
            gs_nCompares = 0;
        }
        ML_StopTiming(ti);
        printf("(%8d %8d cmps) ", nSum/i, nMax);
        VerifySort(rgn, NBIG);
        ML_ReportEventsPerSecond(ti, i, "sorts");
    }
}
*/

void
TestSort()
{
    printf("Testing Sort Routines\n");
    TestQuickSort();
    TestInsertionSort();
    //TestBinaryInsertionSort();
    TestBinarySearch();
    printf("    SortTwo, SortThree, SortFour\n");
    TestSortTwo();
    TestSortThree();
    TestSortFour();
    TestOtherSorts();
    // TestQuickSort2();
}
