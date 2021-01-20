#pragma once
//
//  sort - template functions for fast sorting, searching and shuffling.
//  D. P. Mitchell
//
//  Insertion sort is useful because it is O(N) for almost-sorted arrays.
//  Insertion sort is STABLE.  Uses ML_IsLess()
//
#include "Utilities.h"
#include "IsLess.h"

template <class T>
void
ML_InsertionSort(T rgtItem[], int nItems)
{
    int i, j;
    T tKey;

    for (j = 1; j < nItems; j++) {
        tKey = rgtItem[j];
        for (i = j; i > 0 && ML_IsLess(tKey, rgtItem[i - 1]); i = i - 1)
            rgtItem[i] = rgtItem[i - 1];    // faster than [i+1] = [i]
        rgtItem[i] = tKey;
    }
}
//
//  Efficient, stable, fixed-size sorting routines.  Use ML_IsLess()
//
template <class T>
inline void
ML_SortTwo(T rgtItem[])
{
    T tTemp;

    if (ML_IsLess(rgtItem[1], rgtItem[0])) {
        tTemp = rgtItem[0];
        rgtItem[0] = rgtItem[1];
        rgtItem[1] = tTemp;
    }
}
//
//  2.4814 comparisons on average, worst case 3. best case 2
//
template <class T>
inline void
ML_SortThree(T rgtItem[])
{
    T tTemp;

    if (ML_IsLess(rgtItem[1], rgtItem[0])) {  // insertion sort
        tTemp = rgtItem[0];
        rgtItem[0] = rgtItem[1];
        rgtItem[1] = tTemp;
    }
    if (ML_IsLess(rgtItem[2], rgtItem[1])) {
        tTemp = rgtItem[2];
        rgtItem[2] = rgtItem[1];
        if (ML_IsLess(tTemp, rgtItem[0])) {
            rgtItem[1] = rgtItem[0];
            rgtItem[0] = tTemp;
        } else
            rgtItem[1] = tTemp;
    }
}
//
//  Average 4.5312 compares, worst case 5, best case 4
//  Sort first three and then do binary insertion.
//
template <class T>
inline void
ML_SortFour(T rgtItem[])
{
    T tTemp;

    if (ML_IsLess(rgtItem[1], rgtItem[0])) {
        tTemp = rgtItem[0];
        rgtItem[0] = rgtItem[1];
        rgtItem[1] = tTemp;
    }
    if (ML_IsLess(rgtItem[2], rgtItem[1])) {
        tTemp = rgtItem[2];
        rgtItem[2] = rgtItem[1];
        if (ML_IsLess(tTemp, rgtItem[0])) {
            rgtItem[1] = rgtItem[0];
            rgtItem[0] = tTemp;
        } else
            rgtItem[1] = tTemp;
    }
    if (ML_IsLess(rgtItem[3], rgtItem[1])) {
        tTemp = rgtItem[3];
        rgtItem[3] = rgtItem[2];
        rgtItem[2] = rgtItem[1];
        if (ML_IsLess(tTemp, rgtItem[0])) {
            rgtItem[1] = rgtItem[0];
            rgtItem[0] = tTemp;
        } else {
            rgtItem[1] = tTemp;
        }
    } else {
        if (ML_IsLess(rgtItem[3], rgtItem[2])) {
            tTemp = rgtItem[3];
            rgtItem[3] = rgtItem[2];
            rgtItem[2] = tTemp;
        }
    }
}
//
//  Binary search returns the index of the first element e[i] >= tX.  Uses
//  ML_IsLess().
//
template <class T>
int
ML_BinarySearch(T rgtItem[], int nItems, T &tX)
{
    int iLo, iMid, iHi;

    for (iLo = -1, iHi = nItems; iHi - iLo > 1; ) {
        iMid = (iHi + iLo)/2;
        if (ML_IsLess(tX, rgtItem[iMid]))
            iHi = iMid;
        else
            iLo = iMid;
    }
    return iHi;
}
//
//  Quicksort: uses only ML_IsLess() function.  Quite well optimized.
//
template <class T>
void
ML_QuickSort(T rgtItem[], int nItems)
{
    int i, j;
    T tTemp, tPivot;

    nItems--;   // iLast
    while (nItems >= 32) {
        i = nItems/2;
        //
        //  Sort first, middle and last elements.  This provides median-of-3
        //  partitioning, limits partitioning to only N - 3 remaining items,
        //  and creates sentinals to simplify the inner loop.
        //
        if (ML_IsLess(rgtItem[i], rgtItem[0])) {      // 2.48 compares on average
            tTemp = rgtItem[0];
            rgtItem[0] = rgtItem[i];
            rgtItem[i] = tTemp;
        }
        if (ML_IsLess(rgtItem[nItems], rgtItem[i])) {
            tTemp = rgtItem[nItems];
            rgtItem[nItems] = rgtItem[i];
            if (ML_IsLess(tTemp, rgtItem[0])) {
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
            while(ML_IsLess(rgtItem[++i], tPivot))
                ;
            while(ML_IsLess(tPivot, rgtItem[--j]))
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
            ML_QuickSort(rgtItem ,j + 1);
            rgtItem += i + 1;
            nItems -= i + 1;
        } else {
            ML_QuickSort(rgtItem + i + 1, nItems - i);
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
        if (ML_IsLess(tTemp,rgtItem[j - 1])) {
            do {
                rgtItem[j] = rgtItem[j - 1];
                j = j - 1;
            } while (j > 0 && ML_IsLess(tTemp, rgtItem[j - 1]));    //REVIEW: forloop  faster with MSVC
            rgtItem[j] = tTemp;
        }
    }
}
