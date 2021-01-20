#pragma once
////
//  Priority-queue implementations
//  D. P. Mitchell  07/07/2002
//
#include "IsLess.h"

//
//  Binary Heap: rgt[0] ... rgt[nSize - 1]
//               rgt[0] is the maximum element
//
template <class T>
inline void
ML_DownHeap(T rgt[], int i, int nSize)  // heapify
{
    T tTemp;
    int iChild;

    ML_Assert(i < nSize);
    tTemp = rgt[i];
    iChild = i + i + 2;
    while (iChild < nSize) {
        if (ML_IsLess(rgt[iChild], rgt[iChild - 1]))
            iChild = iChild - 1;    // go down other child
        if (ML_IsLess(tTemp, rgt[iChild])) {
            rgt[i] = rgt[iChild];
            i = iChild;
            iChild = i + i + 2;
        } else
            break;
    }
    if (iChild == nSize && ML_IsLess(tTemp, rgt[nSize - 1])) {
        rgt[i] = rgt[nSize - 1];
        i = nSize - 1;
    }
    rgt[i] = tTemp;
}
//
//  Floyd's trick to reduce compares almost in half.  Swapping a leaf
//  into the root ends up pushing it back down almost to the leaves
//  again.  So push it to the leaf and then UpHeap.  Use only with
//  HeapSort and DeleteMaxHeap!
//
template <class T>
inline void
ML_FloydDownHeap(T rgt[], int i, int nSize)
{
    T tTemp;
    int iChild, iParent, iStart;

    ML_Assert(i < nSize);
    iStart = i;
    tTemp = rgt[i];
    iChild = i + i + 2;
    while (iChild < nSize) {                        // push to leaf
        if (ML_IsLess(rgt[iChild], rgt[iChild - 1]))
            iChild = iChild - 1;
        rgt[i] = rgt[iChild];
        i = iChild;
        iChild = i + i + 2;
    }
    if (iChild == nSize) {
        rgt[i] = rgt[nSize - 1];
        i = nSize - 1;
    }
    iParent = (i - 1) >> 1;
    while (i > iStart) {                            // then upheap
        if (ML_IsLess(rgt[iParent], tTemp)) {
            rgt[i] = rgt[iParent];
            i = iParent;
            iParent = (i - 1) >> 1;
        } else
            break;
    }
    rgt[i] = tTemp;
}
//
//  Floyd's O(N) algorithm
//
template <class T>
inline void
ML_BuildHeap(T rgt[], int nSize)
{
    int i;

    for (i = (nSize >> 1) - 1; i >= 0; --i)
        ML_DownHeap(rgt, i, nSize);
}

template <class T>
inline void
ML_DeleteMaxHeap(T rgt[], int nSize)
{
    if (nSize > 1) {
        rgt[0] = rgt[nSize - 1];
        ML_FloydDownHeap(rgt, 0, nSize - 1);
    }
}

template <class T>
inline void
ML_InsertHeap(T rgt[], T t, int nSize) // UpHeap
{
    int i, iParent;

    i = nSize;      // new size is nSize + 1
    while (i > 0 && ML_IsLess(rgt[iParent = (i - 1) >> 1], t)) {
        rgt[i] = rgt[iParent];
        i = iParent;
    }
    rgt[i] = t;
}

template <class T>
inline void
ML_DeleteHeap(T rgt[], int i, int nSize)
{
    T tTemp;
    int iParent;

    ML_Assert(i < nSize);
    tTemp = rgt[--nSize];
    if (ML_IsLessOrEqual(rgt[i], tTemp)) {
        while (i > 0 && ML_IsLess(rgt[iParent = (i - 1) >> 1], tTemp)) {
            rgt[i] = rgt[iParent];
            i = iParent;
        }
        rgt[i] = tTemp;
    } else {
        rgt[i] = tTemp;
        ML_DownHeap(rgt, i, nSize);
    }
}

template <class T>
inline void
ML_HeapSort(T rgt[], int nSize)
{
    int i;
    T tTemp;

    ML_BuildHeap(rgt, nSize);
    for (i = nSize - 1; i >= 1; --i) {
        tTemp = rgt[0];
        rgt[0] = rgt[i];
        rgt[i] = tTemp;
        ML_FloydDownHeap(rgt, 0, i);
    }
}
