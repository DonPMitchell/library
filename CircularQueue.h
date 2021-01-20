#pragma once
////
//  ML_CircularQueue - circular queue, concurrent producer/consumer
//  D. P. Mitchell  06/02/2001
//

template <class T>
class ML_CircularQueue {
public:
                    ML_CircularQueue() : ptFirst(0) {}
                    ~ML_CircularQueue();

    int             SetLength(int maxItems);
    int             IsFull() const;
    int             IsEmpty() const;
    int             Enqueue(const T &tItem);
    int             Dequeue(T *pt);
private:
                    ML_CircularQueue(const ML_CircularQueue &);   // disallow copying
    ML_CircularQueue&  operator =(const ML_CircularQueue &);

    T * volatile    ptFirst;
    T * volatile    ptIn;
    T * volatile    ptOut;
    T * volatile    ptLimit;
};

template <class T>
inline int
ML_CircularQueue<T>::SetLength(int maxItems)
{
    delete[] ptFirst;
    ptFirst = new T[maxItems + 1];  // size - 1 items in a full queue
    ptIn = ptOut = ptFirst;
    ptLimit = ptFirst + maxItems + 1;
    return ptFirst != 0 && maxItems >= 1;
}

template <class T>
inline int
ML_CircularQueue<T>::IsFull() const
{
    return (ptIn + 1 == ptOut) || (ptIn + 1 == ptLimit && ptOut == ptFirst);
}

template <class T>
inline int
ML_CircularQueue<T>::IsEmpty() const
{
    return ptIn == ptOut;
}

template <class T>
inline int
ML_CircularQueue<T>::Enqueue(const T &tItem)
{
    if (IsFull())
        return 0;
    *ptIn = tItem;
    if (ptIn + 1 == ptLimit)
        ptIn = ptFirst;
    else
        ptIn++;
    return 1;
}

template <class T>
inline int
ML_CircularQueue<T>::Dequeue(T *pt)
{
    if (IsEmpty())
        return 0;
    *pt = *ptOut;
    if (ptOut + 1 == ptLimit)
        ptOut = ptFirst;
    else
        ptOut++;
    return 1;
}

template<class T>
ML_CircularQueue<T>::~ML_CircularQueue()
{
    delete[] ptFirst;
    ptFirst = ptIn = ptOut = ptLimit = 0;
}
