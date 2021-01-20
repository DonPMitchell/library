#pragma once

inline int ML_IsLess(int n1, int n2)                      { return n1 < n2; }
inline int ML_IsLess(unsigned n1, unsigned n2)            { return n1 < n2; }
inline int ML_IsLess(double f1, double f2)                { return f1 < f2; }
inline int ML_IsLess(const char *sz1, const char *sz2)
{
    int nResult, nByte;

    while ((nByte = *(unsigned char *)sz2), !(nResult = *(unsigned char *)sz1 - nByte) && nByte) {
        sz1++;
        sz2++;
    };
    return nResult < 0;
}
