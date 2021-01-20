#include "stdafx.h"
#include "Utilities.h"

//
//  Utilities for timing tasks.  Pin the process to one CPU to stablize for fine measurment.
//
static LARGE_INTEGER    gs_nTicksPerSecond;

void
ML_StartTiming(ML_TimingInfo &ti)
{
    FILETIME ftJunk;

    GetProcessTimes(GetCurrentProcess(), &ftJunk, &ftJunk, &ti.ftStartKernel, &ti.ftStartUser);
    QueryPerformanceCounter(&ti.nStartReal);
    ti.nNullReal = ti.nStartReal;
    ti.ftNullUser = ti.ftStartUser;
    ti.ftNullKernel = ti.ftStartKernel;
}

void
ML_StopTiming(ML_TimingInfo &ti)
{
    FILETIME ftJunk;
    LARGE_INTEGER nElapsedReal, nElapsedKernel, nElapsedUser, nStart, nStop, nNull;

    QueryPerformanceCounter(&ti.nStopReal);
    GetProcessTimes(GetCurrentProcess(), &ftJunk, &ftJunk, &ti.ftStopKernel, &ti.ftStopUser);
    //
    //  Calculate elapsed times in seconds
    //
    if (gs_nTicksPerSecond.LowPart == 0)
        QueryPerformanceFrequency(&gs_nTicksPerSecond);
    nElapsedReal.QuadPart = ti.nStopReal.QuadPart - ti.nNullReal.QuadPart;
    nElapsedReal.QuadPart -= ti.nNullReal.QuadPart - ti.nStartReal.QuadPart;        // deduct null-loop time
    ti.fReal = double(nElapsedReal.QuadPart)/double(gs_nTicksPerSecond.QuadPart);

    nStart.LowPart =  ti.ftStartUser.dwLowDateTime;
    nStart.HighPart = ti.ftStartUser.dwHighDateTime;
    nStop.LowPart =  ti.ftStopUser.dwLowDateTime;
    nStop.HighPart = ti.ftStopUser.dwHighDateTime;
    nNull.LowPart =  ti.ftNullUser.dwLowDateTime;
    nNull.HighPart = ti.ftNullUser.dwHighDateTime;
    nElapsedUser.QuadPart = nStop.QuadPart - nNull.QuadPart;
    nElapsedUser.QuadPart -= nNull.QuadPart - nStart.QuadPart;
    ti.fUser = double(nElapsedUser.QuadPart)/10000000.0;    // 100 nanoseconds per tick

    nStart.LowPart =  ti.ftStartKernel.dwLowDateTime;
    nStart.HighPart = ti.ftStartKernel.dwHighDateTime;
    nStop.LowPart =  ti.ftStopKernel.dwLowDateTime;
    nStop.HighPart = ti.ftStopKernel.dwHighDateTime;
    nNull.LowPart =  ti.ftNullKernel.dwLowDateTime;
    nNull.HighPart = ti.ftNullKernel.dwHighDateTime;
    nElapsedKernel.QuadPart = nStop.QuadPart - nNull.QuadPart;
    nElapsedKernel.QuadPart -= nNull.QuadPart - nStart.QuadPart;
    ti.fKernel = double(nElapsedKernel.QuadPart)/10000000.0;
}

void
ML_NullLoopTiming(ML_TimingInfo &ti)
{
    FILETIME ftJunk;

    QueryPerformanceCounter(&ti.nNullReal);
    GetProcessTimes(GetCurrentProcess(), &ftJunk, &ftJunk, &ti.ftNullKernel, &ti.ftNullUser);
}

void
ML_ReportTiming(ML_TimingInfo ti)
{
    double fUser, fKernel, fReal;

    fUser   = ti.fUser   + 0.00005;     // resolution is 1/64th second
    fKernel = ti.fKernel + 0.00005;
    fReal   = ti.fReal   + 0.00005;     // resolution is very fine

    printf("%9.4f  user %0.4f kernel  %9.4f real\n", fUser, fKernel, fReal);
    //printf("%d.%0.4d user  %d.%0.4d kernel  %d.%0.4d real\n",
    //            int(fUser), int(10000.0*(fUser - floor(fUser))),
    //            int(fKernel), int(10000.0*(fKernel - floor(fKernel))),
    //            int(fReal), int (10000.0*(fReal - floor(fReal))));
}

void
ML_ReportEventsPerSecond(ML_TimingInfo ti, unsigned nEvents, char *szUnit)
{
    double fPerSecond;

    fPerSecond = double(nEvents)/ti.fReal + 0.0005;
    printf("%12.3f %s/sec\n", fPerSecond, szUnit);
    //printf("%8d.%3d %s/sec\n", int(fPerSecond),
    //    int(1000.0*(fPerSecond - floor(fPerSecond))), szUnit);
}
