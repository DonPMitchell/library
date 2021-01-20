//
//  Some system utilties for timing and debugging
//  D.P. Mitchell  2018/09/11.
//
#pragma once
//
//  Assert a condition
//
#if defined(ML_TESTING) || defined(DEBUG) || defined (_DEBUG)
    #define ML_Assert(__expr)                                   \
      do {                                                      \
        if (!(__expr)) {                                        \
            printf("ML_Assert failed in %s\n",__FILE__);        \
            printf(" line %d: (%s)\n", __LINE__, #__expr);      \
             __debugbreak();                                    \
        }                                                       \
      } while (0)
#else
    #define ML_Assert(__expr)       ((void) 0)
#endif
//
//  Timing utilities
//
struct ML_TimingInfo {
    LARGE_INTEGER   nStartReal, nNullReal, nStopReal;           // real-time ticks
    FILETIME        ftStartKernel, ftNullKernel, ftStopKernel;  // CPU times (system and user mode)
    FILETIME        ftStartUser, ftNullUser, ftStopUser;
    double          fUser, fKernel, fReal;                      // seconds
};

extern void ML_StartTiming(ML_TimingInfo &ti);
extern void ML_NullLoopTiming(ML_TimingInfo &ti);
extern void ML_StopTiming(ML_TimingInfo &ti);
extern void ML_ReportTiming(ML_TimingInfo ti);
extern void ML_ReportEventsPerSecond(ML_TimingInfo ti, unsigned nEvents, char *szUnit);
