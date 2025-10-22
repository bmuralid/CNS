#ifndef CNS_INDEX_H_
#define CNS_INDEX_H_

#if (AMREX_SPACEDIM == 2)

#define URHO  0
#define UMX   1
#define UMY   2
#define UEINT 3
#define UY1 4
#define UY2 5
#define NUM_STATE 6

#define NCONS 6
#define NEQNS 4
#define NSP  2

#define QRHO   0
#define QU     1
#define QV     2
#define QPRES  3
#define QTEMP  4
#define QY1    5
#define QY2    6
#define NPRIM  7

#else

#define URHO  0
#define UMX   1
#define UMY   2
#define UMZ   3
#define UEINT 4
#define UY1   5
#define UY2   6
#define NUM_STATE 7

#define NCONS 7
#define NEQNS 4
#define NSP 2

#define QRHO   0
#define QU     1
#define QV     2
#define QW     3
#define QPRES  4
#define QTEMP  5
#define QY1  6
#define QY2  7
#define NPRIM  8

#endif

#define CETA   0
#define CXI    1
#define CLAM   2
#define NCOEF  3

#endif
