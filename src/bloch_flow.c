/*=======================================================================*/
/*                                                                       */
/*  SOURCE_FILE:    BLOCH_FLOW.C                                         */
/*                                                                       */
/*  Copyright 2018: Namgyun Lee                                          */
/*  Author: Namgyun Lee                                                  */
/*  USC (University of Southern California)                              */
/*  namgyunl@usc.edu, ggang56@gmail.com (preferred)                      */
/*=======================================================================*/
/*=======================================================================*/
/*  I N C L U D E S                                                      */
/*=======================================================================*/
/*-----------------------------------------------------------------------*/
/*  Common                                                               */
/*-----------------------------------------------------------------------*/
#define _USE_MATH_DEFINES // To use MATH Constants
#include <time.h>
#include <math.h>

/*-----------------------------------------------------------------------*/
/*  MATLAB                                                               */
/*-----------------------------------------------------------------------*/
/* Macro to define the correct function prototype */
#ifdef _WIN32
#define FORTRAN_COMPLEX_FUNCTIONS_RETURN_VOID 1
#endif

#include "mex.h"
#include "blas.h"

/*=======================================================================*/
/*  G L O B A L   R E F E R E N C E S                                    */
/*=======================================================================*/
/*=======================================================================*/
/*  G L O B A L   D E F I N I T I O N S                                  */
/*=======================================================================*/
/*=======================================================================*/
/*  L O C A L   S Y M B O L   D E F I N I T I O N S                      */
/*=======================================================================*/
/*=======================================================================*/
/*  L O C A L   D A T A   D E F I N I T I O N S                          */
/*=======================================================================*/
/*=======================================================================*/
/*  L O C A L   F U N C T I O N   P R O T O T Y P E S                    */
/*=======================================================================*/
/*=======================================================================*/
/*  F U N C T I O N   P R O T O T Y P E S                                */
/*=======================================================================*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    size_t M;                  /* number of time intervals            */
    size_t N;                  /* number of off-resonance frequencies */
    size_t P;                  /* number of spatial positions         */
    size_t nd;                 /* number of spatial dimensions        */
    size_t idx1;               /* loop variable                       */
    size_t idx2;               /* loop variable                       */
    size_t idx3;               /* loop variable                       */
    size_t ndim;
    size_t dims[3];
    size_t blocksize_M;
    size_t blocksize_P;
    size_t k;

    /* pointers to real and imaginary data of a MATLAB array (mxArray)   */
    double *rb1       = NULL;  /* Rf pulse in G. Can be complex          [G]    (M x 1)           */
    double *ib1       = NULL;
    double *gr        = NULL;  /* 1,2, or 3-dimensional gradient         [G/cm] (M x 1,2,or 3)    */
    double *taus      = NULL;  /* time duration of each b1 and gr point  [sec]  (M x 1)           */
    double *zetas     = NULL;  /* measurement time of each interval      [sec]  (M x 1)           */
    double *df        = NULL;; /* array of off-resonance frequencies     [Hz]   (N x 1)           */
    double *dp        = NULL;; /* array of spatial positions             [cm]   (P x 1,2,or 3)    */
                               /* Width should match width of gr.                                 */
    double *mx0       = NULL;  /* array of starting x magnetization             (P x 1)           */
    double *my0       = NULL;  /* array of starting y magnetization             (P x 1)           */
    double *mz0       = NULL;  /* array of starting z magnetization             (P x 1)           */
    double *ASL_bolus = NULL;  /* ASL bolus signal                       [magnetization/g tissue/sec] (N x 1) */

    double T1;                 /* longitudinal relaxation time           [sec]                    */
    double T2;                 /* transverse relaxation time             [sec]                    */
    int    mode;               /* Bitmask mode:                                                      */
                               /* Bit 0: 0-Simulate from start or M0, 1-Steady State (not supported) */
                               /* Bit 1: 1-Record m at time points. 0-just end time.                 */
    double F;                  /* blood flow                             [mL/g/min]               */
    double lambda;             /* tissue/blood partition coefficient     [mL blood/g tissue]      */
    double m0;                 /* equilibrium magnetization              [magnetization/g tissue] */

    /* Outputs */
    double *m_zeta = NULL;
    double *mx     = NULL;     /* array of the resulting x magnetization  [magnetization/g tissue] (M x P x N) or (P x N) */
    double *my     = NULL;     /* array of the resulting y magnetization  [magnetization/g tissue] (M x P x N) or (P x N) */
    double *mz     = NULL;     /* array of the resulting z magnetization  [magnetization/g tissue] (M x P x N) or (P x N) */

    /* Gradient & position component vectors */
    double *gx = NULL;
    double *gy = NULL;
    double *gz = NULL;
    double *dx = NULL;
    double *dy = NULL;
    double *dz = NULL;

    /* Initialize 3x3 matrices */
    double R[9]  = {0.0}; /* Rotation matrix from axis and angle */
    double E1[9] = {0.0}; /* Relaxation operator 1               */
    double E2[9] = {0.0}; /* Relaxation operator 2               */
    double C1[9] = {0.0}; /* Clearance operator 1                */
    double C2[9] = {0.0}; /* Clearance operator 2                */

    /* Initialize 3 x 1 magnetization vectors */
    double ma[3] = {0.0};
    double mb[3] = {0.0};
    double mc[3] = {0.0};
    double md[3] = {0.0};
    double mt[3] = {0.0};

    /* local scalar variables */
    double T1app;
    double pseudo_M0;;
    double df_;
    double dx_;
    double dy_;
    double dz_;    
    double tau;
    double zeta;
    double b1x;
    double b1y;
    double Bz;
    double B_abs;
    double theta;
    double ux;
    double uy;
    double uz;
    double c;
    double s;
    double one_minus_c;
    double R11;
    double R22;
    double R33;
    double R12a;
    double R12b;
    double R13a;
    double R13b;
    double R23a;
    double R23b;
    double e12;
    double e11;
    double e22;
    double e21;
    double arg;
    double c1;
    double c2;

    /* gyromagnetic ratio */
    double gam = 4257.746778 * 2.0 * M_PI * 1e-2; /* [rad/sec/uT] */

    /* variables for BLAS */
    char *chn     = "N";
    ptrdiff_t n   = 3;
    double done   = 1.0;
    double dzero  = 0.0;
    ptrdiff_t one = 1;

    /*-------------------------------------------------------------------*/
    /* Check for proper number of arguments                              */
    /*-------------------------------------------------------------------*/
    if (nrhs != 16)
    {
        mexErrMsgTxt("16 inputs are required.");
    }
    if (nlhs != 3)
    {
        mexErrMsgTxt("3 outputs are required.");
    }

    /*-------------------------------------------------------------------*/
    /* Set size, dimension related parameters                            */
    /*-------------------------------------------------------------------*/
    M = mxGetNumberOfElements(prhs[0]); /* length for b1, gr, taus, zetas, ASL_bolus */
    N = mxGetNumberOfElements(prhs[6]); /* length for df */
    P = mxGetNumberOfElements(prhs[7]); /* length for dp, mx0, my0, mz0 */

    /*-------------------------------------------------------------------*/
    /* Check the data type of inputs                                     */
    /*-------------------------------------------------------------------*/
    if (!mxIsDouble(prhs[0]))
    {
        mexErrMsgTxt("The 1st input, b1 [G], must be a double-precision real or complex array (M x 1).");
    }
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || (mxGetM(prhs[1]) != M))
    {
        mexErrMsgTxt("The 2nd input, gr [G/cm], must be a double-precision real array (M x 1,2,or 3).");
    }
    if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || (mxGetM(prhs[2]) != M))
    {
        mexErrMsgTxt("The 3rd input, taus [sec], must be a double-precision real array (M x 1).");
    }
    if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || (mxGetM(prhs[3]) != M))
    {
        mexErrMsgTxt("The 4th input, zetas [sec], must be a double-precision real array (M x 1).");
    }
    if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || (mxGetM(prhs[4]) != 1) || (mxGetN(prhs[4]) != 1))
    {
        mexErrMsgTxt("The 5th input, T1 [sec], must be a real scalar (1 x 1).");
    }
    if (!mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) || (mxGetM(prhs[5]) != 1) || (mxGetN(prhs[5]) != 1))
    {
        mexErrMsgTxt("The 6th input, T2 [sec], must be a real scalar (1 x 1).");
    }
    if (!mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]))
    {
        mexErrMsgTxt("The 7th input, df [Hz], must be a double-precision real array (N x 1).");
    }
    if (!mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]))
    {
        mexErrMsgTxt("The 8th input, dp [cm], must be a double-precision real array (P x 1,2,or 3).");
    }
    if (!mxIsDouble(prhs[8]) || mxIsComplex(prhs[8]) || (mxGetM(prhs[8]) != 1) || (mxGetN(prhs[8]) != 1))
    {
        mexErrMsgTxt("The 9th input, mode, must be a real scalar (1 x 1).");
    }
    if (!mxIsDouble(prhs[9]) || mxIsComplex(prhs[9]) || (mxGetM(prhs[9]) != P) || (mxGetN(prhs[9]) != N))
    {
        mexErrMsgTxt("The 10th input, mx0, must be a double-precision real array (P x N).");
    }
    if (!mxIsDouble(prhs[10]) || mxIsComplex(prhs[10]) || (mxGetM(prhs[10]) != P) || (mxGetN(prhs[10]) != N))
    {
        mexErrMsgTxt("The 11th input, my0, must be a double-precision real array (P x N).");
    }
    if (!mxIsDouble(prhs[11]) || mxIsComplex(prhs[11]) || (mxGetM(prhs[11]) != P) || (mxGetN(prhs[11]) != N))
    {
        mexErrMsgTxt("The 12th input, mz0, must be a double-precision real array (P x N).");
    }
    if (!mxIsDouble(prhs[12]) || mxIsComplex(prhs[12]) || (mxGetM(prhs[12]) != 1) || (mxGetN(prhs[12]) != 1))
    {
        mexErrMsgTxt("The 13th input, F [mL/g/min], must be a real scalar (1 x 1).");
    }
    if (!mxIsDouble(prhs[13]) || mxIsComplex(prhs[13]) || (mxGetM(prhs[13]) != 1) || (mxGetN(prhs[13]) != 1))
    {
        mexErrMsgTxt("The 14th input, lambda, must be a real scalar (1 x 1).");
    }
    if (!mxIsDouble(prhs[14]) || mxIsComplex(prhs[14]) || (mxGetM(prhs[14]) != 1) || (mxGetN(prhs[14]) != 1))
    {
        mexErrMsgTxt("The 15th input, m0, must be a real scalar (1 x 1).");
    }
    if (!mxIsDouble(prhs[15]) || mxIsComplex(prhs[15]) || (mxGetM(prhs[15]) != M))
    {
        mexErrMsgTxt("The 16th input, ASL_bolus, must be a double-precision real array (M x 1).");
    }

    /*-------------------------------------------------------------------*/
    /* Set misc. parameters                                              */
    /*-------------------------------------------------------------------*/
    nd = mxGetN(prhs[1]);
    blocksize_M = M * sizeof(double);
    blocksize_P = P * sizeof(double);

    /*-------------------------------------------------------------------*/
    /* Get pointers to the real and imaginary data of MATLAB arrays      */
    /*-------------------------------------------------------------------*/
    rb1       = (double *) mxGetData(prhs[0]);
    gr        = (double *) mxGetData(prhs[1]);
    taus      = (double *) mxGetData(prhs[2]);
    zetas     = (double *) mxGetData(prhs[3]);
    df        = (double *) mxGetData(prhs[6]);
    dp        = (double *) mxGetData(prhs[7]);
    mx0       = (double *) mxGetData(prhs[9]);
    my0       = (double *) mxGetData(prhs[10]);
    mz0       = (double *) mxGetData(prhs[11]);
    ASL_bolus = (double *) mxGetData(prhs[15]);

    if (mxIsComplex(prhs[0]))
    {
        ib1 = (double *) mxGetImagData(prhs[0]);
    } else
    {
        ib1 = (double *) mxMalloc(blocksize_M);
        memset(ib1, '\0', blocksize_M);
    }

    /*-------------------------------------------------------------------*/
    /* Dereference a value directly from a pointer to mxArray            */
    /*-------------------------------------------------------------------*/
    T1     = *((double *) mxGetData(prhs[4]));
    T2     = *((double *) mxGetData(prhs[5]));
    mode   = (int) *((double *) mxGetData(prhs[8]));
    F      = *((double *) mxGetData(prhs[12]));
    lambda = *((double *) mxGetData(prhs[13]));
    m0     = *((double *) mxGetData(prhs[14]));

    /*-------------------------------------------------------------------*/
    /* Create output MATLAB arrays for the return argument               */
    /*-------------------------------------------------------------------*/
    if (mode == 2) /* Record m at time points */
    {
        ndim    = 3;
        dims[0] = M;
        dims[1] = P;
        dims[2] = N;
    } else if (mode == 0)
    {
        ndim    = 2;
        dims[0] = P;
        dims[1] = N;
    }

    /*collection of mc at all zeta's */
    m_zeta = (double *) mxMalloc(3 * blocksize_M); /* M x 3 */

    plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);

    mx = (double *) mxGetData(plhs[0]);
    my = (double *) mxGetData(plhs[1]);
    mz = (double *) mxGetData(plhs[2]);

    /*-------------------------------------------------------------------*/
    /* Set gradient and position component vectors                       */
    /*-------------------------------------------------------------------*/
    if (nd == 1)
    {
        gx = gr;
        gy = (double *) mxMalloc(blocksize_M);
        gz = (double *) mxMalloc(blocksize_M);

        dx = dp;
        dy = (double *) mxMalloc(blocksize_P);
        dz = (double *) mxMalloc(blocksize_P);

        memset(gy, '\0', blocksize_M);
        memset(gz, '\0', blocksize_M);
        memset(dy, '\0', blocksize_P);
        memset(dz, '\0', blocksize_P);
    } else if (nd == 2)
    {
        gx = gr;
        gy = gr + M;
        gz = (double *) mxMalloc(blocksize_M);

        dx = dp;
        dy = dp + P;
        dz = (double *) mxMalloc(blocksize_P);

        memset(gz, '\0', blocksize_M);
        memset(dz, '\0', blocksize_P);
    } else
    {
        gx = gr;
        gy = gr + M;
        gz = gr + 2 * M;

        dx = dp;
        dy = dp + P;
        dz = dp + 2 * P;
    }

    /*-------------------------------------------------------------------*/
    /* Calculate the apparent T1                                         */
    /*-------------------------------------------------------------------*/
    T1app = 1. / (1.0 / T1 + (F / lambda / 60.0)); /* apparent T1 [sec] */

    /*-------------------------------------------------------------------*/
    /* Perform Bloch-Flow simulation                                     */
    /*-------------------------------------------------------------------*/
    /* loop for off-resonance frequencies [Hz] */
    for (idx3 = 0; idx3 < N; idx3++)
    {
        df_ = df[idx3]; /* current off-resonance frequency [Hz] */

        /* loop for spatial positions [cm] */
        for (idx2 = 0; idx2 < P; idx2++)
        {
            dx_ = dx[idx2];
            dy_ = dy[idx2];
            dz_ = dz[idx2];

            k = idx2 + idx3 * P;
            ma[0] = mx0[k];
            ma[1] = my0[k];
            ma[2] = mz0[k];

            /* loop for the RF pulse waveform */
            for (idx1 = 0; idx1 < M; idx1++)
            {
                /*=======================================================*/
                /* Timing diagram                                        */
                /*-------------------------------------------------------*/
                /*    b1(1)       b1(2)      b1(M-1)    b1(M)            */
                /*      |           |          |          |              */
                /*   ma | mb mc  md |          |          |              */
                /*     \|/   |     \|          |          |              */
                /* -----+----x------+----x-----+----x-----+----x-------  */
                /*      t(1)        t(2)       t(M-1)     t(M)           */
                /*      <----------><---------><---------><----------->  */
                /*         tau(1)      tau(2)    tau(M-1)    tau(M)      */
                /*      <--->       <--->      <--->      <--->          */
                /*      zeta(1)     zeta(2)    zeta(M-1)  zeta(M)        */
                /*=======================================================*/
                tau  = taus[idx1];  /* [sec] */
                zeta = zetas[idx1]; /* [sec] */

                /*-------------------------------------------------------*/
                /* Calculate the rotation matrix using                   */
                /* Rodrigues' rotation formula                           */
                /* Free precession is removed and the Bz component is    */
                /* considered during excitation. This makes a huge       */ 
                /* difference from the approach with free precession!    */
                /*-------------------------------------------------------*/
                b1x  = rb1[idx1]; /* [G] */
                b1y  = ib1[idx1]; /* [G] */
                Bz = gx[idx1] * dx_ + gy[idx1] * dy_ + gz[idx1] * dz_ + 2.0 * M_PI * df_ / (gam * 1e2); /* [G] */
                B_abs = sqrt(b1x * b1x + b1y * b1y + Bz * Bz); /* [G] */
                theta = gam * 1e2 * B_abs * tau; /* [rad/sec/uT] * [1e2uT/G] * [G] => [rad/sec] */

                if (B_abs > 0)
                {
                    ux = b1x / B_abs;
                    uy = b1y / B_abs;
                    uz = Bz  / B_abs;
                } else
                {
                    ux = 0.0;
                    uy = 0.0;
                    uz = 0.0;
                }

                /* Rotations are left-handed */
                c = cos(-theta); 
                s = sin(-theta);
                one_minus_c = 1 - c;

                R11  = c + ux * ux * one_minus_c;
                R22  = c + uy * uy * one_minus_c;
                R33  = c + uz * uz * one_minus_c;
                R12a = ux * uy * one_minus_c;
                R12b = uz * s;
                R13a = ux * uz * one_minus_c;
                R13b = uy * s;
                R23a = uy * uz * one_minus_c;
                R23b = ux * s;

                R[0] = R11;
                R[1] = R12a + R12b;
                R[2] = R13a - R13b;
                R[3] = R12a - R12b;
                R[4] = R22;
                R[5] = R23a + R23b;
                R[6] = R13a + R13b;
                R[7] = R23a - R23b;
                R[8] = R33;

                /*-------------------------------------------------------*/
                /* Relaxation operator                                   */
                /* T1 and T2 relaxation over a time interval tau         */
                /* E1 = [exp(-tau/T2)        0             0      ;      */
                /*            0        exp(-tau/T2)        0      ;      */
                /*            0              0       exp(-tau/T1)];      */
                /*-------------------------------------------------------*/
                e12 = exp(-zeta / T2);
                e11 = exp(-zeta / T1);
                E1[0] = e12; E1[4] = e12; E1[8] = e11;

                e22 = exp(-(tau - zeta) / T2);
                e21 = exp(-(tau - zeta) / T1);
                E2[0] = e22; E2[4] = e22; E2[8] = e21;

                /*-------------------------------------------------------*/
                /* Clearance operator                                    */
                /* Describes the clearance of transverse and longitudinal*/
                /* magnetization by venous flow over a time interval     */ 
                /* "tau". F: [mL/g/min] * [min/60 sec] => [mL/g/sec]     */
                /*-------------------------------------------------------*/
                c1 = exp(-F / 60.0 * zeta / lambda);
                C1[0] = c1; C1[4] = c1; C1[8] = c1;

                c2 = exp(-F / 60.0 * (tau - zeta) / lambda);
                C2[0] = c2; C2[4] = c2; C2[8] = c2;

                /*-------------------------------------------------------*/
                /* mb = R * ma;                                          */
                /* y  := alpha * A * x + beta * y                        */
                /* mb := done * R * ma + dzero * mb                      */
                /*-------------------------------------------------------*/
                dgemv(chn, &n, &n, &done, R, &n, ma, &one, &dzero, mb, &one);

                /*-------------------------------------------------------*/
                /* Labeled flow                                          */
                /* Calculate a new pseudo M0 term: M0 + s(t) * T1app,    */
                /* which is the hypothetical equilibrium magnetization   */
                /* for modeling flow during this short interval          */
                /*-------------------------------------------------------*/
                pseudo_M0 = m0 + ASL_bolus[idx1] * T1app;

                /*-------------------------------------------------------*/
                /* mc = C1 * E1 * mb + D1;                               */
                /* mt := done * E1 * mb + dzero * mt                     */
                /* mc := done * C1 * mt + dzero * mc                     */
                /*-------------------------------------------------------*/
                dgemv(chn, &n, &n, &done, E1, &n, mb, &one, &dzero, mt, &one);
                dgemv(chn, &n, &n, &done, C1, &n, mt, &one, &dzero, mc, &one);
                mc[2] = mc[2] + (1.0 - c1 * e11) * pseudo_M0;

                /*-------------------------------------------------------*/
                /* md = C2 * E2 * mc + D2;                               */
                /* mt := done * E2 * mc + dzero * mt                     */
                /* md := done * C2 * mt + dzero * mc                     */
                /*-------------------------------------------------------*/
                dgemv(chn, &n, &n, &done, E2, &n, mc, &one, &dzero, mt, &one);
                dgemv(chn, &n, &n, &done, C2, &n, mt, &one, &dzero, md, &one);
                md[2] = md[2] + (1.0 - c2 * e21) * pseudo_M0;

                /*-------------------------------------------------------*/
                /* m_zeta(idx1,:) = mc;                                  */
                /*-------------------------------------------------------*/
                m_zeta[idx1        ] = mc[0];
                m_zeta[idx1 + M    ] = mc[1];
                m_zeta[idx1 + 2 * M] = mc[2];

                /*-------------------------------------------------------*/
                /* ma = md;                                              */
                /*-------------------------------------------------------*/
                dswap(&n, ma, &one, md, &one);
            } /* end of b1 */

            if (mode == 2) /* Record m at time points */
            {
                k = idx2 * M + idx3 * M * P;
                dcopy(&(ptrdiff_t)M, m_zeta        , &one, mx + k, &one);
                dcopy(&(ptrdiff_t)M, m_zeta + M    , &one, my + k, &one);
                dcopy(&(ptrdiff_t)M, m_zeta + 2 * M, &one, mz + k, &one);
            } else if (mode == 0) /* just end time */
            {
                k = idx2 + idx3 * P;
                mx[k] = m_zeta[M - 1        ];
                my[k] = m_zeta[M - 1 + M    ];
                mz[k] = m_zeta[M - 1 + 2 * M];
            }
        }  /* end of dp */
    } /* end of df */

    /*-------------------------------------------------------------------*/
    /* Free CPU memory                                                   */
    /*-------------------------------------------------------------------*/
    if (!mxIsComplex(prhs[0]))
    {
        mxFree(ib1);
    }

    if (nd == 1)
    {
        mxFree(gy);
        mxFree(gz);
        mxFree(dy);
        mxFree(dz);
    } else if (nd == 2)
    {
        mxFree(gz);
        mxFree(dz);
    }

    mxFree(m_zeta);
}

/*=======================================================================*/
/*  H I S T O R Y                                                        */
/*=======================================================================*/
/* 2018-10-11 Namgyun Lee                                                */
/*      Introducing source-file into BLOCH_FLOW                          */
/* 2018-10-11 Namgyun Lee                                                */
/*      Finished the initial version                                     */
/* 2019-06-11 Namgyun Lee                                                */
/*      Changed the unit from [sec] to [msec]                            */
/* 2019-08-29 Namgyun Lee                                                */
/*      Changed the unit from [msec] to [sec]                            */
/* 2019-09-29 Namgyun Lee                                                */
/*      Added df and dp loops and make the interface identical to        */
/*      Prof. Brian Hargreaves's bloch.c                                 */
