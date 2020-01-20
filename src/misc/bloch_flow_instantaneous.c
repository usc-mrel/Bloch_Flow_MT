/*=======================================================================*/
/*                                                                       */
/*  SOURCE_FILE:    BLOCH_FLOW_INSTANTANEOUS.C                           */
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
    size_t N;  /* length of RFpulses */
    size_t k;  /* loop variable */

    /* pointers to real and imaginary data of a MATLAB array (mxArray) */
    double *rRFpulses = NULL;
    double *iRFpulses = NULL;
    double *TRs       = NULL;
    double *TEs       = NULL;
    double *ASL_bolus = NULL;
    double *m         = NULL;

    double T1;     /* longitudinal relaxation time       [sec]                    */
    double T2;     /* transverse relaxation time         [sec]                    */
    double df;     /* off-resonance                      [Hz]                     */
    double F;      /* blood flow                         [mL/g/min]               */
    double lambda; /* tissue/blood partition coefficient [mL blood/g tissue]      */
    double m0;     /* equilibrium magnetization          [magnetization/g tissue] */

    /* Initialize matrices */
    double Rx_alpha[9] = {0.0}; /* Rotation by alpha about the x-axis */
    double Rz_pphi[9]  = {0.0}; /* Rotation by +phi  about the z-axis */
    double Rz_nphi[9]  = {0.0}; /* Rotation by -phi  about the z-axis */
    double E1[9]       = {0.0}; /* Relaxation operator 1              */
    double E2[9]       = {0.0}; /* Relaxation operator 2              */
    double P1[9]       = {0.0}; /* Free precession 1                  */
    double P2[9]       = {0.0}; /* Free precession 2                  */
    double C1[9]       = {0.0}; /* Clearance operator 1               */
    double C2[9]       = {0.0}; /* Clearance operator 2               */

    /* Initialize vectors */
    double ma[3] = {0.0};
    double mb[3] = {0.0};
    double mc[3] = {0.0};
    double md[3] = {0.0};
    double mt[3] = {0.0};

    /* local scalar variables */
    double TR;
    double TE;
    double rRFpulse;
    double iRFpulse;
    double alpha;
    double phi;
    double bolus;
    double e12;
    double e11;
    double e22;
    double e21;
    double arg;
    double c;
    double s;
    double c1;
    double c2;

    /* variables for BLAS */
    char *chn     = "N";
    ptrdiff_t n   = 3;
    double done   = 1.0;
    double dzero  = 0.0;
    ptrdiff_t one = 1;

    /*-------------------------------------------------------------------*/
    /* Check for proper number of arguments                              */
    /*-------------------------------------------------------------------*/
    if (nrhs != 10)
    {
        mexErrMsgTxt("10 inputs are required.");
    }
    if (nlhs != 1)
    {
        mexErrMsgTxt("1 output is required.");
    }

    /*-------------------------------------------------------------------*/
    /* Check the data type of inputs                                     */
    /*-------------------------------------------------------------------*/
    if (!mxIsDouble(prhs[0]) || !mxIsComplex(prhs[0]))
    {
        mexErrMsgTxt("The 1st input, RFpulses, must be a double-precision complex array of size N.");
    }
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]))
    {
        mexErrMsgTxt("The 2nd input, TRs [sec], must be a double-precision real array of size N.");
    }
    if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]))
    {
        mexErrMsgTxt("The 3rd input, TEs [sec], must be a double-precision real array of size N.");
    }
    if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]))
    {
        mexErrMsgTxt("The 4th input, T1 [sec], must be a scalar, double-precision real array (1 x 1).");
    }
    if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]))
    {
        mexErrMsgTxt("The 5th input, T2 [sec], must be a scalar, double-precision real array (1 x 1).");
    }
    if (!mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]))
    {
        mexErrMsgTxt("The 6th input, df [Hz], must be a scalar, double-precision real array (1 x 1).");
    }
    if (!mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]))
    {
        mexErrMsgTxt("The 7th input, F [mL/g/min], must be a scalar, double-precision real array (1 x 1).");
    }
    if (!mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]))
    {
        mexErrMsgTxt("The 8th input, lambda, must be a scalar, double-precision real array (1 x 1).");
    }
    if (!mxIsDouble(prhs[8]) || mxIsComplex(prhs[8]))
    {
        mexErrMsgTxt("The 9th input, m0, must be a scalar, double-precision real array (1 x 1).");
    }
    if (!mxIsDouble(prhs[9]) || mxIsComplex(prhs[9]))
    {
        mexErrMsgTxt("The 10th input, ASL_bolus, must be a double-precision real array of size N.");
    }

    /*-------------------------------------------------------------------*/
    /* Set size, dimension related parameters                            */
    /*-------------------------------------------------------------------*/
    N = mxGetNumberOfElements(prhs[0]);

    /*-------------------------------------------------------------------*/
    /* Get pointers to the real and imaginary data of MATLAB arrays      */
    /*-------------------------------------------------------------------*/
    rRFpulses = (double *) mxGetData(prhs[0]);
    iRFpulses = (double *) mxGetImagData(prhs[0]);
    TRs       = (double *) mxGetData(prhs[1]);
    TEs       = (double *) mxGetData(prhs[2]);
    ASL_bolus = (double *) mxGetData(prhs[9]);

    /*-------------------------------------------------------------------*/
    /* Dereference a value directly from a pointer to mxArray            */
    /*-------------------------------------------------------------------*/
    T1     = *((double *) mxGetData(prhs[3]));
    T2     = *((double *) mxGetData(prhs[4]));
    df     = *((double *) mxGetData(prhs[5]));
    F      = *((double *) mxGetData(prhs[6]));
    lambda = *((double *) mxGetData(prhs[7]));
    m0     = *((double *) mxGetData(prhs[8]));

    /*-------------------------------------------------------------------*/
    /* Create output MATLAB arrays for the return argument               */
    /*-------------------------------------------------------------------*/
    plhs[0] = mxCreateDoubleMatrix((size_t)3, N, mxREAL);
    m = (double *) mxGetData(plhs[0]);

    /* Set matrices */
    Rx_alpha[0] = 1.0;
    Rz_pphi[8]  = 1.0;
    Rz_nphi[8]  = 1.0;
    P1[8]       = 1.0;
    P2[8]       = 1.0;
    ma[2]       = m0;

    for (k = 0; k < N; k++)
    {
        /*===============================================================*/
        /*    RF(1)         RF(2)         RF(N-1)        RF(N)           */
        /*       |             |             |             |             */
        /*    ma | mb mc    md |             |             |             */
        /*      \|/   |       \|             |             |             */
        /*  -----+----x--------+----x--------+----x--------+----x--------*/ 
        /*       t(1)          t(2)          t(N-1)        t(N)          */
        /*       <--->         <--->         <--->         <--->         */
        /*       TE(1)         TE(2)         TE(N-1)       TE(N)         */
        /*       <------------><------------><------------><------------>*/
        /*            TR(1)         TR(2)        TR(N-1)        TR(N)    */
        /*===============================================================*/
        TR       = TRs[k];
        TE       = TEs[k];
        rRFpulse = rRFpulses[k];
        iRFpulse = iRFpulses[k];
        alpha    = sqrt(rRFpulse * rRFpulse + iRFpulse * iRFpulse); /* RF pulse angle [rad] */
        phi      = atan2(iRFpulse, rRFpulse);                       /* RF pulse angle [rad] */

        /*---------------------------------------------------------------*/
        /* Rotation by alpha about the x-axis (Rotations are left-handed)*/
        /* Rx_alpha = [     1            0             0     ;           */
        /*                  0       cos(alpha)    sin(alpha) ;           */
        /*                  0      -sin(alpha)    cos(alpha)];           */
        /*---------------------------------------------------------------*/
        c = cos(alpha);
        s = sin(alpha);
        Rx_alpha[4] =  c; Rx_alpha[7] = s;
        Rx_alpha[5] = -s; Rx_alpha[8] = c;

        /*---------------------------------------------------------------*/
        /* Rotation by phi about the z-axis (Rotations are left-handed)  */
        /* Rz_pphi = [  cos(phi)   sin(phi)        0     ;               */
        /*             -sin(phi)   cos(phi)        0     ;               */
        /*                 0          0            1    ];               */
        /*---------------------------------------------------------------*/
        c = cos(phi);
        s = sin(phi);
        Rz_pphi[0] =  c; Rz_pphi[3] = s;
        Rz_pphi[1] = -s; Rz_pphi[4] = c;

        /*---------------------------------------------------------------*/
        /* Rotation by -phi about the z-axis (Rotations are left-handed) */
        /* Rz_nphi = [ cos(-phi)   sin(-phi)        0     ;              */
        /*            -sin(-phi)   cos(-phi)        0     ;              */
        /*                0            0            1    ];              */
        /*---------------------------------------------------------------*/
        Rz_nphi[0] = c; Rz_nphi[3] = -s;
        Rz_nphi[1] = s; Rz_nphi[4] =  c;

        /*---------------------------------------------------------------*/
        /* Relaxation operator                                           */
        /* Describes T1 and T2 relaxation over a time interval tau       */
        /* E1 = [exp(-tau/T2)        0             0      ;              */
        /*            0        exp(-tau/T2)        0      ;              */
        /*            0              0       exp(-tau/T1)];              */
        /*---------------------------------------------------------------*/
        e12 = exp(-TE / T2);
        e11 = exp(-TE / T1);
        E1[0] = e12; E1[4] = e12; E1[8] = e11;

        e22 = exp(-(TR - TE) / T2);
        e21 = exp(-(TR - TE) / T1);
        E2[0] = e22; E2[4] = e22; E2[8] = e21;

        /*---------------------------------------------------------------*/
        /* Free precession over a period tau about the z-axis            */
        /* [2*pi rad/cycle] * [Hz] * [sec] => [rad]                      */
        /* P = [ cos(2*pi*df*tau)  sin(2*pi*df*tau)       0 ;            */
        /*      -sin(2*pi*df*tau)  cos(2*pi*df*tau)       0 ;            */
        /*               0                 0              1];            */
        /*---------------------------------------------------------------*/
        arg = 2 * M_PI * df * TE;
        c = cos(arg);
        s = sin(arg);
        P1[0] =  c; P1[3] = s;
        P1[1] = -s; P1[4] = c;

        arg = 2 * M_PI * df * (TR - TE);
        c = cos(arg);
        s = sin(arg);
        P2[0] =  c; P2[3] = s;
        P2[1] = -s; P2[4] = c;

        /*---------------------------------------------------------------*/
        /* Clearance operator                                            */
        /* Describes the clearance of transverse and longitudinal        */
        /* magnetization by venous flow over a time interval "tau"       */
        /* F: [mL/g/min] * [min/60 sec] => [mL/g/sec]                    */
        /*---------------------------------------------------------------*/
        c1 = exp(-F / 60 * TE / lambda);
        C1[0] = c1; C1[4] = c1; C1[8] = c1;

        c2 = exp(-F / 60 * (TR - TE) / lambda);
        C2[0] = c2; C2[4] = c2; C2[8] = c2;

        /*---------------------------------------------------------------*/
        /* Labeled flow                                                  */
        /* Calculate the amount of labeled inflow during the kth interval*/ 
        /* using the hard pulse approximation                            */
        /*---------------------------------------------------------------*/
        ma[2] = ma[2] + ASL_bolus[k] * TR;

        /*---------------------------------------------------------------*/
        /* mb = Rz_nphi * Rx_alpha * Rz_pphi * ma;                       */
        /* y := alpha*A*x + beta*y                                       */
        /*---------------------------------------------------------------*/
        dgemv(chn, &n, &n, &done, Rz_pphi , &n, ma, &one, &dzero, mb, &one);
        dgemv(chn, &n, &n, &done, Rx_alpha, &n, mb, &one, &dzero, mt, &one);
        dgemv(chn, &n, &n, &done, Rz_nphi , &n, mt, &one, &dzero, mb, &one);

        /*---------------------------------------------------------------*/
        /* mc = C1 * P1 * E1 * mb + D1;                                  */
        /*---------------------------------------------------------------*/
        dgemv(chn, &n, &n, &done, E1, &n, mb, &one, &dzero, mc, &one);
        dgemv(chn, &n, &n, &done, P1, &n, mc, &one, &dzero, mt, &one);
        dgemv(chn, &n, &n, &done, C1, &n, mt, &one, &dzero, mc, &one);
        mc[2] = mc[2] + (1.0 - c1 * e11) * m0;

        /*---------------------------------------------------------------*/
        /* md = C2 * P2 * E2 * mc + D2;                                  */
        /*---------------------------------------------------------------*/
        dgemv(chn, &n, &n, &done, E2, &n, mc, &one, &dzero, md, &one);
        dgemv(chn, &n, &n, &done, P2, &n, md, &one, &dzero, mt, &one);
        dgemv(chn, &n, &n, &done, C2, &n, mt, &one, &dzero, md, &one);
        md[2] = md[2] + (1.0 - c2 * e21) * m0;

        /*---------------------------------------------------------------*/
        /* m(:,k) = mc;                                                  */
        /*---------------------------------------------------------------*/
        m[3*k  ] = mc[0];
        m[3*k+1] = mc[1];
        m[3*k+2] = mc[2];

        /*---------------------------------------------------------------*/
        /* ma = md;                                                      */
        /*---------------------------------------------------------------*/
        dswap(&n, ma, &one, md, &one);
    }
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
