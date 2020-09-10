#include <stdio.h>
#include <math.h>
#include "mex.h"

/*SERIES B TERM IN THIRD-BODY PERTURBATION EQUATIONS
 -------------------------------------------------------------------------------------------------
 Function computes series term in alternate expression of a perturbing third-body acceleration in 
 the satellite equations of motion. The process of computing the value of B was implemented in C 
 and compiled as a .mex executable to increase speed when called from Matlab. Function built 
 specifically for use in "third_body.m" function.

 Author: Matthew Buckhout
 Updated: 08/18/2020 

 Inputs:

     - [ARG]        Cosine of Angle between 3rd body and Satellite position vectors     -
     - [TOL]        Tolerance for convergence in B                                      -
     - [resat]      Magnitude of Position Vector, Geocenter to Satellite               [km]
     - [re3]        Magnitude of Position Vector, Geocenter to 3rd body                [km]

 Outputs:

     - [B]          Value of B                                                          -

 Line to call .mex file from Matlab:

     B = series_B_c(ARG,TOL,resat,re3);

 References:
     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 543)

/*======================================= F A C T O R I A L  R O U T I N E ===========================================*/

int factorial(int VAL) /*Function to compute factoriels*/
{
    int n, FACT = 1;

    for (n = 1; n < VAL; n++)
    {
        FACT = FACT * (n + 1);
    }

    return FACT;
}

/*========================================== G A T E W A Y  R O U T I N E ============================================*/

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    int i, j;
    double* ARG,* TOL,* resat,* re3,* B;
    double B1, B2, err, sum;

    ARG = mxGetPr(prhs[0]); /*Argument of Conventional Legendre Polynomials*/
    TOL = mxGetPr(prhs[1]); /*Tolerance for series convergence*/
    resat = mxGetPr(prhs[2]); /*Position Vector, Geocenter to Satellite*/
    re3 = mxGetPr(prhs[3]); /*Position Vector, Geocenter to Third Body*/

    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL); /*Creating Output Arrays*/

    B = mxGetPr(plhs[0]); /*Copying pointers to outputs*/

/*============================================= S U B R O U T I N E =================================================*/
    
    B1 = 0; j = 1; err = 1;
    while (err > *TOL)
    {
        sum = 0;
        for (i = 0; i <= j; i++)
        {
            /*Rodrigues' Formula*/
            sum = sum + (double)pow(((factorial(j)) / (factorial(i) * factorial(j - i))),2) * pow((*ARG - 1), (double)(j - i)) * pow((*ARG + 1), (double)i);
        }

        B2 = B1 + ((1 / (pow(2, (double)j))) * (pow((*resat / *re3), (double)j)) * sum);
        err = abs(B2 - B1);
        B1 = B2;
        j = j + 1;
    
        if (j > 10000)
        {
            B2 = 0;
            break;
        }

    }

    B[0] = B2;

}

