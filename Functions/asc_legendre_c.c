#include <stdio.h>
#include <math.h>
#include "mex.h"

/* ASSOCIATED LEGENDRE POLYNOMIALS
--------------------------------------------------------------------------------------------------------------------- 
Function computes Associated Legendre Polynomials P(l,m) from degree 0 to L and all orders for each degree.
     https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
 
Author: Matthew Buckhout
Updated: 08/06/2020  

Inputs:           

    - [rows]    Number polynomials to compute (all degree-order pairs up to order L)     
    - [L]       Maximum order to compute polynomials for                                 
    - [ARG]     Argument of Legendre Polynomials                                        
    - [P00]     Initial Value, Associated Legendre Polynomial, Degree 0 Order 0          
    - [P10]     Initial Value, Associated Legendre Polynomial, Degree 1 Order 0          
    - [P11]     Initial Value, Associated Legendre Polynomial, Degree 1 Order 1          

Outputs:

    - [P]       Array of computed polynomials, P(l,m)                                               
    - [P1]      Array of computed polynomials, P(l,m+1)                                  

Line to call .mex file from Matlab: 

   B = asc_legendre_c(rows,L,ARG,P00,P10,P11); 
   
Lines for extracting P and P1 from output matrix:

   P = B(:,1);
   P1 = B(:,2);
---------------------------------------------------------------------------------------------------------------------*/

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /*===================== I N P U T   H A N D L I N G   F R A M E W O R K =====================*/
    
    /*1. [nlhs] = Number of Outputs (size of output array)
      2. [plhs] = Array of pointers to Outputs
      3. [nrhs] = Number of Inputs (size of input array)
      4. [prhs] = Array of pointers to Inputs*/

    int inputnum = nrhs;
    int outputnum = nlhs; 
    int i, l, m, step, rows, size, L;
    double *inputpointers[inputnum], *outputpointers;
    double ARG, P00, P10, P11, PLM, PLM1M, PLL, PLP1M;

    /*Getting Input Pointers*/
    for (i = 0; i < inputnum; i++) {
        inputpointers[i] = mxGetPr(prhs[i]); 
    }

    /*Values of Inputs stored in local variables*/
    rows = *inputpointers[0];
    L = *inputpointers[1];
    ARG = *inputpointers[2];
    P00 = *inputpointers[3];
    P10 = *inputpointers[4];
    P11 = *inputpointers[5];

    /*Initializing Arrays for Legendre Polynomials*/
    size = rows + L + 4;
    double VALS1[size];
    double VALS2[size];
    double P[size];
    double P1[size];

    /*Creating Output Array*/
    plhs[0] = mxCreateDoubleMatrix(rows, 2, mxREAL);

    /*Getting Output Pointers*/
    outputpointers = mxGetPr(plhs[0]);
    
    /*=========================== L E G E N D R E   A L G O R I T H M ===========================*/
    step = 0;
    for (l = 1; l < (L + 2); l++) {

        for (m = 0; m < (l + 1); m++) {

            if ((l == 1) && (m == 0)) 
            {
                PLM = P10;
                PLM1M = P00;
            }
            else if ((l == 1) && (m == 1)) 
            {
                PLM = P11;
                PLM1M = 0;
            }
            else if (l == m) 
            {
                PLL = VALS1[step - (l + 1)];
                PLM = -((2 * (l - 1)) + 1) * sqrt(1 - (pow(ARG,2))) * PLL;
                PLM1M = 0;
            }
            else 
            {
                PLM = VALS2[step - l];
                PLM1M = VALS1[step - l];
            }
            PLP1M = (((2 * l) + 1) * ARG * PLM - (l + m) * PLM1M) / (l - m + 1);
            VALS1[step] = PLM; 
            VALS2[step] = PLP1M;
            step = step + 1;
        }

    }

    for (i = 0; i < (size-(l+1)); i++) {

        P[i] = VALS1[i+2];

    }
    for (i = 0; i < (size-l); i++) {

        P1[i] = VALS1[i+3];

    }

    /*==================== O U T P U T   H A N D L I N G   F R A M E W O R K ====================*/

    /*Loading generated outputs into output array*/
    for (i = 0; i < rows; i++) {

        outputpointers[i] = P[i];

    }
    for (i = 0; i < rows; i++) {

        outputpointers[rows + i] = P1[i];

    }

}




