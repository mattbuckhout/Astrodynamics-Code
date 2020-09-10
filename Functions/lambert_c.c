#include <stdio.h>
#include <math.h>
#include "mex.h"
 
/* LAMBERT'S PROBLEM - UNIVERSAL VARIABLE SOLUTION
 ----------------------------------------------------------------------------------------------------------------------
 General solution to Lambert's problem. Computes the velocity vectors of a spacecraft on a trajectory between 2 
 position vectors separated by a specified amount of time. Solution uses f and g functions with universal variables to 
 handle all orbit types. 

 Author: Matthew Buckhout
 Updated: 09/06/2020 

 Inputs:

     - [mu]         Central Body Gravitational Parameter            [km^3/s^2]
     - [R0]         Position Vector (t0)                            [km]
     - [R]          Position Vector (t)                             [km]
     - [dt]         Time between R0 and R (t-t0)                    [sec]
     - [tm]         Transfer Method                                  -
                       tm = +1 (Short Way)
                       tm = -1 (Long Way)

 Outputs:

     - [V0]         Velocity on computed trajectory at R0 (t0)      [km/s]
     - [V]          Velocity on computed trajectory at R (t)        [km/s]

 Lines to call .mex file from Matlab:

      [VS] = lambert_c(mu,R0,R,dt,tm); 
      V0 = VS(:,1); 
      V = VS(:,2);  

 References:
     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 448-454, 461-466)
 ---------------------------------------------------------------------------------------------------------------------*/
 
/*================================================ C 2  A N D  C 3 ===================================================*/
 
double C2C3(double PSI, double *C2, double *C3) /*Function to compute factoriels*/
{

    if (PSI > 1e-8) /*Computing common terms in the Universal Variable formulation*/
    {
       *C2 = (1 - cos(sqrt(PSI)))/PSI;
       *C3 = (sqrt(PSI) - sin(sqrt(PSI)))/sqrt(pow(PSI,3));
    }
    else if (PSI < -1e-8)
    {
       *C2 = (1 - cosh(sqrt(-PSI)))/PSI;
       *C3 = (sinh(sqrt(-PSI)) - sqrt(-PSI))/sqrt(pow(-PSI,3));
    }
    else
    {
       *C2 = (double) 1/2;
       *C3 = (double) 1/6;
    }

}
 
/*========================================== G A T E W A Y  R O U T I N E ============================================*/

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    int inputnum = nrhs;
    int outputnum = nlhs;
    const mwSize* dims_coef;
    const mwSize* dims_bodies;
    double* mu, * R0, * R, * dt, * tm;
    double* vels;

    mu = mxGetPr(prhs[0]); /*Copying pointers to inputs*/
    R0 = mxGetPr(prhs[1]);
    R  = mxGetPr(prhs[2]);
    dt = mxGetPr(prhs[3]);
    tm = mxGetPr(prhs[4]);

                            /*rows*/  /*cols*/
    plhs[0] = mxCreateDoubleMatrix(3, 2, mxREAL); /*Creating Multi-Dim Output Arrays*/

    vels = mxGetPr(plhs[0]); /*Copying pointers to outputs*/

    /*============================================ S U B R O U T I N E ===============================================*/

    int iterations, iterations2;
    double R0x, R0y, R0z, Rx, Ry, Rz, r0, r, V0x, V0y, V0z, Vx, Vy, Vz; 
    double cos_dtheta, A, C2, C3, PSILOW, PSIUP, PSI1, PSI2, dt2, y, CHI, f, g, g_dot;
    
    R0x = *(R0);     /*Initial Position Vector Components*/ 
    R0y = *(R0 + 1);
    R0z = *(R0 + 2); 
    
    Rx = *(R);       /*Final Position Vector Components*/ 
    Ry = *(R + 1);
    Rz = *(R + 2); 
    
    r0 = sqrt(pow(R0x,2) + pow(R0y,2) + pow(R0z,2)); /*Magnitude, Initial Position*/
    r  = sqrt(pow(Rx,2) + pow(Ry,2) + pow(Rz,2)); /*Magnitude, Final Position*/
    
    cos_dtheta = ((R0x*Rx) + (R0y*Ry) + (R0z*Rz))/(r0*r);
    A = (*tm)*sqrt(r*r0*(1 + cos_dtheta));
    
    PSI1 = 0; /*Initial values for iteration*/
    C2 = (double) 1/2;
    C3 = (double) 1/6;
    
    PSIUP = 4*pow(M_PI,2); /*Bounds for iteration*/
    PSILOW = -4*M_PI;
    
    /*BISECTION METHOD*/
    dt2 = (*dt) + 1;
    iterations = 0;
    iterations2 = 0;
    while (abs(dt2 - *dt) > 1e-6)
    {
     
       y = r0 + r + ((A*(PSI1*C3 - 1))/sqrt(C2)); /*Compute initial value of y*/
       
       if ((A > 0) && (y < 0))
       {  
          
          while (y < 0)
          {   
             
             PSILOW = PSILOW + 0.01*M_PI; /*Adjusting lower bound to get positive value for y*/
             PSI1 = (PSIUP + PSILOW)/2;
             
             C2C3(PSI1, &C2, &C3); /*Call C2 and C3 function*/
             y = r0 + r + ((A*(PSI1*C3 - 1))/sqrt(C2)); /*Recalculate y to check condition*/
             
             iterations2++;
             if (iterations2 > 1000)
             {
                break;
             }
             
          }
          
       }
       
       CHI = sqrt(y/C2);
       dt2 = (pow(CHI,3)*C3 + A*sqrt(y))/sqrt(*mu);
       
       if (abs(dt2) <= (*dt)) /*Walking in bounds*/
       {   
          PSILOW = PSI1;
       }
       else
       { 
          PSIUP = PSI1;   
       }
     
       PSI2 = (PSIUP + PSILOW)/2;
     
       C2C3(PSI2, &C2, &C3); /*Call C2 and C3 function*/
     
       PSI1 = PSI2; /*Iterating PSI*/
       iterations++;  
     
       if (iterations > 1000)
       {
          break;
       }
     
    }
             
    f = 1 - (y/r0); /*f and g Functions*/
    g = A*sqrt(y/(*mu));
    g_dot = 1 - (y/r);

    V0x = (Rx - f*(R0x))/g; /*[km/s] Initial Velocity Vector*/
    V0y = (Ry - f*(R0y))/g;
    V0z = (Rz - f*(R0z))/g;
    Vx = (g_dot*Rx - R0x)/g; /*[km/s] Final Velocity Vector*/
    Vy = (g_dot*Ry - R0y)/g;
    Vz = (g_dot*Rz - R0z)/g;
    
    vels[0] = V0x; /*Loading output vectors*/
    vels[1] = V0y;
    vels[2] = V0z;
    vels[3] = Vx;
    vels[4] = Vy;
    vels[5] = Vz;

}