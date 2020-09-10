#include <stdio.h>
#include <math.h>
#include "mex.h"
 
/* HELIOCENTRIC ORBIT PROPAGATION - KEPLER'S PROBLEM, UNIVERSAL VARIABLE APPROACH
 ----------------------------------------------------------------------------------------------------------------------
 Function takes Position & Velocity vectors for a satellite in any inertial frame in an orbit at initial time (t0) and 
 determines Position & Velocity Vectors for the satellite at a future time (t) using the Universal Variable Formulation 
 of Kepler's Problem [CHI,PSI]. Function is valid for all conic orbits and special cases of Equatorial Elliptic, 
 Inclined Circular, and Equatorial Circular orbits. 2-body dynamics defines orbital motion.

 Author: Matthew Buckhout
 Updated: 08/11/2020 

 Inputs:

     - [mu]         Gravitational Parameter of Central Body        [km^3/s^2]
     - [R0]         Position Vector at initial time                [km]
     - [V0]         Velocity Vector at initial time                [km/s]
     - [MJD0]       Initial Time, Modified Julian Date (UT1)       [days]
     - [MJD]        Final Time, Modified Julian Date (UT1)         [days]

 Outputs:

     - [R]          Position Vector at final time                   [km]
     - [V]          Velocity Vector at final time                   [km/s]

 Lines to call .mex file from Matlab:
      
      [RV] = propagator_c(mu,R0,V0,MJD0,MJD);  
      R = RV(:,1);
      V = RV(:,2);

 References:
     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 66-71, 90-102)
 ---------------------------------------------------------------------------------------------------------------------*/
  
/*====================================================  C2 AND C3  ===================================================*/
 
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
 
/*================================================  GATEWAY ROUTINE  =================================================*/

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    int inputnum = nrhs;
    int outputnum = nlhs;
    const mwSize* dims_coef;
    const mwSize* dims_bodies;
    double* mu, * R0, * V0, * MJD0, * MJD;
    double* state_out;

    mu = mxGetPr(prhs[0]); /*Copying pointers to inputs*/
    R0 = mxGetPr(prhs[1]);
    V0 = mxGetPr(prhs[2]);
    MJD0 = mxGetPr(prhs[3]);
    MJD = mxGetPr(prhs[4]);

                            /*rows*/  /*cols*/
    plhs[0] = mxCreateDoubleMatrix(3, 2, mxREAL); /*Creating Multi-Dim Output Arrays*/

    state_out = mxGetPr(plhs[0]); /*Copying pointers to outputs*/

    /*==============================================  SUBROUTINE  ====================================================*/

    int iterations;
    double R0x, R0y, R0z, V0x, V0y, V0z, Rx, Ry, Rz, Vx, Vy, Vz, Ex, Ey, Ez, Hx, Hy, Hz, r0, v0, e, h; 
    double dt, dotRV, eps, alpha, a, p, s, w, CHI0, CHI1, CHI2, PSI, r, f, g, fdot, gdot, C2, C3;
    double check, diff, testval, tolerance;
    
    R0x = *(R0); /*Initial Position Vector Components*/ 
    R0y = *(R0 + 1);
    R0z = *(R0 + 2); 
    
    V0x = *(V0); /*Initial Velocity Vector Components*/ 
    V0y = *(V0 + 1);
    V0z = *(V0 + 2); 
    
    dt = (*MJD - *MJD0)*86400; /*[sec] Change in Time*/
    if (dt == 0)
    {
       
       Rx = R0x;
       Ry = R0y;
       Rz = R0z;
       Vx = V0x;
       Vy = V0y;
       Vz = V0z;
           
    }
    else
    {   
       
       r0 = sqrt(pow(R0x,2) + pow(R0y,2) + pow(R0z,2)); /*Magnitude, Initial Position*/
       v0 = sqrt(pow(V0x,2) + pow(V0y,2) + pow(V0z,2)); /*Magnitude, Initial Velocity*/
    
       dotRV = (R0x*V0x) + (R0y*V0y) + (R0z*V0z); /*Dot Product of R0 and V0*/
       
       Ex = (1/(*mu))*((pow(v0,2) - (*mu/r0))*R0x - (dotRV)*V0x); /*Eccentricity Vector Components*/
       Ey = (1/(*mu))*((pow(v0,2) - (*mu/r0))*R0y - (dotRV)*V0y);
       Ez = (1/(*mu))*((pow(v0,2) - (*mu/r0))*R0z - (dotRV)*V0z);
       
       e = sqrt(pow(Ex,2) + pow(Ey,2) + pow(Ez,2)); /*Eccentricity*/
       
       eps = (pow(v0,2)/2) - (*mu/r0); /*[km^2/s^2]  Specific Mechanical Energy*/
       alpha = -(pow(v0,2)/(*mu)) + (2/r0); /*[1/km]  Parameter (alpha = 1/a)*/
       
       if (alpha != 1) /*Selecting equation for Initial CHI*/
       {
        
          if (e < 1) /*Elliptical and Circular*/
          {
           
             CHI0 = sqrt(*mu)*dt*alpha;
             
          }
          else if (e == 1) /*Parabolic*/
          {
           
             Hx = R0y*V0z - V0y*R0z; /*Angular Momentum Vector Components*/
             Hy = V0x*R0z - R0x*V0z;
             Hz = R0x*V0y - V0x*R0y;
             
             h = sqrt(pow(Hx,2) + pow(Hy,2) + pow(Hz,2)); /*Specific Angular Momentum*/
             
             p = pow(h,2)/(*mu); /*Parabolic Semi-Parameter*/
             s = 2*(acos(3*sqrt(*mu/pow(p,3))*dt)/asin(3*sqrt(*mu/pow(p,3))*dt)); /*Intermediate Angle 1*/
             w = atan(pow(tan(s),(1/3))); /*Intermediate Angle 2*/
             CHI0 = sqrt(p)*2*(cos(2*w)/sin(2*w));
             
          }
          else if (e > 1) /*Hyperbolic*/
          {
           
             a = 1/alpha; /*Semi-Major Axis*/
             CHI0 = (dt/abs(dt))*sqrt(-a)*log((-2*(*mu)*alpha*dt)/(dotRV + (dt/abs(dt))*sqrt(-(*mu)*a)*(1 - (r0*alpha))));
             
          }
       
          /*NEWTON RAPHSON - 1ST ITERATION ---------------------------------------------------------------------------*/ 
          CHI1 = CHI0;
          PSI = pow(CHI0,2)*alpha;
          
          C2C3(PSI, &C2, &C3); /*Call C2 and C3 function*/
          
          r = pow(CHI1,2)*C2 + (dotRV/sqrt(*mu))*CHI1*(1 - PSI*C3) + r0*(1 - PSI*C2);
          CHI2 = CHI1 + ((sqrt(*mu)*dt - pow(CHI1,3)*C3 - (dotRV/sqrt(*mu))*pow(CHI1,2)*C2 - r0*CHI1*(1-PSI*C3))/r);
                    
          /*NEWTON RAPHSON - N+1 ITERATIONS --------------------------------------------------------------------------*/
          diff = 1;
          iterations = 1;
          while (diff > 1e-10)
          {
              
             CHI1 = CHI2; /*Iterating CHI*/      
             PSI = pow(CHI1,2)*alpha; /*Iterating PSI*/
             
             C2C3(PSI, &C2, &C3); /*Call C2 and C3 function*/
             
             r = pow(CHI1,2)*C2 + (dotRV/sqrt(*mu))*CHI1*(1 - PSI*C3) + r0*(1 - PSI*C2);
             CHI2 = CHI1 + ((sqrt(*mu)*dt - pow(CHI1,3)*C3 - (dotRV/sqrt(*mu))*pow(CHI1,2)*C2 - r0*CHI1*(1-PSI*C3))/r);
             
             iterations++;
             if (iterations > 1000)
             {
                break;
             }
             
             diff = abs(CHI2 - CHI1);

          }
          
          f = 1 - (pow(CHI2,2)/r0)*C2; /*f and g Functions*/
          g = dt - (pow(CHI2,3)/sqrt(*mu))*C3;
          fdot = (sqrt(*mu)/(r*r0))*CHI2*(PSI*C3 - 1);
          gdot = 1 - (pow(CHI2,2)/r)*C2;
      
          Rx = f*R0x + g*V0x; /*Final Position Vector Components*/
          Ry = f*R0y + g*V0y;
          Rz = f*R0z + g*V0z;
                    
          Vx = fdot*R0x + gdot*V0x; /*Final Velocity Vector Components*/   
          Vy = fdot*R0y + gdot*V0y;
          Vz = fdot*R0z + gdot*V0z;
 
          check = f*gdot - fdot*g; /*Convergence check for f and g functions*/
          if (abs(check - 1) > 1e-8)
          {
             Rx = Ry = Rz = 0; /*f and g functions failed to converge*/
             Vx = Vy = Vz = 0;
          }         
          
       }  
      
    }
   
    state_out[0] = Rx; /*Loading output vectors*/
    state_out[1] = Ry;
    state_out[2] = Rz;
    state_out[3] = Vx;
    state_out[4] = Vy;
    state_out[5] = Vz;

}
   