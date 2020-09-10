#include <stdio.h>
#include <math.h>
#include "mex.h"

/* PLANETARY EPHEMERIDES - DE430 MODEL (JPL) 
 --------------------------------------------------------------------------------------------------------------------
 Returns position and velocity vectors for any of 11 bodies in the solar system in the International Celestial 
 Reference Frame (ICRF) relative to Solar System Barycenter. Moon position vector given for Geocentric Celestial 
 Reference Frame (GCRF). Ephemerides based on JPL DE430 model.
     https://ssd.jpl.nasa.gov/?planet_eph_export

 Time system used for DE430 model is Barycentric Dynamical Time (TDB). Time intervals in coefficient files are in 
 Julian Days referenced to Terrestrial Dynamical Time (TDT/TT), so the input Modified Julian Date must also be 
 referenced to TT. This implementation uses the DE430t version of coefficient files with TT - TDB included.

 Author: Matthew Buckhout
 Updated: 08/16/2020 

 Inputs:

     - [MJD_TT]        Modified Julian Date, Expressed in TT or TDT                                  [days]
     - [DE430coef]     Matrix of coefficients from JPL DE430 model over period containing MJD_TT      -
     - [bodies]        Vector of identifiers (1-11) specifying which bodies to generate               -
                          ephemerides for 

 Outputs:

     - [ephemerides]   Position & Velocity vectors of selected bodies at MJD_TT                      [km][km/s]
                          [Rx1 Ry1 Rz1 Vx1 Vy1 Vz1; ...
                           Rx1 Ry2 Rz2 Vx2 Vy2 Vz2]                            

 Identifiers:

     1. Mercury
     2. Venus
     3. Earth-Moon Barycenter
     4. Mars
     5. Jupiter
     6. Saturn
     7. Uranus
     8. Neptune
     9. Pluto
    10. Moon (GCRF)
    11. Sun
    12. -
    13. -
    14. -
    15. TT - TDB

 Line to call .mex file from Matlab:

      [ephemerides] = planet_ephemerides_full_c(MJD_TT,transpose(DE430coef),bodies);

 References:
     - JPL Planetary and Lunar Ephemerides, Information and DE430 coefficient files
           https://ssd.jpl.nasa.gov/?planet_eph_export
     - JPL Ephemerides Project (Lukas Bystricky)
           https://people.sc.fsu.edu/~lb13f/projects/space_environment/planet_positions.php
 --------------------------------------------------------------------------------------------------------------------*/

/*========================================== G A T E W A Y  R O U T I N E ============================================*/

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    int inputnum = nrhs;
    int outputnum = nlhs;
    const mwSize* dims_coef;
    const mwSize* dims_bodies;
    double* MJDTT, * coef, * bodies;
    double* ephemerides;

    MJDTT = mxGetPr(prhs[0]); /*Copying pointers to inputs*/
    coef = mxGetPr(prhs[1]);
    bodies = mxGetPr(prhs[2]);

    dims_coef = mxGetDimensions(prhs[1]); /*Getting dimensions of input arrays*/
    dims_bodies = mxGetDimensions(prhs[2]);

                                   /*rows*/  /*cols*/
    /*plhs[0] = mxCreateDoubleMatrix(dims_bodies[1], 3, mxREAL); /*Creating Multi-Dim Output Arrays*/
    plhs[0] = mxCreateDoubleMatrix(dims_bodies[1], 6, mxREAL); /*Creating Multi-Dim Output Arrays*/

    ephemerides = mxGetPr(plhs[0]); /*Copying pointers to outputs*/

    /*============================================ S U B R O U T I N E ===============================================*/

    int startcoef[15] = { 3, 171, 231, 309, 342, 366, 387, 405, 423, 441, 753, 819, 819, 939, 939 };
    int numcoef[15] = { 14, 10, 13, 11, 8, 7, 6, 6, 6, 13, 11, 0, 10, 0, 11 };
    int numsets[15] = { 4, 2, 2, 1, 1, 1, 1, 1, 1, 8, 2, 0, 4, 0, 4 };
      
    int i, j, b, p, n, N, sectionlength, sectionstart, sectionend; 
    int numsections, sectionnum, section, start, end, subinterval_ind;
    int MJDint1a, MJDint2a, MJDsubint1, MJDsubint2, MJDsubint1a, MJDsubint2a, subintspan, startx, starty, startz;
    double MJDspan, tbar, MJDint1, MJDint2, xcoef[14], ycoef[14], zcoef[14], Tn[14], TNM1, TNP1, TN; 
    double dTn[14], dTNM1, dTNP1, dTN, R_ICRFx = 0, R_ICRFy = 0, R_ICRFz = 0, V_ICRFx = 0, V_ICRFy = 0, V_ICRFz = 0;

    sectionlength = *(coef + 1); /*Number of coefficients in each section*/
    numsections = (3 * dims_coef[1]) / (sectionlength + 5); /*Total number of sections in data file*/

    for (i = 0; i < numsections-1; i++) /*Loop to select time interval*/
    {
        sectionnum = i + 1; /*Section number*/
        sectionstart = 3 + (sectionlength + 5) * (sectionnum - 1); /*Index for start of section*/
        sectionend = 3 + ((sectionlength + 5) * (sectionnum)) - 6; /*Index for end of section*/
        MJDint1 = *(coef + sectionstart) - 2400000.5; /*Time Interval for section*/
        MJDint2 = *(coef + sectionstart + 1) - 2400000.5; 

        if (*MJDTT >= MJDint1 && *MJDTT < MJDint2)
        {
            section = sectionnum;
            start = sectionstart;  
            end = sectionend;
            MJDint1a = (int)MJDint1;
            MJDint2a = (int)MJDint2;
        }
    }

    for (b = 0; b < dims_bodies[1]; b++) /*Loop for planetary body*/
    {
        p = *(bodies + b) - 1; /*Identifier for planetary body of interest*/

        subintspan = (MJDint2a - MJDint1a) / numsets[p];
        for (i = 0; i < numsets[p] + 1; i++) /*Loop to select time subinterval*/
        {
            MJDsubint1 = MJDint1a + (i * subintspan);
            MJDsubint2 = MJDint1a + (i + 1) * subintspan;

            if (*MJDTT >= (double)(MJDsubint1) && *MJDTT < (double)(MJDsubint2))
            {
                MJDsubint1a = MJDsubint1;
                MJDsubint2a = MJDsubint2;
                subinterval_ind = i;
            }
        }

        tbar = (*MJDTT - (((double)MJDsubint1a + (double)MJDsubint2a) / 2)) * (2 / ((double)MJDsubint2a - (double)MJDsubint1a));

        startx = start + startcoef[p] + (subinterval_ind * 3 * numcoef[p]) - 1;
        starty = startx + numcoef[p];
        startz = starty + numcoef[p];
        
        j = 0;
        for (i = startx; i < starty; i++) /*Loop for x component*/
        {
            xcoef[j] = *(coef + i);
            j = j + 1;
        }

        j = 0;
        for (i = starty; i < startz; i++) /*Loop for y component*/
        {
            ycoef[j] = *(coef + i);
            j = j + 1;
        }

        j = 0;
        for (i = startz; i < startz + numcoef[p]; i++) /*Loop for z component*/
        {
            zcoef[j] = *(coef + i);
            j = j + 1;
        }
                
        N = numcoef[p]; /*Order for Chebyshev Polynomials*/
        Tn[0] = (double)1; /*Chebyshev Polynomial Order 0*/
        Tn[1] = tbar; /*Chebyshev Polynomial Order 1*/
        dTn[0] = (double)0; /*First Derivative Chebyshev Polynomial Order 0*/
        dTn[1] = (double)1; /*First Derivative Chebyshev Polynomial Order 1*/

        for (n = 0; n < N-1; n++) /*Loop for Chebyshev Polynomials of the First Kind*/
        {
            if (n == 0)
            {
                TNM1 = Tn[0];
                TN = Tn[1];
                dTNM1 = dTn[0];
                dTN = dTn[1];
            }

            TNP1 = (2 * tbar * TN) - TNM1; /*Recurrence Relation: Chebyshev Polynomials*/
            dTNP1 = (2 * tbar * dTN) - dTNM1 + (2 * TN); /*Recurrence Relation: First Derivative of Chebyshev Polynomials*/
            
            TNM1 = TN;
            TN = TNP1;
            Tn[n + 2] = TNP1;
                        
            dTNM1 = dTN;
            dTN = dTNP1;
            dTn[n + 2] = dTNP1;
        }

        R_ICRFx = 0; /*Resetting State Vector components*/
        R_ICRFy = 0;
        R_ICRFz = 0;
        
        V_ICRFx = 0;
        V_ICRFy = 0;
        V_ICRFz = 0;
        
        for (i = 0; i < numcoef[p]; i++) /*Loop to compute position/velocity vector components*/
        {
            R_ICRFx = R_ICRFx + xcoef[i] * Tn[i];
            R_ICRFy = R_ICRFy + ycoef[i] * Tn[i];
            R_ICRFz = R_ICRFz + zcoef[i] * Tn[i];
           
            V_ICRFx = V_ICRFx + xcoef[i] * dTn[i];
            V_ICRFy = V_ICRFy + ycoef[i] * dTn[i];
            V_ICRFz = V_ICRFz + zcoef[i] * dTn[i];
           
        }
        
        MJDspan = (double)(MJDint2a - MJDint1a);
        V_ICRFx = (2 * numsets[p] * V_ICRFx)/(MJDspan * 86400); /*Scaling Velocity by number of subintervals, number of sets, days to seconds*/
        V_ICRFy = (2 * numsets[p] * V_ICRFy)/(MJDspan * 86400);
        V_ICRFz = (2 * numsets[p] * V_ICRFz)/(MJDspan * 86400);
        
        ephemerides[0 + b] = R_ICRFx; /*Loading output vectors*/
        ephemerides[dims_bodies[1] + b] = R_ICRFy;
        ephemerides[(dims_bodies[1] * 2) + b] = R_ICRFz;
        
        ephemerides[(dims_bodies[1] * 3) + b] = V_ICRFx;
        ephemerides[(dims_bodies[1] * 4) + b] = V_ICRFy;
        ephemerides[(dims_bodies[1] * 5) + b] = V_ICRFz;

        
    }

}