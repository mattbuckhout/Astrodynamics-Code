
% SPECIAL PERTURBATIONS - FORCE MODELS AND EQUATIONS OF MOTION
% ---------------------------------------------------------------------------------------------------------------------
% Function takes position vector of satellite in a geocentric orbit at a given time and develops the equations of 
% motion for the satellite, accounting for 2-body acceleration as well as perturbations due to Earth's non-spherical 
% gravitational potential, the Sun, and the Moon. The equations of mottion follow Cowell's Formulation, where the 
% three 2nd order differential equations are re-formed into six first order diff. eqs.
%
% Central body perturbations are developed using the EGM2008 model, containing spherical harmonic coefficients that 
% together with Associated Legendre Polynomials define Earth's Non-Spherical Gravitational Potential. Accelerations are 
% developed using partial derivatives of the potential function. The resulting accelerations account for perturbing 
% effects due to zonal, tesseral, and sectorial harmonics for an L x L field, where L is specified by the user and 
% corresponds to the maximum spherical harmonic degree/order of coefficients from the EGM2008 model.
%
% Third body perturbations are developed using the JPL DE430 model for planetary ephemerides. The DE430 model yields 
% highly accurate position vectors of solar system bodies resulting from numerical simulations of the system. The 
% position vectors of the Sun and Moon are used to compute the perturbing accelerations of each on the satellite.
%
% 
% Author: Matthew Buckhout
% Updated: 08/18/2020 
%
% Inputs:
%
%     - [t]                Array of Time Values, +from Initial Date/Time                                    [sec]
%     - [STATE]            Matrix of Position and Velocity Vectors                                          [km][km/s]
%     - [MJDi]             Initial Modified Julian Date (UT1)                                               [days]
%     - [EOPdata]          Earth Orientation Parameter Table                                                 -
%     - [XTAB]             Selected X coefficients for series terms in IAU2006 Precession-Nutation Model     -
%     - [YTAB]             Selected Y coefficients for series terms in IAU2006 Precession-Nutation Model     -
%     - [sTAB]             Selected s coefficients for series terms in IAU2006 Precession-Nutation Model     -
%     - [EGM2008coef]      Un-Normalaized Spherical Harmonic Coefficients                                    -
%     - [DE430coef]        JPL DE430 Ephemerides Model Coefficients                                          -
%     - [L]                Degree x Order, Field for Spherical Harmonics                                     -
%                             (Aspherical Gravitational Potential)                      
%     - [perturbations]    Vector specifying the perturbing forces to include in simulation: (ex: [1 2])     -
%                                1 - Earth Non-spherical Gravitational Potential
%                                2 - Third Body Accelerations (Sun/Moon)
%
% Outputs:
%
%     - [d_STATE]          Values of Equations of Motion at each value of t                                 [km][km/s]
%
% Functions:
%
%     - EOP
%     - IAU2006_rotations  
%     - function3 
%     - EGM2008
%     - third_body
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 491-564)
% ---------------------------------------------------------------------------------------------------------------------

function [d_STATE] = perturbed_eom(t,STATE,MJDi,EOPdata,XTAB,YTAB,sTAB,EGM2008coef,DE430coef,L,perturbations)

   %Constants
   RE = 6378.1363; %[km] Earth Mean Equatorial Radius
   mu = 3.986004415e5; %[km^3/s^2] Earth Gravitational Parameter
   i = 1; j = 2; k = 3;

   %State Vector (GCRF)
   R_GCRF = [STATE(1); STATE(2); STATE(3)]; %[km]
   V_GCRF = [STATE(4); STATE(5); STATE(6)]; %[km/s]
   r_GCRF = norm(R_GCRF); %[km]
   
   %Time Step - UT1
   MJD_step = MJDi + (t/86400);
   [xp,yp,dUT1,dAT] = EOP(fix(MJD_step),EOPdata);
     
   %State Vector (ITRF)
   length = 'short';
   [ITRF2PEF,PEF2IRE,IRE2GCRF,GCRF2IRE,IRE2PEF,PEF2ITRF] = IAU2006_rotations(MJD_step,dUT1,dAT,xp,yp,length,XTAB,YTAB,sTAB); %Reduction for Rotation Matrices
   R_ITRF = PEF2ITRF*IRE2PEF*GCRF2IRE*R_GCRF; %[km]
   r_ITRF = norm(R_ITRF); %[km] 
      
   % --------------------------------  P E R T U R B I N G  A C C E L E R A T I O N S  --------------------------------
   
   %Aspherical Gravitational Potential
   if (ismember(1,perturbations) == 1)
      ap_nonspher = EGM2008(R_ITRF,L,EGM2008coef); %(ITRF)
   else
      ap_nonspher = [0;0;0];
   end
      
   %Third Bodies
   if (ismember(2,perturbations) == 1)
      ap_thirdbodies = third_body(MJD_step,R_GCRF,DE430coef,EOPdata,[3 10 11]); %(GCRF)
   else
      ap_thirdbodies = [0;0;0];
   end
      
   % ----------------------------------  D I F F E R E N T I A L  E Q U A T I O N S  ----------------------------------
   
   ap_ITRF = ap_nonspher; %Combined Perturbing Accelerations (ITRF) 
   ap_GCRF = ap_thirdbodies; %Combined Perterbing Accelerations (GCRF)
   ap_GCRFtot = (IRE2GCRF*PEF2IRE*ITRF2PEF*ap_ITRF) + ap_GCRF; %[km/s^2] All Perturbing Acceleration to (GCRF)

   %Cowell's Formulation - System of First Order Diff. Eq. (GCRF)
   d_STATE(1) = V_GCRF(i);
   d_STATE(2) = V_GCRF(j);
   d_STATE(3) = V_GCRF(k);
   d_STATE(4) = ((-mu*R_GCRF(i))/(r_GCRF^3)) + ap_GCRFtot(i);
   d_STATE(5) = ((-mu*R_GCRF(j))/(r_GCRF^3)) + ap_GCRFtot(j);
   d_STATE(6) = ((-mu*R_GCRF(k))/(r_GCRF^3)) + ap_GCRFtot(k);
  
end
