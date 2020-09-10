
% 4TH ORDER RUNGE-KUTTA INTEGRATOR - FIXED STEP
% ---------------------------------------------------------------------------------------------------------------------
% Classical 4th Order Runge-Kutta method integrates non-linear differential equations of motion (EOM) for a satellite 
% from an initial state over a given time span and with a fixed time step. At each time step, the 4th Order Runge-Kutta 
% algorithm evaluates the EOM (derivatives of position and velocity) at 4 points near the current time step in order to 
% construct the integrals (position and velocity) at the next time step. This implementation passes the necessary data 
% sets from the outer "propagator_SP" function to the "perturbed_eom" function, which contains the force models, 
% accelerations, and EOM.
%
% Fixed step size most efficient for nearly circular orbits (e < 0.02) where velocity does not change appreciably 
% throughout a single revolution.
%
% Author: Matthew Buckhout
% Updated: 08/12/2020 
%
% Inputs:
%
%     - [functionhandle]   Function Name containing Equations of Motion (@perturbed_eom)                     - 
%     - [tspan]            Time Span for integration                                                        [sec]
%     - [y0]               Initial State Vector [Rx0; Ry0; Rz0; Vx0; Vy0; Vz0]                              [km][km/s]
%     - [h]                Initial Time Step                                                                [sec]
%     - [MJDi]             Initial Modified Julian Date (UT1)                                               [days]
%     - [T]                Orbital Period                                                                   [sec]
%     - [EOPdata]          Earth Orientation Parameter Table                                                 -
%     - [XTAB]             Selected X coefficients for series terms in IAU2006 Precession-Nutation Model     -
%     - [YTAB]             Selected Y coefficients for series terms in IAU2006 Precession-Nutation Model     -
%     - [sTAB]             Selected s coefficients for series terms in IAU2006 Precession-Nutation Model     -
%     - [EGM2008coef]      Un-Normalaized Spherical Harmonic Coefficients                                    -
%     - [DE430coef]        JPL DE430 Ephemerides Model Coefficients                                          -
%     - [L]                Degree x Order, Field for Spherical Harmonics                                     -
%                             (Aspherical Gravitational Potential)                      
%     - [perturbations]    Vector specifying the perturbing forces to include in simulation: (ex: [1 2])     -
%                             1 - Earth Non-spherical Gravitational Potential
%                             2 - Third Body Accelerations (Sun/Moon)
% 
% Outputs:
%
%     - [t]                Array of Time Values, +from Initial Date/Time                                    [sec]
%     - [y]                Matrix of Position and Velocity Vectors                                          [km][km/s]
%
% Format for t and y:
%
%       t       y
%       -----   ------------------------------------------
%         t1      R1x    R1y    R1z    V1x    V1y    V1z
%         t2      R2x    R2y    R2z    V2x    V2y    V2z
%         t3      R3x    R3y    R3z    V3x    V3y    V3z
%         t4      R4x    R4y    R4z    V4x    V4y    V4z
%         t5      R5x    R5y    R5z    V5x    V5y    V5z
%
%                               . . .
%
%         tf      Rfx    Rfy    Rfz    Vfx    Vfy    Vfz
%       -----   ------------------------------------------
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 500-511)
% ---------------------------------------------------------------------------------------------------------------------

function [t,y] = RK4(functionhandle,tspan,y0,h,MJDi,T,EOPdata,XTAB,YTAB,sTAB,EGM2008coef,DE430coef,L,perturbations)

   %Defining Times
   t0 = tspan(1); %[sec] Initial Time
   tf = tspan(2); %[sec] Final Time
   t = [t0:h:tf]'; %[sec] Time vector
   
   %Refining Step size
   h = (tf-t0)/numel(t); %[sec] Refined Time Step
   t = [t0:h:tf]'; %[sec] Refined Time vector
   y = zeros(numel(t),6); %Preallocating
   y(1,:) = y0; %Initial State Vectors
   
   %Count/Print Orbits
   orb = 0;
   cond = 0;
   totalorb = (tf/T);
   if ((fix(totalorb) - totalorb) == 0)
      cond = 1;
   end
      
   % 4th Order Runge-Kutta Integrator
   % ------------------------------------------------------------------------------------------------------------------
   for i=1:1:numel(t)-1

      k1 = functionhandle(t(i),y(i,:),MJDi,EOPdata,XTAB,YTAB,sTAB,EGM2008coef,DE430coef,L,perturbations);
      k2 = functionhandle((t(i) + (h/2)),(y(i,:) + ((h/2)*k1)),MJDi,EOPdata,XTAB,YTAB,sTAB,EGM2008coef,DE430coef,L,perturbations);
      k3 = functionhandle((t(i) + (h/2)),(y(i,:) + ((h/2)*k2)),MJDi,EOPdata,XTAB,YTAB,sTAB,EGM2008coef,DE430coef,L,perturbations);
      k4 = functionhandle((t(i) + h),(y(i,:) + (h*k3)),MJDi,EOPdata,XTAB,YTAB,sTAB,EGM2008coef,DE430coef,L,perturbations);
      
      y(i+1,:) = y(i,:) + (h/6)*(k1 + 2*k2 + 2*k3 + k4); 
   % ------------------------------------------------------------------------------------------------------------------
      
      %Counting Revolutions
      if (t(i) > orb*T)
         orb = orb+1;
         if (cond == 1)
            fprintf('Revolution %d/%d\n',orb,totalorb);
         else
            fprintf('Revolution %d/%d\n',orb,fix(totalorb)+1);
         end
               
      end
      
   end

end
