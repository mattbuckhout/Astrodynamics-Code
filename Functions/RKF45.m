
% 4TH/5TH ORDER RUNGE-KUTTA-FEHLBERG INTEGRATOR - VARIABLE STEP
% ---------------------------------------------------------------------------------------------------------------------
% Function integrates non-linear differential equations of motion for a satellite from an initial state over a given 
% time span using a variable step size. The Runge-Kutta-Fehlberg algorithm performs 4th and 5th order integrations of 
% the equations of motion at each time step and compares the results in order to reduce error. The step size is 
% dynamically changed to maintain good agreement between the 4th and 5th order integrations. In this implementation, 
% the difference between the two integrations of different orders is used to scale the time step, after which the
% integration is repeated with the new time step.
%
% This implementation passes the necessary data sets from the outer "propagator_SP" function to the "perturbed_eom" 
% function, which contains the force models, accelerations, and EOM. 
%
% Variable step size integrator is more appropriate for eccentric orbits (e > 0.02). Accuracy maintained at periapsis 
% where velocity is high and efficiency improved at apoapsis where velocity is low.
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

function [t,y] = RKF45(functionhandle,tspan,y0,h,MJDi,T,EOPdata,XTAB,YTAB,sTAB,EGM2008coef,DE430coef,L,perturbations)

   %Defining Times
   t0 = tspan(1); %[sec] Initial Time
   tf = tspan(2); %[sec] Final Time
   hmin = 60; %[sec] Mininum Time Step
   
   y(1,:) = y0; %Initial State Vectors
   t(1) = t0; %Initial Time Value
      
   %Count/Print Orbits
   orb = 0;
   cond = 0;
   totalorb = (tf/T);
   if ((fix(totalorb) - totalorb) == 0)
      cond = 1;
   end

   % Runge-Kutta-Fehlberg 4/5
   % ------------------------------------------------------------------------------------------------------------------
   
   %Initial 4th/5th Order Runge-Kutta Values
   i = 1;
   k1 = h*functionhandle(t(i),y(i,:),MJDi,EOPdata,XTAB,YTAB,sTAB,EGM2008coef,DE430coef,L,perturbations);
   k2 = h*functionhandle(t(i) + (1/4)*h,y(i,:) + (1/4)*k1,MJDi,EOPdata,XTAB,YTAB,sTAB,EGM2008coef,DE430coef,L,perturbations);
   k3 = h*functionhandle(t(i) + (3/8)*h,y(i,:) + (3/32)*k1 + (9/32)*k2,MJDi,EOPdata,XTAB,YTAB,sTAB,EGM2008coef,DE430coef,L,perturbations);
   k4 = h*functionhandle(t(i) + (12/13)*h,y(i,:) + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3,MJDi,EOPdata,XTAB,YTAB,sTAB,EGM2008coef,DE430coef,L,perturbations);
   k5 = h*functionhandle(t(i) + h,y(i,:) + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4,MJDi,EOPdata,XTAB,YTAB,sTAB,EGM2008coef,DE430coef,L,perturbations);
   k6 = h*functionhandle(t(i) + (1/2)*h,y(i,:) - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5,MJDi,EOPdata,XTAB,YTAB,sTAB,EGM2008coef,DE430coef,L,perturbations);           
   
   while (t(i) < tf) %Loop for each value of time

      d45 = norm((1/360)*k1 - (128/4275)*k3 - (2197/75240)*k4 + (1/50)*k5 + (2/55)*k6); %Diffs 4th and 5th Order
      s = 0.8408*((((1e-8)*h)/d45)^(1/4)); %Scaling parameter for time step
   
      if ((t(i) + h) > tf)
         h = tf - t(i); %Changing time step to match final time
      else
         h = s*h; %Changing time step to match changing orbit speed
      end
      
      k1 = h*functionhandle(t(i),y(i,:),MJDi,EOPdata,XTAB,YTAB,sTAB,EGM2008coef,DE430coef,L,perturbations);
      k2 = h*functionhandle(t(i) + (1/4)*h,y(i,:) + (1/4)*k1,MJDi,EOPdata,XTAB,YTAB,sTAB,EGM2008coef,DE430coef,L,perturbations);
      k3 = h*functionhandle(t(i) + (3/8)*h,y(i,:) + (3/32)*k1 + (9/32)*k2,MJDi,EOPdata,XTAB,YTAB,sTAB,EGM2008coef,DE430coef,L,perturbations);
      k4 = h*functionhandle(t(i) + (12/13)*h,y(i,:) + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3,MJDi,EOPdata,XTAB,YTAB,sTAB,EGM2008coef,DE430coef,L,perturbations);
      k5 = h*functionhandle(t(i) + h,y(i,:) + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4,MJDi,EOPdata,XTAB,YTAB,sTAB,EGM2008coef,DE430coef,L,perturbations);
      k6 = h*functionhandle(t(i) + (1/2)*h,y(i,:) - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5,MJDi,EOPdata,XTAB,YTAB,sTAB,EGM2008coef,DE430coef,L,perturbations);

      y(i+1,:) = y(i,:) + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5; %Integrated State
      t(i+1) = t(i) + h; %Next time step
      i = i+1; 
         
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













