
% ORBIT PROPAGATION - KEPLER'S PROBLEM, UNIVERSAL VARIABLE APPROACH
% -------------------------------------------------------------------------------------------------
% Function takes Position & Velocity vectors for a satellite in any inertial frame in an orbit 
% at initial time (t0) and determines Position & Velocity Vectors for the satellite at a future 
% time (t) using the Universal Variable Formulation of Kepler's Problem [CHI,PSI]. Function is 
% valid for all conic orbits and special cases of Equatorial Elliptic, Inclined Circular, and 
% Equatorial Circular orbits. 2-body dynamics defines orbital motion.
%
% Author: Matthew Buckhout
% Updated: 09/06/2020 
%
% Inputs:
%
%     - [mu]         Gravitational Parameter of Central Body        [km^3/s^2]
%     - [R0]         Position Vector at initial time                [km]
%     - [V0]         Velocity Vector at initial time                [km/s]
%     - [MJD0]       Initial Time, Modified Julian Date (UT1)       [days]
%     - [MJD]        Final Time, Modified Julian Date (UT1)         [days]
%
% Outputs:
%
%     - [R]          Position Vector at final time                   [km]
%     - [V]          Velocity Vector at final time                   [km/s]
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 66-71, 90-102)
% -------------------------------------------------------------------------------------------------

function [R,V] = propagator_HC(mu,R0,V0,MJD0,MJD)

   %Time
   dt = (MJD - MJD0)*86400; %[sec] Change in Time
   if (dt == 0)
      R = R0;
      V = V0;
   else
   
      %Eccentricity
      r0 = norm(R0); %[km] Radius Magnitude
      v0 = norm(V0); %[km/s] Velocity Magnitude
      E = (1/mu)*((v0^2 - (mu/r0))*R0 - (dot(R0,V0))*V0); %Eccentricity Vector
      e = norm(E); %Eccentricity
      
      eps = ((v0^2)/2) - (mu/r0); %[km^2/s^2]  Specific Mechanical Energy
      alpha = -((v0^2)/mu) + (2/r0); %[1/km]  Parameter (alpha = 1/a)
      
      %SELECTING INITIAL CHI
      if (alpha == 1)

         fprintf('First Guess too close to converge\n');

      elseif (e < 1)
         
         %Elliptic and Circular
         CHI0 = sqrt(mu)*dt*alpha;
         
      elseif (e == 1)
         
         %Parabolic
         H = cross(R0,V0); %Angular Momentum Vector      
         h = norm(H); %Specific Angular Momentum Magnitude
         p = (h^2)/mu; %Parabolic Semi-Parameter
         s = 2*acot(3*sqrt(mu/(p^3))*dt); %Intermediate Angle 1
         w = atan((tan(s))^(1/3)); %Intermediate Angle 2
         CHI0 = sqrt(p)*2*cot(2*w);
         
      elseif (e > 1)
         
         %Hyperbolic
         a = 1/alpha; %Semi-Major Axis
         CHI0 = (dt/abs(dt))*sqrt(-a)*log( (-2*mu*alpha*dt)/(dot(R0,V0) + ... 
                (dt/abs(dt))*sqrt(-mu*a)*(1-(r0*alpha))) );
                
      end

      %NEWTON RAPHSON - FIRST ITERATION
      CHI1 = CHI0;
      PSI = (CHI0^2)*alpha;
      %Common Terms in Universal Variable Formulation
      if (PSI > 1e-6)
         c2 = (1-cos(sqrt(PSI)))/PSI;
         c3 = (sqrt(PSI) - sin(sqrt(PSI)))/sqrt(PSI^3);
      elseif (PSI < -1e-6)
         c2 = (1 - cosh(sqrt(-PSI)))/PSI;
         c3 = (sinh(sqrt(-PSI)) - sqrt(-PSI))/sqrt((-PSI)^3);
      else
         c2 = 1/2;
         c3 = 1/6;
      end
      r = (CHI1^2)*c2 + (dot(R0,V0)/sqrt(mu))*CHI1*(1 - PSI*c3) + r0*(1 - PSI*c2);
      CHI2 = CHI1 + ((sqrt(mu)*dt - (CHI1^3)*c3 - (dot(R0,V0)/sqrt(mu))*(CHI1^2)*c2 ... 
             - r0*CHI1*(1 - PSI*c3))/r);
      
      %NEWTON RAPHSON - N+1 ITERATIONS
      diff = 1;
      iterations = 1;
      while (diff > 1e-10)
         
         CHI1 = CHI2;
         
         %Universal Variable 2
         PSI = (CHI1^2)*alpha;
         %Common Terms in Universal Variable Formulation
         if (PSI > 1e-6)
            c2 = (1-cos(sqrt(PSI)))/PSI;
            c3 = (sqrt(PSI) - sin(sqrt(PSI)))/sqrt(PSI^3);
         elseif (PSI < -1e-6)
            c2 = (1 - cosh(sqrt(-PSI)))/PSI;
            c3 = (sinh(sqrt(-PSI)) - sqrt(-PSI))/sqrt((-PSI)^3);
         else
            c2 = 1/2;
            c3 = 1/6;
         end
         r = (CHI1^2)*c2 + (dot(R0,V0)/sqrt(mu))*CHI1*(1 - PSI*c3) + r0*(1 - PSI*c2);
         CHI2 = CHI1 + ((sqrt(mu)*dt - (CHI1^3)*c3 - (dot(R0,V0)/sqrt(mu))*(CHI1^2)*c2 ... 
                - r0*CHI1*(1 - PSI*c3))/r);
             
         iterations = iterations + 1;
         if (iterations > 1000)
            fprintf('propagation_HC - Failed to converge\n');
            break
         end
         
         diff = abs(CHI2 - CHI1);
         
      end
      
      %F AND G FUNCTIONS
      f = 1 - ((CHI2^2)/r0)*c2;
      g = dt - ((CHI2^3)/sqrt(mu))*c3;
      fdot = (sqrt(mu)/(r*r0))*CHI2*(PSI*c3 - 1);
      gdot = 1 - ((CHI2^2)/r)*c2;
      
      %POSITION AND VELOCITY VECTORS (t)
      R = f*R0 + g*V0;
      V = fdot*R0 + gdot*V0;
      
      %CONVERGENCE CHECK
      check = f*gdot - fdot*g;
      if (abs(check-1) > 1e-8);
         fprintf('f and g functions fail\n');
      end
  
   end
end
   