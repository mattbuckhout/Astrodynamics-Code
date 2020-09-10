
% ORBIT PROPAGATION - KEPLER'S PROBLEM, UNIVERSAL VARIABLE APPROACH
% -------------------------------------------------------------------------------------------------
% Function takes Position & Velocity vectors for a satellite in inertial frame in an orbit at 
% initial time (t0) and determines Position & Velocity Vectors for the satellite at a future time 
% (t) using the Universal Variable Formulation of Kepler's Problem [CHI,PSI]. Function is valid for 
% all conic orbits and special cases of Equatorial Elliptic, Inclined Circular, and Equatorial 
% Circular orbits. 2-body dynamics defines orbital motion.
%
% Author: Matthew Buckhout
% Updated: 08/11/2020 
%
% Inputs:
%
%     - [mu]          Gravitational Parameter of Central Body        [km^3/s^2]
%     - [R0]          Position Vector at initial time                [km]
%     - [V0]          Velocity Vector at initial time                [km/s]
%     - [t0]          Initial Time                                   [sec]
%     - [t]           Final Time                                     [sec]
%
% Outputs:
%
%     - [R]          Position Vector at final time                   [km]
%     - [V]          Velocity Vector at final time                   [km/s]
%
% Functions:
%
%     - si_canonical
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 66-71, 90-102)
% -------------------------------------------------------------------------------------------------

% Validation Test Cases:
%  [X]   Elliptical Inclined
%  [ ]   Circular Inclined
%  [ ]   Circular Equatorial
%  [ ]   Elliptical Equatorial
%  [ ]   Parabolic
%  [ ]   Hyperbolic

function [R,V] = propagation_UV(mu,R0,V0,t0,t)

   ER = 6378.1363; %[km] Earth Mean Equatorial Radius 
   mu_SI = mu; %[km^3/s^2] Gravitational Parameter, SI
   mu = 1; %[ER^3/TU^2] Gravitational Parameter, Canonical (unity)
   dt = (t-t0)/(sqrt((ER^3)/mu_SI)); %[TU] Change in Time, Canonical
   [R0,V0] = si_canonical(R0,V0,1); %[ER and ER/TU] Initial Vectors, Canonical
      
   r0 = norm(R0); %[ER] Radius Magnitude
   v0 = norm(V0); %[ER/TU] Velocity Magnitude

   eps = ((v0^2)/2) - (mu/r0); %[km^2/s^2]  Specific Mechanical Energy
   alpha = -((v0^2)/mu) + (2/r0); %[1/ER]  Parameter (alpha = 1/a)

   %SELECTING INITIAL CHI
   if (alpha == 1)

      fprintf('First Guess too close to converge\n');

   elseif (alpha > 0.000001)
      
      %Elliptic and Circular
      CHI0 = sqrt(mu)*dt*alpha;
      
   elseif (abs(alpha) < 0.000001)
      
      %Parabolic
      H = cross(R0,V0); %Angular Momentum Vector      
      h = norm(H); %Specific Angular Momentum Magnitude
      p = (h^2)/mu; %Parabolic Semi-Parameter
      s = 2*acot(3*sqrt(mu/(p^3))*dt); %Intermediate Angle 1
      w = atan((tan(s))^(1/3)); %Intermediate Angle 2
      CHI0 = sqrt(p)*2*cot(2*w);
      
   elseif (alpha < -0.000001)
      
      %Hyperbolic
      a = 1/alpha; %Semi-Major Axis
      CHI0 = (dt/abs(dt))*sqrt(-a)*log( (-2*mu*alpha*dt)/(dot(R0,V0) + ... 
             (dt/abs(dt))*sqrt(-mu*a)*(1-(r0*alpha))) );
             
   end

   %NEWTON RAPHSON - FIRST ITERATION
   CHI1 = CHI0;
   %Universal Variable 2
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
   
   %NEWTON RAPHSON - n+1 ITERATIONS
   iterations = 1;
   while (abs(CHI2 - CHI1) > 1e-10)
      
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
         fprintf('propagation_UV - Failed to converge\n');
         break
      end
      
   end
   
   %F AND G FUNCTIONS
   f = 1 - ((CHI2^2)/r0)*c2;
   g = dt - ((CHI2^3)/sqrt(mu))*c3;
   fdot = (sqrt(mu)/(r*r0))*CHI2*(PSI*c3 - 1);
   gdot = 1 - ((CHI2^2)/r)*c2;
   
   %RADIUS AND VELOCITY VECTORS (t)
   R = f*R0 + g*V0;
   V = fdot*R0 + gdot*V0;
   
   %CONVERGENCE CHECK
   check = f*gdot - fdot*g;
   if (abs(check-1) > 1e-8)
      fprintf('f and g functions fail\n');
   end
   
   %CONVERTING TO SI
   [R,V] = si_canonical(R,V,2); %[km and km/s] Final Vectors, SI
  
end
