
% ORBIT PROPAGATION - KEPLER'S PROBLEM, CLASSICAL ORBITAL ELEMENTS SOLUTION
% -------------------------------------------------------------------------------------------------
% Function takes Position & Velocity vectors for a satellite in inertial frame in an orbit at 
% initial time (t0) and determines Position & Velocity Vectors for the satellite at a future time 
% (t). Function is valid for all conic orbits and special cases of Equatorial Elliptic, Inclined 
% Circular, and Equatorial Circular orbits. 2-body dynamics defines orbital motion.
%
% Author: Matthew Buckhout
% Updated: 08/24/2020 
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
%     -  orbital_elements
%     -  eccentric_anomaly
%     -  kepler_equation
%     -  position_velocity
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 66-71, 87-89)
% -------------------------------------------------------------------------------------------------

% Validation Test Cases:
%  [X]   Elliptical Inclined
%  [ ]   Circular Inclined
%  [X]   Circular Equatorial
%  [X]   Elliptical Equatorial
%  [ ]   Parabolic
%  [ ]   Hyperbolic

function [R,V] = propagation_OE(mu,R0,V0,t0,t)
   
   %ORBITAL ELEMENTS (t0)
   [p,a,e,i,OMEGA,omega,theta0,omega_true,lambda_true0,u0] = orbital_elements(mu,R0,V0);
   
   %CHANGE IN TIME
   dt = t-t0;

   %TRUE ANOMALY (t0)
   if (e == 0)
      EBH0 = lambda_true0;
   else
      [EBH0] = eccentric_anomaly(e,theta0);
   end
   
   %COMPUTING MEAN ANOMALY (t) 
   if (e < 1)
      
      %Elliptic
      p = a*(1-(e^2)); %Semi-Parameter
      n = sqrt(mu/(a^3)); %Mean Motion
      M0 = EBH0 - e*sin(EBH0); %Initial Mean Anomaly
      M = M0 + n*dt; %Final Mean Anomaly
      
   elseif (e > 1)
      
      %Hyperbolic
      p = a*(1-(e^2)); %Semi-Parameter
      nh = sqrt(mu/(-(a^3))); %Hyperbolic Mean Motion
      M0 = e*sinh(EBH0) - EBH0; %Initial Mean Anomaly
      M = M0 + nh*dt; %Final Mean Anomaly
      
   elseif (e == 1)
      
      %Parabolic     
      H0 = cross(R0,V0); %Angular Momentum Vector
      h0 = norm(H0); %Specific Angular Momentum
      p = (h0^2)/mu; %Semi-Parameter
      r = ((p/2)*(1 + (EBH^2))); %Radius Magnitude
      np = 2*sqrt(mu/(p^3)); %Parabolic Mean Motion
      M0 = EBH0 + ((EBH0^3)/3); %Initial Mean Anomaly
      M = M0 + np*dt; %Final Mean Anomaly
      
   end
   
   %COMPUTING ECCENTRIC ANOMALY (t)
   if (e == 0)
      n = sqrt(mu/(a^3)); %Mean Motion
      EBH = EBH0 + n*dt; %Circular Eccentric Anomaly = Mean Anomaly      
   else
      [EBH] = kepler_equation(mu,M,e,dt,p);
   end
            
   %SATELLITE LOCATION ANGLES (t)
   if (e < 1 && e ~= 0)
      
      %Elliptic
      theta = 2*atan(sqrt((1+e)/(1-e))*tan(EBH/2));
      lambda_true = NaN; 
      u = NaN;

   elseif (e > 1)
      
      %Hyperbolic
      theta = 2*atan(sqrt((e+1)/(e-1))*tanh(EBH/2));
      lambda_true = NaN;
      u = NaN;
      
   elseif (e == 1)
      
      %Parabolic
      theta = atan((r*p*EBH)/(r*(p-r)));
      lambda_true = NaN;
      u = NaN;
      
   elseif (e == 0)
      
      %Circular
      u = EBH;
      lambda_true = EBH;
      theta = NaN;
      
   end
   
   %RADIUS AND VELOCITY VECTORS (t)
   [R,V] = position_velocity(mu,p,e,i,OMEGA,omega,theta,omega_true,lambda_true,u);
   
   %UPDATED ORBITAL ELEMENTS AND ADDITIONAL ANGLES (t)
   [a,e,i,OMEGA,omega,theta,omega_true,lambda_true,u] = orbital_elements(mu,R,V);
   
end
