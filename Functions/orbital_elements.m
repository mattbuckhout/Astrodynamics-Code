
% POSITION AND VELOCITY --> ORBITAL ELEMENTS
% -------------------------------------------------------------------------------------------------
% Function takes gravitational parameter for given attracting body, Position vector, and Velocity 
% vector in an inertial frame (GCRF/ICRF) and computes the classical Keplerian orbital elements. 
% Function accounts for special cases of Elliptic Equatorial, Circular Equatorial, and Circular 
% Inclined orbits, returning special parameters to fully describe each if necessary. 
%
% Author: Matthew Buckhout
% Updated: 08/10/2020 
%
% Inputs:
%
%     - [mu]            Gravitational Parameter of Central Body         [km^3/s^2]
%     - [R]             Position Vector                                 [km]
%     - [V]             Velocity Vector                                 [km/s]
%
% Outputs:
%
%     - [p]             Semi-Parameter                                  [km]
%		- [a] 	         Semi-Major Axis                                 [km]
% 		- [e] 	         Eccentricity                                     -
% 		- [i]             Inclination                                     [rad]
% 		- [OMEGA]         Right Ascension of Ascending Node               [rad]
% 		- [omega]         Argument of Periapsis                           [rad]
% 		- [theta]         True Anomaly at Epoch                           [rad]
%     - [omega_true]    True Longitude of Periapsis                     [rad]
%                          (Elliptic equatorial, no OMEGA/omega) 
%     - [lambda_true]   True Longitude                                  [rad]
%                          (Circular equatorial, no OMEGA/theta) 
%     - [u]             Argument of Latitude                            [rad]
%                          (Circular inclined, no OMEGA/omega/theta) 
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 104-121)
% -------------------------------------------------------------------------------------------------

function [p,a,e,i,OMEGA,omega,theta,omega_true,lambda_true,u] = orbital_elements(mu,R,V)

   I = [1 0 0]; %Unit Vector, X direction
   J = [0 1 0]; %Unit Vector, Y direction
   K = [0 0 1]; %Unit Vector, Z direction
   
   %Internal Calculations with row vectors
   if (isrow(R) == 0) 
      R = transpose(R);
      V = transpose(V);
   end

   H = cross(R,V); %Angular Momentum Vector
   r = norm(R); %Radius Magnitude
   v = norm(V); %Velocity Magnitude
   h = norm(H); %Angular Momentum Magnitude

   E = (1/mu)*((v^2 - (mu/r))*R - (dot(R,V))*V); %Eccentricity Vector
   e = norm(E); %Eccentricity
   N = cross(K,H); %Node Line Vector
   n = norm(N); %Node Line Magnitude
   a = (h^2)/(mu*(1 - (e^2))); %Semi-Major Axis
   i = acos(dot(H,K)/h); %Inclination
   
   if ((e < 1) || (e == 0))
      p = a*(1-(e^2)); %Semi-Parameter
   elseif (e > 1)
      p = a*(1-(e^2)); 
   elseif (e == 1)
      p = (h^2)/mu; 
   end
   
   if (e == 0)
      omega_true = NaN;
   else
      if (dot(E,J) >= 0)
         omega_true = acos(dot(E,I)/e); %True Longitude of Periapsis 
      elseif (dot(E,J) < 0)
         omega_true = 2*pi - acos(dot(E,I)/e);
      end
   end
   
   if (n == 0)
      u = NaN;
   else
      if (dot(R,K) >= 0)
         u = acos(dot(N,R)/(n*r)); %Argument of Latitude 
      elseif (dot(R,K) < 0)
         u = 2*pi - acos(dot(N,R)/(n*r));
      end
   end

   if (dot(R,J) >= 0)
      lambda_true = acos(dot(R,I)/r); %True Longitude
   elseif (dot(R,J) < 0)
      lambda_true = 2*pi - acos(dot(R,I)/r);
   end

   if (n == 0)
      OMEGA = NaN;
   else
      if (dot(N,J) >= 0)
         OMEGA = acos(dot(N,I)/n); %Right Ascension of Ascending Node
      else
         OMEGA = (2*pi) - acos(dot(N,I)/n);
      end
   end
   
   if (n == 0)
      omega = NaN;
   else
      if (dot(E,K) >= 0)
         omega = acos(dot(N,E)/(n*e)); %Argument of Periapsis
      else
         omega = (2*pi) - acos(dot(N,E)/(n*e)); 
      end
   end
   
   if (e == 0)
      theta = NaN;
   else
      if (dot(R,V) >= 0)
         theta = acos(dot(E,R)/(e*r)); %True Anomaly at Epoch
      else
         theta = 2*pi - acos(dot(E,R)/(e*r));
      end 
   end
     
end


