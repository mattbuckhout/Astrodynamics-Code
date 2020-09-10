
% ORBITAL ELEMENTS --> POSITION AND VELOCITY
% -------------------------------------------------------------------------------------------------
% Function takes gravitational parameter for given attracting body and the classical Keplerian 
% orbital elements for an orbit and computes the Position and Velocity vector in the inertial frame
% corresponding to the input elements (GCRF/ICRF). Function accounts for special cases of Elliptic 
% Equatorial, Circular Equatorial, and Circular Inclined orbits.  
%
% Author: Matthew Buckhout
% Updated: 08/24/2020 
%
% Inputs:
%
%     - [mu]            Gravitational Parameter of Central Body         [km^3/s^2]
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
% Outputs:
%
%     - [R]             Position Vector                                 [km]
%     - [V]             Velocity Vector                                 [km/s]
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 104-126)
% -------------------------------------------------------------------------------------------------

function [R,V] = position_velocity(mu,p,e,i,OMEGA,omega,theta,omega_true,lambda_true,u)

   % Perifocal Coordinate System Unit Vectors
   P = [1 0 0]; %Unit Vector, Xw direction
   Q = [0 1 0]; %Unit Vector, Yw direction
   W = [0 0 1]; %Unit Vector, Zw direction
   
   %Special Case Logic (Changing elements within function for internal calculation ONLY)
   if (e == 0 && i == 0) %Circular Equatorial
      omega = 0;
      OMEGA = 0;
      theta = lambda_true;
   elseif (e == 0 && i != 0) %Circular Inclined
      omega = 0;
      theta = u;
   elseif (e > 0 && i == 0) %Elliptic Equatorial
      OMEGA = 0;
      omega = omega_true;
   end
     
   r = p/(1 + e*cos(theta)); %Radius Magnitude
      
   Rpf = r*cos(theta)*P + r*sin(theta)*Q; %Radius Vector (Perifocal)
   Vpf = sqrt(mu/p)*(-sin(theta)*P + (e+cos(theta))*Q); %Velocity Vector (Perifocal)

   %PERIFOCAL TO GEOCENTRIC-EQUATORIAL TRANSFORM [PQW-IJK]
   %Rotation Matrix Components
   T11  = cos(OMEGA)*cos(omega) - sin(OMEGA)*sin(omega)*cos(i);
   T12  = -cos(OMEGA)*sin(omega) - sin(OMEGA)*cos(omega)*cos(i);
   T13  = sin(OMEGA)*sin(i);
   T21  = sin(OMEGA)*cos(omega) + cos(OMEGA)*sin(omega)*cos(i);
   T22  = -sin(OMEGA)*sin(omega) + cos(OMEGA)*cos(omega)*cos(i);
   T23  = -cos(OMEGA)*sin(i);
   T31  = sin(omega)*sin(i);
   T32  = cos(omega)*sin(i);
   T33  = cos(i);

   Tmat =   [T11 T12 T13;
             T21 T22 T23;
             T31 T32 T33];

   R = Tmat*transpose(Rpf); %Radius [Rx,Ry,Rz] coordinates for orbit
   V = Tmat*transpose(Vpf); %Velocity [Vx,Vy,Vz] for orbit
   
end
