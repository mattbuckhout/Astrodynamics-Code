
% UNIT CONVERSION - SI TO CANONICAL UNITS
% -------------------------------------------------------------------------------------------------
% Functions converts Position and Velocity vectors from SI units to Canonical Units and vice versa 
% for Geocentric orbits (using Earth Gravitational Parameter). 
% Canonical Units:
%
%     - Gravitational Parameter (mu):  ER^3/TU^2 = 1
%     - Distance Unit (ER):            1 Earth Radius, 6378.1363 km.
%     - Time Unit (TU):                Time for satellite orbiting at 1 ER to travel 1 Radian
%                                         TU = sqrt(ER^3 / mu) sec.
%
% Author: Matthew Buckhout
% Updated: 08/15/2020 
%
% Inputs:
%
%     - [R1]         Position Vector (input units)                   [km] [DU]
%     - [V1]         Velocity Vector (input units)                   [km/s] [ER/TU]
%     - [way]        SI to Canonical (1), Canonical to SI (2)         -
%
% Outputs:
%
%     - [R2]         Position Vector (output units)                  [km] [DU]
%     - [V2]         Velocity Vector (output units)                  [km/s] [ER/TU]
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 232-233)
% -------------------------------------------------------------------------------------------------

function [R2,V2] = si_canonical(R1,V1,way)

   mu = 3.986004415e5; %[km^3/s^2] Earth Gravitational Parameter
   ER = 6378.1363; %[km] Earth Mean Equatorial Radius 
   TU = sqrt((ER^3)/mu); %[sec] Time Unit
   
   %SI to Canonical
   if (way == 1)
      
      R2 = R1/ER; %[ER]
      V2 = (V1/ER)*TU; %[ER/TU]
      
   %Canonical to SI
   elseif (way == 2)
   
      R2 = R1*ER; %[km]
      V2 = (V1/TU)*ER; %[km/s]
   
   end
end