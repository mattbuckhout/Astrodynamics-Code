
% INTERCEPT
% -------------------------------------------------------------------------------------------------
% Determines delta-V vectors to place intercepting spacecraft on orbit to intercept target 
% spacecraft after a given Time of Flight (TOF) and to match orbits with target spacecraft. 
% Function uses solution to Lambert's Problem to generate transfer orbit.
%
% Assumptions:
%
%     - Geocentric orbits only
%     - Two-body dynamics defines motion
%     - All maneuvers employ impulsive burns
%     - Flown transfer orbit segment must not pass into atmosphere or through Earth
%     - Spherical Earth Model
%
% Author: Matthew Buckhout
% Updated: 08/08/2020 
%
% Inputs:
%
%     - [Rtgt0]      Target Spacecraft Position Vector (t0)                   [km]
%     - [Rint0]      Intercepting Spacecraft Position Vector (t0)             [km]
%     - [Vtgt0]      Target Spacecraft Velocity Vector (t0)                   [km/s]
%     - [Vint0]      Intercepting Spacecraft Velocity Vector (t0)             [km/s]
%     - [t0]         Modified Julian Date for Input Vectors (UT1)             [days]
%     - [TOF]        Time of Flight on transfer orbit                         [sec]
%     - [delay]      Time between t0 and Transfer Burn                        [sec]
%     - [tm]         Transfer Method                                           -
%                       tm = +1 (Short Way)
%                       tm = -1 (Long Way)
%                       tm = +0 (Choose Way for Minimum Delta-V)
%
% Outputs:
%
%     - [Rtgt2]      Target & Intercepting S/C Position Vector (t2)           [km]
%     - [Vtgt2]      Target & Intercepting S/C Velocity Vector (t2)           [km/s]
%     - [dV1]        Delta-V to place Intercepting S/C on Transfer Orbit      [km/s]
%     - [dV2]        Delta-V to Match Orbits with Target                      [km/s]
%     - [tm]         Selected Transfer Method                                  -         
%
% Functions:
%
%     - propagation_UV
%     - lambert_UV  
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 448-487)
% -------------------------------------------------------------------------------------------------

function [Rtgt2,Vtgt2,Rtran1,Vtran1,dV1,dV2,tm] = intercept(Rtgt0,Rint0,Vtgt0,Vint0,t0,TOF,delay,tm)

   %Constants
   RE = 6378.1363; %[km]            Earth Mean Equatorial Radius
   mu = 3.986004415e5; %[km^3/s^2]  Earth Gravitational Parameter

   %Time Conversion
   t0sec = t0*86400; %[s] Epoch for Input Vectors
   t1sec = t0sec + delay; %[s] Time of Transfer Burn
   t2sec = t1sec + TOF; %[s] Time of Intercept

   %Interceptor S/C after Delay
   [Rint1,Vint1] = propagation_UV(mu,Rint0,Vint0,t0sec,t1sec); %Target S/C Position and Velocity Vectors (t1)
   %Target S/C after Delay
   [Rtgt1,Vtgt1] = propagation_UV(mu,Rtgt0,Vtgt0,t0sec,t1sec); %Target S/C Position and Velocity Vectors (t1)
   %Target S/C after TOF
   [Rtgt2,Vtgt2] = propagation_UV(mu,Rtgt1,Vtgt1,t1sec,t2sec); %Target S/C Position and Velocity Vectors (t2)
   
   %Selecting Transfer Method for Minimum Delta-V
   if (tm == 0)
      
      tran_n = cross(Rint1,Rtgt2); %Transfer orbit normal vector
      h_n = cross(Rint1,Vint1); %Interceptor S/C orbit normal vector
      
      if (dot(h_n,tran_n) < 0)
         tm = -1; %Long Way
      elseif (dot(h_n,tran_n) > 0)
         tm = 1; %Short Way
      end
      
   end

   %Delta-V Vectors
   [Vtran1,Vtran2] = lambert_UV(Rint1,Rtgt2,TOF,tm); %Intercepting S/C Velocity Vectors on Transfer Orbit (t1),(t2)
   Rtran1 = Rint1; %Intercepting S/C Position (t1)
   
   dV1 = Vtran1 - Vint1; %Delta-V Vector Transfer Burn (t1)
   dV2 = Vtgt2 - Vtran2; %Delta-V Vector Rendevous Burn (t2)

   %Checking Transfer Orbit for Earth Impact/Atmospheric Pass Through
   if (dot(Rint1,Vtran1) < 0) && (dot(Rtgt2,Vtran2) > 0)
      
      eps = ((norm(Vtran1)^2)/2) - (mu/norm(Rint1)); %Specific Mechanical Energy
      a = -mu/(2*eps); %[km] Semi-Major Axis
      h = cross(Rint1,Vtran1); %Specific Angular Momentum
      p = (norm(h)^2)/mu; %[km] Semi-parameter
      e = sqrt((a - p)/a); %Eccentricity
      rp = a*(1 - e); %[km] Periapsis Radius
      
      if (rp <= (RE + 100))
         
         fprintf('\n');
         fprintf('Transfer Periapsis inside atmosphere\n');
         fprintf('Altitude = %0.3f km\n',(rp - RE));
         fprintf('\n');
         
      elseif (rp <= RE)
         
         fprintf('\n');
         fprintf('Transfer Periapsis inside Earth\n');
         fprintf('Altitude = %0.3f km\n',(rp - RE));
         fprintf('\n');
         
      end
      
   end

end