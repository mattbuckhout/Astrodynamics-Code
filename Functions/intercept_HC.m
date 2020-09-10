
% HELIOCENTRIC INTERCEPT - NON-COPLANAR INTERPLANETARY TRAJECTORIES
% ---------------------------------------------------------------------------------------------------------------------
% Determines delta-V vectors to place intercepting spacecraft on orbit to intercept target after a given delay and 
% Time of Flight (TOF). Function uses solution to Lambert's Problem to generate transfer orbit. 
%
% Assumptions:
%
%     - Two-body dynamics defines motion
%     - All maneuvers employ impulsive burns
%
% Author: Matthew Buckhout
% Updated: 09/07/2020 
%
% Inputs:
%
%     - [Rtgt0]      Target Position Vector (t0)                              [km]
%     - [Rint0]      Intercepting Spacecraft Position Vector (t0)             [km]
%     - [Vtgt0]      Target Velocity Vector (t0)                              [km/s]
%     - [Vint0]      Intercepting Spacecraft Velocity Vector (t0)             [km/s]
%     - [MJD0]       Initial Epoch for Input Vectors (UT1)                    [days]
%     - [TOF]        Time of Flight, Transfer Burn to Rendevous               [days]
%     - [delay]      Time between Initial Epoch and Transfer Burn             [days]                         
%     - [tm]         Transfer Method                                           -
%                       tm = +1 (Short Way)
%                       tm = -1 (Long Way)
%                       tm = +0 (Choose Way for Minimum Delta-V)
%
% Outputs:
%
%     - [Rtgt1]      Target Position Vector at t1                             [km]
%     - [Vtgt1]      Target Velocity Vector at t1                             [km/s]
%     - [Rtgt2]      Target & Intercepting S/C Position Vector at t2          [km]
%     - [Vtgt2]      Target & Intercepting S/C Velocity Vector at t2          [km/s]
%     - [Rint1]      Intercepting S/C Position Vector (t1)                    [km]
%     - [Vint1]      Intercepting S/C Velocity Vector (t1)                    [km/s]
%     - [Rtran1]     Intercepting S/C Position Vector on Transfer Orbit (t1)  [km]
%     - [Vtran1]     Intercepting S/C Velocity Vector on Transfer Orbit (t1)  [km/s]
%     - [dV1]        Delta-V to place Intercepting S/C on Transfer Orbit      [km/s]
%     - [dV2]        Delta-V to Match Orbits with Target                      [km/s]
%                    (Relative Heliocentric Velocity)
%     - [tm]         Selected Transfer Method                                  -
%
% Functions:
%
%     - planet_ephemerides_full_c
%     - propagator_c
%     - lambert_c 
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 448-487)
% ---------------------------------------------------------------------------------------------------------------------

function [Rtgt1,Vtgt1,Rtgt2,Vtgt2,Rint1,Vint1,Rtran1,Vtran1,dV1,dV2,tm] = intercept_HC(Rtgt0,Rint0,Vtgt0,Vint0,MJD0,TOF,delay,tm)

   %Constants
   mu = 1.32712428e11; %[km^3/s^2]  Sun Gravitational Parameter
   mode = 2; %Use 2-body dynamics for initial/final positions (1) or ephemerides for accurate positions (2)

   %Time Conversion
   MJD1 = MJD0 + delay; %Time of Transfer Burn
   MJD2 = MJD1 + TOF; %Time of Intercept

   %Target and Intercepting S/C after Delay
   if (MJD0 ~= MJD1)
      
      [RV] = propagator_c(mu,Rint0,Vint0,MJD0,MJD1); %Position and Velocities (HELIOCENTRIC)
      Rint1 = RV(:,1); %[km] (HELIOCENTRIC)
      Vint1 = RV(:,2); %[km/s] (HELIOCENTRIC)
      
      [RV] = propagator_c(mu,Rtgt0,Vtgt0,MJD0,MJD1);
      Rtgt1 = RV(:,1); %[km] (HELIOCENTRIC)
      Vtgt1 = RV(:,2); %[km/s] (HELIOCENTRIC)

   else
   
      Rint1 = Rint0; %[km] (HELIOCENTRIC)
      Rtgt1 = Rtgt0; %[km] (HELIOCENTRIC)
      Vint1 = Vint0; %[km/s] (HELIOCENTRIC)
      Vtgt1 = Vtgt0; %[km/s] (HELIOCENTRIC)
      
   end
   
   %TargeT after TOF
   [RV] = propagator_c(mu,Rtgt1,Vtgt1,MJD1,MJD2); %Target S/C Position and Velocity Vectors (HELIOCENTRIC)
   Rtgt2 = RV(:,1); %[km] (HELIOCENTRIC)
   Vtgt2 = RV(:,2); %[km/s] (HELIOCENTRIC)
   
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
   [VS] = lambert_c(mu,Rint1,Rtgt2,(TOF*24*3600),tm); %Intercepting S/C Velocity Vectors on Transfer Orbit (HELIOCENTRIC)
   Vtran1 = VS(:,1); %[km/s] (HELIOCENTRIC)
   Vtran2 = VS(:,2); %[km/s] (HELIOCENTRIC)
   Rtran1 = Rint1; %[km] Intercepting S/C Position (HELIOCENTRIC)
   
   dV1 = Vtran1 - Vint1; %Delta-V Vector Transfer Burn (HELIOCENTRIC)
   dV2 = Vtgt2 - Vtran2; %Delta-V Vector Rendevous Burn (HELIOCENTRIC)

end