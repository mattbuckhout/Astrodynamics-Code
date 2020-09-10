
% POSITION/VELOCITY VECTORS TO RANGE, AZIMUTH, ELEVATION, + RATES
% -------------------------------------------------------------------------------------------------
% Function takes the state of a satellite (Position and Velocity vectors) in the GCRF frame and 
% finds the range, observation angles, and rates of each in the Topocentric-Horizon (SEZ) frame.
%
% Author: Matthew Buckhout
% Updated: 09/09/2020 
%
% Inputs:
%
%     - [R_GCRF]     Position Vector (GCRF)                          [km]
%     - [V_GCRF]     Velocity Vector (GCRF)                          [km/s]
%     - [Lat]        Geodetic Latitude                               [rad] 
%     - [Long]       Longitude (+East from Greenwich)                [rad] 
%     - [h_ellip]          Site height above reference ellipsoid           [km]
%     - [DATE]       Gregorian Date                                  [yyyy mm dd]
%     - [UTC]        Coordinated Universal Time                      [hh mm ss]   
%
% Outputs:
%
%     - [rho]        Range (site to satellite)                       [km]
%     - [AZ]         Azimuth (+Clockwise from North)                 [rad]
%     - [EL]         Elevation (+Up from Horizon)                    [rad]
%     - [drho]       Range Rate                                      [km/s] 
%     - [dAZ]        Azimuth Rate                                    [rad/s]
%     - [dEL]        Elevation Rate                                  [rad/s]
%
% Functions:
%
%     - julian_date
%     - getdata_EOP_fast
%     - EOP
%     - getdata_IAU2006
%     - IAU2006
%     - site  
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 256-257)
% -------------------------------------------------------------------------------------------------

function [rho,AZ,EL,drho,dAZ,dEL] = RAZEL(R_GCRF,V_GCRF,Lat,Long,h_ellip,DATE,UTC)

   [~,MJD_UTC] = julian_date(DATE,UTC,0); %Modified Juilan Date (UTC)

   [EOPdata] = getdata_EOP_fast(); %Earth Orientation Data Table
   
   [xp,yp,dUT1,dAT] = EOP(MJD_UTC,EOPdata); %Earth Orientation Parameters on Date
   
   MJD = MJD_UTC + (dUT1/86400); %Modified Julian Date (UT1)
   
   [XTAB,YTAB,sTAB] = getdata_IAU2006('short'); %Coefficients for IAU2006 Precession-Nutation Model
   
   [R_ITRF,V_ITRF] = IAU2006(R_GCRF,V_GCRF,MJD,dUT1,dAT,xp,yp,1,'short',XTAB,YTAB,sTAB); %[km] [km/s] Sattelite state transform from GCRF to ITRF
   [Rsite_ITRF] = site(Lat,Long,h_ellip); %[km] Observation Site Position Vector
   RHO_ITRF = R_ITRF - Rsite_ITRF; %[km] Range Vector (ITRF)
   dRHO_ITRF = V_ITRF; %[km/s] Range Rate Vector (ITRF)
   
   RHO_SEZ = [sin(Lat)*cos(Long) sin(Lat)*sin(Long) -cos(Lat); ...
                -sin(Long) cos(Long) 0; ...
                cos(Lat)*cos(Long) cos(Lat)*sin(Long) sin(Lat)]*RHO_ITRF; %[km] Range Vector (SEZ)
   
   dRHO_SEZ = [sin(Lat)*cos(Long) sin(Lat)*sin(Long) -cos(Lat); ...
                -sin(Long) cos(Long) 0; ...
                cos(Lat)*cos(Long) cos(Lat)*sin(Long) sin(Lat)]*dRHO_ITRF; %[km/s] Range Rate Vector (SEZ)
   
   rho = norm(RHO_SEZ); %[km] Range Magnitude
      
   EL = asin(RHO_SEZ(3)/rho); %[rad] Elevation
   
   if (EL ~= (pi/2))
   
      SINB = RHO_SEZ(2)/sqrt((RHO_SEZ(1)^2) + (RHO_SEZ(2)^2));
      COSB = RHO_SEZ(1)/sqrt((RHO_SEZ(1)^2) + (RHO_SEZ(2)^2));
      AZ = abs(atan(SINB/COSB)); %[rad] Azimuth
   
   elseif (EL == (pi/2))

      SINB = dRHO_SEZ(2)/sqrt((dRHO_SEZ(1)^2) + (dRHO_SEZ(2)^2));
      COSB = dRHO_SEZ(1)/sqrt((dRHO_SEZ(1)^2) + (dRHO_SEZ(2)^2));
      AZ = abs(atan(SINB/COSB)); %[rad] Preliminary Azimuth

   end

   %Azimuth Quadrant Checks
   if ((COSB > 0) && (SINB > 0))
      AZ = pi - AZ;
   elseif ((COSB > 0) && (SINB < 0))
      AZ = pi + AZ;
   elseif ((COSB < 0) && (SINB < 0))
      AZ = 2*pi - AZ;
   elseif ((COSB < 0) && (SINB > 0))
      AZ = AZ;
   end
   
   drho = dot(RHO_SEZ,dRHO_SEZ)/rho; %[km/s] Range Rate Magnitude
   dAZ = (dRHO_SEZ(1)*RHO_SEZ(2) - dRHO_SEZ(2)*RHO_SEZ(1))/(RHO_SEZ(1)^2 + RHO_SEZ(2)^2); %[rad/s] Azimuth Rate
   dEL = (dRHO_SEZ(3) - drho*sin(EL))/sqrt(RHO_SEZ(1)^2 + RHO_SEZ(2)^2); %[rad/s] Elevation Rate
   
end
