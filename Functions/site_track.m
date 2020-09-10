
% ORBIT DETERMINATION FROM RANGE, AZIMUTH, ELEVATION, + RATES
% -------------------------------------------------------------------------------------------------
% Takes observations of a satellite in Elevation, Azimuth, Range, and rates of each measured by 
% radar at a given Date and Time (UTC) in the Topocentric-Horizon (SEZ) frame and determines the 
% Position and Velocity vectors of the satellite in the inertial Geocentric Celestial Reference 
% Frame (GCRF).
%
% Author: Matthew Buckhout
% Updated: 09/09/2020 
%
% Inputs:
%
%     - [Lat]        Geodetic Latitude                               [rad] 
%     - [Long]       Longitude (+East from Greenwich)                [rad] 
%     - [h_ellip]          Site height above reference ellipsoid           [km]
%     - [rho]        Range (site to satellite)                       [km]
%     - [AZ]         Azimuth (+Clockwise from North)                 [rad]
%     - [EL]         Elevation (+Up from Horizon)                    [rad]
%     - [drho]       Range Rate                                      [km/s] 
%     - [dAZ]        Azimuth Rate                                    [rad/s]
%     - [dEL]        Elevation Rate                                  [rad/s]
%     - [DATE]       Gregorian Date                                  [yyyy mm dd]
%     - [UTC]        Coordinated Universal Time                      [hh mm ss]   
%
% Outputs:
%
%     - [R_GCRF]     Position Vector (GCRF)                          [km]
%     - [V_GCRF]     Velocity Vector (GCRF)                          [km/s]
%
% Functions:
%
%     - julian_date
%     - getdata_EOP_fast  
%     - getdata_IAU2006
%     - EOP
%     - site
%     - IAU2006 
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 406-411)
% -------------------------------------------------------------------------------------------------

function [R_GCRF,V_GCRF] = site_track(Lat,Long,h_ellip,rho,AZ,EL,drho,dAZ,dEL,DATE,UTC)
   
   [~,MJD_UTC] = julian_date(DATE,UTC,0); %Modified Juilan Date (UTC)

   [EOPdata] = getdata_EOP_fast(); %Earth Orientation Data Table
   
   [xp,yp,dUT1,dAT] = EOP(MJD_UTC,EOPdata); %Earth Orientation Parameters on Date
   
   MJD = MJD_UTC + (dUT1/86400); %Modified Julian Date (UT1)
   
   [XTAB,YTAB,sTAB] = getdata_IAU2006('short'); %Coefficients for IAU2006 Precession-Nutation Model

   [Rsite_ITRF] = site(Lat,Long,h_ellip); %[km] Site Position Vector (ITRF)

   RHO_SEZ = [-rho*cos(EL)*cos(AZ); ...
               rho*cos(EL)*sin(AZ); ...
               rho*sin(EL)]; %[km] Range Vector (SEZ)
               
   RHOdot_SEZ = [(-drho*cos(EL)*cos(AZ) + rho*sin(EL)*cos(AZ)*dEL + rho*cos(EL)*sin(AZ)*dAZ); ... 
                 (drho*cos(EL)*sin(AZ) - rho*sin(EL)*sin(AZ)*dEL + rho*cos(EL)*cos(AZ)*dAZ); ...
                 (drho*sin(EL) + rho*cos(EL)*dEL)]; %[km/s] Range Rate Vector (SEZ)  
                 
   RHO_ITRF = [(sin(Lat)*cos(Long)) (-sin(Long)) (cos(Lat)*cos(Long)); ...
               (sin(Lat)*sin(Long)) (cos(Long)) (cos(Lat)*sin(Long)); ... 
               (-cos(Lat)) 0 sin(Lat)]*RHO_SEZ; %[km] Range Vector (IRTF)

   RHOdot_ITRF = [sin(Lat)*cos(Long) (-sin(Long)) cos(Lat)*cos(Long); ...
                  sin(Lat)*sin(Long) (cos(Long)) cos(Lat)*sin(Long); ... 
                  (-cos(Lat)) 0 sin(Lat)]*RHOdot_SEZ; %[km/s] Range Rate Vector (ITRF)

   R_ITRF = RHO_ITRF + Rsite_ITRF; %[km] Position Vector (ITRF)
   V_ITRF = RHOdot_ITRF; %[km/s] Velocity Vector (ITRF)   
     
   [R_GCRF,V_GCRF] = IAU2006(R_ITRF,V_ITRF,MJD,dUT1,dAT,xp,yp,2,'short',XTAB,YTAB,sTAB); %(GCRF)
   
end
