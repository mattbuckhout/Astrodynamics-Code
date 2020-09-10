
% AZIMUTH/ELEVATION TO TOPOCENTRIC RIGHT ASCENSION/DECLINATION
% -------------------------------------------------------------------------------------------------
% Converts observed Azimuth and Elevation angles at a specified date and time (UTC) to Right 
% Ascension and Declination in the Topocentric Horizon frame.
% 
% Author: Matthew Buckhout
% Updated: 09/09/2020
%
% Inputs:
%                  
%     - [AZ]            Azimuth (+Clockwise from North)                [rad]
%     - [EL]            Elevation (+Up from Horizon)                   [rad]
%     - [Lat]           Geodetic Latitude                              [rad]
%     - [Long]          Longitude (+East from Greenwich)               [rad]
%     - [DATE]          Date of observation                            [yyyy mm dd]                 
%     - [TIME]          Coordinated Universal Time (UTC)               [hh mm ss]                   
%     - [dUT1]          UT1 - UTC                                      [sec]
%     - [dAT]           TAI - UTC                                      [sec]
%  
% Outputs:
%
%     - [RA]            Topocentric Right Ascension 
%                           (+CCW from Topo. Vernal Equinox)           [rad]                                                       
%     - [DEC]           Topocentric Declination                                                      
%                           (+Up from Topo. Equatorial Plane)          [rad]
%
% Functions:
%
%     - convert_time                               
%     - julian_date                                            
%     - local_sidereal     
%
% References: 
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 257-259)
% -------------------------------------------------------------------------------------------------

function [RA,DEC] = AZEL_RADEC(AZ,EL,Lat,Long,DATE,TIME,dUT1,dAT,leap)

   %Local Sidereal Time
   [UT1,~,~,~,~,~,~] = convert_time(DATE,TIME,dUT1,dAT); %[hh:mm:ss]
   [~,MJD_UT1] = julian_date(DATE,UT1,leap); %[days]
   [LST,~] = local_sidereal(Long,MJD_UT1); %[rad] Local Sidereal Time

   DEC = asin(sin(EL)*sin(Lat) + cos(EL)*cos(Lat)*cos(AZ)); %[rad] Topocentric Declination
   
   %Quadrant Checks - LHA
   sin_LHA = -((sin(AZ)*cos(EL))/cos(DEC));
   cos_LHA = (cos(Lat)*sin(EL) - sin(Lat)*cos(AZ)*cos(EL))/cos(DEC);
   tan_LHA = sin_LHA/cos_LHA;
   LHA = atan(tan_LHA);
   if ((cos_LHA > 0) && (sin_LHA > 0) && (tan_LHA > 0))
      LHA = abs(LHA);
   elseif ((cos_LHA < 0) && (sin_LHA > 0) && (tan_LHA < 0))
      LHA = pi - abs(LHA);
   elseif ((cos_LHA < 0) && (sin_LHA < 0) && (tan_LHA > 0))
      LHA = pi + abs(LHA);
   elseif ((cos_LHA > 0) && (sin_LHA < 0) && (tan_LHA < 0))
      LHA = (2*pi) - abs(LHA);
   end
  
   RA = LST - LHA; %[rad] Topocentric Right Ascension

   RA = mod(RA,(2*pi));

end