
% CONVERTING SATELLITE POSITION IN INERTIAL FRAME TO SUB-LATITUDE POINT
% -------------------------------------------------------------------------------------------------
% Takes position for satellite in Geocentric-Equatorial (ECI) frame and finds point directly 
% beneath a satellite on Earth Reference Ellipsoid defined by Geodetic Latitude, Longitude 
% (+East of Greenwich), and Altitude above Mean Sea Level. Used for ground track. Only accounts for
% Earth rotation in transform, no polar motion or precession-nutation included.
%
% Author: Matthew Buckhout
% Updated: 09/09/2020 
%
% Inputs:
%                                                   
%     - [GMST]          Greenwich Mean Sidereal Time                  [rad]
%     - [R_ECI]         Satellite Position Vector,                     [km]
%                       Earth Centered Inertial frame (GCRF)                         
%
% Outputs:
%
%     - [Lat]           Geodetic Latitude                             [deg]
%     - [Long]          Longitude (+East from Greenwich)              [deg]
%     - [h_ellip]       Height above reference ellipsoid              [km]
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 138, 179-197) 
%     - IERS Conventions 2010 (pg. 19, 40)
%           https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html
% -------------------------------------------------------------------------------------------------

function [Lat,Long,h_ellip] = ECI_ground(R_ECI,GMST)
   
   %Earth Reference Ellipsoid - GRS80
   RE = 6378.1363; %[km] Earth Semi-Major Axis 
   fE = 1/298.257223563; % Earth Flattening 
   bE = RE - (fE*RE); %[km] Earth Semi-minor Axis (Polar Radius)
   eE = sqrt(RE^2 - bE^2)/RE; % Earth Reference Ellipsoid Eccentricity
   
   r_delta_sat = sqrt(((R_ECI(1))^2) ... 
                 + ((R_ECI(2))^2)); %[km]      Satellite Equatorial Component
   
   %Geocentric Right Ascension + Quadrant Checks
   alpha1 = (R_ECI(2)/r_delta_sat);
   alpha2 = (R_ECI(1)/r_delta_sat);
   
   alpha = atan((R_ECI(2)/r_delta_sat)/(R_ECI(1)/r_delta_sat));
   
   if ((alpha1 > 0) && (alpha2 > 0))
      alpha = alpha;
   elseif ((alpha1 > 0) && (alpha2 < 0))
      alpha = pi - (-alpha);
   elseif ((alpha1 < 0) && (alpha2 < 0))
      alpha = pi + alpha;
   elseif ((alpha1 < 0) && (alpha2 > 0))
      alpha = (2*pi) - (-alpha); %[rad]    Geocentric Right Ascension
   end
   
   Long = alpha - GMST; %[rad]                  Longitude 

   %Positive East of Greenwich + Values between 0 and 180
   if (Long > 0)
      if (Long > pi)
         Long = (-((2*pi) - Long));
      end
   end
      
   if (Long < 0)
      Long = ((2*pi) + Long);
      if (Long > pi)
         Long = (-((2*pi) - Long));
      end
   end
  
   delta = atan(R_ECI(3)/r_delta_sat); %[rad]  Geocentric Declination
   r = norm(R_ECI); %[km]                      Magnitude of Satellite Position Vector
  
   %Initial Guesses
   Lat0 = delta; %[rad]
   r_delta = r_delta_sat; %[km]
   r_K = R_ECI(3); %[km]
   iterations = 0;
   eps = 1;
   
   %Iterating Geodetic Latitude
   while (eps > 1e-6)
      
      CE = RE/sqrt(1 - ((eE^2)*((sin(Lat0))^2))); 
      Lat = atan((r_K + (CE*(eE^2)*sin(Lat0)))/r_delta);
      eps = abs(Lat - Lat0);
      Lat0 = Lat;
      iterations = iterations + 1;
      
      if (iterations > 1000)
         break
      end
      
   end
      
   h_ellip = (r_delta/(cos(Lat))) - CE; %[km]
   Lat = rad2deg(Lat); %[deg]
   Long = rad2deg(Long); %[deg]
   
end