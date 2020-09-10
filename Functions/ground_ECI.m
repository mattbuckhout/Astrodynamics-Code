
% CONVERTING LATITUDE/LONGITUDE/ALTITUDE TO POSITION VECTOR
% -------------------------------------------------------------------------------------------------
% Takes Geodetic Latitude, Longitude (+East of Greenwich), and Altitude above Mean Sea Level on 
% Earth Reference Ellipsoid and returns the position vector for satellite in Geocentric-Equatorial 
% (ECI) frame. Low fidelity, no reduction to account for polar motion or Precession-Nutation, only
% accounts for the rotation of the Earth when translating into ECI frame.
%
% Author: Matthew Buckhout
% Updated: 09/09/2020 
%
% Inputs:
%                                                   
%     - [GMST]          Greenwich Mean Sidereal Time                   [rad]
%     - [Lat]           Geodetic Latitude                              [deg]
%     - [Long]          Longitude (+East from Greenwich)               [deg]
%     - [h_ellip]       Height above reference ellipsoid               [km]         
%
% Outputs:
%
%     - [R_ECI]         Satellite Position Vector                      [km]
%                       Earth Centered Inertial frame (GCRF)
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 138, 179-197) 
%     - IERS Conventions 2010 (pg. 19, 40)
%           https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html
% -------------------------------------------------------------------------------------------------

function R_ECI = ground_ECI(Lat,Long,h_ellip,GMST)

   %Earth Reference Ellipsoid - GRS80
   RE = 6378.1363; %[km] Earth Semi-Major Axis 
   fE = 1/298.257223563; % Earth Flattening 
   bE = RE - (fE*RE); %[km] Earth Semi-minor Axis (Polar Radius)
   eE = sqrt(RE^2 - bE^2)/RE; % Earth Reference Ellipsoid Eccentricity
   
   %Unit Conversions
   Lat = deg2rad(Lat); %[rad]
   Long = deg2rad(Long); %[rad]
   
   %Auxiliary Quantaties
   CE = RE/sqrt(1 - ((eE^2)*((sin(Lat))^2))); 
   SE = (RE*(1 - (eE^2)))/sqrt(1 - ((eE^2)*((sin(Lat))^2)));

   %Position Vector in ECEF
   R_ECEF = [(CE + h_ellip)*cos(Lat)*cos(Long); ...
                 (CE + h_ellip)*cos(Lat)*sin(Long); ...
                 (SE + h_ellip)*sin(Lat)]; %[km] Position Vector (ECEF)
   
   %Position Vector in ECI   
   r  = norm(R_ECEF); %[km] Magnitude Position Vector (ECEF/ECI)
   rk = (h_ellip + SE)*sin(Lat); %[km] Position Vector (K)
   rj = sqrt(((r^2) - (rk^2))/( (cos(Long + GMST)/sin(Long + GMST))^2 + 1 )); %[km] Position Vector (J)
   ri = rj*(cos(Long + GMST)/sin(Long + GMST)); %[km] Position Vector (I)
   R_ECI = [ri; rj; rk];

end
