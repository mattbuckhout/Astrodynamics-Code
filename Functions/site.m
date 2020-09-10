
% OBSERVATION SITE POSITION VECTOR
% -------------------------------------------------------------------------------------------------
% Function takes Geodetic Latitude, Longitude, and Altitude for an observation site on the surface 
% of Earth and returns the Geocentric position vector of the site expressed in the International 
% Terrestrial Reference Frame (ITRF). Function uses parameters for the GRS80 reference ellipsoid.
%
% Author: Matthew Buckhout
% Updated: 09/09/2020 
%
% Inputs:
%
%     - [Lat]           Geodetic Latitude                               [deg]
%     - [Long]          Longitude (+East from Greenwich)                [deg]
%     - [h_ellip]       Height above reference ellipsoid                [km]
%
% Outputs:
%
%     - [Rsite_ITRF]    Position Vector (ITRF)                          [km]
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 137-144, 406-408)
%     - IERS Conventions 2010 (pg. 19, 40)
%           https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html
% -------------------------------------------------------------------------------------------------

function [Rsite_ITRF] = site(Lat,Long,h_ellip)

   %Earth Reference Ellipsoid - GRS80
   RE = 6378.1363; %[km] Earth Semi-Major Axis 
   fE = 1/298.257223563; % Earth Flattening 
   bE = RE - (fE*RE); %[km] Earth Semi-minor Axis (Polar Radius)
   eE = sqrt(RE^2 - bE^2)/RE; % Earth Reference Ellipsoid Eccentricity

   %Auxiliary Quantaties
   CE = RE/sqrt(1 - ((eE^2)*((sin(Lat))^2))); 
   SE = (RE*(1 - (eE^2)))/sqrt(1 - ((eE^2)*((sin(Lat))^2)));

   %Site Position Vector
   Rsite_ITRF = [(CE + h_ellip)*cos(Lat)*cos(Long); ...
                 (CE + h_ellip)*cos(Lat)*sin(Long); ...
                 (SE + h_ellip)*sin(Lat)]; %[km]  Site Position Vector (ITRF)

end

