
% TWO LINE ELEMENT (TLE) SET TRANSLATOR
% -------------------------------------------------------------------------------------------------
% Takes standard NORAD Two-Line Element (TLE) set and converts to orbital elements, state vectors, 
% Epoch Date, and Time. For epoch info, assumes 21st century (years start with 20--).
% Satellite State extracted from TLE is approximate; the elements contained in the TLE are mean 
% elements in a true-equator, mean-equinox coordinate system calculated to fit many observations 
% using the SGP4 model. 
%
% Author: Matthew Buckhout
% Updated: 08/24/2020 
%
% Inputs:
%
%     - [line1]      String, First line of TLE                        -
%     - [line2]      String, Second line of TLE                       -
%
% Outputs:
%
%     - [p]          Semi-Parameter                                  [km]
%		- [a] 	      Semi-Major Axis                                 [km]
% 		- [e] 	      Eccentricity                                     -
% 		- [i]          Inclination                                     [rad]
% 		- [OMEGA]      Right Ascension of Ascending Node               [rad]
% 		- [omega]      Argument of Periapsis                           [rad]
% 		- [theta]      True Anomaly at Epoch                           [rad]
%     - [M]          Mean Anomaly                                    [rad]
%     - [R]          Position Vector                                 [km]
%     - [V]          Velocity Vector                                 [km/s]
%     - [DATE]       Gregorian Date (UTC)                            [yyyy mm dd]
%     - [TIME]       Time (UTC)                                      [hh mm ss]
%
% Functions:
%
%     - kepler_equation
%     - true_anomaly
%     - position_velocity
%     - DOY_modhms
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 113-116)
%     - Celestrak, Two-Line Element Set Format (Kelso) 
%           https://celestrak.com/columns/v04n05/
% -------------------------------------------------------------------------------------------------

function [p,a,e,i,OMEGA,omega,theta,M,R,V,DATE,TIME] = TLE2state(line1,line2)

%Constants
RE = 6378.1363; %[km] Earth Mean Equatorial Radius 
mu = 3.986004415e5; %[km^3/s^2]

%Mean Motion
n = str2num(line2(53:63)); %[rev/day] Mean Motion
n = (n*(2*pi)/86400); %[rad/s]

%Orbital Elements at Epoch
a       = (mu/(n^2))^(1/3); %[km] Semi-Major Axis
e       = str2num(line2(27:33))/(1e7); % Eccentricity
i       = deg2rad(str2num(line2(9:16))); %[rad] Inclination
omega   = deg2rad(str2num(line2(35:42))); %[rad] Argument of Perigee
OMEGA   = deg2rad(str2num(line2(18:25))); %[rad] Right Ascension of Ascending Node
M       = deg2rad(str2num(line2(44:51))); %[rad] Mean Anomaly at Epoch
[EBH]   = kepler_equation(mu,M,e,0,0); %[rad] Eccentric Anomaly
[theta] = true_anomaly(e,EBH,0,0); %[rad] True Anomaly

if ((e < 1) || (e == 0))
   p = a*(1-(e^2)); %[km] Semi-Parameter
elseif (e > 1)
   p = a*(1-(e^2)); 
elseif (e == 1)
   p = (h^2)/mu; 
end

%Position/Velocity Vectors
[R,V] = position_velocity(mu,p,e,i,OMEGA,omega,theta,0,0,0); %[km and km/s]

%TLE Epoch
yr = str2num(line1(19:20)); %[??yy] Last two digits of year 
yr = yr + 2000; %Assuming 21st century
DOY = str2num(line1(21:32)); %Day of Year
[mo,d,h,m,s] = DOY_modhms(yr,DOY); %Month,Day,Hour,Min,Sec
DATE = [yr mo d]; %Date at TLE epoch
TIME  = [h m s]; %Time UTC at TLE epoch

end



