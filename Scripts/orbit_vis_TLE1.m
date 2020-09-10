
clearvars; clc;

% ORBIT VISUALIZATION FROM TWO-LINE-ELEMENT SET
% ---------------------------------------------------------------------------------------------------------------------
% This script visualizes a Geocentric satellite orbits from Two-Line Element Set (TLE) at a specified Date and Time 
% (UTC). Two plots are generated: a representation of the orbit in 3D and a groundtrack of the orbit over the surface 
% of the Earth. The script can plot 1 revolution behind, and up to 6 revolutions ahead of the specified Date/Time. The 
% satellite state at the TLE epoch is computed from the mean elements in the TLE and propagated to the user specified 
% date using a 2-body propagation routine.
% 
% Due to the use of mean elements in the TLE and the inaccuracy of 2 body propagation over longer periods of time, this
% script can provide a qualitative representation of a satellite orbit close to the epoch for which the TLE was 
% generated, but should not be used to predict the future state of a satellite with great accuracy.
%
% Author: Matthew Buckhout
% Updated: 08/18/2020 
%
% Inputs:
%
%     - [DATE]       Date                                            [yyyy mm dd]
%     - [UTC]        Time (UTC)                                      [hh mm ss]
%     - [numsat]     Vector of Identifiers for TLEs to plot           - 
%                       ex: [1 3 2 4 8 9] 
%     - [now]        Specified Date/Time Selection                    -
%                       1 - Use current Date/Time in your timezone
%                       0 - Use specified Date/Time in UTC
%     - [dst]        Daylight Savings Time                            -
%                       1 - In effect
%                       0 - Not in effect
%
% Functions:
%
%     - getdata_EOP_fast
%     - julian_date
%     - gregorian_date
%     - kepler_equation
%     - true_anomaly
%     - DOYmodhms
%     - position_velocity  
%     - propagation_UV 
%     - EOP
%     - orbit_vis_2
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 113-116)
%     - Celestrak, Two-Line Element Set Format (Kelso) 
%           https://celestrak.com/columns/v04n05/
%     - 3D Earth Example (Ryan Gray - 2004,2006,2013)
%           https://www.mathworks.com/matlabcentral/fileexchange/13823-3d-earth-example
%
% -------------------------------------  U S E R  S P E C I F I E D  I N P U T S  -------------------------------------

DATE   = [2020 09 09]; %[yyyy mm dd] Specified Date
UTC    = [08 52 04]; %[hh mm ss] Specified Time UTC
numorb = 1; %Number of orbits to plot ahead of specified location
now    = 1; %Specified Date/Time Selection
dst    = 1; %Daylight Savings Time

%Satellite State - TLE
line1  = ['1 07530U 74089B   20252.86998931 -.00000027  00000-0  11486-3 0  9990'];
line2  = ['2 07530 101.8087 222.1707 0011908 202.2359 211.6166 12.53644948 96683'];

% -------------------------------------------  I N P U T  H A N D L I N G  --------------------------------------------

% Earth Orientation Parameters 
EOPdata = getdata_EOP_fast;

%Current Date and Time & Daylight Savings Time
if (dst == 1)
   addh = 5; %Central Time
else
   addh = 6; %Central Time
end
if (now == 1)
   
   c = clock;
   DATE = [c(1) c(2) c(3)];
   UTC = [(c(4)+addh) c(5) c(6)];
   
end

%Time Rollover - UTC
[~,MJD_UTC] = julian_date(DATE,UTC,0); %Modified Julian Date (UTC)
[yr,mo,d,h,m,s] = gregorian_date(MJD_UTC); %Gregorian Date/Time (UTC)
DATE = [yr mo d];
UTC  = [h m s];

% -----------------------------------------  C U R R E N T  P O S I T I O N  ------------------------------------------

%Constants
RE = 6378.1363; %[km] Earth Mean Equatorial Radius 
mu = 3.986004415e5; %[km^3/s^2]
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

ap = a*(1 - e) - RE; %[km] Periapsis Altitude
aa = a*(1 + e) - RE; %[km] Apoapsis Altitude

%Epoch Imformation
yr_ep = str2num(line1(19:20)); %[20yy] Last two digits of year 
DOY = str2num(line1(21:32)); % Day of Year
[mo_ep,d_ep,h_ep,m_ep,s_ep] = DOY_modhms(yr_ep,DOY); %Month,Day,Hour,Min,Sec
UTCep  = [h_ep m_ep s_ep];

%Position/Velocity at Epoch
t0 = (d_ep*86400) + (h_ep*3600) + (m_ep*60) + s_ep;
[Rep,Vep] = position_velocity(mu,p,e,i,OMEGA,omega,theta,0,0,0);

%Postion/Velocity at Selected Time
t = (DATE(3)*86400) + (UTC(1)*3600) + (UTC(2)*60) + UTC(3);
[R,V] = propagation_UV(mu,Rep,Vep,t0,t); 

fprintf('\n');

if (now == 1)
   fprintf('TLE Epoch      =  %02d/%02d/20%d | %02d:%02d:%02d UTC\n',mo_ep,d_ep,yr_ep,h_ep,m_ep,fix(s_ep));
   fprintf('Current Epoch  =  %02d/%02d/%d | %02d:%02d:%02d UTC\n',DATE(2),DATE(3),DATE(1),UTC(1),UTC(2),fix(UTC(3)));
   fprintf('Current Epoch  =  %02d/%02d/%d | %02d:%02d:%02d UTC-%d\n',DATE(2),DATE(3),DATE(1),c(4),c(5),fix(c(6)),addh);
elseif (now == 0)
   fprintf('TLE Epoch       =  %02d/%02d/20%d | %02d:%02d:%02d UTC\n',mo_ep,d_ep,yr_ep,h_ep,m_ep,fix(s_ep));
   fprintf('Selected Epoch  =  %02d/%02d/%d | %02d:%02d:%02d UTC\n',DATE(2),DATE(3),DATE(1),UTC(1),UTC(2),fix(UTC(3)));
end   

fprintf('\n');
fprintf('TLE SET\n');
fprintf('---------------------------------------------------------------------\n');
fprintf('%s\n',line1(1:end));
fprintf('%s\n',line2(1:end));
fprintf('---------------------------------------------------------------------\n');
fprintf('\n');

fprintf('ORBITAL ELEMENTS - GEOCENTRIC EQUATORIAL\n');
fprintf('---------------------------------------------\n');
fprintf('Semi-Major Axis          =  %0.3f km\n',a);
fprintf('Eccentricity             =  %0.6f \n',e);
fprintf('Inclination              =  %0.3f degrees\n',rad2deg(i));
fprintf('RA of Ascending Node     =  %0.3f degrees\n',rad2deg(OMEGA));
fprintf('Argument of Periapsis    =  %0.3f degrees\n',rad2deg(omega));
fprintf('True Anomaly at Epoch    =  %0.3f degrees\n',rad2deg(theta));
fprintf('---------------------------------------------\n');
fprintf('\n');

%Time System Parameters at Selected Date
[~,MJD] = julian_date(DATE,UTC,0);
[~,~,dUT1,dAT] = EOP(fix(MJD),EOPdata); 

%Orbit Plots
orbit_vis_2(R,V,DATE,UTC,numorb);
