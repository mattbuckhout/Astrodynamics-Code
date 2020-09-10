
clearvars; clc;

% PRELIMINARY INTERPLANETARY TRAJECTORY
% ---------------------------------------------------------------------------------------------------------------------
% Script determines a heliocentric transfer orbit between two planetary bodies for a specified launch date and arrival 
% date using 2-body dynamics and solutions to Lambert's Problem. The positions of the start body and the target
% body are determined using the JPL DE430 model for planetary ephemerides. The script calculates the delta-V needed
% to leave the starting planet's heliocentric orbit and be placed on the transfer orbit (not including 
% the escape velocity from the starting planet), as well as the velocity in the transfer orbit relative to the target
% body at the arrival date. The propagation and solutions to Lambert's problem are done in a Heliocentric frame, 
% and all vectors are transformed into the standard International Celestial Reference Frame (ICRF) (barycentric) for 
% plotting and output. Inputs can be entered in the "USER SPECIFIED INPUTS" section below the header.
%
% Author: Matthew Buckhout
% Updated: 09/09/2020 
%
% Inputs:
%
%     - [launch_body]      Identifier for Launch Planetary Body                   -
%     - [target_body]      Identifier for Target Planetary Body                   -
%     - [launch_date]      Modified Julian Date of Transfer Orbit Start (UT1)    [days]
%     - [arrival_date]     Modified Julian Date of Arrival (UT1)                 [days]
%     - [tstep]            Time Step for plot of transfer orbit                  [days]
%
% Identifiers for Launch or Target Body:
%
%     1. Mercury
%     2. Venus
%     3. Earth-Moon Barycenter
%     4. Mars
%     5. Jupiter
%     6. Saturn
%     7. Uranus
%     8. Neptune
%     9. Pluto
%
% Data files required:
%
%     DE430t_1550_E.txt
%     DE430t_1650_E.txt
%     DE430t_1750_E.txt
%     DE430t_1850_E.txt
%     DE430t_1950_E.txt
%     DE430t_2050_E.txt
%     DE430t_2150_E.txt
%     DE430t_2250_E.txt
%     DE430t_2350_E.txt
%     DE430t_2450_E.txt
%     DE430t_2550_E.txt
%     IERS_EOP_combined.txt
%
% Functions:
%
%     - getdata_EOP_fast
%     - terrestrial_MJD  
%     - getdata_DE430
%     - planet_ephemerides_full_c
%     - lambert_c
%     - propagator_c
%
% References:
%
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 448-487)
%     - IERS Conventions 2010 (pg. 18)
%           https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html
%
% ----------------------------------------------- USER SPECIFIED INPUTS -----------------------------------------------

launch_body  = 3; %Identifier for Launch Body
target_body  = 4; %Identifier for Target Body

launch_date  = 59049; %[days] Modified Julian Date at Departure
arrival_date = 59256; %[days] Modified Julian Date at Arrival

tstep        = 1; %[days] Time Step for Plot

% ----------------------------------------------- GENERATE EPHEMERIDES ------------------------------------------------

MJD1 = launch_date; %[days] Modified Julian Date at Transfer Burn
MJD2 = arrival_date; %[days] Modified Julian Date at Rendevous
TOF = MJD2 - MJD1; %[days] Time of Flight

EOPdata = getdata_EOP_fast; %Earth Orientation and Time System Data
[MJD1_TT] = terrestrial_MJD(MJD1,EOPdata); %Modified Julian Date expressed in Terrestrial Time (TT)
[MJD2_TT] = terrestrial_MJD(MJD2,EOPdata); %Modified Julian Date expressed in Terrestrial Time (TT)
[DE430coef] = getdata_DE430(MJD1_TT,MJD2_TT); %Coefficients for Planetary Ephemerides 
                                                %(Can use MJD UT1 for slightly less accuracy, difference is small)

[ephemerides] = planet_ephemerides_full_c(MJD1,transpose(DE430coef),[launch_body target_body 11]); %Ephemerides (ICRF)
sun_ephem1 = transpose(ephemerides(3,1:6)); %Sun ephemerides (ICRF)
Rtgt1 = transpose(ephemerides(2,1:3)) - sun_ephem1(1:3); %[km] (HELIOCENTRIC)
Rint1 = transpose(ephemerides(1,1:3)) - sun_ephem1(1:3); %[km] (HELIOCENTRIC)
Vtgt1 = transpose(ephemerides(2,4:6)) - sun_ephem1(4:6); %[km/s] (HELIOCENTRIC)
Vint1 = transpose(ephemerides(1,4:6)) - sun_ephem1(4:6); %[km/s] (HELIOCENTRIC)

[ephemerides] = planet_ephemerides_full_c(MJD2,transpose(DE430coef),[launch_body target_body 11]); %Ephemerides (ICRF)
sun_ephem2 = transpose(ephemerides(3,1:6)); %Sun ephemerides (ICRF)
Rtgt2 = transpose(ephemerides(2,1:3)) - sun_ephem2(1:3); %[km] (HELIOCENTRIC)
Vtgt2 = transpose(ephemerides(2,4:6)) - sun_ephem2(4:6); %[km/s] (HELIOCENTRIC)

% ------------------------------------------------ TRAJECTORY SOLVING -------------------------------------------------

mu = 1.32712428e11; %[km^3/s^2] Sun Gravitational Parameter

%Selecting Transfer Method for Minimum Delta-V
tran_n = cross(Rint1,Rtgt2); %Transfer orbit normal vector
h_n = cross(Rint1,Vint1); %Interceptor S/C orbit normal vector
   
if (dot(h_n,tran_n) < 0) 
   tm = -1; %Long Way
elseif (dot(h_n,tran_n) > 0)
   tm = 1; %Short Way
end
   
[VS] = lambert_c(mu,Rint1,Rtgt2,(TOF*24*3600),tm); %Intercepting S/C Velocity Vectors on Transfer Orbit (HELIOCENTRIC)
Vtran1 = VS(:,1); %[km/s] Velocity on Transfer Orbit at Departure Date (HELIOCENTRIC)
Vtran2 = VS(:,2); %[km/s] Velocity on Transfer Orbit at Arrival Date (HELIOCENTRIC)
Rtran1 = Rint1; %[km] Intercepting S/C Position (HELIOCENTRIC)

dV1 = (Vtran1 - Vint1) + sun_ephem1(4:6); %Delta-V Vector Transfer Burn (ICRF)
dV2 = (Vtgt2 - Vtran2) + sun_ephem2(4:6); %Delta-V Vector Rendevous Burn (ICRF)

j = 1;
for (MJDstep=MJD1:tstep:MJD2) %Loop for Initial, Target, Transfer Orbits

   [ephemerides] = planet_ephemerides_full_c(MJDstep,transpose(DE430coef),[launch_body target_body 11]); %(ICRF)
   
   sun_ephem = ephemerides(3,1:6); %Sun ephemerides (ICRF)
   RLAUNCH(j,:) = ephemerides(1,1:3); %[km] Launch Body Positions (ICRF)
   RTARGET(j,:) = ephemerides(2,1:3); %[km] Target Body Positions (ICRF) 
   RV = propagator_c(mu,Rtran1,Vtran1,MJD1,MJDstep); %[km] Transfer Orbit Positions (HELIOCENTRIC)
   RTRAN(j,1:3) = RV(1:3) + sun_ephem(1:3); %[km] Transfer Orbit Positions (ICRF)
   VTRAN(j,1:3) = RV(4:6) + sun_ephem(4:6); %[km/s] Transfer Orbit Velocities (ICRF)
 
   j = j+1;

end

bodyname = {'Mercury';'Venus';'Earth';'Mars';'Jupiter';'Saturn';'Uranus';'Neptune';'Pluto'};

fprintf('\n');
fprintf('==================== TRANSFER ORBIT INFO ====================\n');
fprintf(' Launch Body: %s\n',bodyname{launch_body});
fprintf(' Target Body: %s\n',bodyname{target_body});
fprintf('\n');
fprintf(' Launch Date     = %d UT1\n',MJD1);
fprintf(' Arrival Date    = %d UT1\n',MJD2);
fprintf(' Time of Flight  = %7.3f days\n',TOF);
fprintf(' Transfer dV     = %7.3f km/s\n',norm(dV1));
fprintf(' Transfer C3     = %7.3f km^2/s^2\n',norm(dV1)^2);
fprintf(' Target Vrel.    = %7.3f km/s\n',norm(dV2));
fprintf('\n');
fprintf('=============================================================\n');
fprintf(' Initial State (ICRF)            Final State (ICRF)\n');
fprintf(' %d UT1                       %d UT1\n',MJD1,MJD2);
fprintf('    Rx = %13.6e km           Rx = %13.6e km\n',RTRAN(1,1),RTRAN(end,1));
fprintf('    Ry = %13.6e km           Ry = %13.6e km\n',RTRAN(1,2),RTRAN(end,2));
fprintf('    Rz = %13.6e km           Rz = %13.6e km\n',RTRAN(1,3),RTRAN(end,3));
fprintf('    Vx = %13.3f km/s         Vx = %13.3f km/s\n',VTRAN(1,1),VTRAN(end,1));
fprintf('    Vy = %13.3f km/s         Vy = %13.3f km/s\n',VTRAN(1,2),VTRAN(end,2));
fprintf('    Vz = %13.3f km/s         Vz = %13.3f km/s\n',VTRAN(1,3),VTRAN(end,3));
fprintf('\n');
fprintf('=============================================================\n');
fprintf(' Transfer Delta-V Vector (ICRF)\n');
fprintf('    Vx = %6.3f km/s \n',dV1(1));             
fprintf('    Vy = %6.3f km/s\n',dV1(2));                
fprintf('    Vz = %6.3f km/s\n',dV1(3));                 
fprintf('\n');
fprintf(' Target Relative Velocity Vector @ Arrival (ICRF)\n'); 
fprintf('    Vx = %6.3f km/s\n',dV2(1));
fprintf('    Vy = %6.3f km/s\n',dV2(2));
fprintf('    Vz = %6.3f km/s\n',dV2(3));
fprintf('=============================================================\n');

% ------------------------------------------------ TRANSFER ORBIT PLOT ------------------------------------------------

AU = 1.49597870700e8; %[km] 1 Astronomical Unit
RLAUNCH_AU = RLAUNCH/AU; %[AU] Launch Body Positions (ICRF)
RTARGET_AU = RTARGET/AU; %[AU] Target Body Positions (ICRF)
RTRAN_AU = RTRAN/AU; %[AU] Transfer Orbit Positions (ICRF)

set(groot,'defaultfigureposition',[700 150 800 700])
figure
p1 = plot3(0,0,0,'ko','MarkerFaceColor','k','MarkerSize',8);
hold on
p2 = plot3(RTARGET_AU(1,1),RTARGET_AU(1,2),RTARGET_AU(1,3),'rd');
hold on
p3 = plot3(RTARGET_AU(end,1),RTARGET_AU(end,2),RTARGET_AU(end,3),'rd','MarkerFaceColor','r');
hold on
p4 = plot3(RLAUNCH_AU(1,1),RLAUNCH_AU(1,2),RLAUNCH_AU(1,3),'bd');
hold on
p5 = plot3(RLAUNCH_AU(end,1),RLAUNCH_AU(end,2),RLAUNCH_AU(end,3),'bd','MarkerFaceColor','b');

hold on
p6 = plot3(RTRAN_AU(:,1),RTRAN_AU(:,2),RTRAN_AU(:,3),'k');
hold on
p7 = plot3(RTARGET_AU(:,1),RTARGET_AU(:,2),RTARGET_AU(:,3),'r');
hold on
p8 = plot3(RLAUNCH_AU(:,1),RLAUNCH_AU(:,2),RLAUNCH_AU(:,3),'b');

axis equal
grid on
set(gca,'FontSize',12);

xlabel('X ICRF (AU)');
ylabel('Y ICRF (AU)');
zlabel('Z ICRF (AU)');

lgd = legend([p1 p4 p2 p5 p3 p6 p7 p8], ...
      sprintf('SSBC'), ...
      sprintf('Launch Body - %d',launch_date), ...
      sprintf('Target Body  - %d',launch_date), ...
      sprintf('Launch Body - %d',arrival_date), ...
      sprintf('Target Body  - %d',arrival_date), ...
      sprintf('Transfer Orbit'), ...
      sprintf('Target Body Orbit'), ...
      sprintf('Launch Body Orbit'));
      
set(lgd,'FontSize',12);     
