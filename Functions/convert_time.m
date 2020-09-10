
% TIME SYSTEM CONVERTER
% -------------------------------------------------------------------------------------------------
% Converts inputs of Date and Time (UTC) to four other time systems: Universal Time (UT1), 
% International Atomic Time (TAI), Terrestrial Time (TT), and Barycentric Dynamical Time (TDB). 
% Function also computes Julian Centuries of each time system.
%
% Author: Matthew Buckhout
% Updated: 08/06/2020 
%
% Inputs  
%
%     - [UTC]              Coordinated Universal Time               [hh mm ss]
%     - [DATE]             Gregorian Date                         [yyyy mm ss]
%     - [dUT1]             UT1 - UTC                                     [sec]
%     - [dTAI]             TAI - UTC                                     [sec]
%
% Outputs
%
%     - [UT1]              Universal Time                           [hh mm ss]
%     - [TAI]              International Atomic Time                [hh mm ss]
%     - [TT]               Terrestrial Time                         [hh mm ss]
%     - [TDB]              Barycentric Dynamical Time               [hh mm ss]
%     - [T_UT1]            Julian Centuries of UT1                          -
%     - [T_TT]             Julian Centuries of TT                           -
%     - [T_TDB]            Julian Centuries of TDB                          -
%
% Functions
%
%     - julian_date        Converts Date/Time to JD and MJD
%     - s2hms              Converts time in seconds to hours, minutes, seconds
%
% References
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 179-197)
% -------------------------------------------------------------------------------------------------

function [UT1,TAI,TT,TDB,T_UT1,T_TT,T_TDB] = convert_time(DATE,UTC,dUT1,dAT)
   
   yr = DATE(1); %Year 
   mo = DATE(2); %Month
   d = DATE(3); %Day
   
   UTC = UTC(1)*3600 + UTC(2)*60 + UTC(3); %[sec] UTC in seconds
   UT1 = UTC + dUT1; %[sec] Universal Time
   TAI = UTC + dAT; %[sec] International Atomic Time
   TT = TAI + 32.184; %[sec] Terrestrial Time
   
   %Universal Time
   [JD_UT1,~] = julian_date([yr mo d],[0 0 UT1],0); %Julian date of UT1
   T_UT1 = (JD_UT1 - 2451545)/36525; %[sec] Julian Centuries of UT1 from J2000
   
   %Terrestrial Time
   [JD_TT] = julian_date([yr mo d],[0 0 TT],0); %Julian date of TT
   T_TT = (JD_TT - 2451545)/36525; %Julian centuries of TT from J2000
   
   %Barycentric Dynamical Time
   ME = ((357.5277233) + (35999.05034*T_TT))*(pi/180); %[rad] Mean Anomaly of Earth's Orbit (approx)
   TDB = TT + (0.001658*sin(ME)) + (0.00001385*sin(2*ME)); %[sec] Barycentric Dynamical Time (approx)
   [JD_TDB] = julian_date([yr mo d],[0 0 TDB],0); %Julian date of TDB
   T_TDB = (JD_TDB - 2451545)/36525; %Julian centuries of TDB from J2000
   
   %Converting Seconds to [H M S] Vector
   [UT1] = s2hms(UT1);
   [TAI] = s2hms(TAI);
   [TT] = s2hms(TT);
   [TDB] = s2hms(TDB);
   
end

   
   
   