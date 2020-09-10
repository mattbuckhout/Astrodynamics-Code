
% CONVERT MODIFIED JULIAN DATE BETWEEN UNIVERSAL TIME (UT1) AND TERRESTRIAL TIME (TT)
% -------------------------------------------------------------------------------------------------
% Returns the Modified Julian Date based in Terrestrial Time (TT) from input Modified Julian Date 
% based in Universal Time (UT1). Designed for quick conversion and use with planetary ephemerides 
% function using DE430 model (tables of coefficients in model are bracketed using MJD expressed in 
% TDT/TT)
%
% Author: Matthew Buckhout
% Updated: 08/16/2020 
%
% Inputs:
%
%     - [MJD_UT1]    Modified Julian Date (UT1)                      [days]
%     - [EOPdata]    Earth Orientation Parameter Table                -
%
% Outputs:
%
%     - [MJD_TT]     Modified Julian Date (TT)                       [days]
%
% Functions:
%
%     - gregorian_date
%     - EOP  
%     - julian_date 
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 194-197)
% -------------------------------------------------------------------------------------------------

function [MJD_TT] = terrestrial_MJD(MJD_UT1,EOPdata)

   [yr,mo,d,h,m,s] = gregorian_date(MJD_UT1); %Date and Time (UT1)

   UT1s = h*3600 + m*60 + s; %[sec] UT1 in seconds
   
   [~,~,dUT1,dAT] = EOP(MJD_UT1,EOPdata); %[sec] Differences between UT1, UTC, TAI

   TTs = ((UT1s - dUT1) + dAT) + 32.184; %[sec] Terrestrial Time in seconds
   
   [~,MJD_TT] = julian_date([yr mo d],[0 0 TTs],0); %Modified Julian Date (TT)

end
