
% CONVERT DAY OF THE YEAR TO MONTH,DAY,HOUR,MINUTE,SECOND
% -------------------------------------------------------------------------------------------------
% Converts day of the year in current year to month, day, hour, minute, second in that year 
% counting from 00:00:00 January 1. Year input can be four digit [yyyy] or last two digits [yy] 
%
% Author: Matthew Buckhout
% Updated: 08/06/2020 
%
% Inputs                   
%
%     - [yr]         Year                       [yyyy]
%     - [DOY]        Day of Year                [days]
%
% Outputs
%
%     - [mo]         Month                          -
%     - [d]          Day                            -
%     - [h]          Hour                           -
%     - [m]          Minute                         -
%     - [s]          Second                         -   
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 179-197)
% -------------------------------------------------------------------------------------------------

function [mo,d,h,m,s] = DOY_modhms(yr,DOY)

   %Constants
   RE = 6378.1363; %[km] Earth Mean Equatorial Radius 
   mu = 3.986004415e5; %[km^3/s^2] Earth Gravitational Parameter
   length_month = [31 28 31 30 31 30 31 31 30 31 30 31]; %Number of days in each month

   %Leap Year Check
   check = (fix(yr/4)/(yr/4)); 
   if (check == 1)
      length_month(2) = 29;
   end

   %Finding Month
   fulldays = fix(DOY);
   mo = 0;
   daycount = 0;
   while (daycount < fulldays)
      
      mo = mo + 1;   
      daycount = daycount + length_month(mo);
      
   end

   d = length_month(mo) - (daycount-fulldays);
   h = (DOY - fix(DOY))*24;
   m = (h - fix(h))*60;
   s = (m - fix(m))*60;
   
   h = fix(h);
   m = fix(m);
   
end







