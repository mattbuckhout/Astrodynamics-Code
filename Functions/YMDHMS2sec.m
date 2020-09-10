
% YEAR, MONTH, DAY, HOUR, MINUTE, SECOND TO SECONDS ELAPSED
% -------------------------------------------------------------------------------------------------
% Takes given Date and Time (UTC) and converts to seconds elapsed since start of that year (January 
% 1, 00:00)
%
% Author: Matthew Buckhout
% Updated: 08/15/2020 
%
% Inputs:
%
%     - [DATE]       Gregorian Date (UTC)                     [yyyy mm dd]
%     - [UTC]        Time (UTC)                               [hh mm ss]
%
% Outputs:
%
%     - [seconds]    Seconds Elapsed in given year            [sec]
% -------------------------------------------------------------------------------------------------

function [seconds] = YMDHMS2sec(DATE,TIME)
   
   yr = DATE(1);
   mo = DATE(2);
   d = DATE(3); 
   h = TIME(1); 
   m = TIME(2);
   s = TIME(3);

   length_month = [31 28 31 30 31 30 31 31 30 31 30 31]; %Number of days in each month

   check = (fix(yr/4)/(yr/4)); %Leap year check
   if (check == 1)
      length_month(2) = 29;
   end

   %Adding all days up to current day
   days = 0;
   for i=1:1:mo-1
      
      days = days + length_month(i);

   end

   seconds = ((days + d)*24*3600) + (h*3600) + (m*60) + s;
      
end