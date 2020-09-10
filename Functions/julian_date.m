
% JULIAN DATE CONVERTER
% -------------------------------------------------------------------------------------------------
% Converts Date and Time into Julian Date and Modified Julian Date with respect to the inputted
% time system (UTC, UT1, etc...), measured from January 1, 4713 BC. Applicable for dates between 
% March 1, 1900 and February 28, 2100. For days during which a leap second is added, "leap" = 1
%
% Author: Matthew Buckhout
% Updated: 08/08/2020 
%
% Inputs:
%
%     - [DATE]       Gregorian Date                               [yyyy mm dd]
%     - [TIME]       Time                                         [hh mm ss]
%     - [leap]       Specifies if leap second added for date       - 
%                    of interest (1) or not (0)
%
% Outputs:
%
%     - [JD]         Julian Date                                  [days]
%     - [MJD]        Modified Julian Date                         [days]
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 186-187)
% -------------------------------------------------------------------------------------------------

function [JD,MJD] = julian_date(DATE,TIME,leap)
   
   yr = DATE(1); %Year
   mo = DATE(2); %Month
   d = DATE(3); %Day
   h = TIME(1); %Hour
   m = TIME(2); %Minute
   s = TIME(3); %Second
   
   %Day does NOT contain leap second
   if (leap == 0)
      
      JD = 367*yr - fix((7*(yr + fix((mo+9)/12)))/4) ... 
           + fix((275*mo)/9) + d + 1721013.5 ... 
           + (((((s/60) + m)/60) + h)/24);
   
   %Day does contain leap second
   elseif (leap == 1)
   
      JD = 367*yr - fix((7*(yr + fix((mo+9)/12)))/4) ... 
           + fix((275*mo)/9) + d + 1721013.5 ... 
           + (((((s/61) + m)/61) + h)/24);
   
   end
   
   %Modified Julian Date
   MJD = JD - 2400000.5;
   
end






