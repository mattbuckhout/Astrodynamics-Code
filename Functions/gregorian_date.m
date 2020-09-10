
% GREGORIAN DATE CONVERTER
% -------------------------------------------------------------------------------------------------
% Converts from Modified Julian Date to Gregorian Date [yyyy mm dd] and Time [hh mm ss]
%
% Author: Matthew Buckhout
% Updated: 08/08/2020 
%
% Inputs:
%
%     - [MJD]        Modified Julian Date       [days]
%
% Outputs:
%
%     - [yr]         Year                       [yyyy]
%     - [mo]         Month                       -
%     - [d]          Day                         -
%     - [h]          Hour                        -
%     - [m]          Minute                      -
%     - [s]          Second                      -
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 203-204)
% -------------------------------------------------------------------------------------------------

function [yr,mo,d,h,m,s] = gregorian_date(MJD)

   length_month = [31 28 31 30 31 30 31 31 30 31 30 31]; %Number of days in each month

   T_1900 = (MJD - 15019)/365.25; %Julian Centuries Epoch 1900 (UT1)

   yr = 1900 + fix(T_1900);
   leapyrs = fix((yr - 1900 - 1)*0.25);
   days = (MJD - 15019) - ((yr - 1900)*365.0 + leapyrs);

   if (days < 1.0)
     
      yr = yr - 1;
      leapyrs = fix((yr - 1900 -1)*0.25);
      days = (MJD - 15019) - ((yr - 1900)*365.0 + leapyrs);
     
   end

   %Leap Year Check
   check = (fix(yr/4)/(yr/4)); 
   if (check == 1)
      length_month(2) = 29;
   end

   DOY = fix(days);

   %Finding Correct Month
   sum1 = 0;
   sum2 = 0;
   ind = 0;
   while (sum1 < DOY) 
      
      sum1 = sum1 + length_month(ind + 1);
      
      if (ind > 0)
         sum2 = sum2 + length_month(ind);
      end
      
      ind = ind + 1;   

   end   

   mo  = ind;
   d   = DOY - sum2;
   tau = (days - DOY)*24;
   h   = fix(tau);
   m   = fix((tau - h)*60);
   s   = (tau - h - (m/60))*3600;

end


