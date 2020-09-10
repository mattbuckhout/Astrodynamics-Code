
% LOCAL AND GREENWICH MEAN SIDEREAL TIME
% -------------------------------------------------------------------------------------------------
% Function returns the Greenwich Mean Sidereal Time (GMST) and Local Sidereal Time (LST) for a 
% given Modified Julian Date (UT1) and Longitude (+East from Greenwich). Longitude input of zero 
% will return GMST for both outputs.
%
% Author: Matthew Buckhout
% Updated: 08/08/2020 
%
% Inputs:
%
%     - [Long]       Longitude (+East from Greenwich)             [rad]
%     - [MJD]        Modified Julian Date (UT1)                   [days]
%
% Outputs:
%
%     - [LST]        Local Sidereal Time                          [rad]
%     - [GMST]       Greenwich Mean Sidereal Time                 [rad]
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 189-192)
% -------------------------------------------------------------------------------------------------

function [LST,GMST] = local_sidereal(Long,MJD)
   
   JD = MJD + 2400000.5;
   T_UT1 = (JD - 2451545)/36525; %          Julian Centuries from J2000
   GMST = 67310.54841 + (((876600*3600) + 8640184.812866)*T_UT1) + (0.093104*(T_UT1^2)) ...
          - ((6.2e-6)*(T_UT1^3)); %[sec]     Greenwich Mean Sidereal Time
   
   %Removing Extra Periods
   if (abs(GMST) >= 86400)
      if (GMST >= 0)
         
         GMST = GMST - (fix(abs(GMST)/86400)*86400);
         GMST = GMST*(15/3600)*(pi/180); %[rad]    Greenwich Mean Sidereal Time
      
      elseif (GMST < 0)
         
         GMST = GMST + (fix(abs(GMST)/86400)*86400);
         GMST = GMST*(15/3600)*(pi/180); %[rad]    Greenwich Mean Sidereal Time
         GMST = (2*pi) + GMST;
         
      end
   end

   if (GMST < 0)
      GMST = (2*pi) - abs(GMST);
   end
      
   LST = GMST + Long; %[rad]         Local Sidereal Time
   
   if (LST < 0)
      LST = (2*pi) - abs(LST);
   end
   
   
end
