
% CONVERT UTC TO TIME ZONE
% -------------------------------------------------------------------------------------------------
% Function takes Date and Time arrays in UTC and converts to Date and Time in specified time zone 
% offset from UTC.
%
% Author: Matthew Buckhout
% Updated: 08/15/2020 
%
% Inputs:
%
%     - [DATE]       Gregorian Date (UTC)                     [yyyy mm dd]
%     - [UTC]        Time (UTC)                               [hh mm ss]
%     - [shift]      Hours Timezone is offset from UTC        [hours]
%
% Outputs:
%
%     - [DATE]       Gregorian Date (Time zone)               [yyyy mm dd]
%     - [TZ]         Time (Time zone)                         [hh mm ss]
%
% Functions:
%
%     - julian_date
%     - gregorian_date  
% -------------------------------------------------------------------------------------------------

function [DATE,TZ] = UTC2TZ(DATE,UTC,shift)

   [~,MJD_UTC] = julian_date(DATE,UTC,0); %Modified Julian Date (UTC)
   
   MJD_TZ = MJD_UTC + shift/24; %Modified Julian Date (TZ)
   
   [yr,mo,d,h,m,s] = gregorian_date(MJD_TZ); %Date and Time (TZ)

   DATE = [yr mo d];
   TZ = [h m s];

end
