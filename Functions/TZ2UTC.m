
% CONVERT TIME ZONE TO UTC
% -------------------------------------------------------------------------------------------------
% Function takes Date and Time arrays in given time zone offset from UTC and converts to Date and 
% Time in UTC.
%
% Author: Matthew Buckhout
% Updated: 08/15/2020 
%
% Inputs:
%
%     - [DATE]       Gregorian Date (Time zone)               [yyyy mm dd]
%     - [TZ]         Time (Time zone)                         [hh mm ss]
%     - [shift]      Hours Timezone is offset from UTC        [hours]
%
% Outputs:
%
%     - [DATE]       Gregorian Date (UTC)                     [yyyy mm dd]
%     - [UTC]        Time (UTC)                               [hh mm ss]
%
% Functions:
%
%     - julian_date
%     - gregorian_date  
% -------------------------------------------------------------------------------------------------

function [DATE,UTC] = TZ2UTC(DATE,TZ,shift)

   [~,MJD_TZ] = julian_date(DATE,TZ,0); %Modified Julian Date (TZ)
   
   MJD_UTC = MJD_TZ - shift/24; %Modified Julian Date (UTC)
   
   [yr,mo,d,h,m,s] = gregorian_date(MJD_UTC); %Date and Time (UTC)

   DATE = [yr mo d];
   UTC = [h m s];

end


