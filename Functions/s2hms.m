
% SECONDS TO HOUR, MINUTE, SECOND
% -------------------------------------------------------------------------------------------------
% Function converts input value of seconds in to Time array of hours, minutes, and seconds
%
% Author: Matthew Buckhout
% Updated: 08/15/2020 
%
% Inputs:
%
%     - [seconds]    Seconds                         [sec]
%
% Outputs:
%
%     - [HMS]        Time Array                      [hh mm ss]
% -------------------------------------------------------------------------------------------------

function [HMS] = s2hms(seconds)
   
   h = fix(seconds/3600);
   m = fix((seconds - (h*3600))/60);
   s = seconds - (h*3600) - (m*60);
   HMS = [h m s];

end
