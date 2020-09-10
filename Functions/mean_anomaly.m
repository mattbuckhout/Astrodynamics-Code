
% MEAN ANOMALY - ALL ORBIT TYPES
% -------------------------------------------------------------------------------------------------
% Calculates the Mean Anomaly for all orbit types depending on inputted Eccentricity.
%
% Author: Matthew Buckhout
% Updated: 08/08/2020 
%
% Inputs:
%
%     - [e]          Eccentricity                                  -
%     - [EBH]        Eccentric Anomaly for given orbit type       [rad]
%
% Outputs:
%
%     - [M]          Mean Anomaly                                 [rad]
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) 
%     - Orbital Mechanics for Engineering Students, 3rd ed. (Curtis)
% -------------------------------------------------------------------------------------------------

function [M] = mean_anomaly(e,EBH)

   %Mean Anomaly
   if ((e < 1) || (e == 0))
      M = EBH - e*sin(EBH); %Elliptical/Circular
   elseif (e > 1) 
      M = e*sinh(EBH) - EBH; %Hyperbolic
   elseif (e == 1) 
      M = EBH + ((EBH^3)/3); %Parabolic
   end

end















