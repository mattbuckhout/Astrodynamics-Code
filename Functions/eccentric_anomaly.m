
% ECCENTRIC ANOMALY - ALL ORBIT TYPES
% --------------------------------------------------------------------------------------------------
% Function generates Eccentric Anomaly, Hyperbolic Eccentric Anomaly, or Parabolic Eccentric Anomaly 
% using given Eccentricity and True Anomaly
%
% Author: Matthew Buckhout
% Updated: 08/06/2020 
%
% Inputs 
%                  
%     - [e]                Eccentricity                                     - 
%     - [theta]            True Anomaly                                  [rad]
%
% Outputs
%
%     - [EBH]              Eccentric Anomaly                             [rad]
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 51-66, 85)
% --------------------------------------------------------------------------------------------------

function [EBH] = eccentric_anomaly(e,theta)

   if (e < 1)
      
      %Eccentric Anomaly
      num = (sin(theta)*sqrt(1 - (e^2)))/(1 + e*cos(theta));
      den = (e + cos(theta))/(1 + e*cos(theta));
      EBH = atan2(num,den);
      
      if (EBH < 0)
         EBH = EBH + 2*pi;
      end

   elseif (e == 1)
      
      %Parabolic Eccentric Anomaly
      EBH = tan(theta/2);
      
   elseif (e > 1)
      
      %Hyperbolic Eccentric Anomaly
      if (theta >= 0 && theta < pi)
         EBH = asinh((sqrt((e^2)-1)*sin(theta))/(1 + e*cos(theta))); 
      elseif (theta >= pi && theta < 2*pi)
         EBH = -asinh((sqrt((e^2)-1)*sin(theta))/(1 + e*cos(theta))); 
      end
      
   end

end


