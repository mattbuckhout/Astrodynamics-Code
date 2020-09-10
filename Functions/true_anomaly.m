
% TRUE ANOMALY - ALL ORBIT TYPES
% --------------------------------------------------------------------------------------------------
% Function calculates True Anomaly from Eccentric Anomaly, Hyperbolic Eccentric Anomaly, or 
% Parabolic Eccentric Anomaly along with the Eccentricity, Semi-Parameter, and corresponding Radius 
% magnitude. 
%
% Author: Matthew Buckhout
% Updated: 08/16/2020 
%
% Inputs 
%                  
%     - [e]             Eccentricity                                     - 
%     - [EBH]           Eccentric Anomaly                               [rad]
%     - [p]             Semi-Parameter (Only required for Parabolic)    [km]
%     - [r]             Radius Magnitude (Only required for Parabolic)  [km]  
%
% Outputs
%
%     - [theta]            True Anomaly                                 [rad]                 
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 51-66, 85)
% --------------------------------------------------------------------------------------------------

function [theta] = true_anomaly(e,EBH,p,r)

   if (e < 1)
      
      theta1 = (sin(EBH)*sqrt(1 - (e^2)))/(1 - e*cos(EBH));
      theta2 = (cos(EBH) - e)/(1 - e*cos(EBH));
      theta = atan2(theta1,theta2);
      
      if (theta < 0)
         theta = theta + (2*pi);
      end
      

   elseif (e == 1)

      theta1 = (p*EBH)/R;
      theta2 = (p - r)/r;
      theta = atan2(theta1,theta2);
      
      if (theta < 0)
         theta = theta + (2*pi);
      end
      
   elseif (e > 1)
      
      theta1 = (-sinh(EBH)*sqrt((e^2) - 1))/(1 - e*cosh(EBH));
      theta2 = (cosh(EBH) - e)/(1 - e*cosh(EBH));
      theta = atan2(theta1,theta2);
      
      if (theta < 0)
         theta = theta + (2*pi);
      end
      
   end

end





