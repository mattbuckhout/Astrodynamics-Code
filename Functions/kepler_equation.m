
% KEPLER'S EQUATION - ALL ORBIT TYPES
% -------------------------------------------------------------------------------------------------
% Function solves Kepler's Equation using Newton-Raphson method and computes either the Eccentric 
% Anomaly, Hyperbolic Eccentric Anomaly, or Parabolic Eccentric Anomaly depending on the orbit type.
%
% Author: Matthew Buckhout
% Updated: 08/08/2020 
%
% Inputs:
%
%     - [mu]         Gravitational Parameter                         [km^3/s^2]
%     - [M]          Mean Anomaly                                    [rad]
%     - [e]          Eccentricity                                     -
%     - [dt]         Change in Time (Only required for Parabolic)    [sec]
%     - [p]          Semi-Parameter (Only required for Parabolic)    [km]
%
% Outputs:
%
%     - [EBH]        Eccentric Anomaly for given orbit type          [rad]
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 51-66, 72-79)
%     - Orbital Mechanics for Engineering Students, 3rd ed. (Curtis) (pg. 147-167)
% -------------------------------------------------------------------------------------------------

function [EBH] = kepler_equation(mu,M,e,dt,p)

   %Eccentric Anomaly
   if (e < 1)
      
      if (M > -pi && M < 0)
         E1 = M - e;
      elseif (M > pi)
         E1 = M - e;
      else 
         E1 = M + e;
      end
      
      E2 = E1 + ((M - E1 + e*sin(E1))/(1 - e*cos(E1))); %First Iteration
      
      diff = 1;
      iterations = 0;
      while (diff > 1e-6) %Converging Eccentric Anomaly
         
         E1 = E2;
         E2 = E1 + ((M - E1 + e*sin(E1))/(1 - e*cos(E1)));
         
         diff = abs(E2) - abs(E1);
         iterations = iterations + 1;
         if (iterations > 100)
            break
         end
         
      end
      EBH = E2;
   
   %Parabolic Eccentric Anomaly
   elseif (e == 1)
      
      np = 2*sqrt(mu/(p^3));
      poly = [(1/3) 0 1 (-np*dt)];
      rts = roots(poly);
      for j=1:1:numel(rts)
   
         output = imag(rts(j)) ~= 0;
         if (output == 0)
            B = rts(j);
         end
      end
      EBH = B;
   
   %Hyperbolic Eccentric Anomaly
   elseif (e > 1)
      
      if (e < 1.6)
         if (M > -pi && M < 0)
            H1 = M - e;
         else
            H1 = M + e;
         end
      else
         if (e < 3.6 && abs(M) > pi)
            H1 = M - (M/abs(M))*e;
         else
            H1 = M/(e-1);
         end
      end
      
      H2 = H1 + ((M - e*sinh(H1) + H1)/(e*cosh(H1) - 1)); %First Iteration
      
      iterations = 0;
      while (abs(H2 - H1) < 0.001)
         
         H1 = H2;
         H2 = H1 + ((M - e*sinh(H1) + H1)/(e*cosh(H1) - 1)); %Converging Eccentric Anomaly
                  
         iterations = iterations + 1;
         if (iterations > 100)
            break
         end
         
      end
      EBH = H2;
            
   end
   
end