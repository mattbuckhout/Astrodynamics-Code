
% LAMBERT'S PROBLEM - UNIVERSAL VARIABLE SOLUTION
% -------------------------------------------------------------------------------------------------
% General solution to Lambert's problem. Computes the velocity vectors of a spacecraft on a 
% trajectory between 2 position vectors separated by a specified amount of time. Solution uses f 
% and g functions with universal variables to handle all orbit types. 
%
%
% Author: Matthew Buckhout
% Updated: 08/08/2020 
%
% Inputs:
%
%     - [mu]         Central Body Gravitational Parameter            [km^3/s^2]
%     - [R0]         Position Vector (t0)                            [km]
%     - [R]          Position Vector (t)                             [km]
%     - [dt]         Time between R0 and R (t-t0)                    [sec]
%     - [tm]         Transfer Method                                  -
%                       tm = +1 (Short Way)
%                       tm = -1 (Long Way)
%
% Outputs:
%
%     - [V0]         Velocity on computed trajectory at R0 (t0)      [km/s]
%     - [V]          Velocity on computed trajectory at R (t)        [km/s]
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 448-454, 461-466)
% -------------------------------------------------------------------------------------------------

function [V0,V] = lambert_UV(mu,R0,R,dt,tm)
   
   r = norm(R); %[km] Final Position Vector Magnitude
   r0 = norm(R0); %[km] Initial Position Vector Magnitude
   
   cos_dtheta = dot(R0,R)/(r0*r);
   
   A = tm*sqrt(r*r0*(1 + cos_dtheta));    
   
   %Initial Values
   PSI1 = 0;
   c2  = 1/2;
   c3  = 1/6;
   
   %Bounds
   PSIup  = 4*(pi^2);
   PSIlow = -4*pi;
      
   %BISECTION METHOD
   dt2 = dt+1;
   iterations = 0;
   iterations2 = 0;
   while (abs(dt2 - dt) > 1e-6)
      
      y = r0 + r + ((A*(PSI1*c3 - 1))/sqrt(c2));
      
      if ((A > 0) && (y < 0))
         
         while (y < 0)
            
            %Adjusting PSIlow for Positive y
            PSIlow = PSIlow + 0.01*pi;
            PSI1 = (PSIup + PSIlow)/2;
            %Common Terms in Universal Variable Formulation
            if (PSI1 > 1e-8)
               c2 = (1-cos(sqrt(PSI1)))/PSI1;
               c3 = (sqrt(PSI1) - sin(sqrt(PSI1)))/sqrt(PSI1^3);
            elseif (PSI1 < -1e-8)
               c2 = (1 - cosh(sqrt(-PSI1)))/PSI1;
               c3 = (sinh(sqrt(-PSI1)) - sqrt(-PSI1))/sqrt((-PSI1)^3);
            else
               c2 = 1/2;
               c3 = 1/6;
            end
            y = r0 + r + ((A*(PSI1*c3 - 1))/sqrt(c2)); 
            iterations2 = iterations2 + 1;
            
            if (iterations2 > 1000)
               break
            end
         
         end
         
      end            
      
      CHI = sqrt(y/c2); 
      dt2 = ((CHI^3)*c3 + A*sqrt(y))/sqrt(mu); 
      
      if (abs(dt2) <= dt)
         PSIlow = PSI1;
      else
         PSIup = PSI1;
      end
      
      PSI2 = (PSIup + PSIlow)/2;
         
      %Common Terms in Universal Variable Formulation
      if (PSI2 > 1e-8)
         c2 = (1-cos(sqrt(PSI2)))/PSI2;
         c3 = (sqrt(PSI2) - sin(sqrt(PSI2)))/sqrt(PSI2^3);
      elseif (PSI2 < -1e-8)
         c2 = (1 - cosh(sqrt(-PSI2)))/PSI2;
         c3 = (sinh(sqrt(-PSI2)) - sqrt(-PSI2))/sqrt((-PSI2)^3);
      else
         c2 = 1/2;
         c3 = 1/6;
      end
      
      PSI1 = PSI2;
      iterations = iterations + 1;
                    
      if (iterations > 1000)
         fprintf('lambert_UV - Failed to converge\n');
         break
      end
      
   end
   
   %f and g Functions
   f = 1-(y/r0);
   g = A*sqrt(y/mu);
   g_dot = 1 - (y/r);
   
   
   %Velocity Vectors
   V0 = (R - f*R0)/g; %[km/s] Initial Velocity Vector
   V  = (g_dot*R - R0)/g; %[km/s] Final Velocity Vector
   
end