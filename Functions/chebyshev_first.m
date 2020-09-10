
% CHEBYSHEV POLYNOMIALS OF THE FIRST KIND (ORDERS 0 TO N)
% -------------------------------------------------------------------------------------------------
% Returns polynomials and first derivatives evaluated for argument (ARG) between -1 and 1 as an 
% array for order 0 to specified order N. Recurrence relations for Chebyshev polynomials: 
%     https://en.wikipedia.org/wiki/Chebyshev_polynomials#First_kind
%
% Author: Matthew Buckhout
% Updated: 08/06/2020  
%
% Inputs
% 
%     - [N]                Maximum order to compute up to                   
%     - [ARG]              Argument of Chebyshev Polynomials                
%
% Outputs
%
%     - [Tn]               Array, Chebyshev Polynomials order 0 to N        
%     - [dTn]              Array, First Derivative of Chebyshev             
%                             Polynomials order 0 to N    
% -------------------------------------------------------------------------------------------------

function [Tn,dTn] = chebyshev_first(N,ARG)
   
   T0 = 1; %Initial values
   T1 = ARG;
   VALS = zeros(N,1);
   VALS(1) = T0;
   VALS(2) = T1;
   
   dT0 = 0; %Initial values
   dT1 = 1; 
   dVALS = zeros(N,1); 
   dVALS(1) = dT0; 
   dVALS(2) = dT1; 
   
   for n=1:1:N-2
      
      if (n == 1)
         
         TNM1 = T0; % T(n-1)
         TN = T1; % T(n)
         dTNM1 = dT0; % dT(n-1)
         dTN = dT1; % dT(n)
         
      end
   
      TNP1 = 2*ARG*TN - TNM1; %  T(n+1) Recurrence relation
      dTNP1 = 2*ARG*dTN - dTNM1 + 2*TN; % dT(n+1) Recurrence relation, first derivative
   
      TNM1 = TN; %Substituting values for iteration
      TN = TNP1;
      dTNM1 = dTN; 
      dTN = dTNP1; 
      VALS(n+2) = TNP1; %Storing outputs
      dVALS(n+2) = dTNP1; 
   
   end
   
   Tn = VALS; %Chebyshev Polynomials
   dTn = dVALS; %First Derivatives
   
end





