
% ASSOCIATED LEGENDRE POLYNOMIALS
% --------------------------------------------------------------------------------------------------------------------- 
% Function computes Associated Legendre Polynomials P(l,m) from degree 0 to L and all orders for each degree.
%     https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
% 
% Author: Matthew Buckhout
% Updated: 08/06/2020 
%
% Inputs
%
%     - [rows]       Number polynomials to compute (all degree-order pairs up to order L)     
%     - [L]          Maximum order to compute polynomials for                                 
%     - [ARG]        Argument of Legendre Polynomials                                        
%     - [P00]        Initial Value, Associated Legendre Polynomial, Degree 0 Order 0          
%     - [P10]        Initial Value, Associated Legendre Polynomial, Degree 1 Order 0          
%     - [P11]        Initial Value, Associated Legendre Polynomial, Degree 1 Order 1          
%
% Outputs
%
%     - [P]          Array of computed polynomials, P(l,m)                                              
%     - [P1]         Array of computed polynomials, P(l,m+1)                                  
%
% ---------------------------------------------------------------------------------------------------------------------     

function [P,P1] = asc_legendre(rows,L,ARG,P00,P10,P11)

   step = 1;
   VALS1 = zeros(rows+4+L,1);
   VALS2 = zeros(rows+4+L,1);
   VALS3 = zeros(rows+4+L,1);
   VALS3(1,1) = 0;
   for l=1:1:L+1
       
      for m=0:1:l

         if (l == 1) && (m == 0)
            PLM = P10;
            PLM1M = P00;
         elseif (l == 1) && (m == 1)
            PLM = P11;
            PLM1M = 0;
         elseif (l == m)
            PLL = VALS1(step-(l+1));
            PLM = -((2*(l-1))+1)*sqrt(1 - (ARG^2))*PLL;
            PLM1M = 0; 
         else
            PLM = VALS2(step-(l));
            PLM1M = VALS1(step-(l));
         end
      
         PLP1M = (((2*l)+1)*ARG*PLM - (l+m)*PLM1M)/(l-m+1);

         VALS1(step,1) = PLM; %P(l,m)
         VALS2(step,1) = PLP1M; %P(l+1,m)

         step = step+1;

      end
      
   end
 
   P = VALS1(3:end-(l+1));
   P1 = VALS1(4:end-(l));

end 

