
% EGM2008 NON-SPHERICAL GRAVITATIONAL MODEL 
% -------------------------------------------------------------------------------------------------
% Function generates perturbing accelerations on a satellite in a Geocentric orbit due to Earth's 
% Non-spherical Gravitational Potential. The non-spherical potential function is developed using 
% spherical harmonic coefficients (S,C) contained within the Earth Gravitational Model (EGM2008) 
% and Associated Legendre Polynomials. Equations for components (x,y,z) of the perturbing 
% acceleration are then developed using partial derivatives of this potential function. The 
% resulting acceleration accounts for perturbing effects due to zonal, tesseral, and sectorial 
% harmonics for an L x L field, where L is specified by the user and corresponds to the maximum 
% spherical harmonic degree/order of coefficients from the EGM2008 model.
%
% Acceleration vector is developed and returned in the body-fixed ITRF frame and as such this 
% function requires the satellite position vector to be inputted in the ITRF frame.
%
%
% Author: Matthew Buckhout
% Updated: 08/11/2020 
%
% Inputs
%
%     - [R_ITRF]           Satellite Position Vector (ITRF)                [km]
%     - [L]                Maximum Degree/Size of Field for Harmonics       -
%     - [EGM2008coef]      Un-Normalaized Spherical Harmonic Coefficients   -
%
% Outputs
%
%     - [ap_aspherical]    Perturbing Acceleration Vector (ITRF)           [km/s^2]
%
% Functions
%
%     - asc_legendre_c         
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 511-524)
%     - Full IERS Conventions 2010
%           https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html
%     - EGM2008 Spherical Harmonic Coefficients (2190 x 2190)
%           https://earth-info.nga.mil/GandG/update/index.php?action=home#tab_wgs84-data
%     - Recurrence Relations for Associated Legendre Polynomials
%           https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
% -------------------------------------------------------------------------------------------------

function [ap_aspherical] = EGM2008(R_ITRF,L,EGM2008coef)

   %Spherical Harmonic Data
   CS = EGM2008coef;
   [rows,~] = size(CS);

   %Constants
   RE = 6378.1363; %[km] Earth Mean Equatorial Radius 
   mu = 3.986004415e5; %[km^3/s^2]

   r_ITRF = norm(R_ITRF); %[km] 

   %Geocentric Latitude Argument of Associated Legendre Polynomials
   ARG = R_ITRF(3)/r_ITRF; %ARG = sin(Geocentric Latitude)
   
   %Longitude of Satellite
   if ((R_ITRF(1) > 0) && (R_ITRF(2) > 0))
      Long = atan(abs(R_ITRF(2))/abs(R_ITRF(1)));
   elseif ((R_ITRF(1) < 0) && (R_ITRF(2) > 0))
      Long = atan(abs(R_ITRF(1))/abs(R_ITRF(2))) + (pi/2);
   elseif ((R_ITRF(1) < 0) && (R_ITRF(2) < 0))
      Long = atan(abs(R_ITRF(2))/abs(R_ITRF(1))) + pi;
   elseif ((R_ITRF(1) > 0) && (R_ITRF(2) < 0))
      Long = atan(abs(R_ITRF(1))/abs(R_ITRF(2))) + (3/2)*pi;
   end
   
   %Initial Values in Legendre Polynomials
   P00 = 1;
   P10 = ARG;
   P11 = -(R_ITRF(2)/r_ITRF)/(sin(Long));
   
   %Algorithm for Associated Legendre Polynomials through Order L (C/mex)
   B = asc_legendre_c(rows,L,ARG,P00,P10,P11); 
   P = B(:,1);
   P1 = B(:,2);

   %Partial Derivatives of Potential Function
   dUdr = (-mu/(r_ITRF^2))*sum(((RE/r_ITRF).^CS(:,1)).*(CS(:,1)+1).*P.*(CS(:,3).*cos(CS(:,2).*Long) + CS(:,4).*sin(CS(:,2).*Long)));
   dUdlat = (mu/r_ITRF)*sum(((RE/r_ITRF).^CS(:,1)).*(P1 - CS(:,2).*(((R_ITRF(3)/R_ITRF(1))*cos(Long)).*P)).*(CS(:,3).*cos(CS(:,2).*Long) + CS(:,4).*sin(CS(:,2).*Long)));
   dUdlong = (mu/r_ITRF)*sum(((RE/r_ITRF).^CS(:,1)).*CS(:,2).*P.*(CS(:,4).*cos(CS(:,2).*Long) - CS(:,3).*sin(CS(:,2).*Long)));

   %Perturbing Acceleration
   ap_aspherical(1) = (((1/r_ITRF)*dUdr) + (R_ITRF(3)/((r_ITRF^2)*sqrt((R_ITRF(1)^2) + (R_ITRF(2)^2))))*dUdlat)*R_ITRF(1) - ((1/((R_ITRF(1)^2) + (R_ITRF(2)^2)))*dUdlong)*R_ITRF(2);
   ap_aspherical(2) = (((1/r_ITRF)*dUdr) + (R_ITRF(3)/((r_ITRF^2)*sqrt((R_ITRF(1)^2) + (R_ITRF(2)^2))))*dUdlat)*R_ITRF(2) + ((1/((R_ITRF(1)^2) + (R_ITRF(2)^2)))*dUdlong)*R_ITRF(1);
   ap_aspherical(3) = ((1/r_ITRF)*dUdr)*R_ITRF(3) - (sqrt((R_ITRF(1)^2) + R_ITRF(2)^2)/(r_ITRF^2))*dUdlat;
   
   ap_aspherical = transpose(ap_aspherical); %Column
   
end




















   
% OLD VERSION, NOT VECTORIZED
%-------------------------------------------------------------------------------------------------   
   
%   %Spherical Harmonic Data
%   C = DATA7; %UnNormalized C coefficients 
%   S = DATA8; %UnNormalized S Coefficients
%   [rows,~] = size(C);
%
%   %Constants
%   RE = 6378.1363; %[km] Earth Mean Equatorial Radius 
%   mu = 3.986004415e5; %[km^3/s^2]
%
%   r_ITRF = norm(R_ITRF); %[km] 
%
%   %Geocentric Latitude Argument of Associated Legendre Polynomials
%   ARG = R_ITRF(3)/r_ITRF; %ARG = sin(Geocentric Latitude)
%   
%   %Longitude of Satellite
%   if ((R_ITRF(1) > 0) && (R_ITRF(2) > 0))
%      Long = atan(abs(R_ITRF(2))/abs(R_ITRF(1)));
%   elseif ((R_ITRF(1) < 0) && (R_ITRF(2) > 0))
%      Long = atan(abs(R_ITRF(1))/abs(R_ITRF(2))) + (pi/2);
%   elseif ((R_ITRF(1) < 0) && (R_ITRF(2) < 0))
%      Long = atan(abs(R_ITRF(2))/abs(R_ITRF(1))) + pi;
%   elseif ((R_ITRF(1) > 0) && (R_ITRF(2) < 0))
%      Long = atan(abs(R_ITRF(1))/abs(R_ITRF(2))) + (3/2)*pi;
%   end
% 
%   %Computing UnNormalized Associated Legendere Polynomials 
%   l = 2;
%   P = transpose(legendre(l,ARG)); %Unnormalized Associated Legendre Polynomials, Degree 2
%   while (l <= L)
%
%      l = l+1;
%      Q = transpose(legendre(l,ARG)); %Unnormalized Associated Legendre Polynomials, Degree l
%      P = [P Q]; %Concatenating results
%      
%   end
%
%   %Summation Terms for Partial Derivatives
%   tabrow = 1; 
%   sumr = 0;
%   sumlat = 0;
%   sumlong = 0;
%   for i=1:1:L-1
%      
%      l = i+1; %l starts at 2
%      
%      for j=1:1:l+1
%         
%         m = j-1; %m starts at 0
%      
%         %Nested Summation WRT Radius
%         sumr = sumr + ((RE/r_ITRF)^l)*(l+1)*P(tabrow)*(C(tabrow)*cos(m*Long) + S(tabrow)*sin(m*Long));
%         %Nested Summation WRT Geocentric Latitude
%         sumlat = sumlat + ((RE/r_ITRF)^l)*(  P(tabrow+1) - m*(((R_ITRF(3)/R_ITRF(1))*cos(Long))*P(tabrow))  )*(  C(tabrow)*cos(m*Long) + S(tabrow)*sin(m*Long)  );
%         %Nested Summation WRT Longitude
%         sumlong = sumlong + ((RE/r_ITRF)^l)*m*P(tabrow)*(  S(tabrow)*cos(m*Long) - C(tabrow)*sin(m*Long)  );
%         
%         tabrow = tabrow+1;
%            
%      end
%
%   end
%
%   %Partial Derivatives of Potential Function
%   dUdr = (-mu/(r_ITRF^2))*sumr;
%   dUdlat = (mu/r_ITRF)*sumlat;
%   dUdlong = (mu/r_ITRF)*sumlong;
%   
%   %Perturbing Acceleration
%   ap_aspherical(1) = (((1/r_ITRF)*dUdr) + (R_ITRF(3)/((r_ITRF^2)*sqrt((R_ITRF(1)^2) + (R_ITRF(2)^2))))*dUdlat)*R_ITRF(1) - ((1/((R_ITRF(1)^2) + (R_ITRF(2)^2)))*dUdlong)*R_ITRF(2);
%   ap_aspherical(2) = (((1/r_ITRF)*dUdr) + (R_ITRF(3)/((r_ITRF^2)*sqrt((R_ITRF(1)^2) + (R_ITRF(2)^2))))*dUdlat)*R_ITRF(2) + ((1/((R_ITRF(1)^2) + (R_ITRF(2)^2)))*dUdlong)*R_ITRF(1);
%   ap_aspherical(3) = ((1/r_ITRF)*dUdr)*R_ITRF(3) - (sqrt((R_ITRF(1)^2) + R_ITRF(2)^2)/(r_ITRF^2))*dUdlat;

   
