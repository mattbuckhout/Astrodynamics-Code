
% TITLE SHORT (100 COLUMN LIMIT)
% -------------------------------------------------------------------------------------------------
% Function computes the combined acceleration on a satellite in Geocentric orbit due to third body 
% perturbations from the moon, sun, and planets. The function determines the accelerations caused 
% by each specified body individually at the input time and satellite position vector, expressed 
% in the GCRF frame, then combines them for a total resultant acceleration vector due to all 
% perturbing third bodies. This acceleration can then be added to the central body acceleration 
% term in the satellite equations of motion.
%
% The array "bodies" specifies which planetary bodies to include in the total perturbing 
% acceleration, must include the value "3" along with the identifiers for other bodies. Identifier 
% 3 refers to the position of the Earth-Moon Barycenter from the Solar System Barycenter, and is 
% needed by this function to transform the ephemerides from the Barycentric ICRF frame to the 
% Geocentric GCRF frame.
%
% Third body perturbations are developed using the JPL DE430 model for planetary ephemerides. The 
% DE430 model yields highly accurate position vectors of solar system bodies resulting from 
% numerical simulations of the system. 
%
%
% Author: Matthew Buckhout
% Updated: 08/18/2020 
%
% Inputs:
%
%     - [MJD_UT1]          Modified Julian Date, Expressed in TT or TDT          [days]
%     - [R_GCRF]           Satellite Position Vector (GCRF)                      [km]
%     - [DE430coef]        JPL DE430 Ephemerides Model Coefficients               -
%     - [EOPdata]          Earth Orientation Parameter Table                      -
%     - [bodies]           Vector of identifiers (1-11) specifying which          -
%                          bodies to generate ephemerides/accelerations 
%                          for (MUST INCLUDE VALUE OF 3)
%
%
% Outputs:
%
%     - [ap_thirdbodies]   Total Acceleration due to third bodies                [km/s^2]
%
% Functions:
%
%     - terrestrial_MJD
%     - planet_ephemerides_c
%     - series_B_c
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 542-544)
% -------------------------------------------------------------------------------------------------

function [ap_thirdbodies] = third_body(MJD_UT1,R_GCRF,DE430coef,EOPdata,bodies)

   ap_thirdbodies = [0; 0; 0]; %Initializing acceleration vector

   %Planetary Ephemerides
   [MJD_TT_step] = terrestrial_MJD(MJD_UT1,EOPdata);
   [ephemerides] = planet_ephemerides_c(MJD_TT_step,transpose(DE430coef),bodies); 

   %Position Vectors
   [~,ind] = find(bodies == 3);
   Rembc_ICRF = transpose(ephemerides(ind,:)); %Earth-Moon Barycenter (ICRF, SS Barycenter)

   %Gravitational Parameters
   mus = [2.2032e4 3.257e5 0 4.305e4 1.268e8 3.794e7 5.794e6 6.809e6 9e2 4902.799 1.32712428e11];

   %Accelerations
   for k=1:1:numel(bodies) %Loop for perturbing 3rd bodies

      if (bodies(k) ~=3) %Ephemerides other than Earth-Moon Barycenter

         if (bodies(k) == 10) %Moon
         
            RE3 = transpose(ephemerides(k,:)); %[km] Position Vector, Geocenter to 3rd body    

         else
         
            RBC3_ICRF = transpose(ephemerides(k,:)); %[km] Position Vector, SS Barycenter to 3rd body
            RE3 = RBC3_ICRF - Rembc_ICRF; %[km] Position Vector, Geocenter to 3rd body

         end

         REsat = R_GCRF; %[km] Position Vector, Geocenter to Satellite
         Rsat3 = RE3 - REsat; %[km] Position Vector, Satellite to 3rd Body
         re3 = norm(RE3); %[km] Magnitudes 
         resat = norm(REsat); %[km]
         rsat3 = norm(Rsat3); %[km]

         ARG = ((rsat3^2) - (resat^2) - (re3^2))/(-2*resat*re3); %[rad] Cosine of Angle between 3rd body and Satellite position vectors
         
         TOL = 1e-8; %Tolerance for convergence in B
         B = series_B_c(ARG,TOL,resat,re3); %Series term B
         BETA = 3*B + 3*(B^2) + (B^3);

         ap_thirdbodies = ap_thirdbodies - (mus(bodies(k))/(re3^3))*(REsat - BETA*Rsat3); %[km/s^2] Acceleration

      end

   end
   
end

