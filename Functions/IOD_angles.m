
% ANGLES ONLY INITIAL ORBIT DETERMINATION
% -------------------------------------------------------------------------------------------------
% Function uses Gauss's Technique for angles only orbit determination; takes observations of an 
% object in Topocentric Right Ascension and Declination at three closely spaced times and 
% determines the geocentric position vectors of the object at each time.
%
% The three position vectors are then used with either the Gibbs or Herrick-Gibbs algorithm to 
% compute the velocity vector at the middle time. The object state is then represented by the 
% Position and Velocity vectors at the middle time.
%
% Best suited for observations with angular separation < 60 degrees 
% High accuracy for observations with angular separation < 10 degrees
%
% Observations recorded in ITRF frame (via Topocentric RADEC)
% Calculations performed in GCRF frame
%
%
% Author: Matthew Buckhout
% Updated: 09/10/2020 
%
% Inputs:
%
%     - [Lat]        Geodetic Latitude                               [deg]
%     - [Long]       Longitude (+East of Greenwich)                  [deg]
%     - [h_ellip]    Height above Reference Ellipsoid                [km]
%     - [dUT1]       UT1 - UTC (middle time)                         [sec]
%     - [dAT]        TAI - UTC (middle time)                         [sec]
%     - [xp]         Position of CEP from IRP (middle time)          ["]
%     - [yp]         Position of CEP from IRP (middle time)          ["]
%     - [ALPHA]      Topocentric Right Ascension                     [rad]
%                       [ alpha1; alpha2; alpha3 ]
%     - [DELTA]      Topocentric Declination                         [rad]
%                       [ delta1; delta2; delta3 ]
%     - [T1]         Gregorian Dates of Observation                   -
%                       [ yyyy1 mm1 dd1; 
%                         yyyy2 mm2 dd2;
%                         yyyy3 mm3 dd3 ]
%     - [T2]         Times UTC of Observation                         -
%                       [ hh1 mm1 ss1;
%                         hh2 mm2 ss2;
%                         hh3 mm3 ss3 ]
%     - [XTAB]       Selected X coefficients for series terms in      -
%                       IAU2006 Precession-Nutation Model
%     - [YTAB]       Selected Y coefficients for series terms in      -
%                       IAU2006 Precession-Nutation Model
%     - [sTAB]       Selected s coefficients for series terms in      -
%                       IAU2006 Precession-Nutation Model
%
% Outputs:
%
%     - [R2]         Position Vector at middle time                  [km]
%     - [V2]         Velocity Vector at middle time                  [km/s]
%
% Functions:
%
%     - site
%     - IAU2006 
%     - orbital_elements
%     - julian_date 
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 413-447)
% -------------------------------------------------------------------------------------------------

function [R2,V2] = IOD_angles(Lat,Long,h_ellip,T1,T2,ALPHA,DELTA,dUT1,dAT,xp,yp,XTAB,YTAB,sTAB)

   %Constants
   RE = 6378.1363; %[km] Earth Mean Equatorial Radius 
   mu = 3.986004415e5; %[km^3/s^2] Earth Gravitational Parameter

   %Preallocation
   UT1 = zeros(3,3);
   JD_UT1 = zeros(1,3);
   L = zeros(3,3);
   RSITE = zeros(3,3);
   R = zeros(3,3);

   %Unit Conversions
   Lat = deg2rad(Lat); %[rad]
   Long = deg2rad(Long); %[rad]
   for i=1:1:3
      t(i) = (T2(i,1)*3600) + (T2(i,2)*60) + T2(i,3); %[sec]
   end

   %Site Vector - ITRF          
   [Rsite_ITRF] = site(Lat,Long,h_ellip); %[km]  
   
   for i=1:1:3
      
      %Line of Sight Vectors
      L(:,i) = [cos(DELTA(i))*cos(ALPHA(i)); ...
                cos(DELTA(i))*sin(ALPHA(i)); ...
                sin(DELTA(i))];

      %Site Vectors - GCRF          
      [~,MJDi] = julian_date(T1(i,:),T2(i,:),0); %(UT1) Modified Julian Date
      [RSITE(:,i),~] = IAU2006(Rsite_ITRF,[0; 0; 0],MJDi,dUT1,dAT,xp,yp,2,'long',XTAB,YTAB,sTAB); %[km]
      
   end

   
   % GAUSS ROUTINE
   fprintf('\n');
   fprintf('Gauss');
   
   %Defining Coefficients
   tau1 = t(1) - t(2); %[sec]
   tau3 = t(3) - t(2); %[sec]
   a1  = tau3/(tau3 - tau1); %[sec]
   a3  = -(tau1/(tau3 - tau1)); %[sec]
   a1u = (tau3*(((tau3 - tau1)^2) - (tau3^2)))/(6*(tau3 - tau1)); %[sec]
   a3u = (-tau1*(((tau3 - tau1)^2) - (tau1^2)))/(6*(tau3 - tau1)); %[sec]
   M = inv(L)*RSITE; %[km]
   d1 = M(2,1)*a1 - M(2,2) + M(2,3)*a3; %[km]
   d2 = M(2,1)*a1u + M(2,3)*a3u; %[sec]
   C = dot(L(:,2),RSITE(:,2)); %[km]

   %Coefficients of 8th order polynomial of r2
   coef = [1; ...
           0; ...
           (-((d1^2) + (2*C*d1) + ((norm(RSITE(:,2)))^2))); ...
           0; ...
           0; ...
           (-2*mu*((C*d2) + (d1*d2))); ...
           0; ...
           0; ...
           (-(mu^2)*(d2^2))];

   %Finding Real, Positive root r2        
   rts = roots(coef);
   for j=1:1:numel(rts)
      
      if (imag(rts(j)) == 0) 
         if (real(rts(j)) >= 0)
            
            r2 = real(rts(j)); %[km]
            
         end
      end
      
   end

   %Calculating Coefficients
   u = mu/(r2^3);
   c1 = a1 + (a1u*u);
   c2 = -1;
   c3 = a3 + (a3u*u);

   %Computing Range Magnitudes
   c_rho = M*[-c1; -c2; -c3];
   rho = [c_rho(1)/c1; c_rho(2)/c2; c_rho(3)/c3]; %[km] Initial Range Estimates

   %Iterating to converge Range
   err1 = err2 = err3 = 1;
   eps = 1e-15;
   rho1 = rho;
   iterations = 0;
   while ((err1 > eps) || (err2 > eps) || (err3 > eps))
      
      for i=1:1:3
         
         R(:,i) = rho1(i)*L(:,i) + RSITE(:,i); %[km] Position Vectors
         r(i) = norm(R(:,i)); %[km] Position Vector Magnitude
         
      end
      
      %Separation Angles
      sep12 = rad2deg(acos(dot(R(:,1),R(:,2))/(norm(R(:,1))*norm(R(:,2))))); %[deg]
      sep23 = rad2deg(acos(dot(R(:,2),R(:,3))/(norm(R(:,2))*norm(R(:,3))))); %[deg]

      if ((sep12 <= 5) && (sep23 <= 5))
         
         % HERRICK-GIBBS ROUTINE 
         if (iterations == 0)
            fprintf(' - Herrick-Gibbs\n');
            fprintf('\n');
         end
         
         %Coplanar Position Vector Test
         Z23 = cross(R(:,2),R(:,3));
         coplanar = abs(90 - rad2deg(acos(dot(Z23,R(:,1))/(norm(Z23)*norm(R(:,1)))))); %[rad]
         if (coplanar < 1e-8) 
         
            dt31 = t(3) - t(1); %[sec]
            dt32 = t(3) - t(2); %[sec]
            dt21 = t(2) - t(1); %[sec]
         
            %Velocity Vector (GCRF)
            V2 = (-dt32*((1/(dt21*dt31)) + (mu/(12*((r(1))^3))))*R(:,1)) ...
                 + ((dt32 - dt21)*((1/(dt21*dt32)) + (mu/(12*((r(2))^3))))*R(:,2)) ...
                 + (dt21*((1/(dt32*dt31)) + (mu/(12*((r(3))^3))))*R(:,3)); %[km/s]
         
         end
         
      elseif ((sep12 > 5) || (sep23 > 5))
         
         
         % GIBBS ROUTINE
         if (iterations == 0)
            fprintf(' - Gibbs\n');
            fprintf('\n');
         end
                  
         %Storing Cross Products
         Z12 = cross(R(:,1),R(:,2));
         Z23 = cross(R(:,2),R(:,3));
         Z31 = cross(R(:,3),R(:,1));
         
         %Coplanar Position Vector Test
         coplanar = abs(90 - rad2deg(acos(dot(Z23,R(:,1))/(norm(Z23)*norm(R(:,1)))))); %[rad]
         if (coplanar < 1e-8) 
            
            %Auxilliary Vectors
            N  = r(1)*Z23 + r(2)*Z31 + r(3)*Z12;
            D  = Z12 + Z23 + Z31;
            S  = (r(2) - r(3))*R(:,1) + (r(3) - r(1))*R(:,2) + (r(1) - r(2))*R(:,3);
            B  = cross(D,R(:,2));
            Lg = sqrt(mu/(norm(N)*norm(D)));
            
            %Coplanar Auxiliary Vector Text
            coplanar_aux = rad2deg(acos(dot(N,D)/(norm(N)*norm(D)))); %[deg]
            if (coplanar_aux < 3)     
            
               %Velocity Vector (GCRF)
               V2 = (Lg/norm(R(:,2)))*B + Lg*S; %[km/s] 
               
            else
               
               fprintf('Auxilliary Vectors not sufficiently coplanar for accurate analysis\n');
                
            end
            
         else
            
            fprintf('Observations not sufficiently coplanar for accurate analysis\n');
         
         end
       
      end
      
      %Orbital Elements
      [p,~,e,~,~,~,theta2,~,~,~] = orbital_elements(mu,R(:,2),V2);
      
      %Computing f and g functions for 3 position vectors
      for i=1:1:3
         
         theta(i) = acos((1/e)*((p/r(i)) - 1)); %True Anomaly 
         dtheta(i) = theta(i) - theta2; %Deviation from central true anomaly
         
         f(i) = 1 - ((r(i)/p)*(1 - cos(dtheta(i))));
         
         g(i) = (r(i)*r(2)*sin(dtheta(i)))/sqrt(mu*p);
         
      end
      
      %Re-Calculating Coefficients
      c1 = g(3)/(f(1)*g(3) - f(3)*g(1));
      c2 = -1;
      c3 = (-g(1))/(f(1)*g(3) - f(3)*g(1));
      
      %Computing NEW Range Magnitudes
      c_rho = M*[-c1; -c2; -c3];
      rho2 = [c_rho(1)/c1; c_rho(2)/c2; c_rho(3)/c3]; %[km] Range Estimates
      
      err1 = abs(rho2(1) - rho1(1));
      err2 = abs(rho2(2) - rho1(2));
      err3 = abs(rho2(3) - rho1(3));
      
      rho1 = rho2;
      iterations = iterations + 1;
      
      if (iterations > 100)
         break
      end
      
   end

   %Position Vector (GCRF)
   R2 = rho(2)*L(:,2) + RSITE(:,2); %[km] 

end