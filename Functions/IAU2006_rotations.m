
% IAU 2006 REDUCTION - ROTATION MATRICES ONLY
% ---------------------------------------------------------------------------------------------------------------------
% Function returns rotation matrices for transformation of vectors between celestial GCRF and body-fixed ITRF frames, 
% accounting for the rotation of the Earth, Polar Motion, and Precession-Nutation. IAU 2006 reduction realized by 
% implementing original IAU 2000 theory in Vallado with updated coefficients from IERS Conventions 2010. 
%
% length = 'short' for X,Y,s files edited down to rows for largest 10 terms in each series for faster computation.
% length = 'long' for full X,Y,s files, slower computation but higher precision.
%
%
% Author: Matthew Buckhout
% Updated: 08/25/2020 
%
% Acronyms:
%
% ITRF - International Terrestrial Reference Frame (Body Fixed)
% GCRF - Geocentric Celestial Reference Frame (Inertial, Earth Centered)
% PEF  - Pseudo Earth Fixed Frame (Body Fixed)
% IRE  - Intermediate Frame of Reference Epoch (Inertial, Earth Centered)
% ICRF - International Celestial Reference Frame (Inertial, Solarsystem Barycenter)
% CIP - Celestial Intermediate Pole (located by angles xp,yp)
% CEO - Celestial Ephemeris Origin (non-rotating origin on equator of CIP)
% TEO - Terrestrial Ephemeris Origin (reference meridian)
%
% Inputs:
%                    
%     - [MJD_UT1]       Modified Julian Date (UT1)                      [days]
%     - [dUT1]          UT1 - UTC                                       [sec]
%     - [dAT]           TAI - UTC                                       [sec]
%     - [xp]            Position of CEP from IRP                        ["]
%     - [yp]            Position of CEP from IRP                        ["]
%     - [way]           ITRF-GCRF (2) or GCRF-ITRF (1)                   -
%     - [length]        Length of X,Y,s data sets                        -
%     - [XTAB]          Selected X coefficients for series terms in      -
%                          IAU2006 Precession-Nutation Model
%     - [YTAB]          Selected Y coefficients for series terms in      -
%                          IAU2006 Precession-Nutation Model
%     - [sTAB]          Selected s coefficients for series terms in      -
%                          IAU2006 Precession-Nutation Model
%
% Outputs:
%
%     - [ITRF2PEF]      Rotation matrix ITRF to PEF                      
%     - [PEF2IRE]       Rotation matrix PEF to IRE                       
%     - [IRE2GCRF]      Rotation matrix IRE to GCRF                      
%     - [GCRF2IRE]      Rotation matrix GCRFA to IRE                     
%     - [IRE2PEF]       Rotation matrix IRE to PEF                       
%     - [PEF2ITRF]      Rotation matrix PEF to ITRF       
%
% Functions:
%
%     - gregorian_date   
%     - convert_time    
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 205-219)
%     - Full IERS Conventions 2010
%           https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html
%     - Earth Orientation Parameter Data from IERS 
%           https://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html
%     - Coefficients of Series Terms for IAU 2000 X,Y,s Series terms from IERS Conventions (2010)
%           http://webtai.bipm.org/iers/convupdt/convupdt_c5.html#Listupdt
% ---------------------------------------------------------------------------------------------------------------------               
 
function [ITRF2PEF,PEF2IRE,IRE2GCRF,GCRF2IRE,IRE2PEF,PEF2ITRF] = IAU2006_rotations(MJD_UT1,dUT1,dAT,xp,yp,length,XTAB,YTAB,sTAB)
   
   %Julian Date (UT1)
   JD_UT1 = MJD_UT1 + 2400000.5;

   %Converting to Date/Time (UTC)
   [yr,mo,d,h,m,s] = gregorian_date(JD_UT1 - (dUT1/86400) - 2400000.5);

   %Obtaining Time Values (input UTC)
   [UT1,~,~,~,~,T_TT,T_TDB] = convert_time([yr mo d],[h m s],dUT1,dAT);

   %Earth Mean Angular Rotation (Vallado 184)
%   omegaE = 1.002737909350795 + (5.9006e-11)*(T_UT1) - (5.9e-15)*(T_UT1^2); %[rev/sidereal day]
%   omegaE = ((2*pi)/(86164.09052))*omegaE; %[rad/s]
%   OMEGAE = [0 0 omegaE]; %[rad/s]
   OMEGAE = [0 0 7.292115e-5]; %[rad/s]

   %Delaunay Arguments
   M_lun  = (pi/180)*(134.96340251 + (1717915923.2178/3600)*(T_TDB) ... 
            + (31.8792/3600)*(T_TDB^2) + (0.051635/3600)*(T_TDB^3) ... 
            - (0.00024470/3600)*(T_TDB^4)); %[rad] Luna Mean Anomaly   

   M_sun  = (pi/180)*(357.52910918 + (129596581.0481/3600)*(T_TDB) ... 
            - (0.5532/3600)*(T_TDB^2) + (0.000136/3600)*(T_TDB^3) ... 
            - (0.00001149/3600)*(T_TDB^4)); %[rad] Sun Mean Anomaly

   F_lun  = (pi/180)*(93.27209062 + (1739527262.8478/3600)*(T_TDB) ... 
            - (12.7512/3600)*(T_TDB^2) - (0.001037/3600)*(T_TDB^3) ... 
            + (0.00000417/3600)*(T_TDB^4)); %[rad] Luna Mean Long. - OM_lun

   D_sun  = (pi/180)*(297.85019547 + (1602961601.2090/3600)*(T_TDB) ... 
            - (6.3706/3600)*(T_TDB^2) + (0.006593/3600)*(T_TDB^3) ... 
            - (0.00003169/3600)*(T_TDB^4)); %[rad] Sun Mean Elongation

   OM_lun = (pi/180)*(125.04455501 - (6962890.5431/3600)*(T_TDB) ... 
            + (7.4722/3600)*(T_TDB^2) + (0.007702/3600)*(T_TDB^3) ...
            - (0.00005939/3600)*(T_TDB^4)); %[rad] Luna Right Ascension of Ascending Node

   %Heliocentric Longitudes of the Planets
   lambda_M_mercury = (4.402608842 + 2608.7903141574*T_TDB); %[rad]
   lambda_M_venus   = (3.176146697 + 1021.3285546211*T_TDB); %[rad]
   lambda_M_earth   = (1.753470314 + 628.3075849991*T_TDB); %[rad]
   lambda_M_mars    = (6.203480913 + 334.0612426700*T_TDB); %[rad]
   lambda_M_jupiter = (0.599546497 + 52.9690962641*T_TDB); %[rad]
   lambda_M_saturn  = (0.874016757 + 21.3299104960*T_TDB); %[rad]
   lambda_M_uranus  = (5.481293872 + 7.4781598567*T_TDB); %[rad]
   lambda_M_neptune = (5.311886287 + 3.8133035638*T_TDB); %[rad]

   %General Precession of Longitude
   pa = (0.02438175*T_TDB) + (0.00000538691*(T_TDB^2)); %[rad]

   %Vector of Fundamental Arguments and Mean Longitudes
   CONSTANTS = [M_lun; M_sun; F_lun; D_sun; OM_lun; ... 
                lambda_M_mercury; lambda_M_venus; lambda_M_earth; ... 
                lambda_M_mars; lambda_M_jupiter; lambda_M_saturn; ... 
                lambda_M_uranus; lambda_M_neptune; pa]; %[rad]

   %Removing Extra Periods          
   CONSTANTS = CONSTANTS - fix(CONSTANTS./(2*pi)).*(2*pi); %[rad]

   %Polar Motion 
   ac = 0.26; %["] Chandler wobble
   aa = 0.12; %["] Annual wobble
   sprime = (0.0015*(((ac^2)/1.2) ... 
            + (aa^2))*T_TT)*(pi/180)*(1/3600); %[rad] Instantaneous Prime Meridian

   %Earth Rotation 
   ERA = 2*pi*(0.7790572732640 + 1.00273781191135448*(JD_UT1 - 2451545.0)); %[rad] Earth Rotation Angle
   ERA = ERA - fix(ERA/(2*pi))*(2*pi); %[rad] Removing extra periods
 
   %Indices for table lookup XTAB, YTAB, sTAB
   switch length
      
      case 'long'
         lowx = [1 1307 1560 1596 1600];
         highx = [1306 1559 1595 1599 1600];
         lowy = [1 963 1240 1270 1275];
         highy = [962 1239 1269 1274 1275];
         lows = [1 34 37 62 66];
         highs = [33 36 61 65 66];
 
      case 'short'
         lowx = [1 11 21 31 35];
         highx = [10 20 30 34 35];
         lowy = [1 11 21 31 36];
         highy = [10 20 30 35 36]; 
         lows = [1 11 14 24 28];
         highs = [10 13 23 27 28];
      
   end 
   
   %Solving for value of X
   times = [1 T_TDB T_TDB^2 T_TDB^3 T_TDB^4];
   diffs = (highx-lowx)+[1 1 1 1 1];
  
   Asc = XTAB(:,2)/(1e6); %Sine Coefficient (Data in micro", converted to ")
   Acc = XTAB(:,3)/(1e6); %Cosine Coefficient (Data in micro", converted to ")
   Tc = vertcat(ones(diffs(1),1),T_TDB*ones(diffs(2),1),(T_TDB^2)*ones(diffs(3),1),(T_TDB^3)*ones(diffs(4),1),(T_TDB^4)*ones(diffs(5),1)); %Column of Times
   apc = XTAB(:,4:17)*CONSTANTS; %Argument
   
   XX = sum(Asc.*Tc.*sin(apc) + Acc.*Tc.*cos(apc)); %Summation
 
   X = -0.016617 ...
       + (2004.191898*T_TDB) ... 
       - (0.4297829*(T_TDB^2)) ... 
       - (0.19861834*(T_TDB^3)) ... 
       - (0.000007578*(T_TDB^4)) ... 
       + (0.0000059285*(T_TDB^5)) + XX; %["] Angle of CIP in ICRF

   %Solving for value of Y
   diffs = (highy-lowy)+[1 1 1 1 1];
  
   Asc = YTAB(:,2)/(1e6); %Sine Coefficient (Data in micro", converted to ")
   Acc = YTAB(:,3)/(1e6); %Cosine Coefficient (Data in micro", converted to ")
   Tc = vertcat(ones(diffs(1),1),T_TDB*ones(diffs(2),1),(T_TDB^2)*ones(diffs(3),1),(T_TDB^3)*ones(diffs(4),1),(T_TDB^4)*ones(diffs(5),1)); %Column of Times
   apc = YTAB(:,4:17)*CONSTANTS; %Argument
   
   YY = sum(Asc.*Tc.*sin(apc) + Acc.*Tc.*cos(apc)); %Summation
 
   Y = -0.006951 ... 
       - (0.025896*(T_TDB)) ... 
       - (22.4072747*(T_TDB^2)) ... 
       + (0.00190059*(T_TDB^3)) ...
       + (0.00111252*(T_TDB^4)) ...
       + (0.0000001358*(T_TDB^5)) + YY; %["] Angle of CIP in ICRF
   
   %Solving for value of s
   diffs = (highs-lows)+[1 1 1 1 1];
  
   Asc = sTAB(:,2)/(1e6); %Sine Coefficient (Data in micro", converted to ")
   Acc = sTAB(:,3)/(1e6); %Cosine Coefficient (Data in micro", converted to ")
   Tc = vertcat(ones(diffs(1),1),T_TDB*ones(diffs(2),1),(T_TDB^2)*ones(diffs(3),1),(T_TDB^3)*ones(diffs(4),1),(T_TDB^4)*ones(diffs(5),1)); %Column of Times
   apc = sTAB(:,4:17)*CONSTANTS; %Argument
   
   ss = sum(Asc.*Tc.*sin(apc) + Acc.*Tc.*cos(apc)); %Summation
   
   %Converting " to rad
   X = (X/3600)*(pi/180); %[rad]
   Y = (Y/3600)*(pi/180); %[rad]

   sRHS = (pi/180)*(1/3600)*(0.000094 ...
          + (0.00380865*(T_TDB)) ...
          - (0.00012268*(T_TDB^2)) ... 
          - (0.07257411*(T_TDB^3)) ... 
          + ss ...
          + (0.00000173*(T_TDB)*sin(OM_lun)) ... 
          + (0.00000357*(T_TDB)*cos(2*OM_lun)) ... 
          + (0.00074352*(T_TDB^2)*sin(OM_lun)) ... 
          + (0.00005691*(T_TDB^2)*sin(2*F_lun - 2*D_sun + 2*OM_lun)) ... 
          + (0.00000984*(T_TDB^2)*sin(2*F_lun + 2*OM_lun)) ... 
          - (0.00000885*(T_TDB^2)*sin(2*OM_lun))); %[rad]

   sLHS = (X*Y)/2; %[rad]
   s = sRHS - sLHS; %[rad] CIO locator         

   %Parameter for rotation
   a = 0.5 + 0.125*((X^2) + (Y^2)); %[rad]

   %Rotation from ITRF to Pseudo-Earth Fixed (PEF)
   xp = xp*(pi/180)*(1/3600); %[rad] converting from "
   yp = yp*(pi/180)*(1/3600); %[rad] converting from "
   W1 = [cos(-sprime) sin(-sprime) 0; -sin(-sprime) cos(-sprime) 0; 0 0 1];
   W2 = [cos(xp) 0 -sin(xp); 0 1 0; sin(xp) 0 cos(xp)]; 
   W3 = [1 0 0; 0 cos(yp) sin(yp); 0 -sin(yp) cos(yp)];
   W  = W1*W2*W3;

   %Rotation from PEF to Intermediate Reference Frame of Epoch (IRE)
   ERAmat = [cos(-ERA) sin(-ERA) 0; -sin(-ERA) cos(-ERA) 0; 0 0 1];

   %Rotation from IRE to Geocentric Celestial Reference Frame (GCRF)
   PN1 = [(1-(a*(X^2))) (-a*X*Y) X; (-a*X*Y) (1-(a*(Y^2))) Y; (-X) (-Y) (1-(a*((X^2)+(Y^2))))];
   PN2 = [cos(s) sin(s) 0; -sin(s) cos(s) 0; 0 0 1];
   PN  = PN1*PN2;

   %Output Rotation Matrices for Both Ways
   ITRF2PEF = W;
   PEF2IRE  = ERAmat;
   IRE2GCRF = PN;
   GCRF2IRE = inv(PN);
   IRE2PEF  = transpose(ERAmat);
   PEF2ITRF = transpose(W);
   
end
