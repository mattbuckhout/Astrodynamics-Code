
clearvars; clc;

% PREDICT LOOK ANGLES FOR SATELLITE FROM GROUND SITE
% ---------------------------------------------------------------------------------------------------------------------
% This script takes a Geocentric satellite orbit from input Two Line Element set (TLE) or state vector (Position and 
% Velocity) and determines look angles (Elevation and Azimuth) for the satellite from a specified ground site over a 
% time period of interest. The script prints a table containing the Date, Time, Range, Elevation, and Azimuth, as well 
% as the situation of the satellite with respect to the Sun and Earth for every time step through the period of 
% interest. Initial state from a TLE is extracted from mean elements fitting a series of observations and predictions 
% will be less accurate than using an actual satellite state measured at a specific epoch. In this iteration the 
% initial state is propagated using 2-body dynamics, and is accurate to around +/-2 degrees in elevation and azimuth 
% within 24 hours of the TLE epoch. Beyond 24 hours accuracy degrades quickly.
%
% Future improvements will include replacing the 2-body propagator with one incorporating General Perturbation 
% techniques to improve accuracy and extend the range of the predictions.
%
% In addition to printing the output table in the command window, two arrays are created in the workspace so the data
% generated can be used in other Matlab applications. First, a matrix "ANGLES_OUT" is generated which contains only 
% Modified Julian Dates (UT1), Range, and look angles. Second, a cell array "OUTPUT" is generated which contains all 
% the information displayed in the command window, including the satellite visibility situation and dates/times in the
% selected timezone.
%
%
% Author: Matthew Buckhout
% Updated: 09/09/2020 
%
% Inputs:
%
%     - [precision]        Length of X,Y,s coefficient sets for IAU2006 Reduction                   -
%                             'short' - Largest 10 terms in each series (lower accuracy, faster)
%                             'long'  - Full coefficient sets (greater accuracy, slower)
%     - [fulltable]        Format for printed output table:                                         - 
%                             0 - Print/store time steps for satellite above local horizon
%                             1 - Print/store all time steps
%     - [input_type]       Form of initial satellite state input:                                   - 
%                             1 - Input State via TLE 
%                             2 - Input state via Vectors                              
%     - [timezone]         Hours offset from UTC for selected timezone                             [hour]
%     - [dt]               Time step for output tables                                             [sec]
%     - [Lat]              Geodetic Latitude                                                       [deg]
%     - [Long]             Longitude (+East from Greenwich)                                        [deg]
%     - [h_ellip]          Ground site height above Reference Ellipsoid                            [km]
%     - [DATE0]            Initial Condition Date (selected timezone)                              [yyyy mm dd]
%     - [DATEs]            Period of Interest Start Date (selected timezone)                       [yyyy mm dd]
%     - [DATEe]            Period of Interest End Date (selected timezone)                         [yyyy mm dd]
%     - [TIME0]            Initial Condition Time (selected timezone)                              [hh mm ss]
%     - [TIMEs]            Period of Interest Start Time (selected timezone)                       [hh mm ss]
%     - [TIMEe]            Period of Interest End Time (selected timezone)                         [hh mm ss]
%
%
% Outputs:
%
%     OUTPUT Format (commmand window/workspace)
%     -----------------------------------------------------------------------------------------
%     MJD (UTC)    Date        Time (UTC-5)   Visibility     Range (km)    AZ (deg)    EL (deg)
%     -----------------------------------------------------------------------------------------
%
%     ANGLES_OUT Format (workspace)
%     ----------------------------------------------
%     MJD (UT1)   Range (km)    AZ (deg)    EL (deg)
%     ----------------------------------------------
%
% Visibility Situations:
%     
%     - Below Horizon : Satellite below local horizon (Elevation < 0)
%     - Illuminated   : Satellite above local horizon, Sun below local horizon, Sun above satellite horizon
%     - Dark          : Satellite above local horizon, Sun below local horizon, Sun below satellite horizon
%     - Day           : Satellite above local horizon, Sun above local horizon
%
% Data files required:
%
%     DE430t_1550_E.txt
%     DE430t_1650_E.txt
%     DE430t_1750_E.txt
%     DE430t_1850_E.txt
%     DE430t_1950_E.txt
%     DE430t_2050_E.txt
%     DE430t_2150_E.txt
%     DE430t_2250_E.txt
%     DE430t_2350_E.txt
%     DE430t_2450_E.txt
%     DE430t_2550_E.txt
%     IAU2006_X_coef.txt
%     IAU2006_Y_coef.txt
%     IAU2006_s_coef.txt
%     IAU2006_X_coef_short.txt
%     IAU2006_Y_coef_short.txt
%     IAU2006_s_coef_short.txt
%     IERS_EOP_combined.txt
%
% Functions:
%
%     - getdata_EOP_fast
%     - TZ2UTC
%     - TLE2state
%     - UTC2TZ
%     - julian_date
%     - EOP
%     - convert_time
%     - getdata_IAU2006
%     - terrestrial_MJD
%     - getdata_DE430
%     - propagation_UV
%     - IAU2006_rotations
%     - site
%
% References:
%
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 829-834)
%
% ----------------------------------------------- USER SPECIFIED INPUTS -----------------------------------------------

precision   = 'long'; %Precision for IAU2006 Reduciton
fulltable   = 1; %Table Format             
input_type  = 1; %Input State Form
timezone    = -5; %[hrs] Hours offset from UTC
dt          = 20; %[sec] Time step for output table

Lat         = 32.6513169; %[deg] Geodetic Latitude of observation site
Long        = -97.0572212; %[deg] Longitude of observation site (+East of Greewich)
h_ellip      = 0.13049; %[km] Height above Reference Ellipsoid

DATE0       = [2020 09 09]; %[yyyy mm dd] Initial Condition Start Date (selected timezone)
DATEs       = [2020 09 10]; %[yyyy mm dd] Period of Interest Start Date (selected timezone)                           
DATEe       = [2020 09 10]; %[yyyy mm dd] Period of Interest End Date (selected timezone)

TIME0       = [08 06 09]; %[hh mm ss] Initial Condition Time (selected timezone)
TIMEs       = [05 30 00]; %[hh mm ss] Period of Interest Start Time (selected timezone)
TIMEe       = [07 45 00]; %[hh mm ss] Period of Interest End Time (selected timezone)

%Satellite Initial State - Input Type 1

line1       = ['1 25544U 98067A   20253.54594774 -.00000354  00000-0  17323-5 0  9991'];
line2       = ['2 25544  51.6459 291.0011 0002027  91.0785  62.1264 15.49214270245093'];

%Satellite Initial State - Input Type 2

R0x         = -400.98; %[km] Initial Position Vector Components (GCRF) 
R0y         = 6345.6; %[km]
R0z         = 2401; %[km]
V0x         = -5.1975; %[km/s] Initial Velocity Vector Components (GCRF)
V0y         = 1.7016; %[km/s]
V0z         = -5.3614; %[km/s]

% ------------------------------------------ INITIAL STATE & TIME FORMATTING ------------------------------------------

%Input Earth Orientation Data File
fprintf('\nLoading Earth Orientation Data... \n');
EOPdata = getdata_EOP_fast; %Earth Orientation Parameters data file

%Convert to UTC
[DATE_UTCs,UTCs] = TZ2UTC(DATEs,TIMEs,timezone);
[DATE_UTCe,UTCe] = TZ2UTC(DATEe,TIMEe,timezone);

if (input_type == 1) %Initial state from TLE
   
   %Initial State and Time from TLE
   [p,a,e,i,OMEGA,omega,theta,M,R0,V0,DATE_UTC0,UTC0] = TLE2state(line1,line2);

   %Convert to Timezone
   [DATE0,TIME0] = UTC2TZ(DATE_UTC0,UTC0,timezone);

elseif (input_type == 2) % Initial state from Vectors
   
   %Convert to UTC
   [DATE_UTC0,UTC0] = TZ2UTC(DATE0,TIME0,timezone);
   R0 = [R0x;R0y;R0z]; %[km] Initial Position Vector (GCRF)
   V0 = [V0x;V0y;V0z]; %[km/s] Initial Velocity Vector (GCRF)
   
end
   
%Modified Julian Date (UTC)
[~,MJD0] = julian_date(DATE_UTC0,UTC0,0);  
[~,MJDs] = julian_date(DATE_UTCs,UTCs,0);  
[~,MJDe] = julian_date(DATE_UTCe,UTCe,0); 
MJD_UTC = transpose(MJDs:(dt/86400):MJDe+dt); %(UTC)  Array of MJD at each time step

if (MJD0 > MJDs)
   clc;
   fprintf('\n***Initial Condition Epoch after Start Period of Interest***\n\n');
else

   %Earth Orientation Parameters
   [~,~,dUT10,dAT0] = EOP((MJD0),EOPdata);
   [~,~,dUT1s,dATs] = EOP((MJDs),EOPdata);                       
   [~,~,dUT1e,dATe] = EOP((MJDe),EOPdata);

   %Times (UT1)
   [UT10,~,~,~,~,~,~] = convert_time(DATE_UTC0,UTC0,dUT10,dAT0);  
   [UT1s,~,~,~,~,~,~] = convert_time(DATE_UTCs,UTCs,dUT1s,dATs);
   [UT1e,~,~,~,~,~,~] = convert_time(DATE_UTCe,UTCe,dUT1e,dATe);

   %Modified Julian Date (UT1)
   [~,MJD0] = julian_date(DATE_UTC0,UT10,0);  
   [~,MJDs] = julian_date(DATE_UTCs,UT1s,0);  
   [~,MJDe] = julian_date(DATE_UTCe,UT1e,0); 
   MJD = transpose(MJDs:(dt/86400):MJDe); %(UT1)  Array of MJD at each time step

   %Input Data Sets
   fprintf('Loading Coefficients... \n');
   [XTAB,YTAB,sTAB] = getdata_IAU2006(precision); 
   MJDs_TT = terrestrial_MJD(MJDs,EOPdata); 
   MJDe_TT = terrestrial_MJD(MJDe,EOPdata); 
   DE430coef = getdata_DE430(MJDs_TT,MJDe_TT); 

   %Constants
   RE = 6378.1363; %[km]               Earth Mean Equatorial Radius 
   mu = 3.986004415e5; %[km^3/s^2]     Earth Gravitational Parameter
   Mmoon = 7.34767309e22; %[kg]        Mass of Moon
   Mearth = 5.972e24; %[kg]            Mass of Earth

   %Preallocating
   OUTPUT = cell([numel(MJD) 11]);

   Latrad = deg2rad(Lat); %[rad]
   Longrad = deg2rad(Long); %[rad]

   [Rs,Vs] = propagation_UV(mu,R0,V0,MJD0*86400,MJDs*86400); %[km] [km/s] Satellite State, Start of period of interest

   clc;
   fprintf('\n');
   if (timezone >= 0)
      fprintf('Initial Condition Epoch:   %d/%02d/%02d | %02d:%02d:%02d (UTC+%d)\n',DATE0(1),DATE0(2),DATE0(3),TIME0(1),TIME0(2),fix(TIME0(3)),timezone);
      fprintf('Start Period of Interest:  %d/%02d/%02d | %02d:%02d:%02d (UTC+%d)\n',DATEs(1),DATEs(2),DATEs(3),TIMEs(1),TIMEs(2),fix(TIMEs(3)),timezone);
      fprintf('End Period of Interest:    %d/%02d/%02d | %02d:%02d:%02d (UTC+%d)\n',DATEe(1),DATEe(2),DATEe(3),TIMEe(1),TIMEe(2),fix(TIMEe(3)),timezone);
   elseif (timezone < 0)
      fprintf('Initial Condition Epoch:   %d/%02d/%02d | %02d:%02d:%02d (UTC%d)\n',DATE0(1),DATE0(2),DATE0(3),TIME0(1),TIME0(2),fix(TIME0(3)),timezone);
      fprintf('Period of Interest:        %d/%02d/%02d | %02d:%02d:%02d (UTC%d)\n',DATEs(1),DATEs(2),DATEs(3),TIMEs(1),TIMEs(2),fix(TIMEs(3)),timezone);
      fprintf('                           %d/%02d/%02d | %02d:%02d:%02d (UTC%d)\n',DATEe(1),DATEe(2),DATEe(3),TIMEe(1),TIMEe(2),fix(TIMEe(3)),timezone);
   end
   fprintf('\n');
   fprintf('Latitude   =  %12.7f degrees\n',Lat);
   fprintf('Longitude  =  %12.7f degrees\n',Long);
   fprintf('\n');

   lineformat = '-----------------------------------------------------------------------------------------\n';
   fprintf(lineformat);
   if (timezone >= 0)
      fprintf('MJD (UTC)    Date        Time (UTC+%d)   Visibility     Range (km)    AZ (deg)    EL (deg)\n',timezone);
   elseif (timezone < 0)
      fprintf('MJD (UTC)    Date        Time (UTC%d)   Visibility     Range (km)    AZ (deg)    EL (deg)\n',timezone);
   end
   fprintf(lineformat);
   format1 = '%0.4f   %d/%02d/%02d     %02d:%02d:%02d    Day             %9.3f     %7.3f     %7.3f\n';
   format2 = '%0.4f   %d/%02d/%02d     %02d:%02d:%02d    Illuminated     %9.3f     %7.3f     %7.3f\n';
   format3 = '%0.4f   %d/%02d/%02d     %02d:%02d:%02d    Dark            %9.3f     %7.3f     %7.3f\n';
   format4 = '%0.4f   %d/%02d/%02d     %02d:%02d:%02d    Below Horizon   %9.3f     %7.3f     %7.3f\n';

   % --------------------------------------------- GENERATE OUTPUT TABLE ----------------------------------------------
   
   %Loop through period of interest
   linecount = 0;
   for i=1:1:numel(MJD) 
      
      [xp,yp,dUT1,dAT] = EOP((MJD(i)),EOPdata); %["] [sec]                             %Earth Orientation Parameters at step
      [R_GCRF,V_GCRF] = propagation_UV(mu,R0,V0,MJD0*86400,MJD(i)*86400); %[km] [km/s]    %Satellite State at step (GCRF)                      
      [ITRF2PEF,PEF2IRE,IRE2GCRF,GCRF2IRE,IRE2PEF,PEF2ITRF] = ... 
                  IAU2006_rotations(MJD(i),dUT1,dAT,xp,yp,precision,XTAB,YTAB,sTAB);      %Rotation Matrices between ITRF and GCRF
      R_ITRF = PEF2ITRF*IRE2PEF*GCRF2IRE*R_GCRF; %[km]                                    %Satellite Position Vector (ITRF)  
      [Rsite_ITRF] = site(Latrad,Longrad,h_ellip); %[km]                                        %Observation Site Position Vector (ITRF)
      Rsite_GCRF = IRE2GCRF*PEF2IRE*ITRF2PEF*Rsite_ITRF; %[km]                            %Observation Site Position Vector (GCRF)
      range_ITRF = R_ITRF - Rsite_ITRF; %[km]                                             %Range Vector (ITRF)
      range_SEZ = [cosd(90-Lat) 0 -sind(90-Lat); 0 1 0; sind(90-Lat) 0 cosd(90-Lat)] ...
                 *[cosd(Long) sind(Long) 0; -sind(Long) cosd(Long) 0; 0 0 1] ...
                 *range_ITRF; %[km]                                                       %Range Vector (SEZ)                                                                     
      
      %Check if Satellite is above local horizon
      if (range_SEZ(3) > 0) 
         
         MJD_TT = terrestrial_MJD(MJD(i),EOPdata); %(TT)                                     %Modified Julian Date at step
         ephemerides = planet_ephemerides_c(MJD_TT,transpose(DE430coef),[3 10 11]); %[km]    %Position Vectors, Earth-Moon Barycenter, Moon, Sun
         Rembc_ICRF = transpose(ephemerides(1,:)); %[km]                                     %Earth-Moon Barycenter (ICRF)
         Rmoon_GCRF = transpose(ephemerides(2,:)); %[km]                                     %Moon Position Vector (GCRF)
         Rsun_ICRF = transpose(ephemerides(3,:)); %[km]                                      %Sun Position Vector (ICRF)
         Rsun_GCRF = Rsun_ICRF - Rembc_ICRF + Rmoon_GCRF*(Mmoon/(Mearth+Mmoon)); %[km]       %Sun Position Vector (GCRF)
         
         %Check if satellite and site illuminated
         if (dot(Rsun_GCRF,Rsite_GCRF) > 0)
            
            condition = 1; %Satellite above horizon, Satellite Illuminated, Site Illuminated
            if (i == 1)
               satset = 0;
            end
                     
         else
            
            %Angle between Sun and Satellite Position Vectors (GCRF)
            zeta = asind(norm(cross(Rsun_GCRF,R_GCRF))/(norm(Rsun_GCRF)*norm(R_GCRF))); %[deg] 
            
            %Perpendicular distance from Sun position vector to satellite
            dist = norm(R_GCRF)*cosd(zeta-90); %[km] 
            
            %Check if satellite is illuminated
            if (dist > RE)
               
               condition = 2; %Satellite above horizon, Satellite Illuminated, Site Dark
               if (i == 1)
                  satset = 0;
               end
               
            else
               
               condition = 3; %Satellite above horizon, Satellite in Shadow, Site Dark
               if (i == 1)
                  satset = 0;
               end
            
            end
            
         end   
         
      else
       
         condition = 4; %Satellite not above site local horizon
         if (i == 1)
            satset = 1;
         end
       
      end
      
      %Angles
      if ((condition == 1) || (condition == 2) || (condition == 3) || (condition == 4))
      
         range_mag = norm(range_SEZ); %[km]
         
         EL = asin(range_SEZ(3)/range_mag); %[rad] Elevation
         
         if (EL ~= (pi/2))
         
            SINB = range_SEZ(2)/sqrt((range_SEZ(1)^2) + (range_SEZ(2)^2));
            COSB = range_SEZ(1)/sqrt((range_SEZ(1)^2) + (range_SEZ(2)^2));
            AZ = abs(atan(SINB/COSB)); %[rad] Azimuth
         
         elseif (EL == (pi/2))

            SINB = drange_SEZ(2)/sqrt((drange_SEZ(1)^2) + (drange_SEZ(2)^2));
            COSB = drange_SEZ(1)/sqrt((drange_SEZ(1)^2) + (drange_SEZ(2)^2));
            AZ = abs(atan(SINB/COSB)); %[rad] Preliminary Azimuth

         end

         %Azimuth Quadrant Checks
         if ((COSB > 0) && (SINB > 0))
            AZ = pi - AZ;
         elseif ((COSB > 0) && (SINB < 0))
            AZ = pi + AZ;
         elseif ((COSB < 0) && (SINB < 0))
            AZ = 2*pi - AZ;
         elseif ((COSB < 0) && (SINB > 0))
            AZ = AZ;
         end
      
      end
      
      %Date and Time
      [yr,mo,d,h,m,s] = gregorian_date(MJD_UTC(i)); %(UTC)
      [DATE,TZ] = UTC2TZ([yr mo d],[h m s],timezone); %(TZ) Selected Timezone
      
      %Output values into Cell Array
      OUTPUT(i,1) = MJD_UTC(i);
      OUTPUT(i,2) = DATE(1);
      OUTPUT(i,3) = DATE(2);
      OUTPUT(i,4) = DATE(3);
      OUTPUT(i,5) = TZ(1);
      OUTPUT(i,6) = TZ(2);
      OUTPUT(i,7) = TZ(3);
      
      %Output values into Matrix 
      ANGLES_OUT(i,1) = MJD(i);
      ANGLES_OUT(i,2) = range_mag;
      ANGLES_OUT(i,3) = AZ;
      ANGLES_OUT(i,4) = EL;
      

      %Outputs
      if (condition == 1)
         if (satset == 1)
            satset = 0;
            if (fulltable == 1)
               if (linecount > 0)
                  fprintf(lineformat);
               end
            end
         end
         OUTPUT(i,1) = MJD_UTC(i);
         OUTPUT(i,2) = DATE(1);
         OUTPUT(i,3) = DATE(2);
         OUTPUT(i,4) = DATE(3);
         OUTPUT(i,5) = TZ(1);
         OUTPUT(i,6) = TZ(2);
         OUTPUT(i,7) = TZ(3);
         OUTPUT(i,9) = range_mag;
         OUTPUT(i,10) = AZ;
         OUTPUT(i,11) = EL;
         OUTPUT(i,8) = 'Day';
         fprintf(format1,MJD_UTC(i),DATE(1),DATE(2),DATE(3),TZ(1),TZ(2),fix(TZ(3)),range_mag,rad2deg(AZ),rad2deg(EL));
         linecount = linecount + 1;

      elseif (condition == 2)
         if (satset == 1)
            satset = 0;
            if (linecount > 0)
               fprintf(lineformat);
            end      
         end
         OUTPUT(i,1) = MJD_UTC(i);
         OUTPUT(i,2) = DATE(1);
         OUTPUT(i,3) = DATE(2);
         OUTPUT(i,4) = DATE(3);
         OUTPUT(i,5) = TZ(1);
         OUTPUT(i,6) = TZ(2);
         OUTPUT(i,7) = TZ(3);
         OUTPUT(i,8) = 'Illuminated';
         OUTPUT(i,9) = range_mag;
         OUTPUT(i,10) = AZ;
         OUTPUT(i,11) = EL;
         fprintf(format2,MJD_UTC(i),DATE(1),DATE(2),DATE(3),TZ(1),TZ(2),fix(TZ(3)),range_mag,rad2deg(AZ),rad2deg(EL));
         linecount = linecount + 1;
      
      elseif (condition == 3)
         if (satset == 1)
            satset = 0;
            if (linecount > 0)
               fprintf(lineformat);
            end      
         end
         OUTPUT(i,1) = MJD_UTC(i);
         OUTPUT(i,2) = DATE(1);
         OUTPUT(i,3) = DATE(2);
         OUTPUT(i,4) = DATE(3);
         OUTPUT(i,5) = TZ(1);
         OUTPUT(i,6) = TZ(2);
         OUTPUT(i,7) = TZ(3);
         OUTPUT(i,9) = range_mag;
         OUTPUT(i,10) = AZ;
         OUTPUT(i,11) = EL;
         OUTPUT(i,8) = 'Dark';
         fprintf(format3,MJD_UTC(i),DATE(1),DATE(2),DATE(3),TZ(1),TZ(2),fix(TZ(3)),range_mag,rad2deg(AZ),rad2deg(EL));
         linecount = linecount + 1;
         
      elseif (condition == 4)
         if (satset == 0)
            satset = 1;
            if (fulltable == 1)
               if (linecount > 0)
                  fprintf(lineformat);
               end
            end      
         end
         OUTPUT(i,1) = MJD_UTC(i);
         OUTPUT(i,2) = DATE(1);
         OUTPUT(i,3) = DATE(2);
         OUTPUT(i,4) = DATE(3);
         OUTPUT(i,5) = TZ(1);
         OUTPUT(i,6) = TZ(2);
         OUTPUT(i,7) = TZ(3);
         OUTPUT(i,8) = 'Below Horizon';
         if (fulltable == 1)
            fprintf(format4,MJD_UTC(i),DATE(1),DATE(2),DATE(3),TZ(1),TZ(2),fix(TZ(3)),range_mag,rad2deg(AZ),rad2deg(EL));
            linecount = linecount + 1;
         end
      end

   end
   fprintf(lineformat);

   clear EOPdata;
end