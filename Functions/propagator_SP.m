
% SPECIAL PERTURBATIONS PROPAGATOR
% ---------------------------------------------------------------------------------------------------------------------
% Function numerically integrates satellite 3-DOF equations of motion for a geocentric orbit, accounting for perturbing 
% accelerations due to conservative forces (Non-spherical Central Body Gravitational Potential & Third Body 
% Accelerations). 
%
% Central body perturbations are developed using the EGM2008 model, containing spherical harmonic coefficients that 
% together with Associated Legendre Polynomials define Earth's Non-Spherical Gravitational Potential. Accelerations are 
% developed using partial derivatives of the potential function. The resulting accelerations account for perturbing 
% effects due to zonal, tesseral, and sectorial harmonics for an L x L field, where L is specified by the user and 
% corresponds to the maximum spherical harmonic degree/order of coefficients from the EGM2008 model.
%
% Third body perturbations are developed using the JPL DE430 model for planetary ephemerides. The DE430 model yields 
% highly accurate position vectors of solar system bodies resulting from numerical simulations of the system. The 
% position vectors of the Sun and Moon are used to compute the perturbing accelerations of each on the satellite.
%
% Several computationally intense processes in force models such as computing the B series term in the third body 
% acceleration equation and computing Associated Legendre Polynomials were implemented as C code and compiled as 
% Matlab Executables (.mex) files for increased speed.
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
%     EGM2008_Norm_Coefficients.txt
%     IAU2006_X_coef_short.txt
%     IAU2006_Y_coef_short.txt
%     IAU2006_s_coef_short.txt
%     IERS_EOP_combined.txt
%
%
% Author: Matthew Buckhout
% Updated: 08/24/2020 
%
% Inputs:
%
%     - [Rin]              Initial Position Vector (GCRF)                           [km]
%     - [Vin]              Initial Velocity Vector (GCRF)                           [km/s]
%     - [DATEi]            Initial Date (for plot label)                            [yyyy mm dd]
%     - [UTCi]             Initial Time (UTC) (for plot label)                      [hh mm ss]
%     - [DATEf]            Final Date (for plot label)                              [yyyy mm dd]
%     - [UTCf]             Final Time (UTC) (for plot label)                        [hh mm ss]  
%     - [now]              Propagate to Current Time (1) or Input Final Time (0)     -
%     - [dst]              Daylight Savings Time in effect for local timezone (1)    -
%     - [L]                Degree x Order, Field for Spherical Harmonics             -
%                             (Non-spherical Gravitational Potential)                      
%     - [SPR]              Time steps per revolution (Reccomend ~100)                -
%     - [perturbations]    Vector specifying the perturbing forces to include        -
%                             in simulation: (ex: [1 2])
%                                1 - Earth Non-spherical Gravitational Potential
%                                2 - Third Body Accelerations (Sun/Moon)
%
% Outputs:
%
%     - [t]                Array of Time Values, +from Initial Date/Time            [sec]
%     - [STATE]            Matrix of Position and Velocity Vectors                  [km][km/s]
%     - [MJDi]             Modified Julian Date, Initial (UT1)                      [days]
%     - [MJDf]             Modified Julian Date, Final (UT1)                        [days]
%
% Format for t and STATE:
%
%       t       STATE
%       -----   ------------------------------------------
%         t1      R1x    R1y    R1z    V1x    V1y    V1z
%         t2      R2x    R2y    R2z    V2x    V2y    V2z
%         t3      R3x    R3y    R3z    V3x    V3y    V3z
%         t4      R4x    R4y    R4z    V4x    V4y    V4z
%         t5      R5x    R5y    R5z    V5x    V5y    V5z
%
%                               . . .
%
%         tf      Rfx    Rfy    Rfz    Vfx    Vfy    Vfz
%       -----   ------------------------------------------
%
% Functions:
%
%     - julian_date
%     - gregorian_date
%     - getdata_EOP_fast
%     - getdata_EGM2008
%     - getdata_IAU2006
%     - YMDHMS2sec
%     - EOP
%     - terrestrial_MJD
%     - orbital_elements
%     - getdata_DE430
%     - perturbed_eom
%     - RK4
%     - RKF45
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 491-564)
%     - Full IERS Conventions 2010
%           https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html
%     - EGM2008 Spherical Harmonic Coefficients (2190 x 2190)
%           https://earth-info.nga.mil/GandG/update/index.php?action=home#tab_wgs84-data
%     - Recurrence Relations for Associated Legendre Polynomials
%           https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
%     - JPL DE430 Planetary Ephemerides Coefficient Files
%           https://ssd.jpl.nasa.gov/?planet_eph_export 
% ---------------------------------------------------------------------------------------------------------------------

function [t,STATE,MJDi,MJDf] = propagator_SP(Rin,Vin,DATEi,UTCi,DATEf,UTCf,now,dst,L,SPR,perturbations)

   %Constants
   RE = 6378.1363; %[km] Earth Mean Equatorial Radius 
   mu = 3.986004415e5; %[km^3/s^2] Earth Gravitational Parameter

   % --------------------------------------------- TIME MANAGEMENT ---------------------------------------------------- 

   %Current Date and Time & Daylight Savings Time
   if (dst == 1)
      addh = 5; %Central Time
   else
      addh = 6; %Central Time
   end
   if (now == 1)
      
      c = clock;
      DATEf = [c(1) c(2) c(3)];
      UTCf = [(c(4)+addh) c(5) c(6)];
      
      %Time Rollover, Target Final Epoch - UTC
      [JD_UTCf,MJD_UTCf] = julian_date(DATEf,UTCf,0); %Modified Julian Date (UTC)
      [yrf,mof,df,hf,mf,sf] = gregorian_date(MJD_UTCf); %Date/Time (UTC)
      DATEf = [yrf mof df]; 
      UTCf  = [hf mf sf];
      
   end

   % -------------------------------------------- LOADING DATA SETS ---------------------------------------------------
   
   clc;
   fprintf('\nLoading Physical Data... \n');
   EOPdata = getdata_EOP_fast; %Earth Orientation Data Table
   EGM2008coef = getdata_EGM2008(L); %Coefficients for EGM2008 Gravitational Model
   [XTAB,YTAB,sTAB] = getdata_IAU2006('short'); %Coefficients for IAU2006 Precession-Nutation Model

   % --------------------------------------------- INPUT HANDLING -----------------------------------------------------

   %Vectors as Columns
   if (isrow(Rin) ~= 0) 
      Rin = transpose(Rin);
      Vin = transpose(Vin);
   end

   %Target Epoch in Seconds
   [secondsi] = YMDHMS2sec(DATEi,UTCi); %[sec] Seconds Elapsed in Current Year to Initial
   [secondsf] = YMDHMS2sec(DATEf,UTCf); %[sec] Seconds Elapsed in Current Year to Final
   tf = secondsf - secondsi; %[sec] Final Time 

   %Converting to Modified Julian Date - Initial Time
   [~,MJDi_UTC] = julian_date(DATEi,UTCi,0); %Modified Julian Date (UTC)
   [~,~,dUT1,dAT] = EOP(fix(MJDi_UTC),EOPdata); %Initial Earth Orientation Data 
   MJDi = MJDi_UTC + (dUT1/86400); %Modified Julian Date (UT1)

   %Converting to Modified Julian Date - Final Time
   [~,MJDf] = julian_date(DATEf,UTCf,0); %Modified Julian Date (UTC)
   [~,~,dUT1,~] = EOP(fix(MJDf),EOPdata); %Initial Earth Orientation Data 
   MJDf = MJDf + (dUT1/86400); %Modified Julian Date (UT1)

   %Terrestrial Time MJD for DE430 Coefficients
   [MJD_TTi] = terrestrial_MJD(MJDi,EOPdata); %Input UT1
   [MJD_TTf] = terrestrial_MJD(MJDf,EOPdata); %Input UT1

   %Initial Orbital Elements
   [p,a,e,i,OMEGA,omega,theta,omega_true,lambda_true,u] = orbital_elements(mu,Rin,Vin);

   %Orbit Period
   T = 2*pi*sqrt((a^3)/mu); %[sec]

   %Coefficients for Planetary Ephemerides
   [DE430coef] = getdata_DE430(MJD_TTi,MJD_TTf); %Data Input
   
   %Printing Simulation Info
   clc;
   fprintf('\n');
   fprintf('======================= PROPAGATOR INFO =====================\n');
   if (isempty(perturbations) == 1)
      fprintf('2 Body Motion - No Perturbations\n');
   elseif (ismember(1,perturbations) == 0 && ismember(2,perturbations) == 0)
      fprintf('2 Body Motion - No Perturbations\n');
   else   
      fprintf('Force Models:\n\n');
      fprintf('   - Central Body Acceleration (2-body term)\n');
   end
   if (ismember(1,perturbations) == 1)
      fprintf('   - EGM2008 %d x %d Non-spherical Gravitational Potential\n',L,L);   
   end
   if (ismember(2,perturbations) == 1)
      fprintf('   - Third Body Sun & Moon w/ DE430 Ephemerides\n');
   end

   % ----------------------------------------- NUMERICAL INTEGRATION --------------------------------------------------

   %Integrator Inputs
   R0_GCRF = [Rin]; %[km] Initial Position
   V0_GCRF = [Vin]; %[km/s] Initial Velocity
   tspan = [0 tf]; %[sec] Time Span for Iteration
   step = fix(T/SPR); %[sec] Time Step
   IC = [R0_GCRF; V0_GCRF]; %Initial Conditions (Column)

   if (e < 0.02)
      
      %4th Order Fixed Step Runge-Kutta
      integrator = 1;
      t0 = tspan(1); %[sec] Initial Time
      tf = tspan(2); %[sec] Final Time
      t = [t0:step:tf]'; %[sec] Time vector
      h = (tf-t0)/numel(t); %[sec] Refined Time Step
      fprintf('\n');
      fprintf('Integrator:\n\n');
      fprintf('   4th Order Runge-Kutta - Fixed Step ');
      fprintf('(%0.2f seconds)\n',h);
      fprintf('=============================================================\n\n');
      fprintf('Computing Trajectory... \n\n');
      [t,STATE] = RK4(@perturbed_eom,tspan,transpose(IC),step,MJDi,T,EOPdata,XTAB,YTAB,sTAB,EGM2008coef, ...
                      DE430coef,L,perturbations); 

   elseif (e >= 0.02)
   
      %4th/5th Order Variable Step Runge-Kutta-Fehlberg
      integrator = 2;
      fprintf('\n');
      fprintf('Integrator:\n\n');
      fprintf('   4/5 Order Runge-Kutta-Fehlberg - Variable Step\n');
      fprintf('=============================================================\n\n');
      fprintf('Computing Trajectory... \n\n');
      [t,STATE] = RKF45(@perturbed_eom,tspan,transpose(IC),step,MJDi,T,EOPdata,XTAB,YTAB,sTAB,EGM2008coef, ...
                        DE430coef,L,perturbations);
 
   end
 
   %Converting from Modified Julian Date - Final Time
   MJDf = MJDi + (t(end)/86400); %(UT1)
   [xpf,ypf,dUT1f,dATf] = EOP(fix(MJDf),EOPdata); %(UT1 - UTC)
   MJDf_UTC = MJDi + ((t(end) - dUT1f)/86400); %(UTC)
   [yrf,mof,df,hf,mf,sf] = gregorian_date(MJDf_UTC); %Gregorian Date/Time (UTC)
   DATEf = [yrf mof df]; %Final Date (UTC)
   UTCf = [hf mf sf]; %Final Time (UTC)
   
   %Re-Printing Simulation Info
   clc;
   fprintf('\n');
   fprintf('======================= PROPAGATOR INFO =====================\n');
   if (isempty(perturbations) == 1)
      fprintf('2 Body Motion - No Perturbations\n');
   elseif (ismember(1,perturbations) == 0 && ismember(2,perturbations) == 0)
      fprintf('2 Body Motion - No Perturbations\n');
   else   
      fprintf('Force Models:\n\n');
      fprintf('   - Central Body Acceleration (2-body term)\n');
   end
   if (ismember(1,perturbations) == 1)
      fprintf('   - EGM2008 %d x %d Non-spherical Gravitational Potential\n',L,L);   
   end
   if (ismember(2,perturbations) == 1)
      fprintf('   - Third Body Sun & Moon w/ DE430 Ephemerides\n');
   end
   fprintf('\n');
   fprintf('Integrator:\n\n');
   if (integrator == 1) 
      fprintf('   4th Order Runge-Kutta - Fixed Step (%0.2f seconds)\n',h);
   elseif (integrator == 2)
      fprintf('   4/5 Order Runge-Kutta-Fehlberg - Variable Step\n');
   end
   fprintf('=============================================================\n');
   fprintf('Simulation complete\n');
   fprintf('%0.2f Revolutions\n\n',tf/T);
   
   %Clearing data sets                                                                     
   clear EGM2008coef;
   clear EOPdata;
   clear XTAB;
   clear YTAB;
   clear sTAB;

end
