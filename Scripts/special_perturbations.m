
clearvars; clc;

% ORBIT PROPAGATION USING SPECIAL PERTURBATION TECHNIQUES
% ---------------------------------------------------------------------------------------------------------------------
% This script runs the "propagator_SP" function to propagate a satellite orbit from a set of initial conditions to a 
% specified Date and Time using Special Perturbation techniques to account for perturbing accelerations due to the 
% Earth's non-spherical mass distribution and the gravitational pull of the Sun and Moon. This script also uses the 
% "trajectory_multi_plot" function to visualize the perturbed orbit from initial state to the final state. Information 
% on the propagator settings, selected integrator, initial and final state, and initial and final orbital elements are 
% displayed in the command window after the trajectory is computed.
%
% Input state and time information can be entered in the "USER SPECIFIED INPUTS" section below the header. For testing 
% and demonstration, example inputs contained in the script for two satellites can be run in place of user inputs. 
%
% More detailed information on the force models used is contained in the "propagator_SP" header.
%
% Author: Matthew Buckhout
% Updated: 08/25/2020
%
% Inputs:
%
%     - [DATEi]            Initial Date                                             [yyyy mm dd]
%     - [UTCi]             Initial Time (UTC)                                       [hh mm ss]
%     - [Rix]              Initial Position Vector, x component (GCRF)              [km]
%     - [Riy]              Initial Position Vector, y component (GCRF)              [km]
%     - [Riz]              Initial Position Vector, z component (GCRF)              [km]
%     - [Vix]              Initial Velocity Vector, x component (GCRF)              [km/s]
%     - [Viy]              Initial Velocity Vector, y component (GCRF)              [km/s]
%     - [Viz]              Initial Velocity Vector, z component (GCRF)              [km/s]
%     - [DATEf]            Final Date                                               [yyyy mm dd]
%     - [UTCf]             Final Time (UTC)                                         [hh mm ss]
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
% Map Images:
%
%     world_color_image.jpg
%     carree_proj_hq.png 
%
% Functions:
%
%     - propagator_SP
%     - trajectory_multi_plot  
%     - julian_date
%     - orbital_elements 
%
% References:
%
%     - JPL HORIZONS Web Interface - Example Initial Conditions
%           https://ssd.jpl.nasa.gov/horizons.cgi
%
% ----------------------------------------------- USER SPECIFIED INPUTS -----------------------------------------------

%Initial State & Epoch
DATEi          = [0000 00 00]; %[yyyy mm dd] 
UTCi           = [00 00 00]; %[hh mm ss]
Rix            = 0; %[km]
Riy            = 0; %[km]
Riz            = 0; %[km]
Vix            = 0; %[km/s]
Viy            = 0; %[km/s]
Viz            = 0; %[km/s]

%Final Epoch
DATEf          = [0000 00 00]; %[yyyy mm dd]
UTCf           = [00 00 00]; %[hh mm ss]

%Settings
now            = 0; % -
dst            = 0; % -
L              = 0; % -
SPR            = 0; % -
perturbations  = [0 0]; % -

% ------------------------------------------- SELECTING INITIAL CONDITIONS --------------------------------------------

clc;
fprintf('\n');
fprintf('========================= SPECIAL PERTURBATIONS ========================\n');
fprintf('Propagate and visualize satellite orbit with graviational perturbations\n\n');
fprintf('Select input method:\n');
fprintf('   1. Use user specified initial conditions in script\n');
fprintf('   2. Use example initial conditions\n');
fprintf('\n');
method = input('Input Method: ');
fprintf('\n');

if (method == 1) %Use initial conditions entered in script

   fprintf('Using pre-loaded inputs\n');

   Rin = [Rix; Riy; Riz]; %[km]
   Vin = [Vix; Viy; Viz]; %[km/s]

elseif (method == 2) %Use example initial conditions
   
   %Example Settings
   now  = 0; 
   dst  = 1; 
   L    = 50; 
   SPR = 100; 
   perturbations = [1 2];
   
   fprintf('Select example satellite:\n');
   fprintf('   1. INTEGRAL\n');
   fprintf('   2. LAGEOS-2\n\n');
   satellite = input('Satellite: ');
       
   if (satellite == 1)

      %INTEGRAL - Example Initial Conditions from JPL HORIZONS
      TDBi = [02 00 00];
      TDB_UT1 = 69.185484;
      UT1_UTC = -0.2429197;
      UTCis = TDBi(1)*3600 - TDB_UT1 - UT1_UTC;
      DATEi = [2020 05 02];
      UTCi = s2hms(UTCis);
      Rin = [-4.639577443051213E+04; -9.038381297426131E+04; 1.135645308636836E+05];
      Vin = [3.915560093189137E-01; -2.541040034639090E-01; -3.153637151000829E-01];

      TDBf = [04 00 00];
      TDB_UT1 = 69.185439;
      UT1_UTC = -0.2445984;
      UTCfs = TDBf(1)*3600 - TDB_UT1 - UT1_UTC;
      DATEf = [2020 05 05];
      UTCf = s2hms(UTCfs);
      Ract = [-2.897316694533984E+04; -9.203320392725250E+04; 9.306987272688538E+04];
      Vact = [5.706971179596264E-01; 1.920381892971841E-01; -8.180831004553567E-01];

   elseif (satellite == 2) 
      
      %LAGEOS-2 - Example Initial Conditions from JPL HORIZONS
      TDBi = [08 00 00];
      TDB_UT1 = 66.184657;
      UT1_UTC = -0.2872785;
      UTCis = TDBi(1)*3600 - TDB_UT1 - UT1_UTC;
      DATEi = [2011 06 12];
      UTCi = s2hms(UTCis);
      Rin = [-5.136350436065671E+03; 1.011098841260996E+04; -3.962127458480597E+03];
      Vin = [-3.995354495114606E+00; -3.512908868979608E-01; 4.184177315369652E+00];
      
      TDBf = [12 00 00];
      TDB_UT1 = 66.184652;
      UT1_UTC = -0.2872785;
      UTCfs = TDBf(1)*3600 - TDB_UT1 - UT1_UTC;
      DATEf = [2011 06 12];
      UTCf = s2hms(UTCfs);
      Ract = [-8.518087487621629E+03; 8.508924437138130E+03; 7.604021293127074E+02];
      Vact = [-2.316436040016994E+00; -2.640366568112984E+00; 4.580317113520911E+00];
               
   else
      
      fprintf('Invalid Input\n\n');
      
   end

else
   
   fprintf('Invalid Input\n\n');

end

% ------------------------------------------------ GENERATING OUTPUTS -------------------------------------------------

if (method == 1 || method == 2)

   if (satellite == 1 || satellite == 2)

      %Generating Perturbed Ephemerides for Satellite
      [t,STATE,MJDi,MJDf] = propagator_SP(Rin,Vin,DATEi,UTCi,DATEf,UTCf,now,dst,L,SPR,perturbations);
      
      [~,MJDi_UTC] = julian_date(DATEi,UTCi,0);
      [~,MJDf_UTC] = julian_date(DATEf,UTCf,0);
      
      Rf = transpose(STATE(end,1:3)); %[km] Final Position Vector
      Vf = transpose(STATE(end,4:6)); %[km/s] Final Velocity Vector
      
      mu = 3.986004415e5; %[km^3/s^2] Earth Gravitational Parameter
      [~,ai,ei,ii,OMEGAi,omegai,thetai,~,~,~] = orbital_elements(mu,Rin,Vin);
      [~,af,ef,iff,OMEGAf,omegaf,thetaf,~,~,~] = orbital_elements(mu,Rf,Vf);

      fprintf('\n========================== RESULTS ==========================\n');
      
      %Display initial state
      fprintf('Initial State (GCRF):\n');
      fprintf('%f | %d-%02d-%02d | %02d:%02d:%02d UTC\n',MJDi_UTC,DATEi(1),DATEi(2),DATEi(3),UTCi(1),UTCi(2),fix(UTCi(3)));
      fprintf('   Rx = %12.5f km\n',Rin(1));   
      fprintf('   Ry = %12.5f km\n',Rin(2));
      fprintf('   Rz = %12.5f km\n',Rin(3));
      fprintf('   Vx = %12.5f km/s\n',Vin(1));
      fprintf('   Vy = %12.5f km/s\n',Vin(2));
      fprintf('   Vz = %12.5f km/s\n',Vin(3));

      %Display final state
      fprintf('\n');
      fprintf('Final State (GCRF):\n');
      fprintf('%f | %d-%02d-%02d | %02d:%02d:%02d UTC\n',MJDf_UTC,DATEf(1),DATEf(2),DATEf(3),UTCf(1),UTCf(2),fix(UTCf(3)));
      fprintf('   Rx = %12.5f km\n',Rf(1));   
      fprintf('   Ry = %12.5f km\n',Rf(2));
      fprintf('   Rz = %12.5f km\n',Rf(3));
      fprintf('   Vx = %12.5f km/s\n',Vf(1));
      fprintf('   Vy = %12.5f km/s\n',Vf(2));
      fprintf('   Vz = %12.5f km/s\n',Vf(3));
      fprintf('\n');
      
      %Display initial and final orbital elements
      fprintf('Orbital Elements               Initial          Final\n');
      fprintf('=============================================================\n'); 
      fprintf('Semi-Major Axis           %12.5f   %12.5f   km\n',ai,af);
      fprintf('Eccentricity              %12.5f   %12.5f   -\n',ei,ef);
      fprintf('Inclination               %12.5f   %12.5f   deg\n',rad2deg(ii),rad2deg(iff));
      fprintf('RAAN                      %12.5f   %12.5f   deg\n',rad2deg(OMEGAi),rad2deg(OMEGAf));
      fprintf('Argument of Periapsis     %12.5f   %12.5f   deg\n',rad2deg(omegai),rad2deg(omegaf));
      fprintf('True Anomaly              %12.5f   %12.5f   deg\n',rad2deg(thetai),rad2deg(thetaf));
      fprintf('=============================================================\n\n'); 
      
      %Plotting 3D Orbit, Groundtrack, Orbital Elements                                                 
      trajectory_multi_plot(t,STATE,MJDi,MJDf,DATEi,DATEf,UTCi,UTCf);

   end
   
end