
% ORBIT VISUALIZATION - 3D ORBIT PLOT
% ---------------------------------------------------------------------------------------------------------------------
% Function takes an initial satellite state vector (R,V) in a Geocentric-Equatorial frame and the corresponding Date
% and Time (UTC) and generates a plot of the orbit through 1 period. The generated plot gives a rough representation 
% of an orbit for qualitative analysis. Determined orbit is defined by 2-body dynamics defines orbital motion and 
% orbits must be circular or 
% elliptical.
%
%
% Author: Matthew Buckhout
% Updated: 08/24/2020 
%
% Inputs:
%
%     - [R]                Position Vector at Date/Time (GCRF/ICRF)        [km]
%     - [V]                Velocity Vector at Date/Time (GCRF/ICRF)        [km/s]
%     - [DATE]             Gregorian Date                                  [yyyy mm dd]
%     - [TIME]             Coordinated Universal Time (UTC)                [hh mm ss]
%     - [centralbody]      Central Attracting Body ('earth' or 'sun')       -
%                             centralbody = 'earth' for Geocentric orbit
%                             centralbody = 'sun' for Heliocentric orbit
%
% Functions:
%
%     - julian_date
%     - getdata_EOP_fast
%     - EOP
%     - orbital_elements
%     - position_velocity
%     - convert_time
%     - local_sidereal
%
% References:
%     -  3D Earth Example [Ryan Gray - 2004,2006,2013]
%        https://www.mathworks.com/matlabcentral/fileexchange/13823-3d-earth-example
% ---------------------------------------------------------------------------------------------------------------------

function orbit_vis_1(R,V,DATE,TIME,centralbody)

   switch centralbody
      case 'earth'
         mu = 3.986004415e5; %[km^3/s^2] Earth Gravitational Parameter
         [~,MJD] = julian_date(DATE,TIME,0); %Modified Julian Date (UTC)
         [EOPdata] = getdata_EOP_fast; %Earth Orientation Parameter Data
         [~,~,dUT1,dAT] = EOP(MJD,EOPdata); %Earth Orientation Paramters
         
      case 'sun'
         mu = 1.32712428e11; %[km^3/s^2] Sun Gravitational Parameter
         
   end
   
   % -------------------------------------------  G E N E R A T E  D A T A  -------------------------------------------
   
   [p,a,e,i,OMEGA,omega,theta,omega_true,lambda_true,u] = orbital_elements(mu,R,V); %Orbital Elements

   if (e < 1) %Circular and Elliptical Orbits
      
      jnum = 100; %Number of state vectors to generate
      for j=1:1:jnum+1 %Loop for 1 period
         
         theta1 = j*(2*pi/jnum) - 2*pi/jnum; %Iterated True Anomaly
    
         [RR,VV] = position_velocity(mu,p,e,i,OMEGA,omega,theta1,omega_true,lambda_true,u); 
         STATE(j,1:3) = RR; %[km] Position Vector
         STATE(j,4:6) = VV; %[km/s] Velocity Vector
         
      end

      [Ra,Va] = position_velocity(mu,p,e,i,OMEGA,omega,pi,omega_true,lambda_true,u); 

      Rp = STATE(1,1:3); %[km] Periapsis Vectors (theta = 0)
      Vp = STATE(1,4:6); %[km/s]

   elseif (e >= 1)

      fprintf('Orbit is Hyperbolic\n\n');
      
   end

   % -------------------------------------------  3 D  O R B I T  P L O T  --------------------------------------------
 
   set(groot,'defaultfigureposition',[700 200 900 750])
   figure('Color','k');

   switch centralbody
      
      case 'earth'
      
         %Plot Textured Ellipsoid Earth Model [modified from Original Code: Ryan Gray 2004,2006,2013]
         %https://www.mathworks.com/matlabcentral/fileexchange/13823-3d-earth-example
         
         %Time and Earth Rotation
         [UT1,~,~,~,~,~,~] = convert_time(DATE,TIME,dUT1,dAT); %Universal Time
         [~,MJD_UT1] = julian_date(DATE,UT1,0); %Modified Julian Date (UT1)
         [~,GMST] = local_sidereal(0,MJD_UT1); %[rad] %Greenwich Mean Standard Time

         %Earth Ellipsoid Model
         RE = 6378.1363; %[km] Earth Mean Equatorial Radius
         bE = 6356.7516005; %[km] Earth Semi-minor Axis (Polar Radius)(Derived)
         eE = 0.081819221456; % Earth Reference Ellipsoid Eccentricity
         npanels = 180; % Number of panels around the equator 
         alpha   = 1; % Transparency
         [x_refellipse, y_refellipse, z_refellipse] = ellipsoid(0, 0, 0, RE, RE, bE, npanels); %Generates Ellipsoid
         globe = surf(x_refellipse, y_refellipse, -z_refellipse,'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
         rotate(globe,[0 0 1],(rad2deg(GMST)+180)); %Rotates Earth to correct GMST (-X axis through Greenwich at start)

         %Texture Map, Projecting World Map onto Ellipsoid
         [world,map] = imread('equirectangular_1.jpg');
         world = im2double(world);
         set(globe, 'FaceColor', 'texturemap', 'CData', world,'FaceAlpha', alpha, 'EdgeColor', 'none');

         %Plot Formatting
         grid on;
         set(gca,'Color','k');
         set(gca,'GridColor',[1 1 1]);
         axis equal;
         
         %Clearing unnecessary data
         clear world;
         clear map;
         hold on
   
         %Plotting Occupied Focus
         plot3(0,0,0,'ko','MarkerFaceColor','k','MarkerSize',10);
         hold on
         %Plot of Orbit 1
         plot3(STATE(:,1),STATE(:,2),STATE(:,3),'LineWidth',2);
         hold on
         %Plot Periapsis
         plot3(Rp(1),Rp(2),Rp(3),'bo','MarkerSize',8);
         hold on
         %Plot Apoapsis
         plot3(Ra(1),Ra(2),Ra(3),'bd','MarkerSize',10);
         hold on
   
      case 'sun'
      
         %Plotting Occupied Focus
         plot3(0,0,0,'yo','MarkerFaceColor','y','MarkerSize',10);
         hold on
         %Plot of Orbit 1
         plot3(STATE(:,1),STATE(:,2),STATE(:,3),'w','LineWidth',2);
         hold on
         %Plot Periapsis
         plot3(Rp(1),Rp(2),Rp(3),'wo','MarkerSize',8);
         hold on
         %Plot Apoapsis
         plot3(Ra(1),Ra(2),Ra(3),'wd','MarkerSize',10);
         hold on
         
   end

   %Bounding Plot
   axis equal
   x1 = xlim;
   y1 = ylim;
   z1 = zlim;
   
   RE = 6378.1363; %[km] Earth Mean Equatorial Radius
   if (abs(x1(1)) <= (RE + 5000))
      x1(1) = -RE - 5000;
   end
   if (abs(x1(2)) <= (RE + 5000))
      x1(2) = RE + 5000;
   end
   if (abs(y1(1)) <= (RE + 5000))
      y1(1) = -RE - 5000;
   end
   if (abs(y1(2)) <= (RE + 5000))
      y1(2) = RE + 5000;
   end   
   if (abs(z1(1)) <= (RE + 5000))
      z1(1) = -RE - 5000;
   end
   if (abs(z1(2)) <= (RE + 5000))
      z1(2) = RE + 5000;
   end  
      
   %Plotting Satellite Location (R)
   plot3(R(1),R(2),R(3),'rd','LineWidth',1.5,'MarkerSize',10);
   hold on   

   %Plot of Equitorial plane
   xx = x1(1):abs(x1(1))+abs(x1(2)):x1(2);
   yy = y1(1):abs(y1(1))+abs(y1(2)):y1(2);
   [xx,yy] = meshgrid(xx,yy);
   zz = zeros(size(xx));
   surf(xx,yy,zz,'FaceColor',[0 0.5 0.5],'FaceAlpha',0.15,'EdgeColor',[0.5 0.5 0.5]);
   
   %Plot Formatting
   set(gca,'Color','k','XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1],'FontSize',12);
   xlim(x1);
   ylim(y1);
   zlim(z1);
   grid on
   
   switch centralbody
      
      case 'earth'

         xlabel('X GCRF (km)');
         ylabel('Y GCRF (km)');
         zlabel('Z GCRF (km)');
 
      case 'sun'

         xlabel('X ICRF (km)');
         ylabel('Y ICRF (km)');
         zlabel('Z ICRF (km)');
   
   end

end