
% ORBIT VISUALIZATION - 3D ORBIT AND GROUND TRACK (GEOCENTRIC EQUATORIAL)
% ---------------------------------------------------------------------------------------------------------------------
% Function takes an initial satellite state vector (R,V) in a Geocentric-Equatorial frame (GCRF) and the corresponding 
% Date and Time (UTC) and generates a 3D plot and groundtrack of the orbit for a specified number of revolutions.
% Groundtrack plot includes 1 revolution before input state, as well as 1-6 revolutions following the input state.
% Determined orbit is defined by 2-body motion. Orbits must be Geocentric and either circular or elliptical.
%
%
% Author: Matthew Buckhout
% Updated: 08/24/2020 
%
% Inputs:
%
% Inputs:
%
%     - [R]                Position Vector at Date/Time (GCRF)             [km]
%     - [V]                Velocity Vector at Date/Time (GCRF)             [km/s]
%     - [DATE]             Gregorian Date                                  [yyyy mm dd]
%     - [TIME]             Coordinated Universal Time (UTC)                [hh mm ss]
%     - [numorb]           Number of revolutions to plot on groundtrack     -
%                             (1-6)
%
% Functions:
%
%     - julian_date
%     - getdata_EOP_fast  
%     - EOP
%     - orbital_elements
%     - position_velocity
%     - gregorian_date
%     - groundtrack
%     - convert_time
%     - local_sidereal 
%
% References:
%     -  3D Earth Example [Ryan Gray - 2004,2006,2013]
%        https://www.mathworks.com/matlabcentral/fileexchange/13823-3d-earth-example
% ---------------------------------------------------------------------------------------------------------------------

function orbit_vis_2(R,V,DATE,TIME,numorb);

   RE = 6378.1363; %[km] Earth Mean Equatorial Radius 
   mu = 3.986004415e5; %[km^3/s^2] Earth Gravitational Parameter

   [~,MJD_UTC] = julian_date(DATE,TIME,0); %Modified Julian Date at input date/time (UTC)
   [EOPdata] = getdata_EOP_fast; %Earth Orientation Parameter Data
   [~,~,dUT1,dAT] = EOP(MJD_UTC,EOPdata); %Earth Orientation Paramters
   
   % ---------------------------------------  G E N E R A T E  3 D  D A T A  ------------------------------------------

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
   
   % -------------------------------  G E N E R A T E  G R O U N D  T R A C K  D A T A  -------------------------------

   [p,a,e,~,~,~,~,~,~,~] = orbital_elements(mu,R,V); %Orbital Elements at Epoch
  
   T = 2*pi*sqrt((a^3)/mu); %[sec] %Orbital Period
   
   MJDprev = MJD_UTC - T/86400; %Modified Julian Date, start of previous revolution (UTC)
   
   [DATEprev(1),DATEprev(2),DATEprev(3), ...
    TIMEprev(1),TIMEprev(2),TIMEprev(3)] = gregorian_date(MJDprev); %Date/Time, start of previous revolution (UTC) 
   
   [LATITUDE,LONGITUDE] = groundtrack(R,V,DATE,TIME,dUT1,dAT,numorb); %Ground Track Data

   [LATprev,LONGprev] = groundtrack(R,V,DATEprev,TIMEprev,dUT1,dAT,1); %Previous Period
      
   %Sizing and Positioning Windows
   set(groot,'defaultfigureposition',[50, 190, 800, 420])
   figure 1
   set(gcf, 'Position',[870, 190, 400, 350]);
   figure 2
   set(gcf, 'Position',[50, 190, 800, 420]); 
   
   % -------------------------------------------  3 D  O R B I T  P L O T  --------------------------------------------
   figure(1,'Color','k');
   
   %Plot Textured Ellipsoid Earth Model [modified from Original Code: Ryan Gray 2004,2006,2013]
   %https://www.mathworks.com/matlabcentral/fileexchange/13823-3d-earth-example
   
   %Time and Earth Rotation
   [UT1,~,~,~,~,~,~] = convert_time(DATE,TIME,dUT1,dAT); %Universal Time
   [~,MJD] = julian_date(DATE,UT1,0); %Modified Julian Date at input date/time (UT1)
   [~,GMST] = local_sidereal(0,MJD); %[rad] %Greenwich Mean Standard Time at input date/time

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
   plot3(STATE(:,1),STATE(:,2),STATE(:,3),'color',[0 0.4470 0.7410],'LineWidth',2);
   hold on
   %Plot Periapsis
   plot3(Rp(1),Rp(2),Rp(3),'o','color',[0 0.4470 0.7410],'MarkerSize',10);
   hold on
   %Plot Apoapsis
   plot3(Ra(1),Ra(2),Ra(3),'d','color',[0 0.4470 0.7410],'MarkerSize',10);
   hold on
   
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
   
   %Axis Labels
   xlabel('X GCRF (km)');
   ylabel('Y GCRF (km)');
   zlabel('Z GCRF (km)');
   
   % --------------------------------------  G R O U N D  T R A C K  P L O T  -----------------------------------------
   figure(2,'Color','k');
   
   %Earth Map
   map = imread('equirectangular_2.png');
   [y,x,z] = size(map);
   image('XData', [-180 180], 'YData', [-90 90],'Cdata',flipud(map));
   hold on

   %Default Color Set
   RGB = [0 0.4470 0.7410; ...
          0.8500 0.3250 0.0980; ...
          0.9290 0.6940 0.1250; ...
          0.4940 0.1840 0.5560; ...
          0.4660 0.6740 0.1880; ...
          0.6350 0.0780 0.1840];
   
   %Plotting Orbit 0
   for m=1:1:numel(LONGprev(1,:))

      LONGq = LONGprev(1,m); %[deg]
      LATq = LATprev(1,m); %[deg]  
      
      %Longitude Handling
      if (LONGq > 180)
         LONGprev(1,m) = LONGq - 360;
         LATprev(1,m) = LATq;   
      elseif (LONGq < -180)
         LONGprev(1,m) = LONGq + 360;
         LATprev(1,m) = LATq;
      else
         LONGprev(1,m) = LONGq;
         LATprev(1,m) = LATq;
      end

   end
   
   %Editing Longitude and Latitude 
   q = 1;
   start(1) = 1;
   for m=1:1:numel(LONGprev)
   
      if (m < numel(LONGprev))
         signdiff = (LONGprev(m)/abs(LONGprev(m))) + (LONGprev(m+1)/abs(LONGprev(m+1)));
      end
      if (abs(LONGprev(m)) < 90)
         signdiff = 1;
      end
      
      if (signdiff == 0)
   
         start(q+1) = m+1;
         stop(q) = m;
         q = q+1;
   
      end
   end

   stop(q) = m;
   
   for q=1:1:numel(start)
   
      plot(LONGprev(start(q):stop(q)),LATprev(start(q):stop(q)),'color', [0 0.7 0.7],'LineWidth',1.5);
      hold on
   
   end

   %Plotting Orbits 1-numorb
   for k=1:1:numorb
 
      for m=1:1:numel(LONGITUDE(:,1))

         LONGq = LONGITUDE(m,k); %[deg]
         LATq = LATITUDE(m,k); %[deg]  
         
         %Longitude Handling
         if (LONGq > 180)
            LONGITUDE(m,k) = LONGq - 360;
            LATITUDE(m,k) = LATq;   
         elseif (LONGq < -180)
            LONGITUDE(m,k) = LONGq + 360;
            LATITUDE(m,k) = LATq;
         else
            LONGITUDE(m,k) = LONGq;
            LATITUDE(m,k) = LATq;
         end

      end
      
      %Editing Longitude and Latitude 
      q = 1;
      start(1) = 1;
      for m=1:1:numel(LONGITUDE(:,1))
      
         if (m < numel(LONGITUDE(:,1)))
            signdiff = (LONGITUDE(m,k)/abs(LONGITUDE(m,k))) + (LONGITUDE(m+1,k)/abs(LONGITUDE(m+1,k)));
         end
         if (abs(LONGITUDE(m,k)) < 90)
            signdiff = 1;
         end
         
         if (signdiff == 0)
      
            start(q+1) = m+1;
            stop(q) = m;
            q = q+1;
      
         end
      end

      stop(q) = m;
      
      for q=1:1:numel(start)
      
         plot(LONGITUDE(start(q):stop(q),k),LATITUDE(start(q):stop(q),k),'color',RGB(k,:),'LineWidth',1.5);
         hold on
      
      end
   end
   
   %Boxing Plot
   plot([-180 -180],[-90 90],'k');
   hold on
   plot([180 180],[-90 90],'k');
   hold on
   plot([-180 180],[-90 -90],'k');
   hold on
   plot([-180 180],[90 90],'k');
   hold on
  
   %Plotting Satellite Location at Selected Time
   plot(LONGITUDE(1,1),LATITUDE(1,1),'rd','LineWidth',1.5,'MarkerSize',10);
   
   %Axes
   axis equal
   box on;
   grid on;
   xlim([-180 180]);
   ylim([-90 90]);
   xticks([-180:15:180]);
   yticks([-90:15:90]);
   set(gca, 'Layer','top')
   set(gca,'Color','k','XColor',[1 1 1],'YColor',[1 1 1],'FontSize',12);
   
   %Legend
   L(1) = plot(nan,nan,'color',[0 0.7 0.7],'LineWidth',1.5);
   ar{1} = sprintf('Orbit %02d',0);
   for k=1:1:numorb
      L(k+1) = plot(nan,nan,'color',RGB(k,:),'LineWidth',1.5);
      ar{k+1} = sprintf('Orbit %02d',k);
   end
   legend(L,ar);

   xlabel('Longitude (deg)');
   ylabel('Latitude (deg)');
      
   %Sizing and Positioning
   figure 1
   set(gcf, 'Position',[1220, 210, 650, 700]);
   figure 2
   set(gcf, 'Position',[30, 210, 1180, 700]); 

end
