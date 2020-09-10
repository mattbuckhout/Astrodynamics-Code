
clearvars; clc;

% ORBIT VISUALIZATION FROM TWO-LINE-ELEMENT SETS
% ---------------------------------------------------------------------------------------------------------------------
% This script visualizes a series of Geocentric satellite orbits in 3D from Two-Line Element Sets (TLEs) at a specified
% Date and Time (UTC). The script can plot up to 10 orbits in a single figure, with each orbit labeled. The satellite
% state at the TLE epoch is computed from the mean elements in the TLE and propagated to the user specified date using
% a 2-body propagation routine. 
%
% Due to the use of mean elements in the TLE and the inaccuracy of 2 body propagation over longer periods of time, this
% script can provide a qualitative representation of a satellite orbit close to the epoch for which the TLE was 
% generated, but should not be used to predict the future state of a satellite with great accuracy.
%
% Author: Matthew Buckhout
% Updated: 08/18/2020 
%
% Inputs:
%
%     - [DATE]       Date                                            [yyyy mm dd]
%     - [UTC]        Time (UTC)                                      [hh mm ss]
%     - [numsat]     Vector of Identifiers for TLEs to plot           - 
%                       ex: [1 3 2 4 8 9] 
%     - [now]        Specified Date/Time Selection                    -
%                       1 - Use current Date/Time in your timezone
%                       0 - Use specified Date/Time in UTC
%     - [dst]        Daylight Savings Time                            -
%                       1 - In effect
%                       0 - Not in effect
%
% Functions:
%
%     - kepler_equation
%     - true_anomaly
%     - DOYmodhms
%     - position_velocity  
%     - propagation_UV 
%     - orbital_elements
%     - julian_date
%     - getdata_EOP_fast
%     - EOP
%     - convert_time
%     - local_sidereal
%
% References:
%     - Fundamentals of Astrodynamics with Applications, 2nd ed. (Vallado) (pg. 113-116)
%     - Celestrak, Two-Line Element Set Format (Kelso) 
%           https://celestrak.com/columns/v04n05/
%     - 3D Earth Example (Ryan Gray - 2004,2006,2013)
%           https://www.mathworks.com/matlabcentral/fileexchange/13823-3d-earth-example
%
% -------------------------------------  U S E R  S P E C I F I E D  I N P U T S  -------------------------------------

%Time Information & Settings
DATE   = [2020 01 29]; %[yyyy mm dd] Specified Date
UTC    = [23 39 35]; %[hh mm ss] Specified Time UTC
now    = 1; %Current time (1) or specified time (0)
dst    = 1; %Daylight Savings Time in effect (1)
numsat = [1 2 3 4 5 6 7 8 9 10]; %Identifiers of satellite orbits to plot

%Satellite 1
name1  = ['HST'];
line11 = ['1 20580U 90037B   20129.00785880  .00000438  00000-0  15673-4 0  9995'];
line12 = ['2 20580  28.4690 227.1744 0002835  20.2350 327.5360 15.09382698449849'];

%Satellite 2
name2  = ['NUSTAR'];
line21 = ['1 38358U 12031A   20128.84725905  .00000924  00000-0  22789-4 0  9994'];
line22 = ['2 38358   6.0254 203.9989 0011562  22.7978 337.2664 14.88295913429956'];

%Satellite 3
name3  = ['FGRST'];
line31 = ['1 33053U 08029A   20129.09992711  .00000470  00000-0  10059-4 0  9995'];
line32 = ['2 33053  25.5838  15.1784 0012014 260.2386  99.6792 15.11334177657203'];

%Satellite 4
name4  = ['IRIS'];
line41 = ['1 39197U 13033A   20129.07574371  .00000070  00000-0  16694-4 0  9994'];
line42 = ['2 39197  97.9581 311.6211 0027732   3.5420 356.5982 14.78721910370132'];

%Satellite 5
name5  = ['SWIFT']; 
line51 = ['1 28485U 04047A   20128.86197786  .00000827  00000-0  26009-4 0  9995'];
line52 = ['2 28485  20.5571 309.0401 0010974 263.8246  96.0938 15.04605572847355'];

%Satellite 6
name6  = ['WISE']; 
line61 = ['1 36119U 09071A   20129.18633385  .00000690  00000-0  26249-4 0  9990'];
line62 = ['2 36119  97.3322 165.9810 0003425 138.0938 222.0565 15.30061148577644'];

%Satellite 7
name7  = ['TESS']; 
line71 = ['1 43435U 18038A   18124.35259418 -.00013193  00000-0  00000+0 0  9999'];
line72 = ['2 43435  28.9141  36.7562 9597568 230.9244 354.5130  0.10934850    10'];

%Satellite 8 
name8  = ['XMM-NEWTON'];
line81 = ['1 25989U 99066A   20129.66380319 -.00000332  00000-0  00000-0 0  9992'];
line82 = ['2 25989  71.1198 322.1430 7182791  93.4109 359.6991  0.50136007 26189'];

%Satellite 9
name9  = ['INTEGRAL'];
line91 = ['1 27540U 02048A   20126.93283063  .00000364  00000-0  00000+0 0  9999'];
line92 = ['2 27540  53.9015 117.4830 8940548 290.1011   2.1921  0.37584322 18232'];

%Satellite 10 
name10  = ['CXO'];
line101 = ['1 25867U 99040B   20129.49755612  .00000633  00000-0  00000+0 0  9994'];
line102 = ['2 25867  63.1186 268.0403 7671163 228.6376   0.1770  0.37801908 20040'];

% -------------------------------------------  I N P U T  H A N D L I N G  --------------------------------------------

%Constants
RE = 6378.1363; %[km] Earth Mean Equatorial Radius 
mu = 3.986004415e5; %[km^3/s^2] Earth Gravitational Parameter

%Current Date and Time & Daylight Savings Time
if (dst == 1)
   addh = 5;
else
   addh = 6;
end
if (now == 1)
   
   c = clock;
   DATE = [c(1) c(2) c(3)];
   UTC = [(c(4)+addh) c(5) c(6)];

end

%Stacking TLEs from all satellites into single matrix
TLE_line1 = [line11; line21; line31; line41; line51; line61; line71; line81; line91; line101];
TLE_line2 = [line12; line22; line32; line42; line52; line62; line72; line82; line92; line102];
names = [name1; name2; name3; name4; name5; name6; name7; name8; name9; name10];

% ----------------------------------------  C U R R E N T  P O S I T I O N S  -----------------------------------------

%Looping through Satellites for State Vectors
for sat=1:1:numel(numsat)
   
   sat_ident = numsat(sat); %Identifier for satellite TLE
   
   line1 = TLE_line1(sat_ident,:); %Satellite TLE line 1
   line2 = TLE_line2(sat_ident,:); %Satellite TLE line 2
   
   %Mean Motion
   n = str2num(line2(53:63)); %[rev/day] Mean Motion
   n = (n*(2*pi)/86400); %[rad/s]

   %Orbital Elements at Epoch
   a       = (mu/(n^2))^(1/3); %[km] Semi-Major Axis
   e       = str2num(line2(27:33))/(1e7); % Eccentricity
   i       = deg2rad(str2num(line2(9:16))); %[rad] Inclination
   omega   = deg2rad(str2num(line2(35:42))); %[rad] Argument of Perigee
   OMEGA   = deg2rad(str2num(line2(18:25))); %[rad] Right Ascension of Ascending Node
   M       = deg2rad(str2num(line2(44:51))); %[rad] Mean Anomaly at Epoch
   [EBH]   = kepler_equation(mu,M,e,0,0); %[rad] Eccentric Anomaly
   [theta] = true_anomaly(e,EBH,0,0); %[rad] True Anomaly

   if ((e < 1) || (e == 0))
      p = a*(1-(e^2)); %[km] Semi-Parameter
   elseif (e > 1)
      p = a*(1-(e^2)); 
   elseif (e == 1)
      p = (h^2)/mu; 
   end

   ap = a*(1 - e) - RE; %[km] Periapsis Altitude
   aa = a*(1 + e) - RE; %[km] Apoapsis Altitude

   %Epoch Imformation
   yr_ep = str2num(line1(19:20)); %[20yy] Last two digits of year 
   DOY = str2num(line1(21:32)); % Day of Year
   [mo_ep,d_ep,h_ep,m_ep,s_ep] = DOY_modhms(yr_ep,DOY); %Month,Day,Hour,Min,Sec
   UTCep  = [h_ep m_ep s_ep];

   %Position/Velocity at Epoch
   t0 = (d_ep*86400) + (h_ep*3600) + (m_ep*60) + s_ep;
   [Rep,Vep] = position_velocity(mu,p,e,i,OMEGA,omega,theta,0,0,0);

   %Postion/Velocity at Selected Time
   t = (DATE(3)*86400) + (UTC(1)*3600) + (UTC(2)*60) + UTC(3);
   [R,V] = propagation_UV(mu,Rep,Vep,t0,t); 

   %Matrix of Satellite State Vectors at Specified Time [Rx Ry Rz Vx Vy Vz]
   STATE(sat,:) = [R' V'];
  
   % ------------------------------------------------  O R B I T S  ---------------------------------------------------
   
   [p,a,e,i,OMEGA,omega,theta,omega_true,lambda_true,u] = orbital_elements(mu,R,V); %Orbital Elements

   if (e < 1) %Circular and Elliptical Orbits
      
      jnum = 100; %Number of state vectors to generate
      for j=1:1:jnum+1 %Loop for 1 period
         
         theta1 = j*(2*pi/jnum) - 2*pi/jnum; %Iterated True Anomaly

         [RR,VV] = position_velocity(mu,p,e,i,OMEGA,omega,theta1,omega_true,lambda_true,u); 
         Rmat(j,1:3) = RR; %[km] Position Vector
         Vmat(j,1:3) = VV; %[km/s] Velocity Vector
         
      end

   elseif (e >= 1)

      fprintf('Orbit is Hyperbolic\n\n');
      
   end
   
   %Position Vectors defining Orbit
   ORBIT(sat,:,:) = Rmat;

end

% --------------------------------------------  3 D  O R B I T  P L O T  ----------------------------------------------

set(groot,'defaultfigureposition',[700 200 900 750])

figure(1,'Color','k');

%Orbit Colors
RGB = [0 0.4470 0.7410; ...
       0.8500 0.3250 0.0980; ...
       0.9290 0.6940 0.1250; ...
       0.4940 0.1840 0.5560; ...
       0.4660 0.6740 0.1880; ...
       0.6350 0.0780 0.1840; ...
       0 0.75 0.75; ...
       0 1 1; ...
       0 1 0; ...
       0 0 1];
       
RGB = [1 1 1; ...
       1 1 1; ...
       1 1 1; ...
       1 1 1; ...
       1 1 1; ...
       1 1 1; ...
       1 1 1; ...
       1 1 1; ...
       1 1 1; ...
       1 1 1];

%Plot Textured Ellipsoid Earth Model [modified from Original Code: Ryan Gray 2004,2013]
yr   = DATE(1);
mo   = DATE(2);
d    = DATE(3);
leap = 0; %Does day contain added leap second

%Earth Orientation Data
[~,MJD] = julian_date(DATE,UTC,0);
EOPdata = getdata_EOP_fast;
[~,~,dUT1,dAT] = EOP(MJD,EOPdata);

%Universal Time
[UT1,~,~,~,~,~,~] = convert_time(DATE,UTC,dUT1,dAT);
%Julian Date
[~,MJD_UT1] = julian_date(DATE,UT1,0);
%Greenwich Mean Standard Time
[~,GMST] = local_sidereal(0,MJD_UT1); %[rad]
 GMST = rad2deg(GMST); %[deg]

%Earth Ellipsoid Model
RE = 6378.1363; %[km]       Earth Mean Equatorial Radius
bE = 6356.7516005; %[km]    Earth Semi-minor Axis (Polar Radius)(Derived)
eE = 0.081819221456; %      Earth Reference Ellipsoid Eccentricity
npanels = 180; %            Number of panels around the equator 
alpha   = 1; %              Transparency
[x_refellipse, y_refellipse, z_refellipse] = ellipsoid(0, 0, 0, RE, RE, bE, npanels);
globe = surf(x_refellipse, y_refellipse, -z_refellipse, ... 
      'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);

rotate(globe,[0 0 1],(GMST+180)); %Rotates Earth to correct GMST (-X axis through Greenwich at start)

%Texture Map
[world,map] = imread('equirectangular_1.jpg');
world = im2double(world);
set(globe, 'FaceColor', 'texturemap', 'CData', world, ... 
         'FaceAlpha', alpha, 'EdgeColor', 'none');

%Plot Formatting
grid on;
set(gca,'Color','k');
set(gca,'GridColor',[1 1 1]);
axis equal;

clear world;
clear map;
hold on
    
%Loop for all other orbits
for sat=1:1:numel(numsat)

   sat_ident = numsat(sat); %Identifier for satellite TLE

   %Plot of Orbit
   plot3(ORBIT(sat,:,1),ORBIT(sat,:,2),ORBIT(sat,:,3),'color',RGB(sat,:),'LineWidth',1.5);
   hold on
   %Plotting Satellite Location
   labelspace = ['  '];
   satname = names(sat_ident,:);
   plot3(STATE(sat,1),STATE(sat,2),STATE(sat,3),'d','MarkerEdgeColor',RGB(sat,:),'MarkerSize',12);
   text(STATE(sat,1),STATE(sat,2),STATE(sat,3),[labelspace satname],'color',[1 1 1],'FontSize',11);
   hold on
  
end

%Bounding Plot
axis equal
x1 = xlim;
y1 = ylim;
z1 = zlim;
   
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

%Plotting Occupied Focus
plot3(0,0,0,'ko','MarkerFaceColor','k','MarkerSize',10);
hold on

%Plot of Equitorial plane
xx = x1(1):abs(x1(1))+abs(x1(2)):x1(2);
yy = y1(1):abs(y1(1))+abs(y1(2)):y1(2);
[xx,yy] = meshgrid(xx,yy);
zz = zeros(size(xx));
surf(xx,yy,zz,'FaceColor',[0 0.5 0.5],'FaceAlpha',0.15,'EdgeColor',[0.5 0.5 0.5]);
hold on

%Plot Formatting
set(gca,'Color','k','XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1],'FontSize',12);
xlim(x1);
ylim(y1);
zlim(z1);
grid on
xlabel('X GCRF (km)');
ylabel('Y GCRF (km)');
zlabel('Z GCRF (km)');
