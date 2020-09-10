      
% PLANETARY EPHEMERIDES - DE430 MODEL (JPL) 
% ---------------------------------------------------------------------------------------------------------------------
% Returns position and velocity vectors for any of 11 bodies in the solar system in the International Celestial 
% Reference Frame (ICRF) relative to Solar System Barycenter. Moon position vector given for Geocentric Celestial 
% Reference Frame (GCRF). Ephemerides based on JPL DE430 model.
%     https://ssd.jpl.nasa.gov/?planet_eph_export
%
% Time system used for DE430 model is Barycentric Dynamical Time (TDB). Time intervals in coefficient files are in 
% Julian Days referenced to Terrestrial Dynamical Time (TDT/TT), so the input Modified Julian Date must also be 
% referenced to TT. This implementation uses the DE430t version of coefficient files with TT - TDB included.
%
% Author: Matthew Buckhout
% Updated: 08/16/2020 
%
% Inputs:
%
%     - [MJD_TT]        Modified Julian Date, Expressed in TT or TDT                                  [days]
%     - [DE430coef]     Matrix of coefficients from JPL DE430 model over period containing MJD_TT      -
%     - [bodies]        Vector of identifiers (1-11) specifying which bodies to generate               -
%                          ephemerides for (input string 'all' for all)
%
% Outputs:
%
%     - [ephemerides]   Position & Velocity vectors of selected bodies at MJD_TT                      [km][km/s]
%                          [Rx1 Ry1 Rz1 Vx1 Vy1 Vz1; ...
%                           Rx1 Ry2 Rz2 Vx2 Vy2 Vz2]                            
%
% Identifiers:
%
%     1. Mercury
%     2. Venus
%     3. Earth-Moon Barycenter
%     4. Mars
%     5. Jupiter
%     6. Saturn
%     7. Uranus
%     8. Neptune
%     9. Pluto
%    10. Moon (GCRF)
%    11. Sun
%    12. -
%    13. -
%    14. -
%    15. TT - TDB
%
% References:
%     - JPL Planetary and Lunar Ephemerides, Information and DE430 coefficient files
%           https://ssd.jpl.nasa.gov/?planet_eph_export
%     - JPL Ephemerides Project (Lukas Bystricky)
%           https://people.sc.fsu.edu/~lb13f/projects/space_environment/planet_positions.php
% ---------------------------------------------------------------------------------------------------------------------

function [ephemerides] = planet_ephemerides(MJD,DE430coef,bodies)

   %Data section indices for DE430t 
   startcoef = [3   171   231   309   342   366   387   405   423   441   753   819   819   939   939];
   numcomp   = [3     3     3     3     3     3     3     3     3     3     3     2     3     3     1];  
   numcoef   = [14   10    13    11     8     7     6     6     6    13    11     0    10     0    11];
   numsets   = [4     2     2     1     1     1     1     1     1     8     2     0     4     0     4];

   %Data file size and section lengths
   [rows,~] = size(DE430coef);
   sectionlength = fix(DE430coef(1,2)/3) + 2;
   numsections = rows/sectionlength;
   
   %Finding location of section for given epoch in coefficients (1142 loops)
   i = 1;
   sect = 1;
   MJDints = zeros(numsections,2);
   while (i <= rows)
      MJDints(sect,1) = DE430coef(i+1,1) - 2400000.5;
      MJDints(sect,2) = DE430coef(i+1,2) - 2400000.5;
      if (MJD > MJDints(sect,1) && (MJD <= MJDints(sect,2)))
         section = sect;
         rowstart = i+1;
         break
      end
      i = i + sectionlength;
      sect = sect+1;    
   end

   %Section of coefficients for given epoch (column vector) (339 loops)
   row2 = 1;
   for row1=rowstart:1:rowstart+sectionlength-2
      inds = (row2*3)-2;
      coefs(inds:inds+2,1) = transpose(DE430coef(row1,:));
      row2 = row2 + 1;
   end
   MJDinterval = [MJDints(section,1) MJDints(section,2)]; %Time interval for selected section
   MJDspan = MJDinterval(2) - MJDinterval(1); %Span of time interval [days]
   
   %Generate all ephemerides or selection
   switch bodies
      case 'all'
         bodies = [1 2 3 4 5 6 7 8 9 10 11];      
      otherwise
         bodies = bodies;
      end   
         
      %Loop through planetary bodies (excluding last 4 sets)
      for b=1:1:numel(bodies)

         %Selecting planetary body
         p = bodies(b);

         %Determining subinterval for time
         MJDsubints = zeros(1,numsets(p));
         for i=1:1:numsets(p)+1
            MJDsubints(i) = MJDinterval(1) + (i-1)*(MJDspan/numsets(p));
            if (i > 1)
               if ((MJD > MJDsubints(i-1)) && (MJD <= MJDsubints(i)))
                  subinterval = [MJDsubints(i-1) MJDsubints(i)];
                  subinterval_ind = i-1;
               end
            end
         end
         
         %Time Scaled to between -1 and 1 for subinterval
         tbar = (MJD - ((subinterval(1) + subinterval(2))/2))*(2/(subinterval(2) - subinterval(1)));
         
         %Start indices for x,y,z coefficients in subinterval
         startx = startcoef(p) + (subinterval_ind - 1)*3*numcoef(p);
         starty = startx + numcoef(p);
         startz = starty + numcoef(p);
         
         %Coefficients for x,y,z in subinterval
         xcoef = coefs(startx:starty-1);
         ycoef = coefs(starty:startz-1);
         zcoef = coefs(startz:startz+numcoef(p)-1);

         %Chebyshev Polynomials of the First Kind
         [Tn,dTn] = chebyshev_first(numcoef(p),tbar);

         %Summation for position vector components
         R_ICRF = [0 0 0];
         V_ICRF = [0 0 0];
         
         for i=1:1:numcoef(p)
            
            R_ICRF(1) = R_ICRF(1) + xcoef(i)*Tn(i); %Position Vector of body in ICRF Barycentered Frame
            R_ICRF(2) = R_ICRF(2) + ycoef(i)*Tn(i);
            R_ICRF(3) = R_ICRF(3) + zcoef(i)*Tn(i);

            V_ICRF(1) = V_ICRF(1) + xcoef(i)*dTn(i); %Velocity Vector of body in ICRF Barycentered Frame
            V_ICRF(2) = V_ICRF(2) + ycoef(i)*dTn(i);
            V_ICRF(3) = V_ICRF(3) + zcoef(i)*dTn(i);
         
         end

         V_ICRF = (2*numsets(p)*V_ICRF)/(MJDspan*86400); %Scaling Velocity by number of subintervals, number of sets, days to seconds
         
         ephemerides(b,:) = [R_ICRF V_ICRF];
         
      end

end
