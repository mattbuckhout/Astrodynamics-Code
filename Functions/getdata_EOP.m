
% READ EARTH ORIENTATION PARAMETERS (EOP) FROM INTERNATIONAL EARTH ROTATION SERVICE (IERS) DATA FILES
% ---------------------------------------------------------------------------------------------------------------------
% Function reads Earth Orientation Parameters from 3 data files and generates a combined table which can be
% interpolated by the EOP.m function to retrieve EOP at date of interest. 
%
% Author: Matthew Buckhout
% Updated: 08/24/2020 
%
% Earth Orientation Parameter Data from IERS 
%     https://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html
%
% IERS filename:  "finals.all(IAU2000)"
% Local filename: "IERS_EOP1.txt"
%
%     Data Format:
%     ---------------------------------------------------------------------------
%        yr mo d      MJD      xp ["]   yp ["]  erxp ["] eryp ["]   UT1-UTC [sec]   
%     ---------------------------------------------------------------------------
%        73  1 2 41684.00 I  0.120733 0.009786  0.136966 0.015902  I 0.8084178 
%        73  1 3 41685.00 I  0.118980 0.011039  0.135656 0.013616  I 0.8056163 
%        73  1 4 41686.00 I  0.117227 0.011039  0.134348 0.013616  I 0.8027895       ...
%        73  1 5 41687.00 I  0.115473 0.009743  0.133044 0.013089  I 0.7998729 
%        73  1 6 41688.00 I  0.113717 0.011236  0.131746 0.009898  I 0.7968144
%     --------------------------------------------------------------------------- 
%
% IERS filename:  "EOP C01 IAU2000(1900-now)"
% Local filename: "IERS_EOP2.txt"
%
%     Data Format:
%     ----------------------------------------------------------------------
%      #an         x(")        x_er(")     y(")        y_er(")    UT1-TAI(s) 
%     ----------------------------------------------------------------------
%       1900.00   -0.008134    0.030848    0.074551    0.023046  -0.0001443  
%       1900.05    0.036878    0.028098    0.022985    0.020404  -0.0000593  
%       1900.10    0.048627    0.030600    0.012210    0.023918   0.0002632      ...
%       1900.15    0.030027    0.021668    0.001851    0.021684  -0.0003384  
%       1900.20    0.028247    0.022871   -0.001633    0.021847   0.0000194   
%     -----------------------------------------------------------------------
%
%     Output Format [EOPdata]:
%     -------------------------------------------------------------------------------------------------
%      MJD(UT1)     yr  mo  d     xp (")      yp (")    UT1-UTC (sec)    UT1-TAI (sec)    TAI-UTC (sec)
%     -------------------------------------------------------------------------------------------------

function [EOPdata] = getdata_EOP()

   % Earth Orientation Parameters File 2 
   fid2 = fopen('IERS_EOP2.txt');
   data2 = zeros(2409,2);
      
   %Loop through number of lines
   dataline = fgetl(fid2); %Reads First line of file
   dataline = fgetl(fid2); %Reads next line of file
   linenum = 1;
   while ischar(dataline)

      data2(linenum,2) = str2double(dataline(60:70)); %%UT1 - TAI [sec]
      yrfrac = str2double(dataline(4:10)); %Year fraction
      [~,MJDtemp2] = julian_date([yrfrac 0 0],[0 0 0],0); %MJD
      data2(linenum,1) = MJDtemp2;
      
      linenum = linenum + 1;
      
      dataline = fgetl(fid2); %Reads next line of file
      
   end
   fclose(fid2);

   % Earth Orientation Parameters File 1 
   fid1 = fopen('IERS_EOP1.txt','r');
   EOPdata = zeros(10813,9);

   %Loop through number of lines
   dataline = fgetl(fid1); %Reads First line of file
   linenum = 1;
   while ischar(dataline)

      MJDtemp = str2double(dataline(7:15)); %Modified Julian Date
   
      EOPdata(linenum,1) = str2double(dataline(7:15)); %MJD
      EOPdata(linenum,2) = str2double(dataline(1:2)); %Year
      EOPdata(linenum,3) = str2double(dataline(3:4)); %Month
      EOPdata(linenum,4) = str2double(dataline(5:6)); %Day
      EOPdata(linenum,5) = str2double(dataline(20:27)); %xp ["]
      EOPdata(linenum,6) = str2double(dataline(39:46)); %yp ["]
      EOPdata(linenum,7) = str2double(dataline(59:68)); %UT1-UTC [sec]
      EOPdata(linenum,8) = interp1(data2(:,1),data2(:,2),MJDtemp); %UT1-TAI [sec] 
      EOPdata(linenum,9) = EOPdata(linenum,7) - EOPdata(linenum,8); %TAI-UTC [sec] 

      if (isnan(EOPdata(linenum,9)) == 1)
         EOPdata(linenum,9) = 37; %TAI-UTC [sec] Value held constant until next leap second added 
         EOPdata(linenum,8) = EOPdata(linenum,7) - EOPdata(linenum,9); %UT1-TAI [sec] 
      end
            
      if (MJDtemp > 51543) %Years past 2000
         
         EOPdata(linenum,2) = EOPdata(linenum,2) + 2000;
         
      elseif (MJDtemp <= 51543) %Years up to 1999 
      
         EOPdata(linenum,2) = EOPdata(linenum,2) + 1900;
      
      end
      
      linenum = linenum + 1;
      
      dataline = fgetl(fid1); %Reads next line of file
      
   end
   fclose(fid1);
   
end

