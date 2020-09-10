
% EARTH ORIENTATION PARAMETERS (EOP) FOR GIVEN DATE
% -------------------------------------------------------------------------------------------------
% Generates EOP at a date of interest for IAU2006 reduction and time conversions. EOPdata variable 
% generated by getdata_EOP.m or getdata_EOP_fast.m functions.
%
% Earth Orientation Parameter Data from IERS 
%     https://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html
%
% Author: Matthew Buckhout
% Updated: 09/07/2020 
%
% Inputs
%
%     - [MJD]              Modified Julian Date (UT1)                    [days]
%     - [EOPdata]          Earth Orientation Parameter Table              -
%
% Outputs
%
%     - [xp]               Position of CEP from IRP (rad)                ["]
%     - [yp]               Position of CEP from IRP (rad)                ["]
%     - [dUT1]             UT1 - UTC                                     [sec]
%     - [dAT]              TAI - UTC                                     [sec]
%
% EOPdata Format:
% -------------------------------------------------------------------------------------------------
%  MJD(UT1)     yr  mo  d     xp (")      yp (")    UT1-UTC (sec)    UT1-TAI (sec)    TAI-UTC (sec)
% -------------------------------------------------------------------------------------------------

function [xp,yp,dUT1,dAT] = EOP(MJD,EOPdata)
   
   %Nearest Daily EOP value
   MJD = round(MJD);
   
   indata = ismember(MJD,EOPdata(:,1));
   
   if (indata == 1)
   
      %Selecting correct EOP Based on Current Date
      [index,~] = ismember(EOPdata(:,1),MJD,'rows');
      ind = find(index == 1);
      xp = EOPdata(ind,5); %["]
      yp = EOPdata(ind,6); %["]
      dUT1 = EOPdata(ind,7); %[sec] UT1-UTC
      dAT = EOPdata(ind,9); %[sec] TAI-UTC
      
      if (isnan(xp) == 1)
         
         fprintf('\nEOP - Specified date outside range of data\n\n');
         xp = 0;
         yp = 0;
         dUT1 = 0;
         dAT = 0;
         
      end
            
   elseif (indata == 0)
      
      fprintf('\nEOP - Specified date outside range of data\n\n');
      xp = 0;
      yp = 0;
      dUT1 = 0;
      dAT = 0;
      
   end
   
end
