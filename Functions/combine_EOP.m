
% GENERATE SINGLE EOP DATA FILE FOR FAST EOP READ
% ---------------------------------------------------------------------------------------------------------------------
% Function can be called once to combine the following data files from IERS into single file:
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
% This function calls the slower getdata_EOP.m function to generate the EOPdata table from 3 data
% files, then writes the table to a single text file "IERS_EOP_combined.txt" that can be called 
% by getdata_EOP_fast.m 
%
% Local filename: "IERS_EOP_combined.txt"
%
%     Data Format:
%     ----------------------------------------------------------------------------------------------------------
%      MJD(UT1)     yr           mo       d     xp (")   yp (")   UT1-UTC (sec)  UT1-TAI (sec)  TAI-UTC (sec)
%     ----------------------------------------------------------------------------------------------------------
%      48622.000000 1992.000000 1.000000 1.000000 0.182987 0.168775 -0.125166      -26.199424     26.074258
%      48623.000000 1992.000000 1.000000 2.000000 0.180614 0.166776 -0.126955      -26.201957     26.075002
%      48624.000000 1992.000000 1.000000 3.000000 0.178183 0.164804 -0.128700      -26.204490     26.075790
%      48625.000000 1992.000000 1.000000 4.000000 0.175821 0.162834 -0.130467      -26.207022     26.076556
%      48626.000000 1992.000000 1.000000 5.000000 0.173528 0.160861 -0.132317      -26.209555     26.077239
%      48627.000000 1992.000000 1.000000 6.000000 0.171271 0.158881 -0.134296      -26.212088     26.077792
%     ----------------------------------------------------------------------------------------------------------
%
% Earth Orientation Parameter Data from IERS 
%     https://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html
%
% Author: Matthew Buckhout
% Updated: 08/24/2020 
% ---------------------------------------------------------------------------------------------------------------------

function [] = combine_EOP()

   EOPdata = getdata_EOP; %Generates EOP table from multiple data files 

   fid = fopen('IERS_EOP_combined.txt','w');
   fprintf(fid, '%f %f %f %f %f %f %f %f %f\n',numel(EOPdata(:,1)),0,0,0,0,0,0,0,0); %Storing # of rows in first row 
   for i=1:1:numel(EOPdata(:,1))

      fprintf(fid, '%f %f %f %f %f %f %f %f %f\n',EOPdata(i,1),EOPdata(i,2),EOPdata(i,3),EOPdata(i,4), ...
              EOPdata(i,5),EOPdata(i,6),EOPdata(i,7),EOPdata(i,8),EOPdata(i,9));

   end
   fclose(fid);

end