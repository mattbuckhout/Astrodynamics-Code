
% READ IAU2006 PRECESSION-NUTATION MODEL COEFFICIENTS
% ---------------------------------------------------------------------------------------------------------------------
% Function reads and stores data files containing coefficients used to buidl X,Y,s series terms in the IAU2006 
% reduction, allowing for transformation between the inertial GCRF and body-fixed ITRF frames.
% 
% Coefficients of X,Y,s Series terms for IAU 2006 Reduction from IERS Conventions (2010)
%     http://webtai.bipm.org/iers/convupdt/convupdt_c5.html#Listupdt
%
% length = 'short' for X,Y,s files edited down to rows for largest 10 terms in each series for faster computation.
% length = 'long' for full X,Y,s files, slower computation but higher precision.
% Data file and stored variable format identical.
%
% Author: Matthew Buckhout
% Updated: 08/25/2020 
%
% Local filenames: "IAU2006_X_coef_short.txt"
%                  "IAU2006_Y_coef_short.txt"
%                  "IAU2006_s_coef_short.txt"
%
%                  "IAU2006_X_coef.txt"
%                  "IAU2006_Y_coef.txt"
%                  "IAU2006_s_coef.txt"
%
%     Data Format:
%     ----------------------------------------------------------------------------------------------------------------
%         i            Axs            Axc   apx1 apx2 apx3 apx4 apx5 apx6 apx7 apx8 apx9 apx10 apx11 apx12 apx13 apx14
%     ----------------------------------------------------------------------------------------------------------------
%         1    -6844318.44        1328.67      0    0    0    0    1    0    0    0    0     0     0     0     0     0
%         2     -523908.04        -544.75      0    0    2   -2    2    0    0    0    0     0     0     0     0     0
%         3      -90552.22         111.23      0    0    2    0    2    0    0    0    0     0     0     0     0     0
%         4       82168.76         -27.64      0    0    0    0    2    0    0    0    0     0     0     0     0     0
%         5       58707.02         470.05      0    1    0    0    0    0    0    0    0     0     0     0     0     0
%         6       28288.28         -34.69      1    0    0    0    0    0    0    0    0     0     0     0     0     0
%         7      -20557.78         -20.84      0    1    2   -2    2    0    0    0    0     0     0     0     0     0
%         8      -15406.85          15.12      0    0    2    0    1    0    0    0    0     0     0     0     0     0
%         9      -11991.74          32.46      1    0    2    0    2    0    0    0    0     0     0     0     0     0
%        10       -8584.95           4.42      0    1   -2    2   -2    0    0    0    0     0     0     0     0     0
%     ----------------------------------------------------------------------------------------------------------------
%
%     ----------------------------------------------------------------------------------------------------------------
%         i            Ays            Ayc   apy1 apy2 apy3 apy4 apy5 apy6 apy7 apy8 apy9 apy10 apy11 apy12 apy13 apy14
%     ----------------------------------------------------------------------------------------------------------------
%         1        1538.18     9205236.26      0    0    0    0    1    0    0    0    0     0     0     0     0     0
%         2        -458.66      573033.42      0    0    2   -2    2    0    0    0    0     0     0     0     0     0
%         3         137.41       97846.69      0    0    2    0    2    0    0    0    0     0     0     0     0     0
%         4         -29.05      -89618.24      0    0    0    0    2    0    0    0    0     0     0     0     0     0
%         5         -17.40       22438.42      0    1    2   -2    2    0    0    0    0     0     0     0     0     0
%         6          31.80       20069.50      0    0    2    0    1    0    0    0    0     0     0     0     0     0
%         7          36.70       12902.66      1    0    2    0    2    0    0    0    0     0     0     0     0     0
%         8         -13.20       -9592.72      0    1   -2    2   -2    0    0    0    0     0     0     0     0     0
%         9        -192.40        7387.02      0    1    0    0    0    0    0    0    0     0     0     0     0     0
%        10           3.92       -6918.22      0    0    2   -2    1    0    0    0    0     0     0     0     0     0
%     ----------------------------------------------------------------------------------------------------------------
%     
%     ----------------------------------------------------------------------------------------------------------------
%         i            As             Asc   aps1 aps2 aps3 aps4 aps5 aps6 aps7 aps8 aps9 aps10 aps11 aps12 aps13 aps14
%     ----------------------------------------------------------------------------------------------------------------
%              1       -2640.73           0.39    0    0    0    0    1    0    0    0    0    0    0    0     0     0
%              2         -63.53           0.02    0    0    0    0    2    0    0    0    0    0    0    0     0     0
%              3         -11.75          -0.01    0    0    2   -2    3    0    0    0    0    0    0    0     0     0
%              4         -11.21          -0.01    0    0    2   -2    1    0    0    0    0    0    0    0     0     0
%              5           4.57           0.00    0    0    2   -2    2    0    0    0    0    0    0    0     0     0
%              6          -2.02           0.00    0    0    2    0    3    0    0    0    0    0    0    0     0     0
%              7          -1.98           0.00    0    0    2    0    1    0    0    0    0    0    0    0     0     0
%              8           1.72           0.00    0    0    0    0    3    0    0    0    0    0    0    0     0     0
%              9           1.41           0.01    0    1    0    0    1    0    0    0    0    0    0    0     0     0
%             10           1.26           0.01    0    1    0    0   -1    0    0    0    0    0    0    0     0     0
%     ----------------------------------------------------------------------------------------------------------------

function [XTAB,YTAB,sTAB] = getdata_IAU2006(length)
   
   switch length
      
      case 'short' %Truncated Data Sets
      
         % Coefficients of Series Terms for IAU2006 X,Y,s from IERS Conventions (2010)
         %IAU 2006 X Coefficients - Short
         fid = fopen('IAU2006_X_coef_short.txt','r');
         [XTAB, count] = fscanf(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [17, 35]);
         fclose(fid);
         XTAB = transpose(XTAB);

         %IAU 2006 Y Coefficients - Short
         fid = fopen('IAU2006_Y_coef_short.txt','r');
         [YTAB, count] = fscanf(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [17, 36]);
         fclose(fid);
         YTAB = transpose(YTAB);

         %IAU 2006 s Coefficients - Short
         fid = fopen('IAU2006_s_coef_short.txt','r');
         [sTAB, count] = fscanf(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [17, 28]);
         fclose(fid);
         sTAB = transpose(sTAB);
      
      case 'long' %Full Data Sets
         
         %IAU 2006 X Coefficients
         fid = fopen('IAU2006_X_coef.txt','r');
         [XTAB, count] = fscanf(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [17, 1600]);
         fclose(fid);
         XTAB = transpose(XTAB);

         %IAU 2006 Y Coefficients
         fid = fopen('IAU2006_Y_coef.txt','r');
         [YTAB, count] = fscanf(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [17, 1275]);
         fclose(fid);
         YTAB = transpose(YTAB);

         %IAU 2006 s Coefficients
         fid = fopen('IAU2006_s_coef.txt','r');
         [sTAB, count] = fscanf(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [17, 66]);
         fclose(fid);
         sTAB = transpose(sTAB);
         
   end

end