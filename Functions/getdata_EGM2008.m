
% READ AND UN-NORMALIZE EGM2008 ASPHERICAL GRAVITATIONAL POTENTIAL COEFFICIENTS
% ---------------------------------------------------------------------------------------------------------------------
% Retrieves normalized spherical harmonic coefficients for the EGM2008 gravitational model from data file and stores 
% the coefficients as a matrix for use in computing perturbing accelerations. Function un-normalizes and stores 
% C and S coefficients for each degree (l) and order (m) combination up to the maximum order specified (L).
%
% Author: Matthew Buckhout
% Updated: 08/24/2020 
%
% EGM2008 Spherical Harmonics (2190 x 2190) (Normalized C and S Coefficients)
%     https://earth-info.nga.mil/GandG/update/index.php?action=home#tab_wgs84-data
%
% NGA filename: "EGM2008 Spherical Harmonics (104MB)" (zip)
% Local filename: "EGM2008_Norm_Coefficients_1000x1000.txt"
% 
%     Data Format:
%     ---------------------------------------------------------------------------------------------------
%         l  m                       C                      S               error C               error S
%     ---------------------------------------------------------------------------------------------------
%         2  0  -4.841651437908150e-04  0.000000000000000e+00 7.481239490000000e-12 0.000000000000000e+00
%         2  1  -2.066155090741760e-10  1.384413891379791e-09 7.063781502000000e-12 7.348347201000000e-12
%         2  2   2.439383573283130e-06 -1.400273703859340e-06 7.230231722000001e-12 7.425816951000001e-12
%         3  0   9.571612070934732e-07  0.000000000000000e+00 5.731430751000000e-12 0.000000000000000e+00
%         3  1   2.030462010478640e-06  2.482004158568721e-07 5.726633183000001e-12 5.976692146000002e-12
%         3  2   9.047878948095281e-07 -6.190054751776180e-07 6.374776928000001e-12 6.401837794000001e-12
%         3  3   7.213217571215681e-07  1.414349261929410e-06 6.029131793000000e-12 6.028311182000001e-12
%         4  0   5.399658666389911e-07  0.000000000000000e+00 4.431111968000000e-12 0.000000000000000e+00
%         4  1  -5.361573893888670e-07 -4.735673465180860e-07 4.568074333000001e-12 4.684043490000000e-12
%         4  2   3.505016239626491e-07  6.624800262758292e-07 5.307840320000001e-12 5.186098530000002e-12
%     ---------------------------------------------------------------------------------------------------
%
%     Output Format [EGM2008coef]:
%     --------------------------------------------------------
%         l  m                       C                      S 
%     --------------------------------------------------------
%         2	 0	 -0.001082626173852223	 0.0000000000000000+00
%         2	 1	 -2.667394752374836e-10	 1.787270648524045e-09
%         2	 2	  1.574615325722917e-06	-9.038727891965667e-07
%         3	 0	  2.532410518567723e-06	 0.0000000000000000+00
%         3	 1	  2.193149631313328e-06	 2.680870894008978e-07
%         3	 2	  3.090439003916488e-07	-2.114306209334826e-07
%         3	 3	  1.005835134088229e-07	 1.972215818357184e-07
%         4	 0	  1.619897599916974e-06	 0.0000000000000000+00
%         4	 1	 -5.086435604395839e-07	-4.492654321438083e-07
%         4	 2	  7.837454574045525e-08	 1.481350372488601e-07 
%     --------------------------------------------------------

function [EGM2008coef] = getdata_EGM2008(L)

   % Rows to read based on maximum Degree and Order (L)
   totalrows = 0;
   for l=2:1:L
      totalrows = totalrows + (l + 1);
   end

   % Normalized Spherical Harmonic Coefficients & Uncertainties (L x L)
   fid = fopen('EGM2008_Norm_Coefficients_1000x1000.txt','r');
   CSnorm = transpose(fscanf(fid,'%d %d %f %f %f %f',[6 totalrows]));
   fclose(fid);

   % UnNormalizing EGM2008 Spherical Harmonic Coefficients (L x L)
   EGM2008coef = zeros(totalrows,4);
   EGM2008coef(:,1) = CSnorm(:,1); %Degree (l)
   EGM2008coef(:,2) = CSnorm(:,2); %Order (m)
   tabrow = 1;
   l = 2;
   m = 0; 
   for tabrow=1:1:totalrows
      
      l = EGM2008coef(tabrow,1); %Degree (l)
      m = EGM2008coef(tabrow,2); %Order (m)
      
      if (m == 0)
         k = 1;
      elseif (m ~= 0)
         k = 2;
      end

      PI = sqrt(factorial(l+m)/(factorial(l-m)*k*(2*l + 1))); %Normalization Factor (Vallado 519)
      EGM2008coef(tabrow,3) = CSnorm(tabrow,3)/PI; %Unnormalized C
      EGM2008coef(tabrow,4) = CSnorm(tabrow,4)/PI; %Unnormalized S

   end

end