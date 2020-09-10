
% SETUP SCRIPT - BUILDING MATLAB EXECUTABLES FROM C CODE
% -------------------------------------------------------------------------------------------------
% Some of the astrodynamics functions in this code repository contain routines that are very 
% computationally intense, requiring complex tasks to be repeated thousands of times to generate 
% the desired set of outputs. For some of these functions, the computationally intense routines
% proved to be prohibitively slow when implemented in the Matlab language, which is optimized for 
% operating on matrices or performing vectorized operations, but is very slow when dealing with 
% element by element operations in nested loops. To improve performance, some routines were 
% translated from Matlab to C in order to speed up operation. These functions can be compiled 
% as Matlab executables (MEX), which can then be called within Matlab functions or scripts as if 
% they were Matlab files, but run with the performance of compiled C code. 
%
% This script must be run once in order to compile all the necessary functions implemented in C 
% code into Matlab executable files. When the MEX files have been built, they will appear
% in the directory in which the source C files reside.
%
% Author: Matthew Buckhout
% Updated: 09/08/2020 
%
% C Source Files
%
%     - asc_legendre_c.c
%     - lambert_c.c
%     - planet_ephemerides_c.c
%     - planet_ephemerides_full_c.c
%     - propagator_c.c
%     - series_B_c.c
%
% MEX Functions Built:
%
%     - asc_legendre_c.mex
%     - lambert_c.mex
%     - planet_ephemerides_c.mex
%     - planet_ephemerides_full_c.mex
%     - propagator_c.mex
%     - series_B_c.mex 
%
% References:
%     - Mathworks Create a C Source MEX File
%           https://www.mathworks.com/help/matlab/matlab_external/standalone-example.html
%     - GNU Octave Guide to Mex-Files
%           https://octave.org/doc/v4.2.0/Mex_002dFiles.html#Mex_002dFiles
% -------------------------------------------------------------------------------------------------

%Compiling C code as Matlab Executables
mex asc_legendre_c.c; 
mex lambert_c.c;
mex planet_ephemerides_c.c;
mex planet_ephemerides_full_c.c;
mex propagator_c.c;
mex series_B_c.c;
