%
% Sample usage of interpolation constant evaluation library.
% Xuefeng LIU
% First version: 2018/09/26
%

% my_intlab_mode_config;

INTERVAL_MODE=1; %The value can be 0 or 1. See readme.md
format long

h=1.0/64;
tri=[0,0;1,0;0,1]; %right triangle with unit legs.
tri=[0,0;1,0;0.5,sqrt(3)/2]; %regular triangle.
tri=[0,0;1,0;cos(0.01*pi/3),sin(0.01*pi/3)]; %narrow triangle.


get_mesh( tri,h,'./mesh/');

norm_idx=1;  %The value can be 0 or 1. See detailed description in each constant_* file.

%% Example 1: the constant for Lagrange interpolation error estimation

%~ c_value = constant_c1_1('./mesh/',norm_idx);

%% Example 2: the constant for Fujino-Morley interpolation error estimation

% c_value = constant_c2_3('./mesh/',norm_idx);

%% Example 3: the constant for Crouzeix-Raviart interpolation in bounding eigenvalues.
c_value = constant_crouzeix_raviart('./mesh/');

%% Example 4: the constant for Fujino-Morley interpolation in bounding eigenvalues.

% c_value = constant_fujino_morley('./mesh/',norm_idx);

display(c_value)
