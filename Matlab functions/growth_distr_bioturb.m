function[fct1]=growth_distr_bioturb(matrix, a, b)
% Author: Muriel Br?ckner, version 23/4/2021
% function to calculate the growth distribution of benthos fct1 based on
% environmental parameters stored in matrix
% matrix: matrix with the environmental parameters
% a minimum threshold 
% b maximum threshold

% maximum of parabolic curve
maxi = (b+a)/2;
% n number of maxima
n = 1;

% calculate parabolic function k
k = n/((maxi-a)*(maxi-b));
% calculate fraction distribution (0-1) based on environmental parameter
fct1 = k*(matrix-a).*(matrix-b);

% set thresholds to 0 to exclude supratidal
fct1(matrix<=a) = 0;
fct1(matrix>=b) = 0;


end
