%% Dynamic benthic species model
% Author: Muriel Bruckner
% Final version 17/1/2021
% Benthic model computing species distribution of several benthic
% organisms, including competition and grazing. The model changes the
% critical bed shear stress and/or the erosion parameter of the mud
% fraction depending on the biomass in the cells. The new spatial
% distribution of the two parameters is exported and used in the
% Delft3D model

%% Initialisation path
clear 
close all
clc

% User - Define directory
directory_head  = 'C:\Bioturbation\CODE HULL\'; % folder with modules and initial files
ID1             = 'estuary';     % name of the simulation file (mdf-file)
name_model      = 'Delft3D model'; % folder of scenario run
directory       = [directory_head, name_model,'\']; % main directory used in computations
cd([directory, 'initial files\']); % directory initial files
% turn this on in case MATLAB is not linked to this folder
% addpath('C:\win64\delft3d_matlab')  


% add paths with functions and scripts of modules
addpath([directory_head,'Matlab benthos modules']);
addpath([directory_head,'Matlab functions']);
%% User defined parameters for dynamic vegetation model

% User - Define parameters of the hydro-morphodynamic computation
mor         = 1;    % 1= include morphology, 0 = exclude morphology 
morf        = 60;   % give manual morfac in case mor = 0
fl_dr       = 0.05;  % Boundary for flooding/drying threshold used in the computations [m]
Restart     = 0;    % =1 if restart from a trim-file when start time is not 0
phyto       = 1;    % to turn on microphytobenthos computations phyto=1
bioturbation = 1;   % to turn on macrozoobenthic computations bioturbation=1

% User - Define time-scales
t_eco_year  = 12; % number ecological time-steps per year
years       = 50; % number of morphological years of entire simulation

%% run dynamic benthic species model
Benthos_model
