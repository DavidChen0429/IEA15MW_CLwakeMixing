%% Data Analysis for Experiments
clear
close all 
addpath('.\Functions');
%clc

%% Data file (Chage this accordingly)
turbineName = '.\Data\NREL5MW\';
caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\';
fileName = '1Turbines_OL_Helix.mat';

Data = load([turbineName caseName fileName]);

% WT1 info 
Cp_store = Data.Cp_store;
Moop1_store = Data.Moop1_store;
PitchAngles = Data.PitchAngles;
Mflap1_store = Data.Mflap1_store;
Medge1_store = Data.Medge1_store;

%% Define Basic Helix Information
D_NREL5MW = 126;
U_inflow = 10;
timeStep = 0.1;

%% Calcuate Power, DEL, PBD
filter = 3000;
% WT1 
PowerTurb1 = calculatePower(filter,Cp_store,D_NREL5MW,U_inflow); % [MW]
DELTurb1 = calculateDEL(filter,Moop1_store, timeStep); % [Nm]
PBDTurb1 = calculatePBD(filter,PitchAngles,Mflap1_store,Medge1_store); % [kNm deg]
fprintf('================================================== \n');
fprintf('The output of Wind turbine 1:\n');
fprintf('    Power Production: %.2f [MW]\n', PowerTurb1);
fprintf('    DEL: %.2e  [Nm]\n', DELTurb1);
fprintf('    PBD: %.2e [kNm deg]\n', PBDTurb1);