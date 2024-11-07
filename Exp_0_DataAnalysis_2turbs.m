%% Data Analysis for Experiments
clear
close all 
addpath('.\Functions');
%clc

%% Data file (Chage this accordingly)
turbineName = '.\Data\NREL5MW\';
caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\';
fileName = '2Turbines_Basic.mat';
basefile = '2Turbines_Basic.mat';

Data = load([turbineName caseName fileName]);
Baseline = load([turbineName caseName basefile]);

% Controlled info
%   WT1
Cp_store = Data.Cp_store;
PitchAngles = Data.PitchAngles;
%       Blade1
Mflap1_store = Data.Mflap1_store;
Medge1_store = Data.Medge1_store;
%       Blade2
Mflap2_store = Data.Mflap2_store;
Medge2_store = Data.Medge2_store;
%       Blade3
Mflap3_store = Data.Mflap3_store;
Medge3_store = Data.Medge3_store;
%   WT2
Cpturb2_store = Data.Cpturb2_store;
PitchAnglesturb2 = Data.PitchAnglesturb2;
%       Blade1
Mflap1turb2_store = Data.Mflap1turb2_store;
Medge1turb2_store = Data.Medge1turb2_store;
%       Blade2
Mflap2turb2_store = Data.Mflap2turb2_store;
Medge2turb2_store = Data.Medge2turb2_store;
%       Blade3
Mflap3turb2_store = Data.Mflap3turb2_store;
Medge3turb2_store = Data.Medge3turb2_store;

%% Define Basic Helix Information
D_NREL5MW = 126;
U_inflow = 10;
timeStep = 0.1;
filter = 1000;
simLength = length(Baseline.Cp_store);

%% Calcuate Power, DEL, PBD
% WT1
PowerTurb1 = calculatePower(filter,Cp_store,D_NREL5MW,U_inflow); % [MW]
DELTurb1 = calculateDEL(filter, ...
    Mflap1_store,Medge1_store, ...
    Mflap2_store,Medge2_store, ...
    Mflap3_store,Medge3_store, ...
    timeStep); % [Nm]
PBDTurb1 = calculatePBD(filter,PitchAngles, ...
    Mflap1_store,Medge1_store, ...
    Mflap2_store,Medge2_store, ...
    Mflap3_store,Medge3_store); % [kNm deg]
fprintf('================================================== \n');
fprintf('Baseline ====================== \n');
fprintf('The output of WT1:\n');
fprintf('    Power Production: %.2f [MW]\n', PowerTurb1);
fprintf('    DEL:\n');
fprintf('        Blade1 Flapwise: %.2e  [Nm]\n', DELTurb1(1))
fprintf('        Blade1 Edgewise: %.2e  [Nm]\n', DELTurb1(2))
fprintf('        Blade2 Flapwise: %.2e  [Nm]\n', DELTurb1(3))
fprintf('        Blade2 Edgewise: %.2e  [Nm]\n', DELTurb1(4))
fprintf('        Blade3 Flapwise: %.2e  [Nm]\n', DELTurb1(5))
fprintf('        Blade3 Edgewise: %.2e  [Nm]\n', DELTurb1(6))
fprintf('    PBD:\n');
fprintf('        Blade1: %.2e [kNm deg]\n', PBDTurb1(1));
fprintf('        Blade1: %.2e [kNm deg]\n', PBDTurb1(2));
fprintf('        Blade1: %.2e [kNm deg]\n', PBDTurb1(3));
% WT2
PowerTurb2 = calculatePower(filter,Cpturb2_store,D_NREL5MW,U_inflow); % [MW]
DELTurb2 = calculateDEL(filter, ...
    Mflap1turb2_store,Medge1turb2_store, ...
    Mflap2turb2_store,Medge2turb2_store, ...
    Mflap3turb2_store,Medge3turb2_store, ...
    timeStep); % [Nm]
PBDTurb2 = calculatePBD(filter,PitchAnglesturb2, ...
    Mflap1turb2_store,Medge1turb2_store, ...
    Mflap2turb2_store,Medge2turb2_store, ...
    Mflap3turb2_store,Medge3turb2_store); % [kNm deg]
fprintf('The output of WT2:\n');
fprintf('    Power Production: %.2f [MW]\n', PowerTurb2);
fprintf('    DEL:\n');
fprintf('        Blade1 Flapwise: %.2e  [Nm]\n', DELTurb2(1))
fprintf('        Blade1 Edgewise: %.2e  [Nm]\n', DELTurb2(2))
fprintf('        Blade2 Flapwise: %.2e  [Nm]\n', DELTurb2(3))
fprintf('        Blade2 Edgewise: %.2e  [Nm]\n', DELTurb2(4))
fprintf('        Blade3 Flapwise: %.2e  [Nm]\n', DELTurb2(5))
fprintf('        Blade3 Edgewise: %.2e  [Nm]\n', DELTurb2(6))
fprintf('    PBD:\n');
fprintf('        Blade1: %.2e [kNm deg]\n', PBDTurb2(1));
fprintf('        Blade1: %.2e [kNm deg]\n', PBDTurb2(2));
fprintf('        Blade1: %.2e [kNm deg]\n', PBDTurb2(3));
