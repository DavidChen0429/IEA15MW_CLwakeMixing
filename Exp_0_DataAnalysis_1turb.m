%% Data Analysis for Experiments
clear
close all 
addpath('.\Functions');
%clc

%% Data file (Chage this accordingly)
turbineName = '.\Data\NREL5MW\';
caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\';
fileName = '1Turbine_Basic.mat';
basefile = '1Turbine_Basic.mat';

Data = load([turbineName caseName fileName]);
Baseline = load([turbineName caseName basefile]);

% Baseline info
Cp_store_bl = Baseline.Cp_store;
PitchAngles_bl = Baseline.PitchAngles;
%   Blade1
Mflap1_store_bl = Baseline.Mflap1_store;
Medge1_store_bl = Baseline.Medge1_store;
%   Blade2
Mflap2_store_bl = Baseline.Mflap2_store;
Medge2_store_bl = Baseline.Medge2_store;
%   Blade3
Mflap3_store_bl = Baseline.Mflap3_store;
Medge3_store_bl = Baseline.Medge3_store;

% Controlled info 
Cp_store = Data.Cp_store;
PitchAngles = Data.PitchAngles;
%   Blade1
Mflap1_store = Data.Mflap1_store;
Medge1_store = Data.Medge1_store;
%   Blade2
Mflap2_store = Data.Mflap2_store;
Medge2_store = Data.Medge2_store;
%   Blade3
Mflap3_store = Data.Mflap3_store;
Medge3_store = Data.Medge3_store;

%% Basic Settings
D_NREL5MW = 126;
U_inflow = 10;
timeStep = 0.1;
filter = 1000;
simLength = length(Baseline.Cp_store);

%% Calcuate Power, DEL, PBD
% Baseline
PowerTurb_bl = calculatePower(filter,Cp_store_bl,D_NREL5MW,U_inflow); % [MW]
DELTurb_bl = calculateDEL(filter, ...
    Mflap1_store_bl,Medge1_store_bl, ...
    Mflap2_store_bl,Medge2_store_bl, ...
    Mflap3_store_bl,Medge3_store_bl, ...
    timeStep); % [Nm]
PBDTurb_bl = calculatePBD(filter,PitchAngles, ...
    Mflap1_store_bl,Medge1_store_bl, ...
    Mflap2_store_bl,Medge2_store_bl, ...
    Mflap3_store_bl,Medge3_store_bl); % [kNm deg]
fprintf('================================================== \n');
fprintf('Baseline ====================== \n');
fprintf('The output of WT1:\n');
fprintf('    Power Production: %.2f [MW]\n', PowerTurb_bl);
fprintf('    DEL:\n');
fprintf('        Blade1 Flapwise: %.2e  [Nm]\n', DELTurb_bl(1))
fprintf('        Blade1 Edgewise: %.2e  [Nm]\n', DELTurb_bl(2))
fprintf('        Blade2 Flapwise: %.2e  [Nm]\n', DELTurb_bl(3))
fprintf('        Blade2 Edgewise: %.2e  [Nm]\n', DELTurb_bl(4))
fprintf('        Blade3 Flapwise: %.2e  [Nm]\n', DELTurb_bl(5))
fprintf('        Blade3 Edgewise: %.2e  [Nm]\n', DELTurb_bl(6))
fprintf('    PBD:\n');
fprintf('        Blade1: %.2e [kNm deg]\n', PBDTurb_bl(1));
fprintf('        Blade1: %.2e [kNm deg]\n', PBDTurb_bl(2));
fprintf('        Blade1: %.2e [kNm deg]\n', PBDTurb_bl(3));

% Controlled
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
fprintf('Controlled ====================== \n');
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
