%% Data Analysis for Experiments
clear
% close all 
addpath('.\Functions');
clc

%% Data file (Chage this accordingly)
turbineName = '.\Data\NREL5MW\';
caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\2Turbines\';
fileName = '2Turbines_CL_Helix_Shear0.2Cheat_mag3.mat';

Data = load([turbineName caseName fileName]);

%% Define Basic Helix Information
D_NREL5MW = 126;
U_inflow = 10;
timeStep = 0.1;
filter = 3500;  % Steady-state value
simLength = length(Data.Cp_store);

%% Get Data
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

%% Calcuate Power, DEL, PBD
% ============== Case
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
fprintf('\n================================================== \n');
fprintf('========== Controlled \n');
fprintf('The output of Upstream WT:\n');
fprintf('    Power Production: %.2f [MW]\n', PowerTurb1);
fprintf('    DEL:\n');
% fprintf('        Blade1 Flapwise: %.2e  [Nm]\n', DELTurb1(1))
% fprintf('        Blade1 Edgewise: %.2e  [Nm]\n', DELTurb1(2))
% fprintf('        Blade2 Flapwise: %.2e  [Nm]\n', DELTurb1(3))
% fprintf('        Blade2 Edgewise: %.2e  [Nm]\n', DELTurb1(4))
% fprintf('        Blade3 Flapwise: %.2e  [Nm]\n', DELTurb1(5))
% fprintf('        Blade3 Edgewise: %.2e  [Nm]\n', DELTurb1(6))
fprintf('        Average Flapwise: %.2e  [Nm]\n', (DELTurb1(1)+DELTurb1(3)+DELTurb1(5))/3)
fprintf('        Average Edgewise: %.2e  [Nm]\n', (DELTurb1(2)+DELTurb1(4)+DELTurb1(6))/3)
fprintf('    PBD:\n');
% fprintf('        Blade1: %.2e [kNm deg]\n', PBDTurb1(1));
% fprintf('        Blade2: %.2e [kNm deg]\n', PBDTurb1(2));
% fprintf('        Blade3: %.2e [kNm deg]\n', PBDTurb1(3));
fprintf('        Average: %.2e [kNm deg]\n', mean(PBDTurb1));
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
fprintf('The output of Downstream WT:\n');
fprintf('    Power Production: %.2f [MW]\n', PowerTurb2);
fprintf('    DEL:\n');
% fprintf('        Blade1 Flapwise: %.2e  [Nm]\n', DELTurb2(1))
% fprintf('        Blade1 Edgewise: %.2e  [Nm]\n', DELTurb2(2))
% fprintf('        Blade2 Flapwise: %.2e  [Nm]\n', DELTurb2(3))
% fprintf('        Blade2 Edgewise: %.2e  [Nm]\n', DELTurb2(4))
% fprintf('        Blade3 Flapwise: %.2e  [Nm]\n', DELTurb2(5))
% fprintf('        Blade3 Edgewise: %.2e  [Nm]\n', DELTurb2(6))
fprintf('        Average Flapwise: %.2e  [Nm]\n', (DELTurb2(1)+DELTurb2(3)+DELTurb2(5))/3)
fprintf('        Average Edgewise: %.2e  [Nm]\n', (DELTurb2(2)+DELTurb2(4)+DELTurb2(6))/3)
fprintf('    PBD:\n');
% fprintf('        Blade1: %.2e [kNm deg]\n', PBDTurb2(1));
% fprintf('        Blade2: %.2e [kNm deg]\n', PBDTurb2(2));
% fprintf('        Blade3: %.2e [kNm deg]\n', PBDTurb2(3));
fprintf('        Average: %.2e [kNm deg]\n', mean(PBDTurb2));