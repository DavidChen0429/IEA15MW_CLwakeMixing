%% Visualize Metrices for Experiments
clear
close all 
addpath('.\Functions');
%clc

%% Data file (Chage this accordingly)
turbineName = '.\Data\NREL5MW\';
caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\2TurbinesNew\';

% Different case
basefile = '2Turbines_Baseline_4D.mat';
UniformOLfile = '2Turbines_OL_Helix_mag3_4D.mat';
UniformCLfile = '2Turbines_CL_Helix_ramp&stop_mag3_4D.mat';
Shear2OLfile = '2Turbines_OL_Helix_Shear0.2_mag3_4D.mat';
Shear2CLfile = '2Turbines_CL_Helix_Shear0.2_mag3_4D.mat';
Shear3OLfile = '2Turbines_OL_Helix_TI6&Shear0.2_mag3_4D.mat';
Shear3CLfile = '2Turbines_CL_Helix_TI6&Shear0.2_mag3_4D.mat';    
TI6OLfile = '2Turbines_OL_Helix_TI6_mag3_4D.mat';
TI6CLfile = '2Turbines_CL_Helix_TI6_mag3_4D.mat';      
BothOLfile = '2Turbines_OL_Helix_TI6&Shear0.2_mag3_4D.mat';
BothCLfile = '2Turbines_CL_Helix_TI6&Shear0.2_mag3_4D.mat';         

tic
Baseline = load([turbineName caseName basefile]);
UniformOL = load([turbineName caseName UniformOLfile]);
UniformCL = load([turbineName caseName UniformCLfile]);
Shear2OL = load([turbineName caseName Shear2OLfile]);
Shear2CL = load([turbineName caseName Shear2CLfile]);
Shear3OL = load([turbineName caseName Shear3OLfile]);
Shear3CL = load([turbineName caseName Shear3CLfile]);
TI6OL = load([turbineName caseName TI6OLfile]);
TI6CL = load([turbineName caseName TI6CLfile]);
BothOL = load([turbineName caseName BothOLfile]);
BothCL = load([turbineName caseName BothCLfile]);
toc

%% Data save to this file 
turbineName = '.\Data\NREL5MW\';
caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\QBladeDeug\';
fileName = 'inputList.mat';

%% Save input data for FASTER SIMULATION
% FFbeta HFbeta Pitch
% Uniform
InputData.Uniform.OL.FF_beta = UniformOL.FF_beta;
InputData.Uniform.OL.HF_beta = UniformOL.HF_beta;
InputData.Uniform.OL.Pitch = UniformOL.PitchAngles;
InputData.Uniform.CL.FF_beta = UniformCL.FF_beta;
InputData.Uniform.CL.HF_beta = UniformCL.HF_beta;
InputData.Uniform.CL.Pitch = UniformCL.PitchAngles;

% Shear
InputData.Shear2.OL.FF_beta = Shear2OL.FF_beta;
InputData.Shear2.OL.HF_beta = Shear2OL.HF_beta;
InputData.Shear2.OL.Pitch = Shear2OL.PitchAngles;
InputData.Shear2.CL.FF_beta = Shear2CL.FF_beta;
InputData.Shear2.CL.HF_beta = Shear2CL.HF_beta;
InputData.Shear2.CL.Pitch = Shear2CL.PitchAngles;

InputData.Shear3.OL.FF_beta = Shear3OL.FF_beta;
InputData.Shear3.OL.HF_beta = Shear3OL.HF_beta;
InputData.Shear3.OL.Pitch = Shear3OL.PitchAngles;
InputData.Shear3.CL.FF_beta = Shear3CL.FF_beta;
InputData.Shear3.CL.HF_beta = Shear3CL.HF_beta;
InputData.Shear3.CL.Pitch = Shear3CL.PitchAngles;

% Turbulence 
InputData.TurbTI6.OL.FF_beta = TI6OL.FF_beta;
InputData.TurbTI6.OL.HF_beta = TI6OL.HF_beta;
InputData.TurbTI6.OL.Pitch = TI6OL.PitchAngles;
InputData.TurbTI6.CL.FF_beta = TI6CL.FF_beta;
InputData.TurbTI6.CL.HF_beta = TI6CL.HF_beta;
InputData.TurbTI6.CL.Pitch = TI6CL.PitchAngles;

% Both
InputData.Both.OL.FF_beta = BothOL.FF_beta;
InputData.Both.OL.HF_beta = BothOL.HF_beta;
InputData.Both.OL.Pitch = BothOL.PitchAngles;
InputData.Both.CL.FF_beta = BothCL.FF_beta;
InputData.Both.CL.HF_beta = BothCL.HF_beta;
InputData.Both.CL.Pitch = BothCL.PitchAngles;

%% Save Data
save([turbineName caseName fileName], 'InputData')