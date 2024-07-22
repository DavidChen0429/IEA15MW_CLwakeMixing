%% This file is created to study the correlation between tilt,yaw to wake center
clear
close all
addpath('.\Functions');

%% Fix Frame 
Fs = 10;  % sampling frequency Hz
Fc = 0.02;  % cutoff frequency Hz
fileName = '600s_Center_FF_baseline.mat';
fileName2 = '600s_Center_FF_yaw.mat';
dataPath = '.\Data\MAT\LiDAR_sampling\';
caseName = 'Uni\Str0.3_U8_1Dd_10Hz_CCW\';
SimData = load([dataPath caseName fileName]);
SimData2 = load([dataPath caseName fileName2]);
% MomentCenter_comparison_Visualization(SimData, SimData2, Fs, Fc)
wakeCenterTraj(SimData.LiDAR_data, SimData2.LiDAR_data, Fs, Fc)
% videoCompare_func(SimData, SimData2)
% MomentCenter_Visualization(SimData, Fs, Fc)

%% Helix Frame
Fs = 10;  % sampling frequency Hz
Fc = 0.02;  % cutoff frequency Hz
fileName = '600s_Center_HF_basecase.mat';
fileName2 = '600s_Center_HF_yaw_l2.mat';
dataPath = '.\Data\MAT\LiDAR_sampling\';
caseName = 'Uni\Str0.3_U8_1Dd_10Hz_CCW\';
SimData = load([dataPath caseName fileName]);
SimData2 = load([dataPath caseName fileName2]);
% MomentCenter_Visualization(SimData, Fs, Fc);
MomentCenter_comparison_Visualization(SimData, SimData2, Fs, Fc)
wakeCenterTraj(SimData.LiDAR_data, SimData2.LiDAR_data, Fs, Fc)
% videoCompare_func(SimData, SimData2)
% ringVisualization(SimData2.LiDAR_data)