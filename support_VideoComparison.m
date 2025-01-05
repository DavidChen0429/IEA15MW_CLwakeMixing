%% This file is created to study the correlation between tilt,yaw to wake center
clear
close all
addpath('.\Functions');

%% Comparison
Fs = 10;  % sampling frequency Hz
Fc = 0.02;  % cutoff frequency Hz
fileName = 'zeroAnimation.mat';  
fileName2 = 'basicAnimation.mat';   
turbineName = '.\Data\NREL5MW\';
caseName = 'Sth\';
SimData = load([turbineName caseName fileName]);
SimData2 = load([turbineName caseName fileName2]);
% BetaCenter_Comparison_Visualization(SimData, SimData2, Fs, Fc)
% wakeCenterTraj(SimData, SimData2, Fs, Fc)
videoCompare_func(SimData, SimData2, 126, ".\Data\demo.avi")
% BetaCenter_Visualization(SimData, Fs, Fc)
% ringVisualization(SimData.LiDAR_data, 126)

%% Only one
Fs = 10;  % sampling frequency Hz
Fc = 0.02;  % cutoff frequency Hz
fileName = 'basicAnimation.mat';  
caseName = 'Sth\';
turbineName = '.\Data\NREL5MW\';
SimData = load([turbineName caseName fileName]);
videoDemo_func(SimData, 126, ".\Data\demoSingle.avi")