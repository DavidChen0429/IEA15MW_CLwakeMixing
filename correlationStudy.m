%% This file is created to study the correlation between tilt,yaw to wake center
clear
close all
addpath('.\Functions');

%% Fix Frame 
Fs = 10;  % sampling frequency Hz
Fc = 0.05;  % cutoff frequency Hz
fileName = 'FF_Uni_basecase.mat';   % Fixed Frame
fileName2 = 'FF_Uni_inflowAngle.mat';   % Fixed Frame
turbineName = '.\Data\NREL5MW\';
caseName = 'Str0.3_U10_1Dd_10Hz_CCW\CRstudy\';
SimData = load([turbineName caseName fileName]);
SimData2 = load([turbineName caseName fileName2]);
% BetaCenter_Comparison_Visualization(SimData, SimData2, Fs, Fc)
% wakeCenterTraj(SimData, SimData2, Fs, Fc)
videoCompare_func(SimData, SimData2, Fs, Fc, 126, ".\Data\inflowAngle_noCen.avi")
% BetaCenter_Visualization(SimData, Fs, Fc)
% ringVisualization(SimData.LiDAR_data, 126)

%% Bias and center's center relationship
close all
Fs = 10;  % sampling frequency Hz
Fc = 0.02;  % cutoff frequency Hz
fileName = '600s_Center_FF_t,y,bias3.mat';
dataPath = '.\Data\MAT\LiDAR_sampling\';
caseName = 'Uni\Str0.3_U8_1Dd_10Hz_CCW\';
SimData = load([dataPath caseName fileName]);
MomentCenter_Visualization(SimData, Fs, Fc)

data = SimData.LiDAR_data;
wakeCenterYo = arrayfun(@(x) x.centerY, data);
wakeCenterZo = arrayfun(@(x) x.centerZ, data);
% wakeCenterYf = lowpassFilter(wakeCenterYo, Fs, Fc);
% wakeCenterZf = lowpassFilter(wakeCenterZo, Fs, Fc);
meanY = slideWindow(wakeCenterYo, 1000, 300);
meanZ = slideWindow(wakeCenterZo, 1000, 300);

figure();
plot(meanY, meanZ, "*");
hold on;
plot(lowpassFilter(wakeCenterYo, Fs, Fc), lowpassFilter(wakeCenterZo, Fs, Fc));
hold off;
xlabel('Y [m]')
ylabel('Z [m]')
xlim([-50 50])
ylim([100 200])
meanValueY = mean(meanY(2:end));
meanValueZ = mean(meanZ(2:end));
fprintf('Y: %.2f, and Z: %.2f\n', meanValueY, meanValueZ);

% rough conclusion
Mtilt_unit = [2.72-(-7.525) 172.66-158.54]/5;   % (2.049 2.824)
Myaw_unit = [2.72-(-19.82) 162.80-158.54]/5;    % (4.508 0.852)


%% Helix Frame
Fs = 10;  % sampling frequency Hz
Fc = 0.05;  % cutoff frequency Hz
fileName = 'HF_Uni_basecase.mat';  
fileName2 = 'HF_Uni_tilt,yaw.mat';   
turbineName = '.\Data\NREL5MW\';
caseName = 'Str0.3_U10_1Dd_10Hz_CCW\CRstudy\';
SimData = load([turbineName caseName fileName]);
SimData2 = load([turbineName caseName fileName2]);
BetaCenter_Comparison_Visualization(SimData, SimData2, Fs, Fc)
wakeCenterTraj(SimData, SimData2, Fs, Fc)
% videoCompare_func(SimData, SimData2, Fs, Fc, 126, ".\Data\inflowAngle.avi")
% BetaCenter_Visualization(SimData, Fs, Fc)
% ringVisualization(SimData.LiDAR_data, 126)

%% Compare different distance and phase info difference
Fs = 10;  % sampling frequency Hz
Fc = 0.05;  % cutoff frequency Hz
fileName = '600s_Center_FF_baseline.mat';
fileName2 = '600s_Center_HF_NTM_A_basecase.mat';
dataPath = '.\Data\MAT\LiDAR_sampling\';
caseName = 'Uni\Str0.3_U8_1Dd_10Hz_CCW\';
case2Name = 'Turbulence\Str0.3_U8_1Dd_10Hz_CCW\';
SimData = load([dataPath caseName fileName]);
SimData2 = load([dataPath case2Name fileName2]);
MomentCenter_comparison_Visualization(SimData, SimData2, Fs, Fc)
wakeCenterTraj(SimData.LiDAR_data, SimData2.LiDAR_data, Fs, Fc)
% videoCompare_func(SimData, SimData2)
% MomentCenter_Visualization(SimData, Fs, Fc)