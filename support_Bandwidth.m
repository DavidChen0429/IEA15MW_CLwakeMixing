%clc
%clear
close all
addpath('.\Functions');

%% Load data
Fs = 10;  % sampling frequency Hz
Fc = 0.02;  % cutoff frequency Hz
filterIndex = 800;
timeStep = 0.1;

turbineName = '.\Data\NREL5MW\';
caseName = 'IDE\';
fileName = 'bandwidthCheck.mat';

Baseline = load([turbineName caseName fileName]);

%% Same St, Different Amplitude
t = (1: (length(Baseline.FF_beta)-filterIndex+1)) * timeStep;
lw = 2;

figure('Name', 'Bandwidth See', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
plot(t, Baseline.HF_helixCenter_filtered(filterIndex:end, 1), 'LineWidth',lw);
hold on
plot(t, Baseline.HF_helixCenter_filtered(filterIndex:end, 2), 'LineWidth',lw);
plot(t, Baseline.HF_beta(filterIndex:end, 1), 'k:','LineWidth',lw);
plot(t, Baseline.HF_beta(filterIndex:end, 2), 'k:','LineWidth',lw);

yline(0, '--', 'LineWidth', 1)
hold off;
xlabel('Time [s]')
ylabel('Magnitude [m]')
title('Step Response')
xlim([0 200])
ylim([-1 5])
legend('z^{e}', 'y^{e}','\beta^e_{tilt}','\beta^e_{yaw}','Location','southeastoutside')
setfigpaper('Width',[30,0.25],'Interpreter','tex','FontSize',20,'linewidth',2)