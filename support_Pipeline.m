%clc
%clear
close all
addpath('.\Functions');

%% Load data
Fs = 10;  % sampling frequency Hz
Fc = 0.02;  % cutoff frequency Hz
filterIndex = 500;
timeStep = 0.1;

turbineName = '.\Data\NREL5MW\';
caseName = 'Pipeline\';
fileName1 = 'St3A3_30_step.mat';
fileName2 = 'St3A3_30.mat';

% St3A3 = load([turbineName caseName fileName1]);
St3A3 = load([turbineName caseName fileName2]);

%% Visualization
t = (1: (length(St3A3.FF_beta)-filterIndex+1)) * timeStep;
lw = 1;

figure('Name', 'Overall Result', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 2, 1)
plot(t, St3A3.FF_beta(filterIndex:end, 1),'LineWidth', lw);
hold on;
plot(t, St3A3.FF_beta(filterIndex:end, 2),'LineWidth', lw);
yline(0, '--', 'LineWidth', 1)
hold off;
xlabel('Time [s]')
ylabel('Magnitude [deg]')
xlim([0 450])
ylim([-5 5])
title('\beta Fixed Frame')
legend('\beta_{tilt}', '\beta_{yaw}', 'Location', 'southeast')

subplot(2, 2, 3);
plot(t, St3A3.HF_beta(filterIndex:end, 1),'LineWidth', lw);
hold on;
plot(t, St3A3.HF_beta(filterIndex:end ,2),'LineWidth', lw);
yline(0, '--', 'LineWidth', 1)
hold off;
xlabel('Time [s]')
ylabel('Magnitude')
xlim([0 450])
ylim([0 3.5])
title('\beta^e Helix Frame')
legend('\beta^e_{tilt}', '\beta^e_{yaw}', 'Location', 'southeast')

subplot(2, 2, 2)
plot(t, St3A3.FF_helixCenter(filterIndex:end, 1) - 92.0026,'LineWidth', lw);
hold on;
plot(t, St3A3.FF_helixCenter(filterIndex:end, 2) + 4.0999,'LineWidth', lw);
plot(t, St3A3.FF_helixCenter_filtered(filterIndex:end, 1) - 92.0026,'LineWidth', lw);
plot(t, St3A3.FF_helixCenter_filtered(filterIndex:end, 2) + 4.0999,'LineWidth', lw);
yline(0, '--', 'LineWidth', 1)
hold off;
xlabel('Time [s]')
ylabel('Magnitude [m]')
xlim([0 450])
ylim([-20 20])
title('Hub Jet Fixed Frame')
legend('z', 'y', 'z_f', 'y_f', 'Location', 'southeast')

subplot(2, 2, 4)
plot(t, St3A3.HF_helixCenter_filtered(filterIndex:end, 1),'LineWidth', lw);
hold on
plot(t, St3A3.HF_helixCenter_filtered(filterIndex:end, 2),'LineWidth', lw);
yline(0, '--', 'LineWidth', 1)
hold off;
xlabel('Time [s]')
ylabel('Magnitude')
xlim([0 450])
ylim([0 10])
title('Hub Jet Helix Frame')
legend('z^{e}', 'y^{e}', 'Location', 'southeast')

setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',1)

%% 2*2 version 
t = (1: (length(St3A3.FF_beta)-filterIndex+1)) * timeStep;
lw = 1;

figure('Name', 'Overall Result', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 1, 1)
plot(t, St3A3.FF_helixCenter(filterIndex:end, 1) - 92.0026,'LineWidth', lw);
hold on;
plot(t, St3A3.FF_helixCenter(filterIndex:end, 2) + 4.0999,'LineWidth', lw);
plot(t, St3A3.FF_helixCenter_filtered(filterIndex:end, 1) - 92.0026,'LineWidth', lw);
plot(t, St3A3.FF_helixCenter_filtered(filterIndex:end, 2) + 4.0999,'LineWidth', lw);
yline(0, '--', 'LineWidth', 1)
hold off;
xlabel('Time [s]')
ylabel('Position [m]')
xlim([0 450])
ylim([-20 20])
title('Fixed Frame')
legend('z', 'y', 'z_f', 'y_f', 'Location', 'southeast')

subplot(2, 1, 2)
plot(t, St3A3.HF_helixCenter_filtered(filterIndex:end, 1),'LineWidth', lw);
hold on
plot(t, St3A3.HF_helixCenter_filtered(filterIndex:end, 2),'LineWidth', lw);
yline(0, '--', 'LineWidth', 1)
hold off;
xlabel('Time [s]')
ylabel('Position [m]')
xlim([0 450])
ylim([-5 15])
title('Helix Frame')
legend('z^{e}', 'y^{e}', 'Location', 'southeast')

setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',1)
