%clc
clear
close all
addpath('.\Functions');

%% Load data
Fs = 10;  % sampling frequency Hz
Fc = 0.02;  % cutoff frequency Hz
timeStep = 0.1;

turbineName = '.\Data\NREL5MW\';
caseName = 'Sth\';
basicFile = 'basic.mat';
turbFile = 'turb.mat';
ovalFile = 'oval.mat';
flowerFile = 'flower.mat';

Baseline = load([turbineName caseName basicFile]);
TurbCase = load([turbineName caseName turbFile]);
OvalCase = load([turbineName caseName ovalFile]);
FlowerCase = load([turbineName caseName flowerFile]);
simLength = length(Baseline.FF_beta);

showSnapshotOption = 'Y';
showOvalOption = 'Y';
showFlowerOption = 'Y';

Font = 20;
lw = 2;
filter = 1000;
color1 = [0, 0.4470, 0.7410];
color2 = [0.8500, 0.3250, 0.0980];

%% 1. Show the Hub Jet
dataLiDAR_A = Baseline.LiDAR_data;
dataLiDAR_B = TurbCase.LiDAR_data;
D = 126;

% reference rotor disc
theta = linspace(0, 2*pi, 20);
y_1Dref = 0 + D/2 * cos(theta);
z_1Dref = 90 + D/2 * sin(theta);

figure('Position', [10, 10, 800, 310]);
counter = 1000;

% Plot the snapshot
subplot(1, 2, 1)
snapshot = dataLiDAR_A(counter);
y = snapshot.y;
z = snapshot.z;
scatter(y, z, 10, snapshot.u_los, 'filled');
hold on
plot(y_1Dref, z_1Dref, "k:", 'LineWidth',1);
hold off;
xlabel('Y [m]')
ylabel('Z [m]')
title('Uniform')
colorbar;
clim([4 10])

subplot(1, 2, 2)
snapshot2 = dataLiDAR_B(counter);
y = snapshot2.y;
z = snapshot2.z;
scatter(y, z, 10, snapshot2.u_los, 'filled');
hold on
plot(y_1Dref, z_1Dref, "k:", 'LineWidth',1);
hold off;
xlabel('Y [m]')
ylabel('Z [m]')
title('Turb')
colorbarHandle = colorbar;
ylabel(colorbarHandle, 'u [m/s]');
clim([4 10])

%% 2. Show Oval 
% Show Trajectory
figure('Name', 'HubJet Trajectory', 'NumberTitle', 'off', 'Position', [100, 100, 600, 600]);
plot(Baseline.FF_helixCenter_filtered(filter:end, 2), Baseline.FF_helixCenter_filtered(filter:end, 1), 'Color',color1, 'LineWidth', lw)
hold on
plot(OvalCase.FF_helixCenter_filtered(filter:end, 2), OvalCase.FF_helixCenter_filtered(filter:end, 1), 'Color',color2, 'LineWidth', lw)
hold off
title('Hub Jet Trajectory')
xlabel('y [m]')
ylabel('z [m]')
xlim([-30 20])
ylim([67 117])
legend('Baisc Helix', 'Oval Helix', 'Location','southeast')
setfigpaper('Width',[15,1],'Interpreter','tex','FontSize',Font,'linewidth',lw)

% Show Input Signal
t = (1:(simLength)) * timeStep;
figure('Name', 'Input Circle vs. Oval', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 1, 1)
plot(t, Baseline.FF_beta(:, 1),'--', 'Color', color1, 'LineWidth', lw)
hold on
plot(t, OvalCase.FF_beta(:, 1), 'Color', color1, 'LineWidth', lw)
hold off
legend('Circle', 'Oval')
xlim([0 t(end)])
xlabel('Time [s]')
ylabel('Magnitude [deg]')
title('\beta_{tilt}')
subplot(2, 1, 2)
plot(t, Baseline.FF_beta(:, 2),'--', 'Color', color2, 'LineWidth', lw)
hold on
plot(t, OvalCase.FF_beta(:, 2), 'Color', color2, 'LineWidth', lw)
hold off
legend('Circle', 'Oval')
xlim([0 t(end)])
xlabel('Time [s]')
ylabel('Magnitude [deg]')
title('\beta_{yaw}')
setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)
        
%% 3. Show Flower
figure('Name', 'Input Circle vs. Flower', 'NumberTitle', 'off', 'Position', [100, 100, 600, 600]);
plot(Baseline.FF_helixCenter_filtered(filter:end, 2), Baseline.FF_helixCenter_filtered(filter:end, 1), 'Color',color1, 'LineWidth', lw)
hold on
plot(FlowerCase.FF_helixCenter_filtered(filter:end, 2), FlowerCase.FF_helixCenter_filtered(filter:end, 1), 'Color',color2, 'LineWidth', lw)
hold off
title('Hub Jet Trajectory')
xlabel('y [m]')
ylabel('z [m]')
xlim([-30 20])
ylim([67 117])
legend('Baisc Helix', 'Flower Helix', 'Location','southeast')
setfigpaper('Width',[15,1],'Interpreter','tex','FontSize',Font,'linewidth',lw)

% Show Input Signal
t = (1:(simLength)) * timeStep;
figure('Name', 'Input Detail', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 1, 1)
plot(t, Baseline.FF_beta(:, 1),'--', 'Color', color1, 'LineWidth', lw)
hold on
plot(t, FlowerCase.FF_beta(:, 1), 'Color', color1, 'LineWidth', lw)
hold off
legend('Circle', 'Flower')
xlim([0 t(end)])
xlabel('Time [s]')
ylabel('Magnitude [deg]')
title('\beta_{tilt}')
subplot(2, 1, 2)
plot(t, Baseline.FF_beta(:, 2),'--', 'Color', color2, 'LineWidth', lw)
hold on
plot(t, FlowerCase.FF_beta(:, 2), 'Color', color2, 'LineWidth', lw)
hold off
legend('Circle', 'Flower')
xlim([0 t(end)])
xlabel('Time [s]')
ylabel('Magnitude [deg]')
title('\beta_{yaw}')
setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)