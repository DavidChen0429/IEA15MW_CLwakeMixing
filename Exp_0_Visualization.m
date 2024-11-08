%% Data Analysis for Experiments
clear
close all 
addpath('.\Functions');
%clc

%% Data file (Chage this accordingly)
turbineName = '.\Data\NREL5MW\';
caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\';
basefile = '2Turbines_OL_Helix_TI5_mag3.mat';
fileName = '2Turbines_CL_Helix_TI5_mag3.mat';

Baseline = load([turbineName caseName basefile]);
Data = load([turbineName caseName fileName]);

%% Basic Settings
D_NREL5MW = 126;
U_inflow = 10;
timeStep = 0.1;
filter = 1000;
simLength = length(Baseline.Cp_store);

%% Visualization
t = (1:(simLength-filter+1)) * timeStep;
lw = 2;
color1 = [0, 0.4470, 0.7410];  % Equivalent to "#0072BD"
color2 = [0.8500, 0.3250, 0.0980];  % Equivalent to "#D95319"

% ============== Overall Result (Input & Output)
% Hub Jet
figure('Name', 'Experiment Input-Output', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 2, 2)
plot(t, Baseline.HF_helixCenter_filtered(filter:end, 1), '--','Color',color1,'LineWidth', lw)
hold on
plot(t, Baseline.HF_helixCenter_filtered(filter:end, 2), '--', 'Color',color2,'LineWidth', lw)
plot(t, Data.HF_helixCenter_filtered(filter:end, 1), 'Color',color1, 'LineWidth', lw)
plot(t, Data.HF_helixCenter_filtered(filter:end, 2), 'Color',color2, 'LineWidth', lw)
plot(t, Data.r(filter:end, 2), 'k:', 'LineWidth', lw)
plot(t, Data.r(filter:end, 2), 'k:', 'LineWidth', lw)
yline(0, '--', 'LineWidth', lw)
hold off
title('Output: Helix Frame')
xlim([0 t(end)])
xlabel('Time [s]')
ylim([-5 15])
ylabel('Magnitude [m]')
legend('z^e_b','y^e_b','z^e','y^e','r_z','r_y','Location','southeast')
subplot(2, 2, 4)
plot(t, Baseline.FF_helixCenter_filtered(filter:end, 1), '--','Color',color1, 'LineWidth', lw)
hold on
plot(t, Baseline.FF_helixCenter_filtered(filter:end, 2), '--','Color',color2, 'LineWidth', lw)
plot(t, Data.FF_helixCenter_filtered(filter:end, 1),'Color',color1,'LineWidth', lw)
plot(t, Data.FF_helixCenter_filtered(filter:end, 2),'Color',color2,'LineWidth', lw)
hold off
title('Output: Fixed Frame')
xlim([0 t(end)])
xlabel('Time [s]')
ylim([-50 150])
ylabel('Position [m]')
legend('z^e_b','y^e_b','z^e','y^e','Location','southeast')
setfigpaper('Width',[40,0.5],'Interpreter','tex','FontSize',15,'linewidth',lw)

% Control Input
% figure('Name', 'Experiment Result: Hub Jet', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 2, 1)
plot(t, Baseline.HF_beta(filter:end, 1), '--','Color',color1, 'LineWidth', lw)
hold on
plot(t, Baseline.HF_beta(filter:end, 2), '--','Color',color2, 'LineWidth', lw)
plot(t, Data.u(filter:end, 1), 'Color',color1, 'LineWidth', lw)
plot(t, Data.u(filter:end, 2), 'Color',color2, 'LineWidth', lw)
hold off
title('Input: Helix Frame')
xlim([0 t(end)])
xlabel('Time [s]')
ylim([-1 10])
ylabel('Magnitude [deg]')
legend('\beta^e_{tilt,b}','\beta^e_{yaw,b}','\beta^e_{tilt}','\beta^e_{yaw}','Location','southeast')
subplot(2, 2, 3)
plot(t, Baseline.FF_beta(filter:end, 1), '--','Color',color1, 'LineWidth', lw)
hold on
plot(t, Baseline.FF_beta(filter:end, 2), '--','Color',color2, 'LineWidth', lw)
plot(t, Data.FF_beta(filter:end, 1), 'Color',color1, 'LineWidth', lw)
plot(t, Data.FF_beta(filter:end, 2), 'Color',color2, 'LineWidth', lw)
hold off
title('Input: Fixed Frame')
xlim([0 t(end)])
xlabel('Time [s]')
ylim([-12.5 12.5])
ylabel('Magnitude [deg]')
legend('\beta_{tilt,b}','\beta_{yaw,b}','\beta_{tilt}','\beta_{yaw}','Location','southeast')
setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',15,'linewidth',lw)

% ============== Hub Jet Trajectory
center_bl = mean(Baseline.FF_helixCenter_filtered(3000:end, :));
center_cl = mean(Data.FF_helixCenter_filtered(3000:end, :));
figure('Name', 'Experiment HubJet Trajectory', 'NumberTitle', 'off', 'Position', [100, 100, 600, 600]);
plot(Baseline.FF_helixCenter_filtered(2000:end, 2), Baseline.FF_helixCenter_filtered(2000:end, 1), 'Color',color1, 'LineWidth', lw)
hold on
plot(Data.FF_helixCenter_filtered(2000:end, 2), Data.FF_helixCenter_filtered(2000:end, 1), 'Color',color2, 'LineWidth', lw)
plot(center_bl(2), center_bl(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color1);
plot(center_cl(2), center_cl(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color2);
plot(0, 90, 'k*', 'MarkerSize', 10);
hold off
title('Hub Jet Trajectory')
xlabel('y [m]')
ylabel('z [m]')
xlim([-20 10])
ylim([77 107])
legend('Baseline', 'Controlled', 'Location','southeast')
setfigpaper('Width',[15,1],'Interpreter','tex','FontSize',15,'linewidth',lw)

%% Video comparison
ringVisualization(Baseline, D_NREL5MW)
videoCompare_func(Baseline,Data,D_NREL5MW,'')