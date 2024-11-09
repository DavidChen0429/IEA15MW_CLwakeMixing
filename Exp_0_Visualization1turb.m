%% Data Analysis for Experiments
clear
close all 
addpath('.\Functions');
%clc

%% Data file (Chage this accordingly)
turbineName = '.\Data\NREL5MW\';
caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\1Turbine\';
basefile = '1Turbine_Baseline.mat';
OLfileName = '1Turbine_OL_Helix_mag3.mat';
CLfileName = '1Turbine_CL_Helix_MIMO_step_mag3.mat';

Baseline = load([turbineName caseName basefile]);
OL = load([turbineName caseName OLfileName]);
CL = load([turbineName caseName CLfileName]);

%% Basic Settings
D_NREL5MW = 126;
U_inflow = 10;
timeStep = 0.1;
filter = 1000;
simLength = length(Baseline.Cp_store);

%% Visualization
t = (1:(simLength-filter+1)) * timeStep;
lw = 2;
color0 = [0.4660 0.6740 0.1880];
color1 = [0, 0.4470, 0.7410];
color2 = [0.8500, 0.3250, 0.0980];

%% ============== Overall Result (Input & Output)
% Hub Jet
figure('Name', 'Experiment Input-Output', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 2, 2)
plot(t, OL.HF_helixCenter_filtered(filter:end, 1), '--','Color',color1,'LineWidth', lw)
hold on
plot(t, OL.HF_helixCenter_filtered(filter:end, 2), '--', 'Color',color2,'LineWidth', lw)
plot(t, CL.HF_helixCenter_filtered(filter:end, 1), 'Color',color1, 'LineWidth', lw)
plot(t, CL.HF_helixCenter_filtered(filter:end, 2), 'Color',color2, 'LineWidth', lw)
plot(t, CL.r(filter:end, 2), 'k:', 'LineWidth', lw)
plot(t, CL.r(filter:end, 2), 'k:', 'LineWidth', lw)
yline(0, '--', 'LineWidth', lw)
hold off
title('Output: Helix Frame')
xlim([0 t(end)])
xlabel('Time [s]')
ylim([-5 15])
ylabel('Magnitude [m]')
legend('z^e_b','y^e_b','z^e','y^e','r_z','r_y','Location','southeast')
subplot(2, 2, 4)
plot(t, OL.FF_helixCenter_filtered(filter:end, 1), '--','Color',color1, 'LineWidth', lw)
hold on
plot(t, OL.FF_helixCenter_filtered(filter:end, 2), '--','Color',color2, 'LineWidth', lw)
plot(t, CL.FF_helixCenter_filtered(filter:end, 1),'Color',color1,'LineWidth', lw)
plot(t, CL.FF_helixCenter_filtered(filter:end, 2),'Color',color2,'LineWidth', lw)
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
plot(t, OL.HF_beta(filter:end, 1), '--','Color',color1, 'LineWidth', lw)
hold on
plot(t, OL.HF_beta(filter:end, 2), '--','Color',color2, 'LineWidth', lw)
plot(t, CL.u(filter:end, 1), 'Color',color1, 'LineWidth', lw)
plot(t, CL.u(filter:end, 2), 'Color',color2, 'LineWidth', lw)
hold off
title('Input: Helix Frame')
xlim([0 t(end)])
xlabel('Time [s]')
ylim([-1 10])
ylabel('Magnitude [deg]')
legend('\beta^e_{tilt,b}','\beta^e_{yaw,b}','\beta^e_{tilt}','\beta^e_{yaw}','Location','southeast')
subplot(2, 2, 3)
plot(t, OL.FF_beta(filter:end, 1), '--','Color',color1, 'LineWidth', lw)
hold on
plot(t, OL.FF_beta(filter:end, 2), '--','Color',color2, 'LineWidth', lw)
plot(t, CL.FF_beta(filter:end, 1), 'Color',color1, 'LineWidth', lw)
plot(t, CL.FF_beta(filter:end, 2), 'Color',color2, 'LineWidth', lw)
hold off
title('Input: Fixed Frame')
xlim([0 t(end)])
xlabel('Time [s]')
ylim([-12.5 12.5])
ylabel('Magnitude [deg]')
legend('\beta_{tilt,b}','\beta_{yaw,b}','\beta_{tilt}','\beta_{yaw}','Location','southeast')
setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',15,'linewidth',lw)

%% ============== Hub Jet Trajectory
center_bl = mean(Baseline.FF_helixCenter_filtered(3000:end, :));
center_ol = mean(OL.FF_helixCenter_filtered(3000:end, :));
center_cl = mean(CL.FF_helixCenter_filtered(3000:end, :));
figure('Name', 'Experiment HubJet Trajectory', 'NumberTitle', 'off', 'Position', [100, 100, 600, 600]);
plot(Baseline.FF_helixCenter_filtered(2000:end, 2), Baseline.FF_helixCenter_filtered(2000:end, 1), 'Color',color0, 'LineWidth', lw)
hold on
plot(OL.FF_helixCenter_filtered(2000:end, 2), OL.FF_helixCenter_filtered(2000:end, 1), 'Color',color1, 'LineWidth', lw)
plot(CL.FF_helixCenter_filtered(2000:end, 2), CL.FF_helixCenter_filtered(2000:end, 1), 'Color',color2, 'LineWidth', lw)
plot(center_bl(2), center_bl(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color0);
plot(center_ol(2), center_ol(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color1);
plot(center_cl(2), center_cl(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color2);
plot(0, 90, 'k*', 'MarkerSize', 10);
hold off
title('Hub Jet Trajectory')
xlabel('y [m]')
ylabel('z [m]')
xlim([-20 10])
ylim([77 107])
legend('Baseline', 'OL', 'CL', 'Location','southeast')
setfigpaper('Width',[15,1],'Interpreter','tex','FontSize',15,'linewidth',lw)

%% Video comparison
% ringVisualization2(Baseline, D_NREL5MW)
% videoCompare_func(Baseline,CL,D_NREL5MW,'')