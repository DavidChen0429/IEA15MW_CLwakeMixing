%% Data Analysis for Experiments
clear
close all 
addpath('.\Functions');
%clc

%% Data file (Chage this accordingly)
turbineName = '.\Data\NREL5MW\';
caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\';
fileName = '1Turbines_CL_Helix_SISO_I.mat';
basefile = '1Turbines_OL_Helix.mat';

Data = load([turbineName caseName fileName]);
Baseline = load([turbineName caseName basefile]);

% Baseline info
Cp_store_bl = Baseline.Cp_store;
Moop1_store_bl = Baseline.Moop1_store;
PitchAngles_bl = Baseline.PitchAngles;
Mflap1_store_bl = Baseline.Mflap1_store;
Medge1_store_bl = Baseline.Medge1_store;

% WT1 info 
Cp_store = Data.Cp_store;
Moop1_store = Data.Moop1_store;
PitchAngles = Data.PitchAngles;
Mflap1_store = Data.Mflap1_store;
Medge1_store = Data.Medge1_store;

%% Basic Settings
D_NREL5MW = 126;
U_inflow = 10;
timeStep = 0.1;
filter = 1000;
simLength = length(Baseline.Cp_store);

%% Calcuate Power, DEL, PBD
% Baseline
PowerTurb_bl = calculatePower(filter,Cp_store_bl,D_NREL5MW,U_inflow); % [MW]
DELTurb_bl = calculateDEL(filter,Moop1_store_bl, timeStep); % [Nm]
PBDTurb_bl = calculatePBD(filter,PitchAngles_bl,Mflap1_store_bl,Medge1_store_bl); % [kNm deg]
fprintf('================================================== \n');
fprintf('The output of Baseline:\n');
fprintf('    Power Production: %.2f [MW]\n', PowerTurb_bl);
fprintf('    DEL: %.2e  [Nm]\n', DELTurb_bl);
fprintf('    PBD: %.2e [kNm deg]\n', PBDTurb_bl);
% Controlled
PowerTurb1 = calculatePower(filter,Cp_store,D_NREL5MW,U_inflow); % [MW]
DELTurb1 = calculateDEL(filter,Moop1_store, timeStep); % [Nm]
PBDTurb1 = calculatePBD(filter,PitchAngles,Mflap1_store,Medge1_store); % [kNm deg]
% fprintf('================================================== \n');
fprintf('The output of Wind turbine 1:\n');
fprintf('    Power Production: %.2f [MW]\n', PowerTurb1);
fprintf('    DEL: %.2e  [Nm]\n', DELTurb1);
fprintf('    PBD: %.2e [kNm deg]\n', PBDTurb1);

%% Visualization
t = (1:(simLength-filter+1)) * timeStep;
lw = 1;
color1 = [0, 0.4470, 0.7410];  % Equivalent to "#0072BD"
color2 = [0.8500, 0.3250, 0.0980];  % Equivalent to "#D95319"

% Hub Jet
figure('Name', 'Experiment Input-Output', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 2, 2)
plot(t, Baseline.HF_helixCenter_filtered(filter:end, 1), '--','Color',color1,'LineWidth', lw)
hold on
plot(t, Baseline.HF_helixCenter_filtered(filter:end, 2), '--', 'Color',color2,'LineWidth', lw)
plot(t, Data.HF_helixCenter_filtered(filter:end, 1), 'Color',color1, 'LineWidth', lw)
plot(t, Data.HF_helixCenter_filtered(filter:end, 2), 'Color',color2, 'LineWidth', lw)
plot(t, Data.r(filter:end, 2), 'k:', 'LineWidth', lw*2)
plot(t, Data.r(filter:end, 2), 'k:', 'LineWidth', lw*2)
yline(0, '--', 'LineWidth', lw*1.5)
hold off
title('Output: Helix Frame')
xlim([0 t(end)])
xlabel('Time [s]')
ylim([-5 15])
ylabel('Magnitude')
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
ylim([-1 6])
ylabel('Magnitude')
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
ylim([-10 10])
ylabel('Magnitude [deg]')
legend('\beta_{tilt,b}','\beta_{yaw,b}','\beta_{tilt}','\beta_{yaw}','Location','southeast')
setfigpaper('Width',[40,0.5],'Interpreter','tex','FontSize',15,'linewidth',lw)