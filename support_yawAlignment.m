%% Data Analysis for Experiments
% clear
close all 
addpath('.\Functions');
%clc

%% Data file (Chage this accordingly)
turbineName = '.\Data\NREL5MW\';
% caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\2TurbinesNew\';
caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\2TurbinesLonger\';

% Different case
basefile = '2Turbines_Baseline_4D.mat';
basefileYaw = '2Turbines_Baseline_4D_1inflowAngle.mat';
OLfileName = '2Turbines_OL_Helix_mag3_4D.mat';
OLfileNameYaw = '2Turbines_OL_Helix_mag3_4D_1inflowAngle.mat';
CLfileName = '2Turbines_CL_Helix_ramp&stop_mag3_4D.mat';
CLfileNameYaw = '2Turbines_CL_Helix_ramp&stop_mag3_4D_1inflowAngle.mat';

Baseline = load([turbineName caseName basefile]);
BaselineYaw = load([turbineName caseName basefileYaw]);
OL = load([turbineName caseName OLfileName]);
OLYaw = load([turbineName caseName OLfileNameYaw]);
CL = load([turbineName caseName CLfileName]);
CLYaw = load([turbineName caseName CLfileNameYaw]);

% Basic settings
D_NREL5MW = 126;
U_inflow = 10;
timeStep = 0.1;
filter = 3000;

%% Visualization
trajOption = 'Y';
overallDetailOption = 'Y';

simLength = length(Baseline.Power_store);
t = (1:(simLength-filter+1)) * timeStep;
lw = 2;
Font = 20;
color0 = [0.4660 0.6740 0.1880];
color1 = [0, 0.4470, 0.7410];
color2 = [0.8500, 0.3250, 0.0980];
color3 = [0.3010 0.7450 0.9330];

if strcmp(trajOption, 'Y')
    center_bl = mean(Baseline.FF_helixCenter_filtered(filter:end, :));
    center_blyaw = mean(BaselineYaw.FF_helixCenter_filtered(filter:end, :));
    center_ol = mean(OL.FF_helixCenter_filtered(filter:end, :));
    center_olyaw = mean(OLYaw.FF_helixCenter_filtered(filter:end, :));
    center_cl = mean(CL.FF_helixCenter_filtered(filter:end, :));
    center_clyaw = mean(CLYaw.FF_helixCenter_filtered(filter:end, :));
    
    figure('Name', 'HubJet Trajectory', 'NumberTitle', 'off', 'Position', [100, 100, 1050, 300]);
    subplot(1, 3, 1)
    plot(Baseline.FF_helixCenter_filtered(filter:end, 2), Baseline.FF_helixCenter_filtered(filter:end, 1), 'Color',color0, 'LineWidth', lw)
    hold on
    plot(BaselineYaw.FF_helixCenter_filtered(filter:end, 2), BaselineYaw.FF_helixCenter_filtered(filter:end, 1), 'Color',color1, 'LineWidth', lw)
    plot(center_bl(2), center_bl(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color0);
    plot(center_blyaw(2), center_blyaw(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color1);
    plot(0, 90, 'k*', 'MarkerSize', 10);
    hold off
    title('Baseline')
    xlabel('y [m]')
    ylabel('z [m]')
    xlim([-30 20])
    ylim([67 117])
    legend('Ideal', 'Skewed', 'Location','southeast')

    subplot(1, 3, 2)
    plot(OL.FF_helixCenter_filtered(filter:end, 2), OL.FF_helixCenter_filtered(filter:end, 1), 'Color',color0, 'LineWidth', lw)
    hold on
    plot(OLYaw.FF_helixCenter_filtered(filter:end, 2), OLYaw.FF_helixCenter_filtered(filter:end, 1), 'Color',color2, 'LineWidth', lw)
    plot(center_ol(2), center_ol(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color0);
    plot(center_olyaw(2), center_olyaw(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color2);
    plot(0, 90, 'k*', 'MarkerSize', 10);
    hold off
    title('Open-loop')
    xlabel('y [m]')
    ylabel('z [m]')
    xlim([-30 20])
    ylim([67 117])
    legend('Ideal', 'Skewed', 'Location','southeast')

    subplot(1, 3, 3)
    plot(CL.FF_helixCenter_filtered(filter:end, 2), CL.FF_helixCenter_filtered(filter:end, 1), 'Color',color0, 'LineWidth', lw)
    hold on
    plot(CLYaw.FF_helixCenter_filtered(filter:end, 2), CLYaw.FF_helixCenter_filtered(filter:end, 1), 'Color',color3, 'LineWidth', lw)
    plot(center_cl(2), center_cl(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color0);
    plot(center_clyaw(2), center_clyaw(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color2);
    plot(0, 90, 'k*', 'MarkerSize', 10);
    hold off
    title('Closed-loop')
    xlabel('y [m]')
    ylabel('z [m]')
    xlim([-30 20])
    ylim([67 117])
    legend('Ideal', 'Skewed', 'Location','southeast')
    setfigpaper('Width',[40,0.3],'Interpreter','tex','FontSize',Font,'linewidth',lw)
end

if strcmp(overallDetailOption, 'Y')
    % input HF
    figure('Name', 'Input Detail', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(2, 1, 1)
    plot(t, OL.HF_beta(filter:end, 1),'Color',color0,'LineWidth', lw)
    hold on
    plot(t, OLYaw.HF_beta(filter:end, 1), 'Color',color1,'LineWidth', lw)
    plot(t, CLYaw.HF_beta(filter:end, 1), 'Color',color2,'LineWidth', lw)
    hold off
    title('\beta^e_{tilt}')
    xlim([0 t(end)])
    xlabel('Time [s]')
    ylabel('Magnitude [m]')
    legend('Ideal','OL','CL','Location','southeast')
    subplot(2, 1, 2)
    plot(t, OL.HF_beta(filter:end, 2),'Color',color0,'LineWidth', lw)
    hold on
    plot(t, OLYaw.HF_beta(filter:end, 2), 'Color',color1,'LineWidth', lw)
    plot(t, CLYaw.HF_beta(filter:end, 2), 'Color',color2,'LineWidth', lw)
    hold off
    title('\beta^e_{yaw}')
    xlim([0 t(end)])
    xlabel('Time [s]')
    ylabel('Magnitude [m]')
    legend('Ideal','OL','CL','Location','southeast')
    setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)

    % output HF
    figure('Name', 'Output Detail', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(2, 1, 2)
    plot(t, OL.HF_helixCenter_filtered(filter:end, 1),'Color',color0,'LineWidth', lw)
    hold on
    plot(t, OLYaw.HF_helixCenter_filtered(filter:end, 1), 'Color',color1,'LineWidth', lw)
    plot(t, CLYaw.HF_helixCenter_filtered(filter:end, 1), 'Color',color2,'LineWidth', lw)
    yline(0, '--', 'LineWidth', lw)
    hold off
    title('y^e')
    xlim([0 t(end)])
    xlabel('Time [s]')
    ylabel('Magnitude [m]')
    legend('Ideal','OL','CL','Location','southeast')
    subplot(2, 1, 1)
    plot(t, OL.HF_helixCenter_filtered(filter:end, 2),'Color',color0,'LineWidth', lw)
    hold on
    plot(t, OLYaw.HF_helixCenter_filtered(filter:end, 2), 'Color',color1,'LineWidth', lw)
    plot(t, CLYaw.HF_helixCenter_filtered(filter:end, 2), 'Color',color2,'LineWidth', lw)
    yline(0, '--', 'LineWidth', lw)
    hold off
    title('z^e')
    xlim([0 t(end)])
    xlabel('Time [s]')
    ylabel('Magnitude [m]')
    legend('Ideal','OL','CL','Location','southeast')
    setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)
end