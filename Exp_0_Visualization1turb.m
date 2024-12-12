%% Data Analysis for Experiments
clear
close all 
addpath('.\Functions');
%clc

%% Data file (Chage this accordingly)
turbineName = '.\Data\NREL5MW\';
caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\1Turbine\';
basefile = '1Turbine_Baseline.mat';
OLfileName = '1Turbine_Baseline.mat';
CLfileName = 'SISO_PI_2chls_suc_ramp&stop_mag3.mat';
% SISO_PI_2chls_suc_ramp&stop_mag3
% SISO_PI_2chls_fail_ramp&stop_mag3

Baseline = load([turbineName caseName basefile]);
OL = load([turbineName caseName OLfileName]);
CL = load([turbineName caseName CLfileName]);

%% Basic Settings
D_NREL5MW = 126;
U_inflow = 10;
timeStep = 0.1;
filter0 = 1000;
filter1 = 500;
filter = 1000;
simLength = length(Baseline.Cp_store);
DeadtimeDelay = 112;

%% Visualization Option 
overallDetailOption = 'N';
coorFrame = 'HF';
trajOption = 'N';
errorOption = 'Y';
videoOption = 'N';

%% Visualization
t = (1:(simLength-filter+1)) * timeStep;
lw = 2;
Font = 20;
color0 = [0.4660 0.6740 0.1880];
color1 = [0, 0.4470, 0.7410];
color2 = [0.8500, 0.3250, 0.0980];

%% ============== Overeall Detailed Visualization
if strcmp(overallDetailOption, 'Y')
    t0 = (1:(simLength-filter0+1)) * timeStep;
    if strcmp(coorFrame,'HF')
        % Input
        figure('Name', 'Input HF', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
        subplot(2, 1, 2)
        plot(t0, OL.HF_beta(filter0:end, 1),'Color',color0,'LineWidth', lw)
        hold on
        plot(t0, CL.u(filter0:end, 1), 'Color',color1, 'LineWidth', lw)
        hold off
        title('\beta^e_{yaw}')
        xlim([0 t0(end)])
        xlabel('Time [s]')
        ylim([-1 5])
        ylabel('Magnitude [deg]')
        legend('OL','CL','Location','southeast')
        subplot(2, 1, 1)
        hold on
        plot(t0, OL.HF_beta(filter0:end, 2),'Color',color0,'LineWidth', lw)
        plot(t0, CL.u(filter0:end, 2),'Color',color2, 'LineWidth', lw)
        hold off
        title('\beta^e_{tilt}')
        xlim([0 t0(end)])
        xlabel('Time [s]')
        ylim([-1 5])
        ylabel('Magnitude [deg]')
        legend('OL','CL','Location','southeast')
        setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)
        
        % Output
        r1 = ones(simLength-filter0+1, 1)*8.5606;
        figure('Name', 'Output Detail HF', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
        subplot(2, 1, 2)
        plot(t0, OL.HF_helixCenter_filtered(filter0:end, 1),'Color',color0,'LineWidth', lw)
        hold on
        plot(t0, CL.HF_helixCenter_filtered(filter0:end, 1),'Color',color1,'LineWidth', lw)
        plot(t0, r1, 'k:', 'LineWidth', lw)
        yline(0, '--', 'LineWidth', lw)
        hold off
        title('y^e')
        xlim([0 t0(end)])
        xlabel('Time [s]')
        ylim([-10 20])
        ylabel('Position [m]')
        legend('OL','CL','r','Location','southeast')
        subplot(2, 1, 1)
        r2 = ones(simLength-filter0+1, 1)*9.0661;
        plot(t0, OL.HF_helixCenter_filtered(filter0:end, 2),'Color',color0,'LineWidth', lw)
        hold on
        plot(t0, CL.HF_helixCenter_filtered(filter0:end, 2),'Color',color2,'LineWidth', lw)
        plot(t0, r2, 'k:', 'LineWidth', lw)
        yline(0, '--', 'LineWidth', lw)
        hold off
        title('z^e')
        xlim([0 t0(end)])
        xlabel('Time [s]')
        ylim([-10 20])
        ylabel('Position [m]')
        legend('OL','CL','r','Location','southeast')
        setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)
    elseif strcmp(coorFrame, 'FF')
        % Input
        figure('Name', 'Input FF', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
        subplot(2, 1, 2)
        plot(t0, OL.FF_beta(filter0:end, 1),'Color',color0,'LineWidth', lw)
        hold on
        plot(t0, CL.FF_beta(filter0:end, 1), 'Color',color1, 'LineWidth', lw)
        hold off
        title('\beta_{yaw}')
        xlim([0 t0(end)])
        xlabel('Time [s]')
        ylabel('Magnitude [deg]')
        legend('OL','CL','Location','southeast')
        subplot(2, 1, 1)
        hold on
        plot(t0, OL.FF_beta(filter0:end, 2),'Color',color0,'LineWidth', lw)
        plot(t0, CL.FF_beta(filter0:end, 2),'Color',color2, 'LineWidth', lw)
        hold off
        title('\beta_{tilt}')
        xlim([0 t0(end)])
        xlabel('Time [s]')
        ylabel('Magnitude [deg]')
        legend('OL','CL','Location','southeast')
        setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)

        % Output
        figure('Name', 'Output Detail FF', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
        subplot(2, 1, 2)
        plot(t0, OL.FF_helixCenter_filtered(filter0:end, 1),'Color',color0,'LineWidth', lw)
        hold on
        plot(t0, CL.FF_helixCenter_filtered(filter0:end, 1),'Color',color1,'LineWidth', lw)
        hold off
        title('z')
        xlim([0 t0(end)])
        xlabel('Time [s]')
        ylabel('Position [m]')
        legend('OL','CL','Location','southeast')
        subplot(2, 1, 1)
        plot(t0, OL.FF_helixCenter_filtered(filter0:end, 2),'Color',color0,'LineWidth', lw)
        hold on
        plot(t0, CL.FF_helixCenter_filtered(filter0:end, 2),'Color',color2,'LineWidth', lw)
        hold off
        title('y')
        xlim([0 t0(end)])
        xlabel('Time [s]')
        ylabel('Position [m]')
        legend('OL','CL','Location','southeast')
        setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)
    end
end

%% ============== Hub Jet Trajectory
if strcmp(trajOption, 'Y')
    center_bl = mean(Baseline.FF_helixCenter_filtered(filter1:end, :));
    center_ol = mean(OL.FF_helixCenter_filtered(filter1:end, :));
    center_cl = mean(CL.FF_helixCenter_filtered(filter1:end, :));
    figure('Name', 'Experiment HubJet Trajectory', 'NumberTitle', 'off', 'Position', [100, 100, 600, 600]);
    plot(Baseline.FF_helixCenter_filtered(filter1:end, 2), Baseline.FF_helixCenter_filtered(filter1:end, 1), 'Color',color0, 'LineWidth', lw)
    hold on
    plot(OL.FF_helixCenter_filtered(filter1:end, 2), OL.FF_helixCenter_filtered(filter1:end, 1), 'Color',color1, 'LineWidth', lw)
    plot(CL.FF_helixCenter_filtered(filter1:end, 2), CL.FF_helixCenter_filtered(filter1:end, 1), 'Color',color2, 'LineWidth', lw)
    plot(0, 90, 'k*', 'MarkerSize', 10);
    plot(center_bl(2), center_bl(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color0);
    plot(center_ol(2), center_ol(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color1);
    plot(center_cl(2), center_cl(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color2);
    hold off
    title('Hub Jet Trajectory')
    xlabel('y [m]')
    ylabel('z [m]')
    xlim([-30 20])
    ylim([67 117])
    legend('Baseline', 'OL', 'CL', 'Hub', 'Location','southeast')
    setfigpaper('Width',[15,1],'Interpreter','tex','FontSize',15,'linewidth',lw)
end

%% Error comparison
if strcmp(errorOption, 'Y')
    figure('Name', 'Experiment Error', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(2, 1, 1)
    plot(t, (CL.r(filter:end, 1)-CL.y(filter:end, 1)), 'Color',color1,'LineWidth',lw)
    hold on
    plot(t, (CL.ym(filter:end, 1)-CL.ytilda(filter:end, 1)),'--', 'Color',color2,'LineWidth',lw)
    yline(0, 'k--','LineWidth',lw)
    xline(20, '--k', 'Activate CL Ctrl', 'LabelOrientation', 'horizontal', 'LineWidth', 1);
    hold off
    xlabel('Time [s]')
    ylabel('Magnitude [m]')
    xlim([0 500])
    ylim([-20 15])
    legend('e_p', 'e_f', 'Location','southeast')
    title('Error Comparsion z^e')
    subplot(2, 1, 2)
    plot(t, (CL.r(filter:end, 2)-CL.y(filter:end, 2)), 'Color',color1,'LineWidth',lw)
    hold on
    plot(t, (CL.ym(filter:end, 2)-CL.ytilda(filter:end, 2)),'--', 'Color',color2,'LineWidth',lw)
    yline(0, 'k--','LineWidth',lw)
    xline(20, '--k', 'Activate CL Ctrl', 'LabelOrientation', 'horizontal', 'LineWidth', 1);
    hold off
    xlabel('Time [s]')
    ylabel('Magnitude [m]')
    xlim([0 500])
    ylim([-20 15])
    legend('e_p', 'e_f', 'Location','southeast')
    title('Error Comparsion y^e')
    setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',15,'linewidth',lw)
end

%% Video comparison
if strcmp(videoOption, 'Y')
    ringVisualization2(Baseline, D_NREL5MW)
    videoCompare_func(Baseline,CL,D_NREL5MW,'')
end