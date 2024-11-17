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
DeadtimeDelay = 112;

%% Visualization Option 
overallOption = 'N';
trajOption = 'Y';
errorOption = 'Y';
videoOption = 'N';

%% Visualization
t = (1:(simLength-filter+1)) * timeStep;
lw = 2;
color0 = [0.4660 0.6740 0.1880];
color1 = [0, 0.4470, 0.7410];
color2 = [0.8500, 0.3250, 0.0980];

%% ============== Overall Result (Input & Output)
% Hub Jet
if strcmp(overallOption, 'Y')
    figure('Name', 'Experiment Input-Output', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(2, 2, 2)
    plot(t, OL.HF_helixCenter_filtered(filter:end, 1), '--','Color',color1,'LineWidth', lw)
    hold on
    plot(t, OL.HF_helixCenter_filtered(filter:end, 2), '--', 'Color',color2,'LineWidth', lw)
    plot(t, CL.HF_helixCenter_filtered(filter:end, 1), 'Color',color1, 'LineWidth', lw)
    plot(t, CL.HF_helixCenter_filtered(filter:end, 2), 'Color',color2, 'LineWidth', lw)
    plot(t, delayseq(CL.r(filter:end, 1), DeadtimeDelay), 'k:', 'LineWidth', lw)
    plot(t, delayseq(CL.r(filter:end, 2), DeadtimeDelay), 'k:', 'LineWidth', lw)
    yline(0, '--', 'LineWidth', lw)
    hold off
    title('Output: Helix Frame')
    xlim([0 t(end)])
    xlabel('Time [s]')
    ylim([-5 15])
    ylabel('Magnitude [m]')
    legend('z^e_{ol}','y^e_{ol}','z^e_{cl}','y^e_{cl}','r_z','r_y','Location','southeastoutside')
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
    legend('z^e_{ol}','y^e_{ol}','z^e_{cl}','y^e_{cl}','Location','southeastoutside')
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
    legend('\beta^e_{tilt,ol}','\beta^e_{yaw,ol}','\beta^e_{tilt,cl}','\beta^e_{yaw,cl}','Location','southeastoutside')
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
    ylim([-15 15])
    ylabel('Magnitude [deg]')
    legend('\beta_{tilt,ol}','\beta_{yaw,ol}','\beta_{tilt,cl}','\beta_{yaw,cl}','Location','southeastoutside')
    setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',15,'linewidth',lw)
end

%% ============== Hub Jet Trajectory
if strcmp(trajOption, 'Y')
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
    xlim([-30 20])
    ylim([67 117])
    legend('Baseline', 'OL', 'CL', 'Location','southeast')
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