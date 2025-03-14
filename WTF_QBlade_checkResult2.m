%% Data Analysis for Experiments
clear
close all 
addpath('.\Functions');
%clc

%% Data file (Chage this accordingly)
turbineName = '.\Data\NREL5MW\';
caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\2TurbinesNew\';

% Different case
windCase = 'Uni&Turb'; % Uni&Shear, Uni&Turb
if strcmp(windCase, 'Uni&Shear')
    OLfileNameUni = '2Turbines_OL_Helix_mag3_4D.mat';
    OLfileNameShear = '2Turbines_OL_Helix_Shear0.2_mag3_4D.mat';
    OL = load([turbineName caseName OLfileNameUni]);
    CL = load([turbineName caseName OLfileNameShear]);
elseif strcmp(windCase, 'Uni&Turb')
    OLfileNameUni = '2Turbines_OL_Helix_mag3_4D.mat';
    OLfileNameTurb = '2Turbines_OL_Helix_TI6_mag3_4D.mat';
    OL = load([turbineName caseName OLfileNameUni]);
    CL = load([turbineName caseName OLfileNameTurb]);
end

%% Overall Settings
overallOption = 'N';
overallDetailOption = 'Y';
trajOption = 'Y';
videoOption = 'N';
powerAnalysis = 'N';
DELAnalysis = 'N';
PBDAnalysis = 'Y';
powerDELAnalysis = 'Y';
windAnalysis = 'Y';

% Basic Settings
D_NREL5MW = 126;
U_inflow = 10;
timeStep = 0.1;
filter = 4000;
DeadtimeDelay = 112; % change to 112 when showing whole process

% Visualization
simLength = length(OL.Power_store);
t = (1:(simLength-filter+1)) * timeStep;
lw = 2;
Font = 20;
color0 = [0.4660 0.6740 0.1880];
color1 = [0, 0.4470, 0.7410];
color2 = [0.8500, 0.3250, 0.0980];

% ============== Overall Result (Input & Output)
% Hub Jet
if strcmp(overallOption, 'Y')
    filter = 1000;
    t = (1:(simLength-filter+1)) * timeStep;
    figure('Name', 'Experiment Input-Output', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(2, 2, 2)
    plot(t, OL.HF_helixCenter_filtered(filter:end, 1), '--','Color',color1,'LineWidth', lw)
    hold on
    plot(t, OL.HF_helixCenter_filtered(filter:end, 2), '--', 'Color',color2,'LineWidth', lw)
    plot(t, CL.HF_helixCenter_filtered(filter:end, 1), 'Color',color1, 'LineWidth', lw)
    plot(t, CL.HF_helixCenter_filtered(filter:end, 2), 'Color',color2, 'LineWidth', lw)
    yline(0, '--', 'LineWidth', lw)
    hold off
    title('Output: Helix Frame')
    xlim([0 t(end)])
    xlabel('Time [s]')
    ylim([-1 20])
    ylabel('Magnitude [m]')
    legend('z^e_b','y^e_b','z^e','y^e','Location','southeastoutside')
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
    legend('z^e_b','y^e_b','z^e','y^e','Location','southeastoutside')
%     setfigpaper('Width',[40,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)
    
    % Control Input
    % figure('Name', 'Experiment Result: Hub Jet', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(2, 2, 1)
    plot(t, OL.HF_beta(filter:end, 1), '--','Color',color1, 'LineWidth', lw)
    hold on
    plot(t, OL.HF_beta(filter:end, 2), '--','Color',color2, 'LineWidth', lw)
    plot(t, CL.HF_beta(filter:end, 1), 'Color',color1, 'LineWidth', lw)
    plot(t, CL.HF_beta(filter:end, 2), 'Color',color2, 'LineWidth', lw)
    hold off
    title('Input: Helix Frame')
    xlim([0 t(end)])
    xlabel('Time [s]')
    ylim([-1 10])
    ylabel('Magnitude [deg]')
    legend('\beta^e_{tilt,b}','\beta^e_{yaw,b}','\beta^e_{tilt}','\beta^e_{yaw}','Location','southeastoutside')
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
    legend('\beta_{tilt,b}','\beta_{yaw,b}','\beta_{tilt}','\beta_{yaw}','Location','southeastoutside')
%     setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)
    
    filter = 3500;
    t = (1:(simLength-filter+1)) * timeStep;
end

% ============== Wind Field Visualization
if strcmp(windAnalysis, 'Y')
    [meanU1, TI1] = calculateTI(OL);
    [meanU2, TI2] = calculateTI(CL);
    figure('Name', 'Wind', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(2, 1, 1)
    plot(t, meanU1(filter:end), 'Color',color1,'LineWidth', lw)
    hold on
    plot(t, meanU2(filter:end),'Color',color2,'LineWidth', lw)
    hold off
    title('meanU')
    xlim([0 t(end)])
    xlabel('Time [s]')
%     ylim([-1 5])
    ylabel('Speed [m/s]')
    legend('Uni','Case','Location','southeast')

    subplot(2, 1, 2)
    plot(t, TI1(filter:end),'Color',color1,'LineWidth', lw)
    hold on
    plot(t, TI2(filter:end),'Color',color2,'LineWidth', lw)
    hold off
    title('TI')
    xlim([0 t(end)])
    xlabel('Time [s]')
%     ylim([-1 5])
    ylabel('Percentage')
    legend('Uni','Case','Location','southeast')    
end

% ============== Overeall Detailed Visualization
if strcmp(overallDetailOption, 'Y')
    % Input Hub Jet
    figure('Name', 'Input HF', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(2, 1, 1)
    plot(t, OL.HF_beta(filter:end, 1), '--','Color',color1,'LineWidth', lw)
    hold on
    plot(t, CL.HF_beta(filter:end, 1), 'Color',color1, 'LineWidth', lw)
    hold off
    title('\beta^e_{tilt}')
    xlim([0 t(end)])
    xlabel('Time [s]')
    ylim([-1 5])
    ylabel('Magnitude [deg]')
    legend('OL','CL','Location','southeast')
    subplot(2, 1, 2)
    hold on
    plot(t, OL.HF_beta(filter:end, 2), '--','Color',color2,'LineWidth', lw)
    plot(t, CL.HF_beta(filter:end, 2),'Color',color2, 'LineWidth', lw)
    hold off
    title('\beta^e_{yaw}')
    xlim([0 t(end)])
    xlabel('Time [s]')
    ylim([-1 5])
    ylabel('Magnitude [deg]')
    legend('Uni', 'Case','Location','southeast')
%     setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)

    % Input Fixed Frame
    figure('Name', 'Input FF', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(2, 1, 1)
    plot(t, OL.FF_beta(filter:end, 1), '--','Color',color1,'LineWidth', lw)
    hold on
    plot(t, CL.FF_beta(filter:end, 1), 'Color',color1, 'LineWidth', lw)
    hold off
    title('\beta_{tilt}')
    xlim([0 t(end)])
    xlabel('Time [s]')
%     ylim([-1 5])
    ylabel('Magnitude [deg]')
    legend('OL','CL','Location','southeast')
    subplot(2, 1, 2)
    hold on
    plot(t, OL.FF_beta(filter:end, 2), '--','Color',color2,'LineWidth', lw)
    plot(t, CL.FF_beta(filter:end, 2),'Color',color2, 'LineWidth', lw)
    hold off
    title('\beta_{yaw}')
    xlim([0 t(end)])
    xlabel('Time [s]')
%     ylim([-1 5])
    ylabel('Magnitude [deg]')
    legend('Uni', 'Case','Location','southeast')
%     setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)

    % hub Jet Helix Frame
    figure('Name', 'Output Detail', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(2, 1, 1)
    plot(t, OL.HF_helixCenter_filtered(filter:end, 1), '--','Color',color1,'LineWidth', lw)
    hold on
    plot(t, CL.HF_helixCenter_filtered(filter:end, 1), 'Color',color1,'LineWidth', lw)
    yline(0, '--', 'LineWidth', lw)
    hold off
    title('z^e')
    xlim([0 t(end)])
    xlabel('Time [s]')
    ylim([-1 15])
    ylabel('Magnitude [m]')
    legend('OL','CL','Location','southeast')
    subplot(2, 1, 2)
    plot(t, OL.HF_helixCenter_filtered(filter:end, 2), '--','Color',color2,'LineWidth', lw)
    hold on
    plot(t, CL.HF_helixCenter_filtered(filter:end, 2), 'Color',color2,'LineWidth', lw)
    yline(0, '--', 'LineWidth', lw)
    hold off
    title('y^e')
    xlim([0 t(end)])
    xlabel('Time [s]')
    ylim([-1 15])
    ylabel('Position [m]')
    legend('Uni', 'Case','Location','southeast')
%     setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)
end

% ============== Hub Jet Trajectory
if strcmp(trajOption, 'Y')
    center_ol = mean(OL.FF_helixCenter_filtered(filter:end, :));
    center_cl = mean(CL.FF_helixCenter_filtered(filter:end, :));
    figure('Name', 'HubJet Trajectory', 'NumberTitle', 'off', 'Position', [100, 100, 600, 600]);
    plot(OL.FF_helixCenter_filtered(filter:end, 2), OL.FF_helixCenter_filtered(filter:end, 1), 'Color',color1, 'LineWidth', lw)
    hold on
    plot(CL.FF_helixCenter_filtered(filter:end, 2), CL.FF_helixCenter_filtered(filter:end, 1), 'Color',color2, 'LineWidth', lw)
    plot(center_ol(2), center_ol(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color1);
    plot(center_cl(2), center_cl(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color2);
    plot(0, 90, 'k*', 'MarkerSize', 10);
    hold off
    title('Hub Jet Trajectory')
    xlabel('y [m]')
    ylabel('z [m]')
    xlim([-30 20])
    ylim([67 117])
    legend('Uni', 'Case', 'Location','southeast')
%     setfigpaper('Width',[15,1],'Interpreter','tex','FontSize',Font,'linewidth',lw)
end

% ============== Video comparison
if strcmp(videoOption, 'Y')
%     ringVisualization2(Baseline, D_NREL5MW)
    videoCompare_func(OL,CL,D_NREL5MW,'')
end

% ============== Power and Fatigue Analysis
% Upstream WT1
% =========== Uniform
OL_result.WT1.power = calculatePower(filter,OL.Power_store,D_NREL5MW,U_inflow); % [MW]
OL_result.WT1.DEL = calculateDEL(filter, ...
    OL.Mflap1_store,OL.Medge1_store, ...
    OL.Mflap2_store,OL.Medge2_store, ...
    OL.Mflap3_store,OL.Medge3_store, ...
    timeStep); % [Nm]
OL_result.WT1.PBD = calculatePBD(filter,OL.PitchAngles, ...
    OL.Mflap1_store,OL.Medge1_store, ...
    OL.Mflap2_store,OL.Medge2_store, ...
    OL.Mflap3_store,OL.Medge3_store); % [kNm deg]
% Downstream WT2
OL_result.WT2.power = calculatePower(filter,OL.Powerturb2_store,D_NREL5MW,U_inflow); % [MW]
OL_result.WT2.DEL = calculateDEL(filter, ...
    OL.Mflap1turb2_store,OL.Medge1turb2_store, ...
    OL.Mflap2turb2_store,OL.Medge2turb2_store, ...
    OL.Mflap3turb2_store,OL.Medge3turb2_store, ...
    timeStep); % [Nm]
OL_result.WT2.PBD = calculatePBD(filter,OL.PitchAnglesturb2, ...
    OL.Mflap1turb2_store,OL.Medge1turb2_store, ...
    OL.Mflap2turb2_store,OL.Medge2turb2_store, ...
    OL.Mflap3turb2_store,OL.Medge3turb2_store); % [kNm deg]
OL_result.WT1.DEL_flapwise = mean(OL_result.WT1.DEL(1)+OL_result.WT1.DEL(3)+OL_result.WT1.DEL(5));
OL_result.WT1.DEL_edgewise = mean(OL_result.WT1.DEL(2)+OL_result.WT1.DEL(4)+OL_result.WT1.DEL(6));
OL_result.WT2.DEL_flapwise = mean(OL_result.WT2.DEL(1)+OL_result.WT2.DEL(3)+OL_result.WT2.DEL(5));
OL_result.WT2.DEL_edgewise = mean(OL_result.WT2.DEL(2)+OL_result.WT2.DEL(4)+OL_result.WT2.DEL(6));
OL_result.All.power = OL_result.WT1.power + OL_result.WT2.power;
OL_result.All.DEL_flapwise = OL_result.WT1.DEL_flapwise + OL_result.WT2.DEL_flapwise;
OL_result.All.DEL_edgewise = OL_result.WT1.DEL_edgewise + OL_result.WT2.DEL_edgewise;
OL_result.All.PBD = OL_result.WT1.PBD + OL_result.WT2.PBD;

% =========== Case
CL_result.WT1.power = calculatePower(filter,CL.Power_store,D_NREL5MW,U_inflow); % [MW]
CL_result.WT1.DEL = calculateDEL(filter, ...
    CL.Mflap1_store,CL.Medge1_store, ...
    CL.Mflap2_store,CL.Medge2_store, ...
    CL.Mflap3_store,CL.Medge3_store, ...
    timeStep); % [Nm]
CL_result.WT1.PBD = calculatePBD(filter,CL.PitchAngles, ...
    CL.Mflap1_store,CL.Medge1_store, ...
    CL.Mflap2_store,CL.Medge2_store, ...
    CL.Mflap3_store,CL.Medge3_store); % [kNm deg]
% Downstream WT2
CL_result.WT2.power = calculatePower(filter,CL.Powerturb2_store,D_NREL5MW,U_inflow); % [MW]
CL_result.WT2.DEL = calculateDEL(filter, ...
    CL.Mflap1turb2_store,CL.Medge1turb2_store, ...
    CL.Mflap2turb2_store,CL.Medge2turb2_store, ...
    CL.Mflap3turb2_store,CL.Medge3turb2_store, ...
    timeStep); % [Nm]
CL_result.WT2.PBD = calculatePBD(filter,CL.PitchAnglesturb2, ...
    CL.Mflap1turb2_store,CL.Medge1turb2_store, ...
    CL.Mflap2turb2_store,CL.Medge2turb2_store, ...
    CL.Mflap3turb2_store,CL.Medge3turb2_store); % [kNm deg]
CL_result.WT1.DEL_flapwise = mean(CL_result.WT1.DEL(1)+CL_result.WT1.DEL(3)+CL_result.WT1.DEL(5));
CL_result.WT1.DEL_edgewise = mean(CL_result.WT1.DEL(2)+CL_result.WT1.DEL(4)+CL_result.WT1.DEL(6));
CL_result.WT2.DEL_flapwise = mean(CL_result.WT2.DEL(1)+CL_result.WT2.DEL(3)+CL_result.WT2.DEL(5));
CL_result.WT2.DEL_edgewise = mean(CL_result.WT2.DEL(2)+CL_result.WT2.DEL(4)+CL_result.WT2.DEL(6));
CL_result.All.power = CL_result.WT1.power + CL_result.WT2.power;
CL_result.All.DEL_flapwise = CL_result.WT1.DEL_flapwise + CL_result.WT2.DEL_flapwise;
CL_result.All.DEL_edgewise = CL_result.WT1.DEL_edgewise + CL_result.WT2.DEL_edgewise;
CL_result.All.PBD = CL_result.WT1.PBD + CL_result.WT2.PBD;

% === Power 
if strcmp(powerAnalysis, 'Y')
    x = [1 2 3];
    deltaPower = [(OL_result.WT1.power-OL_result.WT1.power)/(OL_result.WT1.power) (OL_result.WT2.power-OL_result.WT2.power)/(OL_result.WT2.power) (OL_result.All.power-OL_result.All.power)/(OL_result.All.power);
                  (CL_result.WT1.power-OL_result.WT1.power)/(OL_result.WT1.power) (CL_result.WT2.power-OL_result.WT2.power)/(OL_result.WT2.power) (CL_result.All.power-OL_result.All.power)/(OL_result.All.power)];
    figure('Name', 'Power', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    bar(x,deltaPower*100); % convert to 100%
    xticks(x); 
    xticklabels({'WT1', 'WT2', 'WT1+WT2'}); 
    ylabel('\Delta Power [%]')
    legend('OL', 'CL', 'Location','northeast')
%     setfigpaper('Width',[20,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)
end

% === DEL 
if strcmp(DELAnalysis, 'Y')
    x = [1 2 3];
    deltaDEL = [(CL_result.WT1.DEL_flapwise-OL_result.WT1.DEL_flapwise)/(OL_result.WT1.DEL_flapwise) (CL_result.WT2.DEL_flapwise-OL_result.WT2.DEL_flapwise)/(OL_result.WT2.DEL_flapwise) (CL_result.All.DEL_flapwise-OL_result.All.DEL_flapwise)/(OL_result.All.DEL_flapwise);
                (CL_result.WT1.DEL_edgewise-OL_result.WT1.DEL_edgewise)/(OL_result.WT1.DEL_edgewise) (CL_result.WT2.DEL_edgewise-OL_result.WT2.DEL_edgewise)/(OL_result.WT2.DEL_edgewise) (CL_result.All.DEL_edgewise-OL_result.All.DEL_edgewise)/(OL_result.All.DEL_edgewise)];
    figure('Name', 'DEL', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    bar(x,deltaDEL*100);
    xticks(x); 
    xticklabels({'WT1', 'WT2', 'WT1+WT2'}); 
    ylabel('\Delta DEL [%]')
    legend('BR flapwise', 'BR Edgewise', 'Location','southeast')
%     setfigpaper('Width',[20,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)
end

% === Power and DEL (Result Make Sense or Not)
if strcmp(powerDELAnalysis, 'Y')
    filter2 = 2000;
    x = [1 2 3];
    deltaPower = [(OL_result.WT1.power-OL_result.WT1.power)/(OL_result.WT1.power) (OL_result.WT2.power-OL_result.WT2.power)/(OL_result.WT2.power) (OL_result.All.power-OL_result.All.power)/(OL_result.All.power);
                  (CL_result.WT1.power-OL_result.WT1.power)/(OL_result.WT1.power) (CL_result.WT2.power-OL_result.WT2.power)/(OL_result.WT2.power) (CL_result.All.power-OL_result.All.power)/(OL_result.All.power)];
    deltaDELf = [(OL_result.WT1.DEL_flapwise-OL_result.WT1.DEL_flapwise)/(OL_result.WT1.DEL_flapwise) (OL_result.WT2.DEL_flapwise-OL_result.WT2.DEL_flapwise)/(OL_result.WT2.DEL_flapwise) (OL_result.All.DEL_flapwise-OL_result.All.DEL_flapwise)/(OL_result.All.DEL_flapwise);
                 (CL_result.WT1.DEL_flapwise-OL_result.WT1.DEL_flapwise)/(OL_result.WT1.DEL_flapwise) (CL_result.WT2.DEL_flapwise-OL_result.WT2.DEL_flapwise)/(OL_result.WT2.DEL_flapwise) (CL_result.All.DEL_flapwise-OL_result.All.DEL_flapwise)/(OL_result.All.DEL_flapwise)];
    deltaDELe = [(OL_result.WT1.DEL_edgewise-OL_result.WT1.DEL_edgewise)/(OL_result.WT1.DEL_edgewise) (OL_result.WT2.DEL_edgewise-OL_result.WT2.DEL_edgewise)/(OL_result.WT2.DEL_edgewise) (OL_result.All.DEL_edgewise-OL_result.All.DEL_edgewise)/(OL_result.All.DEL_edgewise);
                 (CL_result.WT1.DEL_edgewise-OL_result.WT1.DEL_edgewise)/(OL_result.WT1.DEL_edgewise) (CL_result.WT2.DEL_edgewise-OL_result.WT2.DEL_edgewise)/(OL_result.WT2.DEL_edgewise) (CL_result.All.DEL_edgewise-OL_result.All.DEL_edgewise)/(OL_result.All.DEL_edgewise)];
    figure('Name', 'Power & DEL', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(1, 3, 1)
    bar(x,deltaPower*100);
    xticks(x); 
    xticklabels({'WT1', 'WT2', 'WT1+WT2'}); 
    ylabel('\Delta Power [%]')
    legend('Uni', 'Case', 'Location','northeast')
    title('Power')

    subplot(1, 3, 2)
    bar(x,deltaDELf*100);
    xticks(x); 
    xticklabels({'WT1', 'WT2', 'WT1+WT2'}); 
    ylabel('\Delta DEL [%]')
    legend('Uni', 'Case', 'Location','northeast')
    title('DEL Flapwise')
    
    subplot(1, 3, 3)
    bar(x,deltaDELe*100);
    xticks(x); 
    xticklabels({'WT1', 'WT2', 'WT1+WT2'}); 
    ylabel('\Delta DEL [%]')
    legend('Uni', 'Case', 'Location','northeast')
    title('DEL Edgewise')
    setfigpaper('Width',[40,0.3],'Interpreter','tex','FontSize',Font,'linewidth',lw)
end

% === PBD
% PBD is only compared between the open-loop and closed-loop because
% baseline does not have PBD since Helix is not activated
if strcmp(PBDAnalysis, 'Y')
    deltaPBD = mean(CL_result.WT1.PBD) - mean(OL_result.WT1.PBD);
    disp(deltaPBD)
end