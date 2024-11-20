%% Data Analysis for Experiments
clear
close all 
addpath('.\Functions');
%clc

%% Data file (Chage this accordingly)
turbineName = '.\Data\NREL5MW\';
caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\2Turbines\';

% Different case
windCase = 'Uniform'; % Uniform, Shear, Turb, Both, Shear2
overallOption = 'N';
overallDetailOption = 'Y';
trajOption = 'N';
videoOption = 'N';
powerAnalysis = 'N';
DELAnalysis = 'N';
PBDAnalysis = 'N';

% Basic Settings
D_NREL5MW = 126;
U_inflow = 10;
timeStep = 0.1;
filter = 3500;
DeadtimeDelay = 112; % change to 112 when showing whole process

if strcmp(windCase, 'Uniform')
    basefile = '2Turbines_Baseline.mat';
    OLfileName = '2Turbines_OL_Helix_mag3.mat';
    CLfileName = '2Turbines_CL_Helix_ramp&stop_mag3.mat';
elseif strcmp(windCase, 'Shear')
    basefile = '2Turbines_Baseline_Shear0.2.mat';
    OLfileName = '2Turbines_OL_Helix_Shear0.2_mag3.mat';
    CLfileName = '2Turbines_CL_Helix_Shear0.2_mag3.mat';
elseif strcmp(windCase, 'Turb')
    basefile = '2Turbines_Baseline_TI6.mat';
    OLfileName = '2Turbines_OL_Helix_TI6_mag3.mat';
    CLfileName = '2Turbines_CL_Helix_TI6_mag3.mat';    
elseif strcmp(windCase, 'Both')
    basefile = '2Turbines_Baseline_TI6&Shear0.2.mat';
    OLfileName = '2Turbines_OL_Helix_TI6&Shear0.2_mag3.mat';
    CLfileName = '2Turbines_CL_Helix_TI6&Shear0.2_mag3.mat';   
elseif strcmp(windCase, 'Shear2')
    basefile = '2Turbines_Baseline_Shear0.2Cheat.mat';
    OLfileName = '2Turbines_OL_Helix_Shear0.2_mag3.mat';
    CLfileName = '2Turbines_CL_Helix_Shear0.2_mag3.mat';    
end

Baseline = load([turbineName caseName basefile]);
OL = load([turbineName caseName OLfileName]);
CL = load([turbineName caseName CLfileName]);

%% Visualization
simLength = length(Baseline.Cp_store);
t = (1:(simLength-filter+1)) * timeStep;
lw = 2;
Font = 20;
color0 = [0.4660 0.6740 0.1880];
color1 = [0, 0.4470, 0.7410];
color2 = [0.8500, 0.3250, 0.0980];

%% ============== Overall Result (Input & Output)
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
    plot(t, delayseq(CL.r(filter:end, 1), DeadtimeDelay), 'k:', 'LineWidth', lw)
    plot(t, delayseq(CL.r(filter:end, 2), DeadtimeDelay), 'k:', 'LineWidth', lw)
    yline(0, '--', 'LineWidth', lw)
    hold off
    title('Output: Helix Frame')
    xlim([0 t(end)])
    xlabel('Time [s]')
    ylim([-1 20])
    ylabel('Magnitude [m]')
    legend('z^e_b','y^e_b','z^e','y^e','r_z','r_y','Location','southeastoutside')
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
    plot(t, CL.u(filter:end, 1), 'Color',color1, 'LineWidth', lw)
    plot(t, CL.u(filter:end, 2), 'Color',color2, 'LineWidth', lw)
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

%% ============== Overeall Detailed Visualization
if strcmp(overallDetailOption, 'Y')
    % Control Input
    figure('Name', 'Input Detail', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(2, 1, 1)
    plot(t, OL.HF_beta(filter:end, 1), '--','Color',color1,'LineWidth', lw)
    hold on
    plot(t, CL.u(filter:end, 1), 'Color',color1, 'LineWidth', lw)
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
    plot(t, CL.u(filter:end, 2),'Color',color2, 'LineWidth', lw)
    hold off
    title('\beta^e_{yaw}')
    xlim([0 t(end)])
    xlabel('Time [s]')
    ylim([-1 5])
    ylabel('Magnitude [deg]')
    legend('OL','CL','Location','southeast')
%     setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)

    % hub Jet
    figure('Name', 'Output Detail', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(2, 1, 1)
    plot(t, OL.HF_helixCenter_filtered(filter:end, 1), '--','Color',color1,'LineWidth', lw)
    hold on
    plot(t, CL.HF_helixCenter_filtered(filter:end, 1), 'Color',color1,'LineWidth', lw)
    plot(t, delayseq(CL.r(filter:end, 1), 0), 'k:', 'LineWidth', lw)
    yline(0, '--', 'LineWidth', lw)
    hold off
    title('z^e')
    xlim([0 t(end)])
    xlabel('Time [s]')
    ylim([-1 15])
    ylabel('Magnitude [m]')
    legend('OL','CL','r','Location','southeast')
    subplot(2, 1, 2)
    plot(t, OL.HF_helixCenter_filtered(filter:end, 2), '--','Color',color2,'LineWidth', lw)
    hold on
    plot(t, CL.HF_helixCenter_filtered(filter:end, 2), 'Color',color2,'LineWidth', lw)
    plot(t, delayseq(CL.r(filter:end, 2), 0), 'k:', 'LineWidth', lw)
    yline(0, '--', 'LineWidth', lw)
    hold off
    title('y^e')
    xlim([0 t(end)])
    xlabel('Time [s]')
    ylim([-1 15])
    ylabel('Position [m]')
    legend('OL','CL','r','Location','southeast')
%     setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)
end

%% ============== Hub Jet Trajectory
if strcmp(trajOption, 'Y')
    center_bl = mean(Baseline.FF_helixCenter_filtered(3000:end, :));
    center_ol = mean(OL.FF_helixCenter_filtered(3000:end, :));
    center_cl = mean(CL.FF_helixCenter_filtered(3000:end, :));
    figure('Name', 'HubJet Trajectory', 'NumberTitle', 'off', 'Position', [100, 100, 600, 600]);
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
%     setfigpaper('Width',[15,1],'Interpreter','tex','FontSize',Font,'linewidth',lw)
end

%% Video comparison
if strcmp(videoOption, 'Y')
%     ringVisualization2(Baseline, D_NREL5MW)
    videoCompare_func(OL,CL,D_NREL5MW,'')
end

%% Power and Fatigue Analysis
% Note!!! Open-loop is used as benchmark rather than baseline
% ============== Baseline
% Upstream WT1
BL_result.WT1.power = calculatePower(filter,Baseline.Cp_store,D_NREL5MW,U_inflow); % [MW]
BL_result.WT1.DEL = calculateDEL(filter, ...
    Baseline.Mflap1_store,Baseline.Medge1_store, ...
    Baseline.Mflap2_store,Baseline.Medge2_store, ...
    Baseline.Mflap3_store,Baseline.Medge3_store, ...
    timeStep); % [Nm]
BL_result.WT1.PBD = calculatePBD(filter,Baseline.PitchAngles, ...
    Baseline.Mflap1_store,Baseline.Medge1_store, ...
    Baseline.Mflap2_store,Baseline.Medge2_store, ...
    Baseline.Mflap3_store,Baseline.Medge3_store); % [kNm deg]
% Downstream WT2
BL_result.WT2.power = calculatePower(filter,Baseline.Cpturb2_store,D_NREL5MW,U_inflow); % [MW]
BL_result.WT2.DEL = calculateDEL(filter, ...
    Baseline.Mflap1turb2_store,Baseline.Medge1turb2_store, ...
    Baseline.Mflap2turb2_store,Baseline.Medge2turb2_store, ...
    Baseline.Mflap3turb2_store,Baseline.Medge3turb2_store, ...
    timeStep); % [Nm]
BL_result.WT2.PBD = calculatePBD(filter,Baseline.PitchAnglesturb2, ...
    Baseline.Mflap1turb2_store,Baseline.Medge1turb2_store, ...
    Baseline.Mflap2turb2_store,Baseline.Medge2turb2_store, ...
    Baseline.Mflap3turb2_store,Baseline.Medge3turb2_store); % [kNm deg]
BL_result.WT1.DEL_flapwise = BL_result.WT1.DEL(1)+BL_result.WT1.DEL(3)+BL_result.WT1.DEL(5);
BL_result.WT1.DEL_edgewise = BL_result.WT1.DEL(2)+BL_result.WT1.DEL(4)+BL_result.WT1.DEL(6);
BL_result.WT2.DEL_flapwise = BL_result.WT2.DEL(1)+BL_result.WT2.DEL(3)+BL_result.WT2.DEL(5);
BL_result.WT2.DEL_edgewise = BL_result.WT2.DEL(2)+BL_result.WT2.DEL(4)+BL_result.WT2.DEL(6);
BL_result.All.power = BL_result.WT1.power + BL_result.WT2.power;
BL_result.All.DEL_flapwise = BL_result.WT1.DEL_flapwise + BL_result.WT2.DEL_flapwise;
BL_result.All.DEL_edgewise = BL_result.WT1.DEL_edgewise + BL_result.WT2.DEL_edgewise;
BL_result.All.PBD = BL_result.WT1.PBD + BL_result.WT2.PBD;

% ============== Open-Loop
OL_result.WT1.power = calculatePower(filter,OL.Cp_store,D_NREL5MW,U_inflow); % [MW]
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
OL_result.WT2.power = calculatePower(filter,OL.Cpturb2_store,D_NREL5MW,U_inflow); % [MW]
OL_result.WT2.DEL = calculateDEL(filter, ...
    OL.Mflap1turb2_store,OL.Medge1turb2_store, ...
    OL.Mflap2turb2_store,OL.Medge2turb2_store, ...
    OL.Mflap3turb2_store,OL.Medge3turb2_store, ...
    timeStep); % [Nm]
OL_result.WT2.PBD = calculatePBD(filter,OL.PitchAnglesturb2, ...
    OL.Mflap1turb2_store,OL.Medge1turb2_store, ...
    OL.Mflap2turb2_store,OL.Medge2turb2_store, ...
    OL.Mflap3turb2_store,OL.Medge3turb2_store); % [kNm deg]
OL_result.WT1.DEL_flapwise = OL_result.WT1.DEL(1)+OL_result.WT1.DEL(3)+OL_result.WT1.DEL(5);
OL_result.WT1.DEL_edgewise = OL_result.WT1.DEL(2)+OL_result.WT1.DEL(4)+OL_result.WT1.DEL(6);
OL_result.WT2.DEL_flapwise = OL_result.WT2.DEL(1)+OL_result.WT2.DEL(3)+OL_result.WT2.DEL(5);
OL_result.WT2.DEL_edgewise = OL_result.WT2.DEL(2)+OL_result.WT2.DEL(4)+OL_result.WT2.DEL(6);
OL_result.All.power = OL_result.WT1.power + OL_result.WT2.power;
OL_result.All.DEL_flapwise = OL_result.WT1.DEL_flapwise + OL_result.WT2.DEL_flapwise;
OL_result.All.DEL_edgewise = OL_result.WT1.DEL_edgewise + OL_result.WT2.DEL_edgewise;
OL_result.All.PBD = OL_result.WT1.PBD + OL_result.WT2.PBD;

% ============== Closed-Loop
CL_result.WT1.power = calculatePower(filter,CL.Cp_store,D_NREL5MW,U_inflow); % [MW]
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
CL_result.WT2.power = calculatePower(filter,CL.Cpturb2_store,D_NREL5MW,U_inflow); % [MW]
CL_result.WT2.DEL = calculateDEL(filter, ...
    CL.Mflap1turb2_store,CL.Medge1turb2_store, ...
    CL.Mflap2turb2_store,CL.Medge2turb2_store, ...
    CL.Mflap3turb2_store,CL.Medge3turb2_store, ...
    timeStep); % [Nm]
CL_result.WT2.PBD = calculatePBD(filter,CL.PitchAnglesturb2, ...
    CL.Mflap1turb2_store,CL.Medge1turb2_store, ...
    CL.Mflap2turb2_store,CL.Medge2turb2_store, ...
    CL.Mflap3turb2_store,CL.Medge3turb2_store); % [kNm deg]
CL_result.WT1.DEL_flapwise = CL_result.WT1.DEL(1)+CL_result.WT1.DEL(3)+CL_result.WT1.DEL(5);
CL_result.WT1.DEL_edgewise = CL_result.WT1.DEL(2)+CL_result.WT1.DEL(4)+CL_result.WT1.DEL(6);
CL_result.WT2.DEL_flapwise = CL_result.WT2.DEL(1)+CL_result.WT2.DEL(3)+CL_result.WT2.DEL(5);
CL_result.WT2.DEL_edgewise = CL_result.WT2.DEL(2)+CL_result.WT2.DEL(4)+CL_result.WT2.DEL(6);
CL_result.All.power = CL_result.WT1.power + CL_result.WT2.power;
CL_result.All.DEL_flapwise = CL_result.WT1.DEL_flapwise + CL_result.WT2.DEL_flapwise;
CL_result.All.DEL_edgewise = CL_result.WT1.DEL_edgewise + CL_result.WT2.DEL_edgewise;
CL_result.All.PBD = CL_result.WT1.PBD + CL_result.WT2.PBD;

% === Power 
if strcmp(powerAnalysis, 'Y')
    x = [1 2 3];
    deltaPower = [(OL_result.WT1.power-BL_result.WT1.power)/(BL_result.WT1.power) (OL_result.WT2.power-BL_result.WT2.power)/(BL_result.WT2.power) (OL_result.All.power-BL_result.All.power)/(BL_result.All.power);
                  (CL_result.WT1.power-BL_result.WT1.power)/(BL_result.WT1.power) (CL_result.WT2.power-BL_result.WT2.power)/(BL_result.WT2.power) (CL_result.All.power-BL_result.All.power)/(BL_result.All.power)];
    figure('Name', 'Power', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    bar(x,deltaPower*100); % convert to 100%
    xticks(x); 
    xticklabels({'WT1', 'WT2', 'WT1+WT2'}); 
    ylabel('\Delta Power [%]')
    legend('OL', 'CL', 'Location','northeast')
%     setfigpaper('Width',[20,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)
end

% === DEL 
% if strcmp(DELAnalysis, 'Y')
%     x = [1 2 3];
%     deltaDEL = [(OL_result.WT1.DEL_flapwise-BL_result.WT1.DEL_flapwise)/(BL_result.WT1.DEL_flapwise) (OL_result.WT2.DEL_flapwise-BL_result.WT2.DEL_flapwise)/(BL_result.WT2.DEL_flapwise) (OL_result.All.DEL_flapwise-BL_result.All.DEL_flapwise)/(BL_result.All.DEL_flapwise);
%                 (OL_result.WT1.DEL_edgewise-BL_result.WT1.DEL_edgewise)/(BL_result.WT1.DEL_edgewise) (OL_result.WT2.DEL_edgewise-BL_result.WT2.DEL_edgewise)/(BL_result.WT2.DEL_edgewise) (OL_result.All.DEL_edgewise-BL_result.All.DEL_edgewise)/(BL_result.All.DEL_edgewise);
%                 (CL_result.WT1.DEL_flapwise-BL_result.WT1.DEL_flapwise)/(BL_result.WT1.DEL_flapwise) (CL_result.WT2.DEL_flapwise-BL_result.WT2.DEL_flapwise)/(BL_result.WT2.DEL_flapwise) (CL_result.All.DEL_flapwise-BL_result.All.DEL_flapwise)/(BL_result.All.DEL_flapwise);
%                 (CL_result.WT1.DEL_edgewise-BL_result.WT1.DEL_edgewise)/(BL_result.WT1.DEL_edgewise) (CL_result.WT2.DEL_edgewise-BL_result.WT2.DEL_edgewise)/(BL_result.WT2.DEL_edgewise) (CL_result.All.DEL_edgewise-BL_result.All.DEL_edgewise)/(BL_result.All.DEL_edgewise)];
%     figure('Name', 'DEL', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
%     bar(x,deltaDEL*100);
%     xticks(x); 
%     xticklabels({'WT1', 'WT2', 'WT1+WT2'}); 
%     ylabel('\Delta DEL [%]')
% %     ylim([-0.5 3])
%     legend('OL Flapwise', 'OL Edgewise','CL Flapwise', 'CL Edgewise', 'Location','northeast')
%     setfigpaper('Width',[20,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)
% end
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

% === PBD
% PBD is only compared between the open-loop and closed-loop because
% baseline does not have PBD since Helix is not activated
if strcmp(PBDAnalysis, 'Y')
    deltaPBD = mean(CL_result.WT1.PBD) - mean(OL_result.WT1.PBD);
    disp(deltaPBD)
end