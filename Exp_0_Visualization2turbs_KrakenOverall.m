%% Data Analysis for Experiments
% clear
close all 
addpath('.\Functions');
%clc

%% Data file (Chage this accordingly)
turbineName = '.\Data\NREL5MW\';
caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\2TurbinesLonger\';

ControlOption = 'Kraken'; % Kraken, Helix
% Different case
if strcmp(ControlOption, 'Helix')
    basefile = '2Turbines_OL_Helix_mag3_4D.mat';
    OLfileName = '2Turbines_OL_Helix_ShearReCenter_mag3_4D.mat';
    CLfileName = '2Turbines_OL_Helix_TI6_mag3_4D.mat';
    FKfileName = '2Turbines_OL_Helix_TI6&Shear0.2_mag3_4D.mat';
elseif strcmp(ControlOption, 'Kraken')
    basefile = '2Turbines_OL_Helix_mag3_4D.mat';
    OLfileName = '2Turbines_CL_Helix_ShearReCenter3_mag3_4D.mat';
    CLfileName = '2Turbines_CL_Helix_TI6_mag3_4D.mat';
    FKfileName = '2Turbines_CL_Helix_BothReCenter2_mag3_4D.mat';
end 

Baseline = load([turbineName caseName basefile]);
OL = load([turbineName caseName OLfileName]);
CL = load([turbineName caseName CLfileName]);
FL = load([turbineName caseName FKfileName]);

%% Overall Settings
overallOption = 'N';
flowAnalysis = 'N';
rareDataAnalysis = 'N';
overallDetailOption = 'Y';
coorFrame = 'HF';
trajOption = 'Y';
videoOption = 'N';
powerAnalysis = 'N';
DELAnalysis = 'N';
PBDAnalysis = 'N';
powerDELAnalysis = 'Y';
storyTellingBasic = 'Y';

% Basic Settings
D_NREL5MW = 126;
U_inflow = 10;
timeStep = 0.1;
filter = 3000;
filter0 = 5000;    % Overall
filter2 = 1;    % Wind info
filter3 = 1;    % Rare data
DeadtimeDelay = 112; % change to 112 when showing whole process
Str = 0.3;                          % Strouhal number              
Freq = Str*U_inflow/D_NREL5MW;      % From Str, in Hz

%% Visualization
simLength = length(Baseline.Power_store);
t = (1:(simLength-filter+1)) * timeStep;
lw = 2;
Font = 20;
color0 = [0.4660 0.6740 0.1880];
color1 = [0, 0.4470, 0.7410];
color2 = [0.8500, 0.3250, 0.0980];
color3 = [0.3010 0.7450 0.9330];

% ============== Rare Data Visualization
% Focus on the second wind turbine
if strcmp(rareDataAnalysis, 'Y')
    figure('Name', 'Rare Data --- Power', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    t3 = (1:(simLength-filter3+1)) * timeStep;
    % Power
    subplot(2, 1, 1)
    plot(t3, OL.Powerturb2_store(filter3:end)/1e3,'Color',color1,'LineWidth', lw)
    hold on
    plot(t3, CL.Powerturb2_store(filter3:end)/1e3,'Color',color2,'LineWidth', lw)
    plot(t3, Baseline.Powerturb2_store(filter3:end)/1e3,'Color',color0,'LineWidth', lw)
    hold off
    legend('Shear','Turbulence','Uniform')
    xlim([0 t3(end)])
    xlabel('Time [s]')
    ylabel('Power [MW]')
    title('Power Production')
    subplot(2, 1, 2)
    plot(t3, OL.TorqueStoreTurb2(filter3:end),'Color',color1,'LineWidth', lw)
    hold on
    plot(t3, CL.TorqueStoreTurb2(filter3:end),'Color',color2,'LineWidth', lw)
    plot(t3, Baseline.TorqueStoreTurb2(filter3:end),'Color',color0,'LineWidth', lw)
    hold off
    legend('Shear','Turbulence','Uniform')
    xlim([0 t3(end)])
    xlabel('Time [s]')
    ylabel('Torque [Nm]')
    title('Generator Torque')
    
    % Fatigue WT2 (Time Domain)
    figure('Name', 'Rare Data --- Fatigue', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(3, 2, 1)
    plot(t3, OL.Moop1turb2_store(filter3:end),'Color',color1,'LineWidth', lw)
    hold on
    plot(t3, CL.Moop1turb2_store(filter3:end),'Color',color2,'LineWidth', lw)
    plot(t3, Baseline.Moop1turb2_store(filter3:end),'Color',color0,'LineWidth', lw)
    hold off
    legend('Shear','Turbulence','Uniform')
    xlim([0 t3(end)])
    xlabel('Time [s]')
    ylabel('Moop - Blade1 [NM]')
    subplot(3, 2, 3)
    plot(t3, OL.Moop2turb2_store(filter3:end),'Color',color1,'LineWidth', lw)
    hold on
    plot(t3, CL.Moop2turb2_store(filter3:end),'Color',color2,'LineWidth', lw)
    plot(t3, Baseline.Moop2turb2_store(filter3:end),'Color',color0,'LineWidth', lw)
    hold off
    legend('Shear','Turbulence','Uniform')
    xlim([0 t3(end)])
    xlabel('Time [s]')
    ylabel('Moop - Blade2 [NM]')
    subplot(3, 2, 5)
    plot(t3, OL.Moop3turb2_store(filter3:end),'Color',color1,'LineWidth', lw)
    hold on
    plot(t3, CL.Moop3turb2_store(filter3:end),'Color',color2,'LineWidth', lw)
    plot(t3, Baseline.Moop3turb2_store(filter3:end),'Color',color0,'LineWidth', lw)
    hold off
    legend('Shear','Turbulence','Uniform')
    xlim([0 t3(end)])
    xlabel('Time [s]')
    ylabel('Moop - Blade3 [NM]')

    subplot(3, 2, 2)
    plot(t3, OL.Mip1turb2_store(filter3:end),'Color',color1,'LineWidth', lw)
    hold on
    plot(t3, CL.Mip1turb2_store(filter3:end),'Color',color2,'LineWidth', lw)
    plot(t3, Baseline.Mip1turb2_store(filter3:end),'Color',color0,'LineWidth', lw)
    hold off
    legend('Shear','Turbulence','Uniform')
    xlim([0 t3(end)])
    xlabel('Time [s]')
    ylabel('Mip- Blade1 [NM]')
    subplot(3, 2, 4)
    plot(t3, OL.Mip2turb2_store(filter3:end),'Color',color1,'LineWidth', lw)
    hold on
    plot(t3, CL.Mip2turb2_store(filter3:end),'Color',color2,'LineWidth', lw)
    plot(t3, Baseline.Mip2turb2_store(filter3:end),'Color',color0,'LineWidth', lw)
    hold off
    legend('Shear','Turbulence','Uniform')
    xlim([0 t3(end)])
    xlabel('Time [s]')
    ylabel('Mip - Blade2 [NM]')
    subplot(3, 2, 6)
    plot(t3, OL.Mip3turb2_store(filter3:end),'Color',color1,'LineWidth', lw)
    hold on
    plot(t3, CL.Mip3turb2_store(filter3:end),'Color',color2,'LineWidth', lw)
    plot(t3, Baseline.Mip3turb2_store(filter3:end),'Color',color0,'LineWidth', lw)
    hold off
    legend('Shear','Turbulence','Uniform')
    xlim([0 t3(end)])
    xlabel('Time [s]')
    ylabel('Mip - Blade3 [NM]')

    % Fatigue (PSD)
    figure('Name', 'Rare Data --- Fatigue PSD', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    [MoopM1,Foop1] = pwelch(OL.Moop1turb2_store(filter3:end),[],[],[],1/timeStep);
    [MoopM2,Foop2] = pwelch(OL.Moop2turb2_store(filter3:end),[],[],[],1/timeStep);
    [MoopM3,Foop3] = pwelch(OL.Moop3turb2_store(filter3:end),[],[],[],1/timeStep);
    [MipM1,Fip1] = pwelch(OL.Mip1turb2_store(filter3:end),[],[],[],1/timeStep);
    [MipM2,Fip2] = pwelch(OL.Mip2turb2_store(filter3:end),[],[],[],1/timeStep);
    [MipM3,Fip3] = pwelch(OL.Mip3turb2_store(filter3:end),[],[],[],1/timeStep);
    
    subplot(1, 2, 1)
    semilogx(Foop1, mag2db(MoopM1),'LineWidth',lw)
    hold on
    semilogx(Foop2, mag2db(MoopM2),'LineWidth',lw)
    semilogx(Foop3, mag2db(MoopM3),'LineWidth',lw)
    xline(Freq, '--', 'LineWidth', lw)
    xline(8/60, '--', 'LineWidth', lw)
    hold off
    grid on
    xlabel('Frequency [Hz]');
    ylabel('Amplitude [dB]');
    legend('Moop1', 'Moop2','Moop3','f_e','1P','Location','southeast')
    title('Moop PSD')
    subplot(1, 2, 2)
    semilogx(Fip1, mag2db(MipM1),'LineWidth',lw)
    hold on
    semilogx(Fip2, mag2db(MipM2),'LineWidth',lw)
    semilogx(Fip3, mag2db(MipM3),'LineWidth',lw)
    xline(Freq, '--', 'LineWidth', lw)
    xline(8/60, '--', 'LineWidth', lw)
    hold off
    grid on
    xlabel('Frequency [Hz]');
    ylabel('Amplitude [dB]');
    legend('Mip1', 'Mip2','Mip3','f_e','1P','Location','southeast')
    title('Mip PSD')
end

% ============== Overeall Detailed Visualization
if strcmp(overallDetailOption, 'Y')
    t0 = (1:(simLength-filter0+1)) * timeStep;
    if strcmp(coorFrame,'HF')
        % Output
        r1 = ones(simLength-filter0+1, 1)*8.5606;
        figure('Name', 'Output Detail HF', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
        subplot(2, 1, 2)
        plot(t0, Baseline.HF_helixCenter_filtered(filter0:end, 1),'Color',color0,'LineWidth', lw)
        hold on
        plot(t0, OL.HF_helixCenter_filtered(filter0:end, 1),'Color',color1,'LineWidth', lw)
        plot(t0, CL.HF_helixCenter_filtered(filter0:end, 1),'Color',color2,'LineWidth', lw)
        plot(t0, FL.HF_helixCenter_filtered(filter0:end, 1),'Color',color3,'LineWidth', lw)
%         plot(t0, r1, 'k:', 'LineWidth', lw)
        yline(0, 'k-', 'LineWidth', lw)
        hold off
        title('y^e')
        xlim([0 t0(end)])
        xlabel('Time [s]')
        ylim([-1 15])
        ylabel('Position [m]')
        legend('Uniform','Shear','Turbulence','S&T','Location','southeast')
        subplot(2, 1, 1)
        r2 = ones(simLength-filter0+1, 1)*9.0661;
        plot(t0, Baseline.HF_helixCenter_filtered(filter0:end, 2),'Color',color0,'LineWidth', lw)
        hold on
        plot(t0, OL.HF_helixCenter_filtered(filter0:end, 2),'Color',color1,'LineWidth', lw)
        plot(t0, CL.HF_helixCenter_filtered(filter0:end, 2),'Color',color2,'LineWidth', lw)
        plot(t0, FL.HF_helixCenter_filtered(filter0:end, 2),'Color',color3,'LineWidth', lw)
%         plot(t0, r2, 'k:', 'LineWidth', lw)
        yline(0, 'k-', 'LineWidth', lw)
        hold off
        title('z^e')
        xlim([0 t0(end)])
        xlabel('Time [s]')
        ylim([-1 15])
        ylabel('Position [m]')
        legend('Uniform','Shear','Turbulence','S&T','Location','southeast')
        setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)
    elseif strcmp(coorFrame, 'FF')
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

% ============== Hub Jet Trajectory
if strcmp(trajOption, 'Y')
    center_bl = mean(Baseline.FF_helixCenter_filtered(filter:end, :));
    center_ol = mean(OL.FF_helixCenter_filtered(filter:end, :))+[0 0.5];
    center_cl = mean(CL.FF_helixCenter_filtered(filter:end, :));
    center_fl = mean(FL.FF_helixCenter_filtered(filter:end, :));
    
    figure('Name', 'HubJet Trajectory', 'NumberTitle', 'off', 'Position', [100, 100, 1050, 300]);
    subplot(1, 3, 1)
    plot(Baseline.FF_helixCenter_filtered(filter:end, 2), Baseline.FF_helixCenter_filtered(filter:end, 1), 'Color',color0, 'LineWidth', lw)
    hold on
    plot(OL.FF_helixCenter_filtered(5000:end, 2), OL.FF_helixCenter_filtered(5000:end, 1), 'Color',color1, 'LineWidth', lw)
    plot(0, 90, 'k*', 'MarkerSize', 10);
    plot(center_bl(2), center_bl(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color0);
    plot(center_ol(2), center_ol(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color1);
    hold off
    xlabel('y [m]')
    ylabel('z [m]')
    xlim([-30 20])
    ylim([60 110])
    legend('Uniform', 'Shear','Hub', 'Location','southeast')

    subplot(1, 3, 2)
    plot(Baseline.FF_helixCenter_filtered(filter:end, 2), Baseline.FF_helixCenter_filtered(filter:end, 1), 'Color',color0, 'LineWidth', lw)
    hold on
    plot(CL.FF_helixCenter_filtered(filter:end, 2), CL.FF_helixCenter_filtered(filter:end, 1), 'Color',color2, 'LineWidth', lw)
    plot(0, 90, 'k*', 'MarkerSize', 10);
    plot(center_bl(2), center_bl(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color0);
    plot(center_cl(2), center_cl(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color2);
    hold off
    xlabel('y [m]')
    ylabel('z [m]')
    xlim([-30 20])
    ylim([60 110])
    legend('Uniform', 'Turbulence','Hub', 'Location','southeast')

    subplot(1, 3, 3)
    plot(Baseline.FF_helixCenter_filtered(filter:end, 2), Baseline.FF_helixCenter_filtered(filter:end, 1), 'Color',color0, 'LineWidth', lw)
    hold on
    plot(FL.FF_helixCenter_filtered(5000:end, 2), FL.FF_helixCenter_filtered(5000:end, 1), 'Color',color3, 'LineWidth', lw)
    plot(0, 90, 'k*', 'MarkerSize', 10);
    plot(center_bl(2), center_bl(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color0);
    plot(center_fl(2), center_fl(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color2);
    hold off
    xlabel('y [m]')
    ylabel('z [m]')
    xlim([-30 20])
    ylim([60 110])
    legend('Uniform', 'S&T','Hub', 'Location','southeast')
    setfigpaper('Width',[40,0.3],'Interpreter','tex','FontSize',Font,'linewidth',lw)
end

% ============== Video comparison
if strcmp(videoOption, 'Y')
%     ringVisualization2(Baseline, D_NREL5MW)
    videoCompare_func(OL,CL,D_NREL5MW,'')
end

% ============== Power and Fatigue Analysis
% =========== Uniform
% Upstream WT1
BL_result.WT1.power = calculatePower(filter,Baseline.Power_store,D_NREL5MW,U_inflow); % [MW]
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
BL_result.WT2.power = calculatePower(filter,Baseline.Powerturb2_store,D_NREL5MW,U_inflow); % [MW]
BL_result.WT2.DEL = calculateDEL(filter, ...
    Baseline.Mflap1turb2_store,Baseline.Medge1turb2_store, ...
    Baseline.Mflap2turb2_store,Baseline.Medge2turb2_store, ...
    Baseline.Mflap3turb2_store,Baseline.Medge3turb2_store, ...
    timeStep); % [Nm]
BL_result.WT2.PBD = calculatePBD(filter,Baseline.PitchAnglesturb2, ...
    Baseline.Mflap1turb2_store,Baseline.Medge1turb2_store, ...
    Baseline.Mflap2turb2_store,Baseline.Medge2turb2_store, ...
    Baseline.Mflap3turb2_store,Baseline.Medge3turb2_store); % [kNm deg]
BL_result.WT1.DEL_flapwise = mean(BL_result.WT1.DEL(1)+BL_result.WT1.DEL(3)+BL_result.WT1.DEL(5));
BL_result.WT1.DEL_edgewise = mean(BL_result.WT1.DEL(2)+BL_result.WT1.DEL(4)+BL_result.WT1.DEL(6));
BL_result.WT2.DEL_flapwise = mean(BL_result.WT2.DEL(1)+BL_result.WT2.DEL(3)+BL_result.WT2.DEL(5));
BL_result.WT2.DEL_edgewise = mean(BL_result.WT2.DEL(2)+BL_result.WT2.DEL(4)+BL_result.WT2.DEL(6));
BL_result.All.power = BL_result.WT1.power + BL_result.WT2.power;
BL_result.All.DEL_flapwise = BL_result.WT1.DEL_flapwise + BL_result.WT2.DEL_flapwise;
BL_result.All.DEL_edgewise = BL_result.WT1.DEL_edgewise + BL_result.WT2.DEL_edgewise;
BL_result.All.PBD = BL_result.WT1.PBD + BL_result.WT2.PBD;

% =========== Shear
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

% =========== Turbulence
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

% =========== Shear & Turbulence
FL_result.WT1.power = calculatePower(filter,FL.Power_store,D_NREL5MW,U_inflow); % [MW]
FL_result.WT1.DEL = calculateDEL(filter, ...
    FL.Mflap1_store,FL.Medge1_store, ...
    FL.Mflap2_store,FL.Medge2_store, ...
    FL.Mflap3_store,FL.Medge3_store, ...
    timeStep); % [Nm]
FL_result.WT1.PBD = calculatePBD(filter,FL.PitchAngles, ...
    FL.Mflap1_store,FL.Medge1_store, ...
    FL.Mflap2_store,FL.Medge2_store, ...
    FL.Mflap3_store,FL.Medge3_store); % [kNm deg]
% Downstream WT2
FL_result.WT2.power = calculatePower(filter,FL.Powerturb2_store,D_NREL5MW,U_inflow); % [MW]
FL_result.WT2.DEL = calculateDEL(filter, ...
    FL.Mflap1turb2_store,FL.Medge1turb2_store, ...
    FL.Mflap2turb2_store,FL.Medge2turb2_store, ...
    FL.Mflap3turb2_store,FL.Medge3turb2_store, ...
    timeStep); % [Nm]
FL_result.WT2.PBD = calculatePBD(filter,FL.PitchAnglesturb2, ...
    FL.Mflap1turb2_store,FL.Medge1turb2_store, ...
    FL.Mflap2turb2_store,FL.Medge2turb2_store, ...
    FL.Mflap3turb2_store,FL.Medge3turb2_store); % [kNm deg]
FL_result.WT1.DEL_flapwise = mean(FL_result.WT1.DEL(1)+FL_result.WT1.DEL(3)+FL_result.WT1.DEL(5));
FL_result.WT1.DEL_edgewise = mean(FL_result.WT1.DEL(2)+FL_result.WT1.DEL(4)+FL_result.WT1.DEL(6));
FL_result.WT2.DEL_flapwise = mean(FL_result.WT2.DEL(1)+FL_result.WT2.DEL(3)+FL_result.WT2.DEL(5));
FL_result.WT2.DEL_edgewise = mean(FL_result.WT2.DEL(2)+FL_result.WT2.DEL(4)+FL_result.WT2.DEL(6));
FL_result.All.power = FL_result.WT1.power + FL_result.WT2.power;
FL_result.All.DEL_flapwise = FL_result.WT1.DEL_flapwise + FL_result.WT2.DEL_flapwise;
FL_result.All.DEL_edgewise = FL_result.WT1.DEL_edgewise + FL_result.WT2.DEL_edgewise;
FL_result.All.PBD = FL_result.WT1.PBD + FL_result.WT2.PBD;

% === Power and DEL (Result Make Sense or Not)
if strcmp(powerDELAnalysis, 'Y')
    customColors = [0, 0.4470, 0.7410;  % Blue
                    0.8500, 0.3250, 0.0980;  % Orange
                    0.3010, 0.7450, 0.9330]; % Yellow
    filter2 = 2000;
    x = [1 2 3];
    deltaPower = [(OL_result.WT1.power-BL_result.WT1.power)/(BL_result.WT1.power) (OL_result.WT2.power-BL_result.WT2.power)/(BL_result.WT2.power) (OL_result.All.power-BL_result.All.power)/(BL_result.All.power);
                  (CL_result.WT1.power-BL_result.WT1.power)/(BL_result.WT1.power) (CL_result.WT2.power-BL_result.WT2.power)/(BL_result.WT2.power) (CL_result.All.power-BL_result.All.power)/(BL_result.All.power);
                  (FL_result.WT1.power-BL_result.WT1.power)/(BL_result.WT1.power) (FL_result.WT2.power-BL_result.WT2.power)/(BL_result.WT2.power) (FL_result.All.power-BL_result.All.power)/(BL_result.All.power)];
    deltaDELf = [(OL_result.WT1.DEL_flapwise-BL_result.WT1.DEL_flapwise)/(BL_result.WT1.DEL_flapwise) (OL_result.WT2.DEL_flapwise-BL_result.WT2.DEL_flapwise)/(BL_result.WT2.DEL_flapwise) (OL_result.All.DEL_flapwise-BL_result.All.DEL_flapwise)/(BL_result.All.DEL_flapwise);
                 (CL_result.WT1.DEL_flapwise-BL_result.WT1.DEL_flapwise)/(BL_result.WT1.DEL_flapwise) (CL_result.WT2.DEL_flapwise-BL_result.WT2.DEL_flapwise)/(BL_result.WT2.DEL_flapwise) (CL_result.All.DEL_flapwise-BL_result.All.DEL_flapwise)/(BL_result.All.DEL_flapwise);
                 (FL_result.WT1.DEL_flapwise-BL_result.WT1.DEL_flapwise)/(BL_result.WT1.DEL_flapwise) (FL_result.WT2.DEL_flapwise-BL_result.WT2.DEL_flapwise)/(BL_result.WT2.DEL_flapwise) (FL_result.All.DEL_flapwise-BL_result.All.DEL_flapwise)/(BL_result.All.DEL_flapwise)];
    deltaDELe = [(OL_result.WT1.DEL_edgewise-BL_result.WT1.DEL_edgewise)/(BL_result.WT1.DEL_edgewise) (OL_result.WT2.DEL_edgewise-BL_result.WT2.DEL_edgewise)/(BL_result.WT2.DEL_edgewise) (OL_result.All.DEL_edgewise-BL_result.All.DEL_edgewise)/(BL_result.All.DEL_edgewise);
                 (CL_result.WT1.DEL_edgewise-BL_result.WT1.DEL_edgewise)/(BL_result.WT1.DEL_edgewise) (CL_result.WT2.DEL_edgewise-BL_result.WT2.DEL_edgewise)/(BL_result.WT2.DEL_edgewise) (CL_result.All.DEL_edgewise-BL_result.All.DEL_edgewise)/(BL_result.All.DEL_edgewise);
                 (FL_result.WT1.DEL_edgewise-BL_result.WT1.DEL_edgewise)/(BL_result.WT1.DEL_edgewise) (FL_result.WT2.DEL_edgewise-BL_result.WT2.DEL_edgewise)/(BL_result.WT2.DEL_edgewise) (FL_result.All.DEL_edgewise-BL_result.All.DEL_edgewise)/(BL_result.All.DEL_edgewise)];
    x_labels = {'WT1', 'WT2', 'All'}; 
    figure('Name', 'Power & DEL', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(1, 3, 1)
    b = bar(x,deltaPower'*100);
    for i = 1:length(b)
        b(i).FaceColor = 'flat';
        b(i).CData = repmat(customColors(i, :), size(deltaPower, 1), 1); % Apply custom colors
    end
    xticks(x); 
    xticklabels(x_labels); 
    ylabel('\Delta Power [%]')
%     ylim([-15 10])
    legend('S', 'T','S&T', 'Location','southeast')
    title('Power')

    subplot(1, 3, 2)
    b = bar(x,deltaDELf'*100);
    for i = 1:length(b)
        b(i).FaceColor = 'flat';
        b(i).CData = repmat(customColors(i, :), size(deltaPower, 1), 1); % Apply custom colors
    end
    xticks(x); 
    xticklabels(x_labels); 
    ylabel('\Delta DEL [%]')
    legend('S', 'T','S&T', 'Location','southeast')
    title('DEL Flapwise')
    
    subplot(1, 3, 3)
    b = bar(x,deltaDELe'*100);
    for i = 1:length(b)
        b(i).FaceColor = 'flat';
        b(i).CData = repmat(customColors(i, :), size(deltaPower, 1), 1); % Apply custom colors
    end
    xticks(x); 
    xticklabels(x_labels); 
    ylabel('\Delta DEL [%]')
    legend('S', 'T','S&T', 'Location','southeast')
    title('DEL Edgewise')
    setfigpaper('Width',[40,0.3],'Interpreter','tex','FontSize',Font,'linewidth',lw)
end

% === PBD
% PBD is only compared between the open-loop and closed-loop because
% baseline does not have PBD since Helix is not activated
if strcmp(PBDAnalysis, 'Y')
    disp(mean(OL_result.WT1.PBD) - mean(BL_result.WT1.PBD))
    disp(mean(CL_result.WT1.PBD) - mean(BL_result.WT1.PBD))
    disp(mean(FL_result.WT1.PBD) - mean(BL_result.WT1.PBD))

    disp((mean(OL_result.WT1.PBD) - mean(BL_result.WT1.PBD))/mean(BL_result.WT1.PBD))
    disp((mean(CL_result.WT1.PBD) - mean(BL_result.WT1.PBD))/mean(BL_result.WT1.PBD))
    disp((mean(FL_result.WT1.PBD) - mean(BL_result.WT1.PBD))/mean(BL_result.WT1.PBD))
end

% ============== Story telling
if strcmp(storyTellingBasic, 'Y')
    % Uniform
    VBLpower = [BL_result.WT1.power BL_result.WT2.power BL_result.All.power];
    VBLDELflap = [BL_result.WT1.DEL_flapwise BL_result.WT2.DEL_flapwise BL_result.All.DEL_flapwise];
    VBLDELedge = [BL_result.WT1.DEL_edgewise BL_result.WT2.DEL_edgewise BL_result.All.DEL_edgewise];
    % Shear
    OLpower = [OL_result.WT1.power OL_result.WT2.power OL_result.All.power];
    OLDELflap = [OL_result.WT1.DEL_flapwise OL_result.WT2.DEL_flapwise OL_result.All.DEL_flapwise];
    OLDELedge = [OL_result.WT1.DEL_edgewise OL_result.WT2.DEL_edgewise OL_result.All.DEL_edgewise];
    % Turbulence
    CLpower = [CL_result.WT1.power CL_result.WT2.power CL_result.All.power];
    CLDELflap = [CL_result.WT1.DEL_flapwise CL_result.WT2.DEL_flapwise CL_result.All.DEL_flapwise];
    CLDELedge = [CL_result.WT1.DEL_edgewise CL_result.WT2.DEL_edgewise CL_result.All.DEL_edgewise];
    % Shear&Turbulence
    FLpower = [FL_result.WT1.power FL_result.WT2.power FL_result.All.power];
    FLDELflap = [FL_result.WT1.DEL_flapwise FL_result.WT2.DEL_flapwise FL_result.All.DEL_flapwise];
    FLDELedge = [FL_result.WT1.DEL_edgewise FL_result.WT2.DEL_edgewise FL_result.All.DEL_edgewise];
    
    x_labels = {'WT1', 'WT2', 'All'};   

    % == Shear Power
    figure('Name', 'onlyPower Shear vs ref', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    bar([VBLpower;OLpower]', 'EdgeColor', 'k', 'LineWidth', 1.5);
    hold on
    bar([VBLpower;VBLpower]', 'FaceColor', 'none', 'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
    hold off
    set(gca, 'XTickLabel', x_labels);
    legend({'Uniform', 'Shear'}, 'Location', 'northwest');
    ylabel('Power [MW]');
    setfigpaper('Width',[30,0.4],'Interpreter','tex','FontSize',Font,'linewidth',lw);

    % == Shear
    figure('Name', 'Shear vs ref', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    % Power
    subplot(1, 3, 1)
    b1 = bar(OLpower, 'EdgeColor', 'k', 'FaceColor', color1, 'LineWidth', 1.5);
    hold on;
    b2 = bar(VBLpower, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
    set(gca, 'XTickLabel', x_labels);
    legend([b1, b2], {'Shear', 'Uni'}, 'Location', 'northwest');
    ylabel('Power [MW]');
    hold off;
    
    % DEL Flapwise
    subplot(1, 3, 2)
    b1 = bar(OLDELflap, 'EdgeColor', 'k', 'FaceColor', color1, 'LineWidth', 1.5);
    hold on;
    b2 = bar(VBLDELflap, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
    set(gca, 'XTickLabel', x_labels);
    legend([b1, b2], {'Shear', 'Uni'}, 'Location', 'northwest');
    ylabel('DEL Flapwise [Nm]');
    hold off;
    
    % DEL Edgewise
    subplot(1, 3, 3)
    b1 = bar(OLDELedge, 'EdgeColor', 'k', 'FaceColor', color1, 'LineWidth', 1.5);
    hold on;
    b2 = bar(VBLDELedge, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
    set(gca, 'XTickLabel', x_labels);
    legend([b1, b2], {'Shear', 'Uni'}, 'Location', 'northwest');
    ylabel('DEL Edgewise [Nm]');
    hold off;
    setfigpaper('Width',[40,0.3],'Interpreter','tex','FontSize',Font,'linewidth',lw);

    % == Turbulence
    figure('Name', 'Turb vs ref', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    % Power
    subplot(1, 3, 1)
    b1 = bar(CLpower, 'EdgeColor', 'k', 'FaceColor', color2, 'LineWidth', 1.5);
    hold on;
    b2 = bar(VBLpower, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
    set(gca, 'XTickLabel', x_labels);
    legend([b1, b2], {'Turb', 'Uni'}, 'Location', 'northwest');
    ylabel('Power [MW]');
    hold off;
    
    % DEL Flapwise
    subplot(1, 3, 2)
    b1 = bar(CLDELflap, 'EdgeColor', 'k', 'FaceColor', color2, 'LineWidth', 1.5);
    hold on;
    b2 = bar(VBLDELflap, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
    set(gca, 'XTickLabel', x_labels);
    legend([b1, b2], {'Turb', 'Uni'}, 'Location', 'northwest');
    ylabel('DEL Flapwise [Nm]');
    hold off;
    
    % DEL Edgewise
    subplot(1, 3, 3)
    b1 = bar(CLDELedge, 'EdgeColor', 'k', 'FaceColor', color2, 'LineWidth', 1.5);
    hold on;
    b2 = bar(VBLDELedge, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
    set(gca, 'XTickLabel', x_labels);
    legend([b1, b2], {'Turb', 'Uni'}, 'Location', 'northwest');
    ylabel('DEL Edgewise [Nm]');
    hold off;
    setfigpaper('Width',[40,0.3],'Interpreter','tex','FontSize',Font,'linewidth',lw);

    % == Shear & Turbulence
    figure('Name', 'S&T vs ref', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    % Power
    subplot(1, 3, 1)
    b1 = bar(FLpower, 'EdgeColor', 'k', 'FaceColor', color3, 'LineWidth', 1.5);
    hold on;
    b2 = bar(VBLpower, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
    set(gca, 'XTickLabel', x_labels);
    legend([b1, b2], {'S&T', 'Uni'}, 'Location', 'northwest');
    ylabel('Power [MW]');
    hold off;
    
    % DEL Flapwise
    subplot(1, 3, 2)
    b1 = bar(FLDELflap, 'EdgeColor', 'k', 'FaceColor', color3, 'LineWidth', 1.5);
    hold on;
    b2 = bar(VBLDELflap, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
    set(gca, 'XTickLabel', x_labels);
    legend([b1, b2], {'S&T', 'Uni'}, 'Location', 'northwest');
    ylabel('DEL Flapwise [Nm]');
    hold off;
    
    % DEL Edgewise
    subplot(1, 3, 3)
    b1 = bar(FLDELedge, 'EdgeColor', 'k', 'FaceColor', color3, 'LineWidth', 1.5);
    hold on;
    b2 = bar(VBLDELedge, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
    set(gca, 'XTickLabel', x_labels);
    legend([b1, b2], {'S&T', 'Uni'}, 'Location', 'northwest');
    ylabel('DEL Edgewise [Nm]');
    hold off;
    setfigpaper('Width',[40,0.3],'Interpreter','tex','FontSize',Font,'linewidth',lw);

    % == Three Together
    colors = [color1; color2; color3];
    figure('Name', '3 vs ref', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    % Power
    subplot(1, 3, 1)
    b1 = bar([OLpower;CLpower;FLpower]', 'EdgeColor', 'k', 'LineWidth', 1.5);
    hold on;
    b2 = bar(VBLpower, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
    for i = 1:length(b1)
        b1(i).FaceColor = colors(i, :);
    end
    set(gca, 'XTickLabel', x_labels);
    legend([b1, b2], {'S', 'T', 'S&T', 'Uni'}, 'Location', 'northwest');
    ylabel('Power [MW]');
    hold off;
    
    % DEL Flapwise
    subplot(1, 3, 2)
    b1 = bar([OLDELflap;CLDELflap;FLDELflap]', 'EdgeColor', 'k', 'LineWidth', 1.5);
    hold on;
    b2 = bar(VBLDELflap, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
    for i = 1:length(b1)
        b1(i).FaceColor = colors(i, :);
    end
    set(gca, 'XTickLabel', x_labels);
    legend([b1, b2], {'S', 'T', 'S&T', 'Uni'}, 'Location', 'northwest');
    ylabel('DEL Flapwise [Nm]');
    hold off;
    
    % DEL Edgewise
    subplot(1, 3, 3)
    b1 = bar([OLDELedge;CLDELedge;FLDELedge]', 'EdgeColor', 'k', 'LineWidth', 1.5);
    hold on;
    b2 = bar(VBLDELedge, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
    for i = 1:length(b1)
        b1(i).FaceColor = colors(i, :);
    end
    set(gca, 'XTickLabel', x_labels);
    legend([b1, b2], {'S', 'T', 'S&T', 'Uni'}, 'Location', 'northwest');
    ylabel('DEL Edgewise [Nm]');
    hold off;
    setfigpaper('Width',[40,0.3],'Interpreter','tex','FontSize',Font,'linewidth',lw);


%     % == Open-Loop vs. Closed-Loop
%     figure('Name', 'OL&CL vs ref', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
%     % Power
%     subplot(1, 3, 1)
%     b1 = bar([OLpower; CLpower]', 'grouped');
%     hold on
%     b2 = bar([VBLpower; VBLpower]', 'grouped');
%     hold off
%     b2(1).FaceColor = 'none'; 
%     b2(1).EdgeColor = 'k';
%     b2(1).LineStyle = '--';
%     b2(1).LineWidth = 1.5;
%     b2(2).FaceColor = 'none'; 
%     b2(2).EdgeColor = 'k';
%     b2(2).LineStyle = '--';
%     b2(2).LineWidth = 1.5;
%     set(gca, 'XTickLabel', x_labels);
%     legend([b1(1), b1(2), b2], {'OL', 'CL', 'Ref'}, 'Location', 'northwest');
%     ylabel('Power [MW]');
% 
%     % DEL Flapwise
%     subplot(1, 3, 2)
%     b1 = bar([OLDELflap; CLDELflap]', 'grouped');
%     hold on
%     b2 = bar([VBLDELflap; VBLDELflap]', 'grouped');
%     hold off
%     b2(1).FaceColor = 'none'; 
%     b2(1).EdgeColor = 'k';
%     b2(1).LineStyle = '--';
%     b2(1).LineWidth = 1.5;
%     b2(2).FaceColor = 'none'; 
%     b2(2).EdgeColor = 'k';
%     b2(2).LineStyle = '--';
%     b2(2).LineWidth = 1.5;
%     set(gca, 'XTickLabel', x_labels);
%     legend([b1(1), b1(2), b2], {'OL', 'CL', 'Ref'}, 'Location', 'northwest');
%     ylabel('DEL Flapwise [Nm]');
% 
%     % DEL Edgewise
%     subplot(1, 3, 3)
%     b1 = bar([OLDELedge; CLDELedge]', 'grouped');
%     hold on
%     b2 = bar([VBLDELedge; VBLDELedge]', 'grouped');
%     hold off
%     b2(1).FaceColor = 'none'; 
%     b2(1).EdgeColor = 'k';
%     b2(1).LineStyle = '--';
%     b2(1).LineWidth = 1.5;
%     b2(2).FaceColor = 'none'; 
%     b2(2).EdgeColor = 'k';
%     b2(2).LineStyle = '--';
%     b2(2).LineWidth = 1.5;
%     set(gca, 'XTickLabel', x_labels);
%     legend([b1(1), b1(2), b2], {'OL', 'CL', 'Ref'}, 'Location', 'northwest');
%     ylabel('DEL Edgewise [Nm]');
%     setfigpaper('Width',[40,0.3],'Interpreter','tex','FontSize',Font,'linewidth',lw);
end