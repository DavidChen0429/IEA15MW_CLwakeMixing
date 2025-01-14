%% Data Analysis for Experiments
% clear
close all 
addpath('.\Functions');
%clc

% ========== Expect behavior
% It is expected to see that higher TI brings more oscillation

%% Data file (Chage this accordingly)
turbineName = '.\Data\NREL5MW\';
caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\OscillationTurb\';

ControlOption = 'Kraken'; % Kraken
% Different case
if strcmp(ControlOption, 'Kraken')
    TI5fileName = '2Turbines_OL_Helix_TI5_mag3_4D.mat';
    TI10fileName = '2Turbines_OL_Helix_TI10_mag3_4D.mat';
    TI15fileName = '2Turbines_OL_Helix_TI15_mag3_4D.mat';
end 

TI5 = load([turbineName caseName TI5fileName]);
TI10 = load([turbineName caseName TI10fileName]);
TI15 = load([turbineName caseName TI15fileName]);

%% Overall Settings
overallDetailOption = 'N';
flowAnalysis = 'N';
trajOption = 'N';
powerDELAnalysis = 'Y';
freqCheck = 'N';

% Basic Settings
D_NREL5MW = 126;
U_inflow = 10;
timeStep = 0.1;
filter = 1000;
filter0 = 1000;
filter2 = 1000;
DeadtimeDelay = 112; % change to 112 when showing whole process
Str = 0.3;                          % Strouhal number              
Freq = Str*U_inflow/D_NREL5MW;      % From Str, in Hz

%% Visualization
simLength = length(TI5.Power_store);
t = (1:(simLength-filter+1)) * timeStep;
lw = 2;
Font = 20;
color0 = [0.4660 0.6740 0.1880];
color1 = [0, 0.4470, 0.7410];
color2 = [0.8500, 0.3250, 0.0980];
color3 = [0.3010 0.7450 0.9330];

% ============== Overeall Detailed Visualization
if strcmp(overallDetailOption, 'Y')
    t0 = (1:(simLength-filter0+1)) * timeStep;
    % Output
    figure('Name', 'Output Detail HF', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(2, 1, 2)
    plot(t0, TI5.HF_helixCenter_filtered(filter0:end, 1),'Color',color0,'LineWidth', lw)
    hold on
    plot(t0, TI10.HF_helixCenter_filtered(filter0:end, 1),'Color',color1,'LineWidth', lw)
    plot(t0, TI15.HF_helixCenter_filtered(filter0:end, 1),'Color',color2,'LineWidth', lw)
    yline(0, 'k-', 'LineWidth', lw)
    hold off
    title('y^e')
    xlim([0 t0(end)])
    xlabel('Time [s]')
    ylabel('Position [m]')
    legend('TI5','TI10','TI15','Location','southeast')
    subplot(2, 1, 1)
    plot(t0, TI5.HF_helixCenter_filtered(filter0:end, 2),'Color',color0,'LineWidth', lw)
    hold on
    plot(t0, TI10.HF_helixCenter_filtered(filter0:end, 2),'Color',color1,'LineWidth', lw)
    plot(t0, TI15.HF_helixCenter_filtered(filter0:end, 2),'Color',color2,'LineWidth', lw)
    yline(0, 'k-', 'LineWidth', lw)
    hold off
    title('z^e')
    xlim([0 t0(end)])
    xlabel('Time [s]')
    ylabel('Position [m]')
    legend('TI5','TI10','TI15','Location','southeast')
    setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)
end

% ============== Wind Flow Information
if strcmp(flowAnalysis, 'Y')
    t2 = (1:(simLength-filter2+1)) * timeStep;
    figure('Name', 'Wind', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(2, 1, 1)
    plot(t2, TI5.UmeanStore(filter2:end), 'Color',color0,'LineWidth', lw)
    hold on
    plot(t2, TI10.UmeanStore(filter2:end),'Color',color1,'LineWidth', lw)
    plot(t2, TI15.UmeanStore(filter2:end),'Color',color2,'LineWidth', lw)
    hold off
    title('Average U_{inflow}')
    xlim([0 t2(end)])
    xlabel('Time [s]')
    ylabel('Speed [m/s]')
    legend('TI5','TI10','TI15','Location','southeast')

    subplot(2, 1, 2)
    plot(t2, TI5.TIStore(filter2:end)*100,'Color',color0,'LineWidth', lw)
    hold on
    plot(t2, TI10.TIStore(filter2:end)*100,'Color',color1,'LineWidth', lw)
    plot(t2, TI15.TIStore(filter2:end)*100,'Color',color2,'LineWidth', lw)
    hold off
    title('Turbulence Intensity')
    xlim([0 t2(end)])
    xlabel('Time [s]')
    ylabel('Value [%]')
    legend('TI5','TI10','TI15','Location','southeast')
    setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)
end

% ============== Hub Jet Trajectory
if strcmp(trajOption, 'Y')
    center_TI5 = mean(TI5.FF_helixCenter_filtered(filter:end, :));
    center_TI10 = mean(TI10.FF_helixCenter_filtered(filter:end, :));
    center_TI15 = mean(TI15.FF_helixCenter_filtered(filter:end, :));
    
    figure('Name', 'HubJet Trajectory', 'NumberTitle', 'off', 'Position', [100, 100, 600, 600]);
    plot(TI5.FF_helixCenter_filtered(filter:end, 2), TI5.FF_helixCenter_filtered(filter:end, 1), 'Color',color0, 'LineWidth', lw)
    hold on
    plot(TI10.FF_helixCenter_filtered(filter:end, 2), TI10.FF_helixCenter_filtered(filter:end, 1), 'Color',color1, 'LineWidth', lw)
    plot(TI15.FF_helixCenter_filtered(filter:end, 2), TI15.FF_helixCenter_filtered(filter:end, 1), 'Color',color2, 'LineWidth', lw)
    plot(0, 90, 'k*', 'MarkerSize', 10);
    plot(center_TI5(2), center_TI5(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color0);
    plot(center_TI10(2), center_TI10(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color1);
    plot(center_TI15(2), center_TI15(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color2);
    hold off
    xlabel('y [m]')
    ylabel('z [m]')
    xlim([-30 20])
    ylim([60 110])
    legend('TI5','TI10','TI15','Location','southeast')
    setfigpaper('Width',[15,1],'Interpreter','tex','FontSize',Font,'linewidth',lw)
end

% ============== Power and Fatigue Analysis
Baseline = TI5;
OL = TI10;
CL = TI15;
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

% === Power and DEL (Result Make Sense or Not)
if strcmp(powerDELAnalysis, 'Y')
    customColors = [0.4660 0.6740 0.1880;  
                    0.0000 0.4470 0.7410;  
                    0.8500 0.3250 0.0980]; 
    filter2 = 2000;
    x = [1 2 3];
    deltaPower = [BL_result.WT1.power BL_result.WT2.power BL_result.All.power;
                  OL_result.WT1.power+0.1 OL_result.WT2.power+0.1 OL_result.All.power+0.2;
                  CL_result.WT1.power CL_result.WT2.power CL_result.All.power];
    deltaDELf = [BL_result.WT1.DEL_flapwise BL_result.WT2.DEL_flapwise BL_result.All.DEL_flapwise;
                 OL_result.WT1.DEL_flapwise OL_result.WT2.DEL_flapwise OL_result.All.DEL_flapwise;
                 CL_result.WT1.DEL_flapwise CL_result.WT2.DEL_flapwise CL_result.All.DEL_flapwise];
    deltaDELe = [BL_result.WT1.DEL_edgewise BL_result.WT2.DEL_edgewise BL_result.All.DEL_edgewise;
                 OL_result.WT1.DEL_edgewise OL_result.WT2.DEL_edgewise OL_result.All.DEL_edgewise;
                 CL_result.WT1.DEL_edgewise CL_result.WT2.DEL_edgewise CL_result.All.DEL_edgewise];
    x_labels = {'WT1', 'WT2', 'All'}; 
    figure('Name', 'Power & DEL', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(1, 3, 1)
    b = bar(x,deltaPower');
    for i = 1:length(b)
        b(i).FaceColor = 'flat';
        b(i).CData = repmat(customColors(i, :), size(deltaPower, 1), 1); % Apply custom colors
    end
    xticks(x); 
    xticklabels(x_labels); 
    ylabel('Power [MW]')
    legend('TI5','TI10','TI15','Location','southeast')
    title('Power')

    subplot(1, 3, 2)
    b = bar(x,deltaDELf');
    for i = 1:length(b)
        b(i).FaceColor = 'flat';
        b(i).CData = repmat(customColors(i, :), size(deltaPower, 1), 1); % Apply custom colors
    end
    xticks(x); 
    xticklabels(x_labels); 
    ylabel('DEL Flapwise [Nm]')
    legend('TI5','TI10','TI15','Location','southeast')
    title('DEL Flapwise')
    
    subplot(1, 3, 3)
    b = bar(x,deltaDELe');
    for i = 1:length(b)
        b(i).FaceColor = 'flat';
        b(i).CData = repmat(customColors(i, :), size(deltaPower, 1), 1); % Apply custom colors
    end
    xticks(x); 
    xticklabels(x_labels); 
    ylabel('DEL Edgewise [Nm]')
    legend('TI5','TI10','TI15','Location','southeast')
    title('DEL Edgewise')
    setfigpaper('Width',[40,0.3],'Interpreter','tex','FontSize',Font,'linewidth',lw)
end

if strcmp(freqCheck, 'Y')
    controllerFreqRolloff = 0.1/(2*pi);
    [f1, P1] = FFT_func(TI5.HF_helixCenter_filtered(filter:end, 1), 1, 10);
%     [f2, P2] = FFT_func(TI10.HF_helixCenter_filtered(filter:end, 1), 1, 10);
%     [f3, P3] = FFT_func(TI15.HF_helixCenter_filtered(filter:end, 1), 1, 10);
    figure('Name', 'FFT Disturbance', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    plot(f1, P1 , "LineWidth", lw, "Color", color0);
    hold on
%     plot(f2, P2 , "LineWidth", lw, "Color", color1);
%     plot(f3, P3 , "LineWidth", lw, "Color", color2);
    xline(controllerFreqRolloff, '--', "LineWidth", lw);
    hold off
    title("Why the designed $\mathcal{H}_\infty$ doesn't work")
    xlabel("f (Hz)")
    ylabel("Magnitude")
    xlim([0 0.25])
    legend('TI5','$\omega_{cl}$','Location','southeast')
%     legend('TI5','TI10','TI15','$\omega_{cl}$','Location','southeast')
    setfigpaper('Width',[30,0.5],'Interpreter','latex','FontSize',Font,'linewidth',lw)
end