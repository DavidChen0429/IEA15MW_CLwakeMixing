%% Data Analysis for Experiments
% clear
close all 
addpath('.\Functions');
clc

%% Data file (Chage this accordingly)
turbineName = '.\Data\NREL5MW\';
caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\2TurbinesLonger\';

% Different case
windCase = 'BothRC2';
% Uniform
% Turb 
% ShearRC3, Shear, ShearRC, ShearRC2
% BothRC2, Both, BothRC,
% SkewedRC2, SkewedRC

veryBaseOLfile = '2Turbines_OL_Helix_mag3_4D.mat';
if strcmp(windCase, 'Uniform')
    basefile = '2Turbines_Baseline_4D.mat';
    OLfileName = '2Turbines_OL_Helix_mag3_4D.mat';
    CLfileName = '2Turbines_CL_Helix_ramp&stop_mag3_4D.mat';
    filter = 3000;
elseif strcmp(windCase, 'Shear')
    basefile = '2Turbines_Baseline_Shear0.2_4D.mat';
    OLfileName = '2Turbines_OL_Helix_Shear0.2_mag3_4D.mat';
    CLfileName = '2Turbines_CL_Helix_Shear0.2_mag3_4D.mat';
    filter = 3000;
elseif strcmp(windCase, 'ShearRC')
    basefile = '2Turbines_Baseline_ShearReCenter_4D.mat';
    OLfileName = '2Turbines_OL_Helix_ShearReCenter_mag3_4D.mat';
    CLfileName = '2Turbines_CL_Helix_ShearReCenter_mag3_4D.mat';
    filter = 3000;
elseif strcmp(windCase, 'ShearRC2')
    basefile = '2Turbines_Baseline_ShearReCenter_4D.mat';
    OLfileName = '2Turbines_OL_Helix_ShearReCenter2_mag3.3_4D.mat';
    CLfileName = '2Turbines_CL_Helix_ShearReCenter2_mag3_4D.mat';
    filter = 3500;
elseif strcmp(windCase, 'ShearRC3')
    basefile = '2Turbines_Baseline_ShearReCenter_4D.mat';
    OLfileName = '2Turbines_OL_Helix_ShearReCenter_mag3_4D.mat';
    CLfileName = '2Turbines_CL_Helix_ShearReCenter3_mag3_4D.mat';
    filter = 3000;
elseif strcmp(windCase, 'Turb')
    basefile = '2Turbines_Baseline_TI6_4D.mat';
    OLfileName = '2Turbines_OL_Helix_TI6_mag3_4D.mat';
    CLfileName = '2Turbines_CL_Helix_TI6_mag3_4D.mat';    
    filter = 3000;
elseif strcmp(windCase, 'Both')
    basefile = '2Turbines_Baseline_TI6&Shear0.2_4D.mat';
    OLfileName = '2Turbines_OL_Helix_TI6&Shear0.2_mag3_4D.mat';
    CLfileName = '2Turbines_CL_Helix_TI6&Shear0.2_mag3_4D.mat';
    filter = 3000;
elseif strcmp(windCase, 'SkewedRC')
    basefile = '2Turbines_Baseline_Skewed1ReCenter_4D.mat';
    OLfileName = '2Turbines_OL_Helix_Skewed1ReCenter_mag3_4D.mat';
    CLfileName = '2Turbines_CL_Helix_Skewed1ReCenter_mag3_4D.mat';
    filter = 3000;
elseif strcmp(windCase, 'SkewedRC2')
    basefile = '2Turbines_Baseline_Skewed1ReCenter_4D.mat';
    OLfileName = '2Turbines_OL_Helix_Skewed1ReCenter_mag3_4D.mat';
    CLfileName = '2Turbines_CL_Helix_Skewed1ReCenter2_mag3_4D.mat'; 
    filter = 3000;
elseif strcmp(windCase, 'BothRC')
    basefile = '2Turbines_Baseline_BothReCenter_4D.mat';
    OLfileName = '2Turbines_OL_Helix_BothReCenter_mag3_4D.mat';
    CLfileName = '2Turbines_CL_Helix_BothReCenter_mag3_4D.mat'; 
    filter = 3000;
elseif strcmp(windCase, 'BothRC2')
    basefile = '2Turbines_Baseline_BothReCenter_4D.mat';
    OLfileName = '2Turbines_OL_Helix_BothReCenter_mag3_4D.mat';
    CLfileName = '2Turbines_CL_Helix_BothReCenter2_mag3_4D.mat'; 
    filter = 5000;
end

veryBase = load([turbineName caseName veryBaseOLfile]);
Baseline = load([turbineName caseName basefile]);
OL = load([turbineName caseName OLfileName]);
CL = load([turbineName caseName CLfileName]);
disp(mean(CL.HF_beta(filter:end, :)))

%% Overall Settings
errorOption = 'N';
flowAnalysis = 'Y';
rareDataAnalysis = 'N';
overallDetailOption = 'N';
coorFrame = 'HF';
trajOption = 'Y';
videoOption = 'N';
powerAnalysis = 'N';
DELAnalysis = 'N';
PBDAnalysis = 'N';
powerDELAnalysis = 'N';
numericalAnalysis = 'Y';
storyTellingBasic = 'Y';
showCompnent = 'N';


% Basic Settings
D_NREL5MW = 126;
U_inflow = 10;
timeStep = 0.1;
% filter = 3500; % 3000 for steady-state, 3500 ShearRC3, 5000 BothRC2
filter0 = filter; % Overall result
filter2 = 1; % Wind info
filter3 = 1000; % Rare data
DeadtimeDelay = 112; % change to 112 when showing whole process
Str = 0.3;                          % Strouhal number              
Freq = Str*U_inflow/D_NREL5MW;      % From Str, in Hz
x_labels = {'WT1', 'WT2', 'All'};   

%% Visualization
simLength = length(Baseline.Power_store);
t = (1:(simLength-filter+1)) * timeStep;
lw = 2;
Font = 20;
color0 = [0.4660 0.6740 0.1880];
color1 = [0, 0.4470, 0.7410];
color2 = [0.8500, 0.3250, 0.0980];

% ============== Error for Tracking 
if strcmp(errorOption, 'Y')
    figure('Name', 'Experiment Error', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(2, 1, 1)
    plot(t, (CL.r(filter:end, 1)-CL.y(filter:end, 1)), 'Color',color1,'LineWidth',lw)
    hold on
    plot(t, (CL.ym(filter:end, 1)-CL.ytilda(filter:end, 1)),'--', 'Color',color2,'LineWidth',lw)
    yline(0, 'k--','LineWidth',lw)
    xline(20, '--k', 'Activate Ctrl', 'LabelOrientation', 'horizontal', 'LineWidth', 1);
    hold off
    xlabel('Time [s]')
    ylabel('Magnitude [m]')
    xlim([0 500])
    ylim([-20 15])
    legend('e_p', 'e_f', 'Location','southeast')
    title('z^e')
    subplot(2, 1, 2)
    plot(t, (CL.r(filter:end, 2)-CL.y(filter:end, 2)), 'Color',color1,'LineWidth',lw)
    hold on
    plot(t, (CL.ym(filter:end, 2)-CL.ytilda(filter:end, 2)),'--', 'Color',color2,'LineWidth',lw)
    yline(0, 'k--','LineWidth',lw)
    xline(20, '--k', 'Activate Ctrl', 'LabelOrientation', 'horizontal', 'LineWidth', 1);
    hold off
    xlabel('Time [s]')
    ylabel('Magnitude [m]')
    xlim([0 500])
    ylim([-20 15])
    legend('e_p', 'e_f', 'Location','southeast')
    title('y^e')
    setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)
end

% ============== Wind Flow Information
if strcmp(flowAnalysis, 'Y')
    t2 = (1:(simLength-filter2+1)) * timeStep;
    figure('Name', 'Wind', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(2, 1, 1)
    plot(t2, OL.UmeanStore(filter2:end), 'Color',color1,'LineWidth', lw)
    hold on
    plot(t2, CL.UmeanStore(filter2:end),'Color',color2,'LineWidth', lw)
    hold off
    title('Average U_{inflow}')
    xlim([0 t2(end)])
    xlabel('Time [s]')
    ylabel('Speed [m/s]')
    legend('OL','CL','Location','southeast')

    subplot(2, 1, 2)
    plot(t2, OL.TIStore(filter2:end)*100,'Color',color1,'LineWidth', lw)
    hold on
    plot(t2, CL.TIStore(filter2:end)*100,'Color',color2,'LineWidth', lw)
    hold off
    title('Turbulence Intensity')
    xlim([0 t2(end)])
    xlabel('Time [s]')
    ylabel('Value [%]')
    legend('OL','CL','Location','southeast')
    setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)

    filter2 = 5000;
    t2 = (1:(simLength-filter2+1)) * timeStep;
    plot(t2, OL.UmeanStore(filter2:end), 'Color',color1,'LineWidth', lw)
    hold on
    plot(t2, CL.UmeanStore(filter2:end),'Color',color2,'LineWidth', lw)
    hold off
    title('Average U_{inflow}')
    xlim([0 t2(end)])
    xlabel('Time [s]')
    ylabel('Speed [m/s]')
    legend('OL','CL','Location','southeast')
    setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)

    fprintf('\n==================== Flow Comparison \n');
    fprintf('========== Wind Speed \n');
    fprintf('   OL: %.2e m/s\n', mean(OL.UmeanStore(filter2:end)));
    fprintf('   CL: %.2e m/s\n', mean(CL.UmeanStore(filter2:end)));
    fprintf('========== Turbulence Intensity \n');
    fprintf('   OL: %.2e %%\n', mean(OL.TIStore(filter2:end)*100));
    fprintf('   CL: %.2e %%\n', mean(CL.TIStore(filter2:end)*100));   
end

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
    hold off
    legend('OL','CL')
    xlim([0 t3(end)])
    xlabel('Time [s]')
    ylabel('Power [MW]')
    title('Power Production')
    subplot(2, 1, 2)
    plot(t3, OL.TorqueStoreTurb2(filter3:end),'Color',color1,'LineWidth', lw)
    hold on
    plot(t3, CL.TorqueStoreTurb2(filter3:end),'Color',color2,'LineWidth', lw)
    hold off
    legend('OL','CL')
    xlim([0 t3(end)])
    xlabel('Time [s]')
    ylabel('Torque [Nm]')
    title('Generator Torque')
    setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)

    figure('Name', 'Rare Data --- Pitch', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    % Pitch
    plot(t3, OL.PitchAngles(filter3:end, 1))
    hold on
    plot(t3, CL.PitchAngles(filter3:end, 1))
    hold off
    legend('OL','CL')
    xlim([0 t3(end)])
    xlabel('Time [s]')
    ylabel('Pitch [deg]')
    title('Pitch Signal')
    legend('OL', 'CL')
    setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)
    
    % Fatigue (Time Domain)
    figure('Name', 'Rare Data --- Fatigue', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(3, 2, 1)
    plot(t3, OL.Moop1turb2_store(filter3:end),'Color',color1,'LineWidth', lw)
    hold on
    plot(t3, CL.Moop1turb2_store(filter3:end),'Color',color2,'LineWidth', lw)
    hold off
    legend('OL','CL')
    xlim([0 t3(end)])
    xlabel('Time [s]')
    ylabel('Moop - Blade1 [NM]')
    subplot(3, 2, 3)
    plot(t3, OL.Moop2turb2_store(filter3:end),'Color',color1,'LineWidth', lw)
    hold on
    plot(t3, CL.Moop2turb2_store(filter3:end),'Color',color2,'LineWidth', lw)
    hold off
    legend('OL','CL')
    xlim([0 t3(end)])
    xlabel('Time [s]')
    ylabel('Moop - Blade2 [NM]')
    subplot(3, 2, 5)
    plot(t3, OL.Moop3turb2_store(filter3:end),'Color',color1,'LineWidth', lw)
    hold on
    plot(t3, CL.Moop3turb2_store(filter3:end),'Color',color2,'LineWidth', lw)
    hold off
    legend('OL','CL')
    xlim([0 t3(end)])
    xlabel('Time [s]')
    ylabel('Moop - Blade3 [NM]')

    subplot(3, 2, 2)
    plot(t3, OL.Mip1turb2_store(filter3:end),'Color',color1,'LineWidth', lw)
    hold on
    plot(t3, CL.Mip1turb2_store(filter3:end),'Color',color2,'LineWidth', lw)
    hold off
    legend('OL','CL')
    xlim([0 t3(end)])
    xlabel('Time [s]')
    ylabel('Mip- Blade1 [NM]')
    subplot(3, 2, 4)
    plot(t3, OL.Mip1turb2_store(filter3:end),'Color',color1,'LineWidth', lw)
    hold on
    plot(t3, CL.Mip1turb2_store(filter3:end),'Color',color2,'LineWidth', lw)
    hold off
    legend('OL','CL')
    xlim([0 t3(end)])
    xlabel('Time [s]')
    ylabel('Mip - Blade2 [NM]')
    subplot(3, 2, 6)
    plot(t3, OL.Mip1turb2_store(filter3:end),'Color',color1,'LineWidth', lw)
    hold on
    plot(t3, CL.Mip1turb2_store(filter3:end),'Color',color2,'LineWidth', lw)
    hold off
    legend('OL','CL')
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
    xline(8.5/60, '--', 'LineWidth', lw)
    hold off
    grid on
    xlabel('Frequency [Hz]');
    ylabel('Amplitude [dB]');
    legend('Moop1','Moop2','Moop3','f_e','1P','Location','southeast')
    title('Moop PSD')
    subplot(1, 2, 2)
    semilogx(Fip1, mag2db(MipM1),'LineWidth',lw)
    hold on
    semilogx(Fip2, mag2db(MipM2),'LineWidth',lw)
    semilogx(Fip3, mag2db(MipM3),'LineWidth',lw)
    xline(Freq, '--', 'LineWidth', lw)
    xline(8.5/60, '--', 'LineWidth', lw)
    hold off
    grid on
    xlabel('Frequency [Hz]');
    ylabel('Amplitude [dB]');
    legend('Mip1','Mip2','Mip3','f_e','1P','Location','southeast')
    title('Mip PSD')
end

% ============== Show Different Component
if strcmp(showCompnent, 'Y')
    figure('Name', 'show Compnent', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(2, 1, 1)
    t3 = (1000:9000) * timeStep;
    plot(t3, CL.ym(filter3:end, 1),'Color',color1,'LineWidth', lw)
    hold on
    plot(t3, CL.y(filter3:end, 1),'Color',color0,'LineWidth', lw)
    plot(t3, CL.ytilda(filter3:end, 1),'Color',color0,'LineWidth', lw, 'LineStyle', '--')
    hold off
    xlabel('Time [sec]')
    ylabel('Position [m]')
    title('Hub Jet $z^e$')
    legend('$y_m$', '$y$', '$\tilde{y}$', 'Interpreter', 'latex')

    subplot(2, 1, 2)
    plot(t3, CL.ym(filter3:end, 2),'Color',color2,'LineWidth', lw)
    hold on
    plot(t3, CL.y(filter3:end, 2),'Color',color0,'LineWidth', lw)
    plot(t3, CL.ytilda(filter3:end, 2),'Color',color0,'LineWidth', lw, 'LineStyle', '--')
    hold off
    xlabel('Time [sec]')
    ylabel('Position [m]')
    title('Hub Jet $y^e$')
    legend('$y_m$', '$y$', '$\tilde{y}$', 'Interpreter', 'latex')
    setfigpaper('Width',[30,0.5],'Interpreter','latex','FontSize',Font,'linewidth',lw)
end

% ============== Overeall Detailed Visualization
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
        ylim([-1 8])
        ylabel('[deg]')
        legend('OL','CL','Location','southeast')
        subplot(2, 1, 1)
        hold on
        plot(t0, OL.HF_beta(filter0:end, 2),'Color',color0,'LineWidth', lw)
        plot(t0, CL.u(filter0:end, 2),'Color',color2, 'LineWidth', lw)
        hold off
        title('\beta^e_{tilt}')
        xlim([0 t0(end)])
        xlabel('Time [s]')
        ylim([-1 8])
        ylabel('[deg]')
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
        subplot(2, 1, 1)
        plot(t0, OL.FF_helixCenter_filtered(filter0:end, 1),'Color',color0,'LineWidth', lw)
        hold on
        plot(t0, CL.FF_helixCenter_filtered(filter0:end, 1),'Color',color2,'LineWidth', lw)
        hold off
        title('z')
        xlim([0 t0(end)])
        xlabel('Time [s]')
        ylabel('Position [m]')
        legend('OL','CL','Location','southeast')
        subplot(2, 1, 2)
        plot(t0, OL.FF_helixCenter_filtered(filter0:end, 2),'Color',color0,'LineWidth', lw)
        hold on
        plot(t0, CL.FF_helixCenter_filtered(filter0:end, 2),'Color',color1,'LineWidth', lw)
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
    if strcmp(windCase,'ShearRC3')
        CL.FF_helixCenter_filtered(:, 2) = CL.FF_helixCenter_filtered(:, 2) + 0.5;
    elseif strcmp(windCase, 'BothRC2')
        CL.FF_helixCenter_filtered(:, 2) = CL.FF_helixCenter_filtered(:, 2) + 0.7;
    end
%     filter = 4000
    center_blb = mean(veryBase.FF_helixCenter_filtered(filter:end, :));
    center_bl = mean(Baseline.FF_helixCenter_filtered(filter:end, :));
    center_ol = mean(OL.FF_helixCenter_filtered(filter:end, :));
    center_cl = mean(CL.FF_helixCenter_filtered(filter:end, :));
    figure('Name', 'HubJet Trajectory', 'NumberTitle', 'off', 'Position', [100, 100, 600, 600]);
%     plot(Baseline.FF_helixCenter_filtered(filter:end, 2), Baseline.FF_helixCenter_filtered(filter:end, 1), 'Color',color0, 'LineWidth', lw)
%     hold on
    plot(OL.FF_helixCenter_filtered(filter:end, 2), OL.FF_helixCenter_filtered(filter:end, 1), 'Color',color1, 'LineWidth', lw)
    hold on
    plot(CL.FF_helixCenter_filtered(filter:end, 2), CL.FF_helixCenter_filtered(filter:end, 1), 'Color',color2, 'LineWidth', lw)
    plot(0, 90, 'k*', 'MarkerSize', 10);
    plot(center_blb(2), center_blb(1), 'k+', 'MarkerSize', 10);
%     plot(center_bl(2), center_bl(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color0);
    plot(center_ol(2), center_ol(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color1);
    plot(center_cl(2), center_cl(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color2);
    hold off
    xlabel('y [m]')
    ylabel('z [m]')
    xlim([-30 20])
    ylim([60 110])
    legend('OL','CL','Hub','Uni OL', 'Location','southeast')
    setfigpaper('Width',[15,1],'Interpreter','tex','FontSize',Font,'linewidth',lw)
end

% ============== Video comparison
if strcmp(videoOption, 'Y')
%     ringVisualization2(Baseline, D_NREL5MW)
    videoCompare_func(OL,CL,D_NREL5MW,'')
end

% ============== Power and Fatigue Analysis
% =========== Very Baseline --- Uniform OL
% Upstream WT1
VBL_result.WT1.power = calculatePower(filter,veryBase.Power_store,D_NREL5MW,U_inflow); % [MW]
% BL_result.WT1.power2 = calculatePower2(filter,Baseline.Power_store);
VBL_result.WT1.DEL = calculateDEL(filter, ...
    veryBase.Mflap1_store,veryBase.Medge1_store, ...
    veryBase.Mflap2_store,veryBase.Medge2_store, ...
    veryBase.Mflap3_store,veryBase.Medge3_store, ...
    timeStep); % [Nm]
VBL_result.WT1.PBD = calculatePBD(filter,veryBase.PitchAngles, ...
    veryBase.Mflap1_store,veryBase.Medge1_store, ...
    veryBase.Mflap2_store,veryBase.Medge2_store, ...
    veryBase.Mflap3_store,veryBase.Medge3_store); % [kNm deg]
% Downstream WT2
VBL_result.WT2.power = calculatePower(filter,veryBase.Powerturb2_store,D_NREL5MW,U_inflow); % [MW]
% BL_result.WT2.power2 = calculatePower2(filter,Baseline.Powerturb2_store);
VBL_result.WT2.DEL = calculateDEL(filter, ...
    veryBase.Mflap1turb2_store,veryBase.Medge1turb2_store, ...
    veryBase.Mflap2turb2_store,veryBase.Medge2turb2_store, ...
    veryBase.Mflap3turb2_store,veryBase.Medge3turb2_store, ...
    timeStep); % [Nm]
VBL_result.WT2.PBD = calculatePBD(filter,veryBase.PitchAnglesturb2, ...
    veryBase.Mflap1turb2_store,veryBase.Medge1turb2_store, ...
    veryBase.Mflap2turb2_store,veryBase.Medge2turb2_store, ...
    veryBase.Mflap3turb2_store,veryBase.Medge3turb2_store); % [kNm deg]
VBL_result.WT1.DEL_flapwise = mean(VBL_result.WT1.DEL(1)+VBL_result.WT1.DEL(3)+VBL_result.WT1.DEL(5));
VBL_result.WT1.DEL_edgewise = mean(VBL_result.WT1.DEL(2)+VBL_result.WT1.DEL(4)+VBL_result.WT1.DEL(6));
VBL_result.WT2.DEL_flapwise = mean(VBL_result.WT2.DEL(1)+VBL_result.WT2.DEL(3)+VBL_result.WT2.DEL(5));
VBL_result.WT2.DEL_edgewise = mean(VBL_result.WT2.DEL(2)+VBL_result.WT2.DEL(4)+VBL_result.WT2.DEL(6));
VBL_result.All.power = VBL_result.WT1.power + VBL_result.WT2.power;
% BL_result.All.power2 = BL_result.WT1.power2 + BL_result.WT2.power2;
VBL_result.All.DEL_flapwise = VBL_result.WT1.DEL_flapwise + VBL_result.WT2.DEL_flapwise;
VBL_result.All.DEL_edgewise = VBL_result.WT1.DEL_edgewise + VBL_result.WT2.DEL_edgewise;
VBL_result.All.PBD = VBL_result.WT1.PBD + VBL_result.WT2.PBD;

% =========== Baseline
% Upstream WT1
BL_result.WT1.power = calculatePower(filter,Baseline.Power_store,D_NREL5MW,U_inflow); % [MW]
% BL_result.WT1.power2 = calculatePower2(filter,Baseline.Power_store);
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
% BL_result.WT2.power2 = calculatePower2(filter,Baseline.Powerturb2_store);
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
% BL_result.All.power2 = BL_result.WT1.power2 + BL_result.WT2.power2;
BL_result.All.DEL_flapwise = BL_result.WT1.DEL_flapwise + BL_result.WT2.DEL_flapwise;
BL_result.All.DEL_edgewise = BL_result.WT1.DEL_edgewise + BL_result.WT2.DEL_edgewise;
BL_result.All.PBD = BL_result.WT1.PBD + BL_result.WT2.PBD;

% =========== Open-Loop
OL_result.WT1.power = calculatePower(filter,OL.Power_store,D_NREL5MW,U_inflow); % [MW]
% OL_result.WT1.power2 = calculatePower2(filter,OL.Power_store);
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
% OL_result.WT2.power2 = calculatePower2(filter,OL.Powerturb2_store);
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
% OL_result.All.power2 = OL_result.WT1.power2 + OL_result.WT2.power2;
OL_result.All.DEL_flapwise = OL_result.WT1.DEL_flapwise + OL_result.WT2.DEL_flapwise;
OL_result.All.DEL_edgewise = OL_result.WT1.DEL_edgewise + OL_result.WT2.DEL_edgewise;
OL_result.All.PBD = OL_result.WT1.PBD + OL_result.WT2.PBD;

% =========== Closed-Loop
CL_result.WT1.power = calculatePower(filter,CL.Power_store,D_NREL5MW,U_inflow); % [MW]
% CL_result.WT1.power2 = calculatePower2(filter,CL.Power_store);
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
% CL_result.WT2.power2 = calculatePower2(filter,CL.Powerturb2_store);
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
% CL_result.All.power2 = CL_result.WT1.power2 + CL_result.WT2.power2;
CL_result.All.DEL_flapwise = CL_result.WT1.DEL_flapwise + CL_result.WT2.DEL_flapwise;
CL_result.All.DEL_edgewise = CL_result.WT1.DEL_edgewise + CL_result.WT2.DEL_edgewise;
CL_result.All.PBD = CL_result.WT1.PBD + CL_result.WT2.PBD;

% === Power and DEL (Result Make Sense or Not)
if strcmp(powerDELAnalysis, 'Y')
    filter2 = 2000;
    x = [1 2 3];
    deltaPower = [(OL_result.WT1.power-BL_result.WT1.power)/(BL_result.WT1.power) (OL_result.WT2.power-BL_result.WT2.power)/(BL_result.WT2.power) (OL_result.All.power-BL_result.All.power)/(BL_result.All.power);
                  (CL_result.WT1.power-BL_result.WT1.power)/(BL_result.WT1.power) (CL_result.WT2.power-BL_result.WT2.power)/(BL_result.WT2.power) (CL_result.All.power-BL_result.All.power)/(BL_result.All.power)];
    deltaDELf = [(OL_result.WT1.DEL_flapwise-BL_result.WT1.DEL_flapwise)/(BL_result.WT1.DEL_flapwise) (OL_result.WT2.DEL_flapwise-BL_result.WT2.DEL_flapwise)/(BL_result.WT2.DEL_flapwise) (OL_result.All.DEL_flapwise-BL_result.All.DEL_flapwise)/(BL_result.All.DEL_flapwise);
                 (CL_result.WT1.DEL_flapwise-BL_result.WT1.DEL_flapwise)/(BL_result.WT1.DEL_flapwise) (CL_result.WT2.DEL_flapwise-BL_result.WT2.DEL_flapwise)/(BL_result.WT2.DEL_flapwise) (CL_result.All.DEL_flapwise-BL_result.All.DEL_flapwise)/(BL_result.All.DEL_flapwise)];
    deltaDELe = [(OL_result.WT1.DEL_edgewise-BL_result.WT1.DEL_edgewise)/(BL_result.WT1.DEL_edgewise) (OL_result.WT2.DEL_edgewise-BL_result.WT2.DEL_edgewise)/(BL_result.WT2.DEL_edgewise) (OL_result.All.DEL_edgewise-BL_result.All.DEL_edgewise)/(BL_result.All.DEL_edgewise);
                 (CL_result.WT1.DEL_edgewise-BL_result.WT1.DEL_edgewise)/(BL_result.WT1.DEL_edgewise) (CL_result.WT2.DEL_edgewise-BL_result.WT2.DEL_edgewise)/(BL_result.WT2.DEL_edgewise) (CL_result.All.DEL_edgewise-BL_result.All.DEL_edgewise)/(BL_result.All.DEL_edgewise)];

    figure('Name', 'Power & DEL', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(1, 3, 1)
    bar(x,deltaPower*100);
    xticks(x); 
    xticklabels({'WT1', 'WT2', 'WT1+WT2'}); 
    ylabel('\Delta Power [%]')
    legend('OL', 'CL', 'Location','northeast')
    title('Power')

    subplot(1, 3, 2)
    bar(x,deltaDELf*100);
    xticks(x); 
    xticklabels({'WT1', 'WT2', 'WT1+WT2'}); 
    ylabel('\Delta DEL [%]')
    legend('OL', 'CL', 'Location','northeast')
    title('DEL Flapwise')
    
    subplot(1, 3, 3)
    bar(x,deltaDELe*100);
    xticks(x); 
    xticklabels({'WT1', 'WT2', 'WT1+WT2'}); 
    ylabel('\Delta DEL [%]')
    legend('OL', 'CL', 'Location','northeast')
    title('DEL Edgewise')
    setfigpaper('Width',[40,0.3],'Interpreter','tex','FontSize',Font,'linewidth',lw)

%     % Look into Power
%     deltaPower2WT1 = [(OL_result.WT1.power2-BL_result.WT1.power2)./(BL_result.WT1.power2);
%                       (CL_result.WT1.power2-BL_result.WT1.power2)./(BL_result.WT1.power2)];
%     deltaPower2WT2 = [(OL_result.WT2.power2-BL_result.WT2.power2)./(BL_result.WT2.power2);
%                       (CL_result.WT2.power2-BL_result.WT2.power2)./(BL_result.WT2.power2)];
%     deltaPower2All = [(OL_result.All.power2-BL_result.All.power2)./(BL_result.All.power2);
%                       (CL_result.All.power2-BL_result.All.power2)./(BL_result.All.power2)];
%     
%     figure('Name', 'Just Power', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
%     subplot(1, 3, 1)
%     bar([1, 2], deltaPower2WT1*100);
%     xticks(x); 
%     xticklabels({'OL', 'CL'}); 
%     ylabel('\Delta Power [%]')
%     title('WT1 Power')
% 
%     subplot(1, 3, 2)
%     bar([1, 2], deltaPower2WT2*100);
%     xticks(x); 
%     xticklabels({'OL', 'CL'}); 
%     ylabel('\Delta Power [%]')
%     title('WT2 Power')
% 
%     subplot(1, 3, 3)
%     bar([1, 2], deltaPower2All*100);
%     xticks(x); 
%     xticklabels({'OL', 'CL'}); 
%     ylabel('\Delta Power [%]')
%     title('All Power')
%     setfigpaper('Width',[40,0.3],'Interpreter','tex','FontSize',Font,'linewidth',lw)
end

% === PBD
% PBD is only compared between the open-loop and closed-loop because
% baseline does not have PBD since Helix is not activated
if strcmp(PBDAnalysis, 'Y')
    deltaPBD = mean(CL_result.WT1.PBD) - mean(OL_result.WT1.PBD);
    disp(deltaPBD)
end

% === Numerical Comparsion between OL and CL
if strcmp(numericalAnalysis, 'Y')
    % Power
    PowWT1vs = (CL_result.WT1.power - OL_result.WT1.power)/OL_result.WT1.power;
    PowWT2vs = (CL_result.WT2.power - OL_result.WT2.power)/OL_result.WT2.power;
    PowAllvs = (CL_result.All.power - OL_result.All.power)/OL_result.All.power;
    % DEL
    DELeWT1vs = (CL_result.WT1.DEL_edgewise - OL_result.WT1.DEL_edgewise)/OL_result.WT1.DEL_edgewise;
    DELeWT2vs = (CL_result.WT2.DEL_edgewise - OL_result.WT2.DEL_edgewise)/OL_result.WT2.DEL_edgewise;
    DELeAllvs = (CL_result.All.DEL_edgewise - OL_result.All.DEL_edgewise)/OL_result.All.DEL_edgewise;
    DELfWT1vs = (CL_result.WT1.DEL_flapwise - OL_result.WT1.DEL_flapwise)/OL_result.WT1.DEL_flapwise;
    DELfWT2vs = (CL_result.WT2.DEL_flapwise - OL_result.WT2.DEL_flapwise)/OL_result.WT2.DEL_flapwise;
    DELfAllvs = (CL_result.All.DEL_flapwise - OL_result.All.DEL_flapwise)/OL_result.All.DEL_flapwise;
    % PBD
    PBDWT1vs = (mean(CL_result.WT1.PBD) - mean(OL_result.WT1.PBD))/mean(OL_result.WT1.PBD);
    deltaPBD = mean(CL_result.WT1.PBD) - mean(OL_result.WT1.PBD);

    % print result
    fprintf('\n==================== Power, DEL, PBD Comparison \n');
    fprintf('========== Power \n');
    fprintf('   WT1: %.2e %%\n', PowWT1vs*100);
    fprintf('   WT2: %.2e %%\n', PowWT2vs*100);
    fprintf('   ALL: %.2e %%\n', PowAllvs*100);
    fprintf('========== DEL Flapwise \n');
    fprintf('   WT1: %.2e %%\n', DELfWT1vs*100);
    fprintf('   WT2: %.2e %%\n', DELfWT2vs*100);
    fprintf('   ALL: %.2e %%\n', DELfAllvs*100);
    fprintf('========== DEL Edgewise \n');
    fprintf('   WT1: %.2e %%\n', DELeWT1vs*100);
    fprintf('   WT2: %.2e %%\n', DELeWT2vs*100);
    fprintf('   ALL: %.2e %%\n', DELeAllvs*100);
    fprintf('========== PBD \n');
    fprintf('   WT1: %.2e %%\n', PBDWT1vs*100);
    fprintf('   WT1: %.2e [kNm deg]\n', deltaPBD);
    
end

if strcmp(storyTellingBasic, 'Y')
    % Very Basic (Uniform Open-Loop Helix)
    VBLpower = [VBL_result.WT1.power VBL_result.WT2.power VBL_result.All.power];
    VBLDELflap = [VBL_result.WT1.DEL_flapwise VBL_result.WT2.DEL_flapwise VBL_result.All.DEL_flapwise];
    VBLDELedge = [VBL_result.WT1.DEL_edgewise VBL_result.WT2.DEL_edgewise VBL_result.All.DEL_edgewise];
    % Pre-Control (Open-Loop Helix)
    OLpower = [OL_result.WT1.power OL_result.WT2.power OL_result.All.power];
    OLDELflap = [OL_result.WT1.DEL_flapwise OL_result.WT2.DEL_flapwise OL_result.All.DEL_flapwise];
    OLDELedge = [OL_result.WT1.DEL_edgewise OL_result.WT2.DEL_edgewise OL_result.All.DEL_edgewise];
    % Post-Control (Closed-Loop Helix)
    CLpower = [CL_result.WT1.power CL_result.WT2.power CL_result.All.power];
    CLDELflap = [CL_result.WT1.DEL_flapwise CL_result.WT2.DEL_flapwise CL_result.All.DEL_flapwise];
    CLDELedge = [CL_result.WT1.DEL_edgewise CL_result.WT2.DEL_edgewise CL_result.All.DEL_edgewise];

%     % == Open-Loop System
%     figure('Name', 'OL vs ref', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
%     % Power
%     subplot(1, 3, 1)
%     b2 = bar(VBLpower, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
%     hold on;
%     b1 = bar(OLpower, 'EdgeColor', 'k', 'LineWidth', 1.5);
%     set(gca, 'XTickLabel', x_labels);
%     legend([b1, b2], {'Uniform', 'Conditions'}, 'Location', 'northwest');
%     ylabel('Power [MW]');
%     hold off;
%     
%     % DEL Flapwise
%     subplot(1, 3, 2)
%     b2 = bar(VBLDELflap, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
%     hold on;
%     b1 = bar(OLDELflap, 'EdgeColor', 'k', 'LineWidth', 1.5);
%     set(gca, 'XTickLabel', x_labels);
%     legend([b1, b2], {'Uniform', 'Conditions'}, 'Location', 'northwest');
%     ylabel('DEL Flapwise [Nm]');
%     hold off;
%     
%     % DEL Edgewise
%     subplot(1, 3, 3)
%     b2 = bar(VBLDELedge, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
%     hold on;
%     b1 = bar(OLDELedge, 'EdgeColor', 'k', 'LineWidth', 1.5);
%     set(gca, 'XTickLabel', x_labels);
%     legend([b1, b2], {'Uniform', 'Conditions'}, 'Location', 'northwest');
%     ylabel('DEL Edgewise [Nm]');
%     hold off;
%     setfigpaper('Width',[40,0.3],'Interpreter','tex','FontSize',Font,'linewidth',lw);
% 
%     % == Closed-Loop System
%     figure('Name', 'CL vs ref', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
%     % Power
%     subplot(1, 3, 1)
%     b2 = bar(VBLpower, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
%     hold on;
%     b1 = bar(CLpower, 'EdgeColor', 'k', 'LineWidth', 1.5);
%     set(gca, 'XTickLabel', x_labels);
%     legend([b1, b2], {'Uniform', 'Conditions'}, 'Location', 'northwest');
%     ylabel('Power [MW]');
%     hold off;
%     
%     % DEL Flapwise
%     subplot(1, 3, 2)
%     b2 = bar(VBLDELflap, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
%     hold on;
%     b1 = bar(CLDELflap, 'EdgeColor', 'k', 'LineWidth', 1.5);
%     set(gca, 'XTickLabel', x_labels);
%     legend([b1, b2], {'Uniform', 'Conditions'}, 'Location', 'northwest');
%     ylabel('DEL Flapwise [Nm]');
%     hold off;
%     
%     % DEL Edgewise
%     subplot(1, 3, 3)
%     b2 = bar(VBLDELedge, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
%     hold on;
%     b1 = bar(CLDELedge, 'EdgeColor', 'k', 'LineWidth', 1.5);
%     set(gca, 'XTickLabel', x_labels);
%     legend([b1, b2], {'Uniform', 'Conditions'}, 'Location', 'northwest');
%     ylabel('DEL Edgewise [Nm]');
%     hold off;
%     setfigpaper('Width',[40,0.3],'Interpreter','tex','FontSize',Font,'linewidth',lw);

    % == Open-Loop vs. Closed-Loop
    figure('Name', 'OL&CL vs ref', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    % Power
    subplot(1, 3, 1)
    b1 = bar([OLpower; CLpower]', 'grouped');
    hold on
    b2 = bar([VBLpower; VBLpower]', 'grouped');
    hold off
    b2(1).FaceColor = 'none'; 
    b2(1).EdgeColor = 'k';
    b2(1).LineStyle = '--';
    b2(1).LineWidth = 1.5;
    b2(2).FaceColor = 'none'; 
    b2(2).EdgeColor = 'k';
    b2(2).LineStyle = '--';
    b2(2).LineWidth = 1.5;
    set(gca, 'XTickLabel', x_labels);
    legend([b1(1), b1(2), b2], {'OL', 'CL', 'Uni'}, 'Location', 'northwest');
    ylabel('Power [MW]');

    % DEL Flapwise
    subplot(1, 3, 2)
    b1 = bar([OLDELflap; CLDELflap]', 'grouped');
    hold on
    b2 = bar([VBLDELflap; VBLDELflap]', 'grouped');
    hold off
    b2(1).FaceColor = 'none'; 
    b2(1).EdgeColor = 'k';
    b2(1).LineStyle = '--';
    b2(1).LineWidth = 1.5;
    b2(2).FaceColor = 'none'; 
    b2(2).EdgeColor = 'k';
    b2(2).LineStyle = '--';
    b2(2).LineWidth = 1.5;
    set(gca, 'XTickLabel', x_labels);
    legend([b1(1), b1(2), b2], {'OL', 'CL', 'Uni'}, 'Location', 'northwest');
    ylabel('DEL Flapwise [Nm]');

    % DEL Edgewise
    subplot(1, 3, 3)
    b1 = bar([OLDELedge; CLDELedge]', 'grouped');
    hold on
    b2 = bar([VBLDELedge; VBLDELedge]', 'grouped');
    hold off
    b2(1).FaceColor = 'none'; 
    b2(1).EdgeColor = 'k';
    b2(1).LineStyle = '--';
    b2(1).LineWidth = 1.5;
    b2(2).FaceColor = 'none'; 
    b2(2).EdgeColor = 'k';
    b2(2).LineStyle = '--';
    b2(2).LineWidth = 1.5;
    set(gca, 'XTickLabel', x_labels);
    legend([b1(1), b1(2), b2], {'OL', 'CL', 'Uni'}, 'Location', 'northwest');
    ylabel('DEL Edgewise [Nm]');
    setfigpaper('Width',[40,0.3],'Interpreter','tex','FontSize',Font,'linewidth',lw);
end