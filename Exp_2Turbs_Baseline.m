%% Helix of NREL5MW in script
clear
close all 
addpath('.\Functions');
%clc

%% Define paths
UserPath = 'C:\Users\DAVID CHEN\Desktop\TU_Delft\Thesis\IEA15MW_CLwakeMixing\'; 
QBladePath = 'C:\Users\DAVID CHEN\Desktop\TU_Delft\Thesis\QBladeEE_2.0.6.4\'; 
SourcePath = [UserPath 'Source\'];
DllPath = [QBladePath 'QBladeEE_2.0.6.dll'];
simFile = [SourcePath 'NREL5MW_2turbines_4D.sim'];
addpath('.\Functions');

loadlibrary(DllPath,'QBladeLibInclude.h','alias','QBladeDLL') 
m = libfunctions('QBladeDLL');
if isempty(m)
    fprintf('Error')
end

%% Data file (Chage this accordingly)
simTime = 9000;     % in timestep, actual time is simTime*timestep(Q-blade define)
timeStep = 0.1;    % same with tdhe Q-blade setting
simLen = simTime * timeStep; % seconds
saveOption = 'Y';

turbineName = '.\Data\NREL5MW\';
caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\2TurbinesLonger\';
fileName = '2Turbines_Baseline_4D_1inflowAngle.mat';
QprName = '2Turbines_Baseline_4D.qpr';

%% Load project and Initialize simulation
%this is setup using relative path and depends on the location of this file
calllib('QBladeDLL','createInstance',2,64)  % 64 for ring
calllib('QBladeDLL','setLibraryPath',DllPath)   % set lib path
calllib('QBladeDLL','loadSimDefinition',simFile)
calllib('QBladeDLL','initializeSimulation')

% Variables we care
valuestr = 'Rotational Speed [rpm]';
valuestr2 = 'Gen. HSS Torque [Nm]';
valuestr3 = 'Tip Speed Ratio [-]';
Azimu1 = 'Azimuthal Position Blade 1 [deg]';
Azimu2 = 'Azimuthal Position Blade 2 [deg]';
Azimu3 = 'Azimuthal Position Blade 3 [deg]';
Pit1 = 'Pitch Angle Blade 1 [deg]';
Pit2 = 'Pitch Angle Blade 2 [deg]';
Pit3 = 'Pitch Angle Blade 3 [deg]';
PowerVar = 'Aerodynamic Power [kW]';
CpVar = 'Power Coefficient [-]';
Moop1Var = 'Aero. OOP RootBend. Mom. Blade 1 [Nm]';
Mip1Var = 'Aero. IP RootBend. Mom. Blade 1 [Nm]';
Moop2Var = 'Aero. OOP RootBend. Mom. Blade 2 [Nm]';
Mip2Var = 'Aero. IP RootBend. Mom. Blade 2 [Nm]';
Moop3Var = 'Aero. OOP RootBend. Mom. Blade 3 [Nm]';
Mip3Var = 'Aero. IP RootBend. Mom. Blade 3 [Nm]';

%% Load internal model
buf_sys = load('Model\RightTransform_Azimuth96\ModelOrder4_noise1p_opposite_decoupled.mat');
decoupled_sys = buf_sys.OLi;

% Construct delayed model
DeadtimeDelay = 112;
z = tf('z', timeStep);
G = tf(buf_sys.OLi);
decoupled_delayed_sys = z^(-DeadtimeDelay) .* G;
decoupled_delayed_sys = ss(decoupled_delayed_sys);

%% Set Turbulent Wind
U_inflow = 10;        % Inflow wind speed, same with the Q-blade setting
D_NREL5MW = 126;     % Rotor diameter
Hub_NREL5MW = 90;   % Hub height

%% Defining Torque Control Setting
% This need to be changed when inflow windspeed is varied
K = 2.24;
N = 97;          % Gearbox ratio

%% Defining Helix Control Setting
Str = 0.3;                          % Strouhal number
Helix_amplitude = 3;                % Helix amplitude                
Freq = Str*U_inflow/D_NREL5MW;      % From Str, in Hz
omega_e = Freq*2*pi;
AzimuthOffset = 96; % 6 for pi/2 shift ;96 for pi shift (right relationship)

t = linspace(1, simLen, simTime);
sigTilt_e = Helix_amplitude*ones(simTime, 1);  % basic
sigYaw_e = Helix_amplitude*ones(simTime, 1);   % basic

% Reference is not used, but for comparison with CLctrl
Trigger = 0;
r = referenceGenerator(simTime,Trigger,(simTime*4)/5,'zero',2,0);

%% Defining LiDAR sampling 
% When you change this, don't forget to change the name of data.mat
LiDAR_x = 1*D_NREL5MW;   % Definition of x is pointing downwind 1*D_NREL5MW
LiDAR_y = 0;
LiDAR_z = Hub_NREL5MW;   % Wind height
LiDAR_num_sample = 80;   
MeasurementPos = 4*D_NREL5MW-50;

%% Simulation
% pre-define array to speed up code
TSR_store = zeros(simTime, 1);
Power_store = zeros(simTime, 1);
Cp_store = zeros(simTime, 1);
% blade 1
Moop1_store = zeros(simTime, 1);
Mip1_store = zeros(simTime, 1);
Mflap1_store = zeros(simTime, 1);
Medge1_store = zeros(simTime, 1);
% blade 2
Moop2_store = zeros(simTime, 1);
Mip2_store = zeros(simTime, 1);
Mflap2_store = zeros(simTime, 1);
Medge2_store = zeros(simTime, 1);
% blade 3
Moop3_store = zeros(simTime, 1);
Mip3_store = zeros(simTime, 1);
Mflap3_store = zeros(simTime, 1);
Medge3_store = zeros(simTime, 1);

TSRturb2_store = zeros(simTime, 1);
Powerturb2_store = zeros(simTime, 1);
Cpturb2_store = zeros(simTime, 1);
% blade 1
Moop1turb2_store = zeros(simTime, 1);
Mip1turb2_store = zeros(simTime, 1);
Mflap1turb2_store = zeros(simTime, 1);
Medge1turb2_store = zeros(simTime, 1);
% blade 2
Moop2turb2_store = zeros(simTime, 1);
Mip2turb2_store = zeros(simTime, 1);
Mflap2turb2_store = zeros(simTime, 1);
Medge2turb2_store = zeros(simTime, 1);
% blade 3
Moop3turb2_store = zeros(simTime, 1);
Mip3turb2_store = zeros(simTime, 1);
Mflap3turb2_store = zeros(simTime, 1);
Medge3turb2_store = zeros(simTime, 1);

TorqueStoreTurb1 = zeros(simTime, 1);
TorqueStoreTurb2 = zeros(simTime, 1);

UmeanStore = zeros(simTime, 1);     % flow info
TIStore = zeros(simTime, 1);        % flow info

FF_beta = zeros(simTime, 2);
HF_beta = zeros(simTime, 2);
FF_helixCenter_filtered = zeros(simTime, 2);
HF_helixCenter_filtered = zeros(simTime, 2);
PitchAngles = zeros(simTime, 3);
PitchAnglesturb2 = zeros(simTime, 3);
FF_helixCenter = zeros(simTime, 2);
HF_helixCenter = zeros(simTime, 2);
templateStruct = struct('x', [], 'y', [], 'z', [], 'u_x', [], 'u_y', [], 'u_z', [], 'u_norm', [], 'u_los', []);
LiDAR_data(simTime, 1) = templateStruct;

% Sliding window
ws_filter = 100;
ws_centering = ceil(1/(Freq * timeStep));

%% Real-time LPF
[b_fir, n] = FIR_LPF(1/timeStep, 0.05);
filterState1 = zeros(n, 1);
filterState2 = zeros(n, 1);
filterState3 = zeros(n, 1);
filterState4 = zeros(n, 1);

%% Simulation
% start simulation
tic
f = waitbar(0,'Initializing Simulation');
for i = 1:1:simTime
    calllib('QBladeDLL','advanceTurbineSimulation')
    
    % Get current value
    Azimuth1 = calllib('QBladeDLL','getCustomData_at_num', Azimu1, 0, 0);
    Azimuth2 = calllib('QBladeDLL','getCustomData_at_num', Azimu2, 0, 0);
    Azimuth3 = calllib('QBladeDLL','getCustomData_at_num', Azimu3, 0, 0);

    omega = calllib('QBladeDLL','getCustomData_at_num',valuestr, 0, 0);
    genTorqueQB = calllib('QBladeDLL','getCustomData_at_num',valuestr2, 0, 0);
    TSR = calllib('QBladeDLL','getCustomData_at_num',valuestr3, 0, 0);
    Pitch1 = calllib('QBladeDLL','getCustomData_at_num', Pit1, 0, 0);
    Pitch2 = calllib('QBladeDLL','getCustomData_at_num', Pit2, 0, 0);
    Pitch3 = calllib('QBladeDLL','getCustomData_at_num', Pit3, 0, 0);
    Power = calllib('QBladeDLL','getCustomData_at_num', PowerVar, 0, 0);
    Cp = calllib('QBladeDLL','getCustomData_at_num', CpVar, 0, 0);
    Moop1 = calllib('QBladeDLL','getCustomData_at_num', Moop1Var, 0, 0);
    Mip1 = calllib('QBladeDLL','getCustomData_at_num', Mip1Var, 0, 0);
    Moop2 = calllib('QBladeDLL','getCustomData_at_num', Moop2Var, 0, 0);
    Mip2 = calllib('QBladeDLL','getCustomData_at_num', Mip2Var, 0, 0);
    Moop3 = calllib('QBladeDLL','getCustomData_at_num', Moop3Var, 0, 0);
    Mip3 = calllib('QBladeDLL','getCustomData_at_num', Mip3Var, 0, 0);

    omega_turb2 = calllib('QBladeDLL','getCustomData_at_num',valuestr, 0, 1);
    genTorqueQB_turb2 = calllib('QBladeDLL','getCustomData_at_num',valuestr2, 0, 1);
    TSR_turb2 = calllib('QBladeDLL','getCustomData_at_num',valuestr3, 0, 1);
    Pitch1_turb2 = calllib('QBladeDLL','getCustomData_at_num', Pit1, 0, 1);
    Pitch2_turb2 = calllib('QBladeDLL','getCustomData_at_num', Pit2, 0, 1);
    Pitch3_turb2 = calllib('QBladeDLL','getCustomData_at_num', Pit3, 0, 1);
    Power_turb2 = calllib('QBladeDLL','getCustomData_at_num', PowerVar, 0, 1);
    Cp_turb2 = calllib('QBladeDLL','getCustomData_at_num', CpVar, 0, 1);
    Moop1_turb2 = calllib('QBladeDLL','getCustomData_at_num', Moop1Var, 0, 1);
    Mip1_turb2 = calllib('QBladeDLL','getCustomData_at_num', Mip1Var, 0, 1);
    Moop2_turb2 = calllib('QBladeDLL','getCustomData_at_num', Moop2Var, 0, 1);
    Mip2_turb2 = calllib('QBladeDLL','getCustomData_at_num', Mip2Var, 0, 1);
    Moop3_turb2 = calllib('QBladeDLL','getCustomData_at_num', Moop3Var, 0, 1);
    Mip3_turb2 = calllib('QBladeDLL','getCustomData_at_num', Mip3Var, 0, 1);

    % Define transform matrix 
    invMBC = [1 cosd(Azimuth1+AzimuthOffset) sind(Azimuth1+AzimuthOffset);
              1 cosd(Azimuth2+AzimuthOffset) sind(Azimuth2+AzimuthOffset);
              1 cosd(Azimuth3+AzimuthOffset) sind(Azimuth3+AzimuthOffset)];
    invR_helix = [cos(omega_e*t(i)) sin(omega_e*t(i)); 
                  -sin(omega_e*t(i)) cos(omega_e*t(i))];

    % ==================== LiDAR data sampling (Circle) 
    windspeed = Circle_LiDAR_Parallel(LiDAR_x, LiDAR_y, LiDAR_z, D_NREL5MW, LiDAR_num_sample); 
    windspeed2 = Circle_LiDAR_Parallel(MeasurementPos, LiDAR_y, LiDAR_z, D_NREL5MW, LiDAR_num_sample); 
    wakeCenter = HelixCenter(windspeed, U_inflow, D_NREL5MW);
    FF_helixCenter(i, :) = [wakeCenter(1) wakeCenter(2)]; % Z(tilt), Y(yaw)
    % Get the helix center from the helix frame
    % LPF the single element
    [FF_helixCenter_filtered(i, 1), filterState1] = filter(b_fir, 1, FF_helixCenter(i, 1), filterState1);
    [FF_helixCenter_filtered(i, 2), filterState2] = filter(b_fir, 1, FF_helixCenter(i, 2), filterState2);
    % Get the mean
    meanZ = Hub_NREL5MW;
    meanY = 0;
    if i > ws_centering
        meanZ = mean(FF_helixCenter(i-ws_centering:i, 1));
        meanY = mean(FF_helixCenter(i-ws_centering:i, 2));
    end
    % Low pass filter
    % Centering
%     centerZ = wakeCenter(1) - meanZ;  % 91.2632
%     centerY = wakeCenter(2) - meanY;  % -4.9713
    centerZ = wakeCenter(1) - 92.4522;  % data derived from the basecase
    centerY = wakeCenter(2) + 4.1492;   % data derived from the basecase
    center_e = invR_helix * [centerZ; centerY];
    [HF_helixCenter_filtered(i, 1), filterState3] = filter(b_fir, 1, center_e(1), filterState3);
    [HF_helixCenter_filtered(i, 2), filterState4] = filter(b_fir, 1, center_e(2), filterState4);
    % Sign change because of opposite model
    HF_helixCenter_filtered(i, :) = HF_helixCenter_filtered(i, :) * [-1 0; 0 1];

    % ==================== Control
    % I. Torque control to maintain optimal TSR of 9 
    omega_g = omega*N;                      % rotor to generator
    genTorque = K.*(omega_g*(2*pi/60))^2;
    omega_g_turb2 = omega_turb2*N;                      % rotor to generator
    genTorque_turb2 = K.*(omega_g_turb2*(2*pi/60))^2;

    % II. Wake mixing
    % 1. Get tilt and yaw signals
    beta_tilt_e = 0;
    beta_yaw_e = 0;
    % 2. Inverse MBC 
    % 3. Blade pitch signal
    betaTiltYaw = invR_helix * [beta_tilt_e; 
                                beta_yaw_e];    
    betaBlade_Helix = invMBC * [0; 
                                betaTiltYaw(1); 
                                betaTiltYaw(2)];    

    % Send control signal to qblade
    calllib('QBladeDLL','setControlVars_at_num',[genTorque 0 ...
        0 0 0],0)
    calllib('QBladeDLL','setControlVars_at_num',[genTorque_turb2 0 ...
        0 0 0],1)
    
    % ==================== Compute Flow Information
    Umean = mean(windspeed2.u_los);
    TI = std(windspeed2.u_los)/Umean;
    UmeanStore(i) = Umean;
    TIStore(i) = TI;

    % ==================== Store values 
%     omega_store(i,:) = omega;
%     genTorqueQB_store(i,:) = genTorqueQB;
%     genTorque_store(i,:) = genTorque;
    TSR_store(i) = TSR;
    Power_store(i) = Power;
    Cp_store(i) = Cp;
    % blade 1
    Moop1_store(i) = Moop1;
    Mip1_store(i) = Mip1;
    Mflap1_store(i) = Moop1*cosd(Pitch1) + Mip1*sind(Pitch1);
    Medge1_store(i) = -Moop1*sind(Pitch1) + Mip1*cosd(Pitch1);
    % blade 2
    Moop2_store(i) = Moop2;
    Mip2_store(i) = Mip2;
    Mflap2_store(i) = Moop2*cosd(Pitch2) + Mip2*sind(Pitch2);
    Medge2_store(i) = -Moop2*sind(Pitch2) + Mip2*cosd(Pitch2);
    % blade 3
    Moop3_store(i) = Moop3;
    Mip3_store(i) = Mip3;
    Mflap3_store(i) = Moop3*cosd(Pitch3) + Mip3*sind(Pitch3);
    Medge3_store(i) = -Moop3*sind(Pitch3) + Mip3*cosd(Pitch3);

    TSRturb2_store(i) = TSR_turb2;
    Powerturb2_store(i) = Power_turb2;
    Cpturb2_store(i) = Cp_turb2;
    % blade 1
    Moop1turb2_store(i) = Moop1_turb2;
    Mip1turb2_store(i) = Mip1_turb2;
    Mflap1turb2_store(i) = Moop1_turb2*cosd(Pitch1_turb2) + Mip1_turb2*sind(Pitch1_turb2);
    Medge1turb2_store(i) = -Moop1_turb2*sind(Pitch1_turb2) + Mip1_turb2*cosd(Pitch1_turb2);
    % blade 2
    Moop2turb2_store(i) = Moop2_turb2;
    Mip2turb2_store(i) = Mip2_turb2;
    Mflap2turb2_store(i) = Moop2_turb2*cosd(Pitch2_turb2) + Mip2_turb2*sind(Pitch2_turb2);
    Medge2turb2_store(i) = -Moop2_turb2*sind(Pitch2_turb2) + Mip2_turb2*cosd(Pitch2_turb2);
    % blade 3
    Moop3turb2_store(i) = Moop3_turb2;
    Mip3turb2_store(i) = Mip3_turb2;
    Mflap3turb2_store(i) = Moop3_turb2*cosd(Pitch3_turb2) + Mip3_turb2*sind(Pitch3_turb2);
    Medge3turb2_store(i) = -Moop3_turb2*sind(Pitch3_turb2) + Mip3_turb2*cosd(Pitch3_turb2);
    
    FF_beta(i,:) = [betaTiltYaw(1) betaTiltYaw(2)];
    HF_beta(i,:) = [beta_tilt_e beta_yaw_e];
%     AzimuthAngles(i,:) = [Azimuth1 Azimuth2 Azimuth3];
    PitchAngles(i,:) = [Pitch1 Pitch2 Pitch3];
    PitchAnglesturb2(i,:) = [Pitch1_turb2 Pitch2_turb2 Pitch3_turb2];
    FF_helixCenter(i, :) = [wakeCenter(1) wakeCenter(2)]; % Z(tilt), Y(yaw)
    HF_helixCenter(i, :) = [center_e(1) center_e(2)];   % Ze(tilt), Ye(yaw) 
%     LiDAR_data(i) = windspeed;

    TorqueStoreTurb1(i) = genTorqueQB;
    TorqueStoreTurb2(i) = genTorqueQB_turb2;

    waitbar(i/simTime, f, sprintf('Simulation Running: %.1f%%', (i/simTime)*100));

end
close(f)
if strcmp(saveOption, 'Y')
%     calllib('QBladeDLL','storeProject', [turbineName caseName QprName]) 
    save([turbineName caseName fileName], 'FF_helixCenter', ...
                                      'FF_helixCenter_filtered', ...
                                      'HF_helixCenter', ...
                                      'HF_helixCenter_filtered', ...
                                      'FF_beta', ...
                                      'HF_beta', ...
                                      'Power_store', ...
                                      'Powerturb2_store', ...
                                      'Cp_store', ...
                                      'Cpturb2_store', ...
                                      'Moop1_store', ...
                                      'Mip1_store', ...
                                      'Mflap1_store', ...
                                      'Medge1_store', ...
                                      'Moop2_store', ...
                                      'Mip2_store', ...
                                      'Mflap2_store', ...
                                      'Medge2_store', ...
                                      'Moop3_store', ...
                                      'Mip3_store', ...
                                      'Mflap3_store', ...
                                      'Medge3_store', ...
                                      'Moop1turb2_store', ...
                                      'Mip1turb2_store', ...
                                      'Mflap1turb2_store', ...
                                      'Medge1turb2_store', ...
                                      'Moop2turb2_store', ...
                                      'Mip2turb2_store', ...
                                      'Mflap2turb2_store', ...
                                      'Medge2turb2_store', ...
                                      'Moop3turb2_store', ...
                                      'Mip3turb2_store', ...
                                      'Mflap3turb2_store', ...
                                      'Medge3turb2_store', ...
                                      'PitchAngles', ...
                                      'PitchAnglesturb2', ...
                                      'UmeanStore', ...
                                      'TIStore', ...
                                      'TorqueStoreTurb1', ...
                                      'TorqueStoreTurb2');
end
calllib('QBladeDLL','closeInstance')
toc 

%% Visualization
trigger_time = Trigger * timeStep;

% Overall input and output
figure('Name', 'Overall Result', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 2, 1)
plot((1:length(FF_beta)) * timeStep, FF_beta(:, 1));
hold on;
plot((1:length(FF_beta)) * timeStep, FF_beta(:, 2));
xline(trigger_time, '--k', 'Activate CL Ctrl', 'LabelOrientation', 'horizontal', 'LineWidth', 1);
hold off;
xlabel('Time [s]')
title('\beta FF')
legend('\beta_{tilt}', '\beta_{yaw}')
subplot(2, 2, 3);
plot((1:length(HF_beta)) * timeStep, HF_beta(:, 1));
hold on;
plot((1:length(HF_beta)) * timeStep, HF_beta(: ,2));
xline(trigger_time, '--k', 'Activate CL Ctrl', 'LabelOrientation', 'horizontal', 'LineWidth', 1);
hold off;
xlabel('Time [s]')
title('\beta_e HF')
legend('\beta^e_{tilt}', '\beta^e_{yaw}')
subplot(2, 2, 2)
plot((1:length(HF_beta)) * timeStep, FF_helixCenter(:, 1));
hold on;
plot((1:length(HF_beta)) * timeStep, FF_helixCenter(:, 2));
plot((1:length(HF_beta)) * timeStep, FF_helixCenter_filtered(:, 1));
plot((1:length(HF_beta)) * timeStep, FF_helixCenter_filtered(:, 2));
xline(trigger_time, '--k', 'Activate CL Ctrl', 'LabelOrientation', 'horizontal', 'LineWidth', 1);
hold off;
xlabel('Time [s]')
title('Center FF')
legend('z', 'y', 'z_f', 'y_f')
subplot(2, 2, 4)
% plot((1:length(HF_beta)) * timeStep, HF_helixCenter(:, 1));
% hold on;
% plot((1:length(HF_beta)) * timeStep, HF_helixCenter(:, 2));
plot((1:length(HF_beta)) * timeStep, HF_helixCenter_filtered(:, 1));
hold on
plot((1:length(HF_beta)) * timeStep, HF_helixCenter_filtered(:, 2));
plot((1:length(r)) * timeStep, delayseq(r(:, 1), DeadtimeDelay),'m--','LineWidth', 0.5)
plot((1:length(r)) * timeStep, delayseq(r(:, 2), DeadtimeDelay),'k--','LineWidth', 0.5)
xline(trigger_time, '--k', 'Activate CL Ctrl', 'LabelOrientation', 'horizontal', 'LineWidth', 1);
yline(0, '--', 'LineWidth', 1)
hold off;
xlabel('Time [s]')
title('Center HF')
% legend('z_e', 'y_e', 'z_{e,f}', 'y_{e,f}')
legend('z_{e,f}', 'y_{e,f}')

% ringVisualization2(LiDAR_data, D_NREL5MW)

%% Unload Library 
% unloadlibrary 'QBladeDLL'