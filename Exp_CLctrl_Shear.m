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
simFile = [SourcePath 'NREL5MW_1turbine_turbulence.sim'];
addpath('.\Functions');

loadlibrary(DllPath,'QBladeLibInclude.h','alias','QBladeDLL') 
m = libfunctions('QBladeDLL');
if isempty(m)
    fprintf('Error')
end

%% Data file 
turbineName = '.\Data\NREL5MW\';
caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\';
fileName = '1Turbines_CL_Helix_Shear.mat';
QprName = '1Turbines_CL_Helix_Shear.qpr';

%% Load project and Initialize simulation
%this is setup using relative path and depends on the location of this file
calllib('QBladeDLL','createInstance',2,64)  % 64 for ring
calllib('QBladeDLL','setLibraryPath',DllPath)   % set lib path
calllib('QBladeDLL','loadSimDefinition',simFile)
calllib('QBladeDLL','initializeSimulation')
simTime = 6000;     % in timestep, actual time is simTime*timestep(Q-blade define)
timeStep = 0.1;    % same with the Q-blade setting
simLen = simTime * timeStep; % seconds

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
Wind_Height = Hub_NREL5MW;
dimension = D_NREL5MW;     % span dim*dim meters
grid_point = 50;     % sqaure grid
Turb_time = 10;      % Simulation length of the windfield in seconds
Turb_dt = timeStep;  % Temporal resolution of the windfield
Turb_class = 'C';    % A, B, C
Turb_type = 'NTM';   % NTM, ETM, etc   
seed = 43;
vertInf = 0;         % Vertical inflow angle in degrees
horInf = 0;          % Horizontal inflow angle in degrees
% calllib('QBladeDLL', 'addTurbulentWind', ...
%     U_inflow,Hub_NREL5MW,Hub_NREL5MW,dimension,grid_point, ...
%     Turb_time,Turb_dt,Turb_class,Turb_type,seed,vertInf,horInf,1)

%% Defining Torque Control Setting
% This need to be changed when inflow windspeed is varied
K = 2.24;
N = 97;          % Gearbox ratio

%% Defining Helix Control Setting
Str = 0.3;                          % Strouhal number
Helix_amplitude = 0;                % Helix amplitude                
Freq = Str*U_inflow/D_NREL5MW;      % From Str, in Hz
omega_e = Freq*2*pi;
AzimuthOffset = 96; % 6 (2\pi) & 96; History -35

t = linspace(1, simLen, simTime);
% sigTilt_e = Helix_amplitude * ones(simTime, 1);  % basic
% sigYaw_e = 0 * ones(simTime, 1);                 % basic

% % Step input to test basic properties
% steps = [0*ones(1, simTime/5) Helix_amplitude*ones(1, simTime*2/5) 0*ones(1, simTime*2/5)];
% steps = [0*ones(1, simTime/10) Helix_amplitude*ones(1, simTime/10) 0*ones(1, simTime/10) Helix_amplitude*ones(1, simTime/10) 2*ones(1, simTime/10) -2*ones(1, simTime/10) 0*ones(1, simTime/10) Helix_amplitude*ones(1, simTime/10) -2*ones(1, simTime/10) 0*ones(1, simTime/10)];
sigTilt_e = 0 * ones(simTime, 1);   % steps
sigYaw_e = 0 * ones(simTime, 1);    % steps

% figure;
% plot(t, sigTilt_e);
% hold on
% plot(t, sigYaw_e);
% hold off
% legend('\beta_{tilt,e}', '\beta_{yaw,e}')

%% Define CL Ctrl setting
Trigger = ceil(simTime/5);      % Time that CL ctrl is triggered
e = zeros(simTime, 2);      % error
% integral_error = 0;         % error for integrator
u = zeros(simTime, 2);      % control input
y = zeros(simTime, 2);      % internal model output
ym = zeros(simTime, 2);     % WT measurement
ytilda = zeros(simTime, 2); % delayed sys output
ybuf_fir = zeros(simTime, 2);
bufy_error = zeros(simTime, 2);
yc = zeros(simTime, 2);     % combined output

% state space variables (add 1 due to the loop simulation)
xM = zeros(simTime+1, size(decoupled_sys.A, 1));
xMd = zeros(simTime+1, size(decoupled_delayed_sys.A, 1));

% Controller Design
W_s = tf([1, 1.6], [100, 1]);  % Emphasizes performance and disturbance rejection
W_t = tf([1, 1], [5, 1]);  % Emphasizes robustness and noise rejection
W_s_d = c2d(W_s, timeStep, 'tustin');
W_t_d = c2d(W_t, timeStep, 'tustin');
P = augw(decoupled_sys, W_s_d, [], W_t_d);  % augw creates the weighted augmented plant
ncont = 2; 
nmeas = 2; 
[K_hinf,CL,gamma] = hinfsyn(P,nmeas,ncont);
A_K = K_hinf.A;
B_K = K_hinf.B;
C_K = K_hinf.C;
D_K = K_hinf.D;

% H-inf variables
xk = zeros(simTime+ 1, length(A_K));
uk = y; % property of H inf
yk = zeros(simTime, length(C_K(:, 1)));

% Create reference
r = zeros(simTime, 2);     
% 1. Steps
reference_magnitude = [8.25 8.42];
r(Trigger:end, 1) = reference_magnitude(1)*ones(simTime+1-Trigger, 1);   % z_e
r(Trigger:end, 2) = reference_magnitude(2)*ones(simTime+1-Trigger, 1);   % y_e
% 2. Ramp
% reference_slope = [0.0025 0.0025]; % Define the slope of the ramp signal
% for tt = Trigger:simTime
%     r(tt, 1) = reference_slope(1) * (tt - Trigger);   % z_e ramp signal
%     r(tt, 2) = reference_slope(2) * (tt - Trigger);   % y_e ramp signal
% end
% 3. Ramp and Stop
% reference_slope = [0.0025 0.0015]; % Define the slope of the ramp signal
% endTime = (simTime*3)/5;
% for tt = Trigger:endTime
%     r(tt, 1) = reference_slope(1) * (tt - Trigger);   % z_e ramp signal
%     r(tt, 2) = reference_slope(2) * (tt - Trigger);   % y_e ramp signal
% end
% r(endTime:end, 1) = 5*ones(simTime+1-endTime, 1);
% r(endTime:end, 2) = 3*ones(simTime+1-endTime, 1);
% 4. Random step
% steps = cat(2, ...
%     0*ones(1, Trigger), 1*ones(1, (simTime-Trigger)/5), ...
%     -2*ones(1, (simTime-Trigger)/5), 2*ones(1, (simTime-Trigger)/5), ...
%     -1*ones(1, (simTime-Trigger)/5), 0*ones(1, (simTime-Trigger)/5));
% r(:, 1) = steps;
% r(:, 2) = steps;

% figure()
% plot(r(:, 1))
% hold on
% plot(r(:, 2))
% hold off

%% Defining LiDAR sampling 
% When you change this, don't forget to change the name of data.mat
LiDAR_x = 1*D_NREL5MW;   % Definition of x is pointing downwind
LiDAR_y = 0;
LiDAR_z = Hub_NREL5MW;   % Wind height
LiDAR_num_sample = 80;   % 5(ring) to speed up sampling, only 4 valid points

%% Simulation
% pre-define array to speed up code
TSR_store = zeros(simTime, 1);
Power_store = zeros(simTime, 1);
Moop1_store = zeros(simTime, 1);
Mip1_store = zeros(simTime, 1);
Mflap1_store = zeros(simTime, 1);
Medge1_store = zeros(simTime, 1);
Cp_store = zeros(simTime, 1);
FF_beta = zeros(simTime, 2);
HF_beta = zeros(simTime, 2);
FF_helixCenter_filtered = zeros(simTime, 2);
HF_helixCenter_filtered = zeros(simTime, 2);
PitchAngles = zeros(simTime, 3);
FF_helixCenter = zeros(simTime, 2);
HF_helixCenter = zeros(simTime, 2);
templateStruct = struct('x', [], 'y', [], 'z', [], 'u_x', [], 'u_y', [], 'u_z', [], 'u_norm', [], 'u_los', []);
LiDAR_data(simTime, 1) = templateStruct;

% Sliding window
ws_filter = 100;
ws_centering = ceil(1/(Freq * timeStep));

%% Real-time LPF
Fs = 1/timeStep;
Fc = 0.05;
Wn = Fc / (Fs / 2);

% Finite Impulse Response LPF (small phase lag in real-time)
n = 50; % Filter order
b_fir = fir1(n, Wn, 'low');
filterState1 = zeros(n, 1);
filterState2 = zeros(n, 1);
filterState3 = zeros(n, 1);
filterState4 = zeros(n, 1);

%% Adaptive filter for Smith Predictor
filter_order_adpFIR = 80;
omega_adpFIR = pi / (8 * DeadtimeDelay);
Wn_adpFIR = omega_adpFIR / (Fs / 2);
SP_adpFIR = fir1(filter_order_adpFIR, Wn_adpFIR, 'low');
filterState_adpFIR1 = zeros(filter_order_adpFIR, 1);
filterState_adpFIR2 = zeros(filter_order_adpFIR, 1);
% freqz(SP_adpFIR, 1);

%% Simulation
% start simulation
tic
f = waitbar(0,'Initializing Simulation');
for i = 1:1:simTime
    calllib('QBladeDLL','advanceTurbineSimulation')
    
    % Get current value
    omega = calllib('QBladeDLL','getCustomData_at_num',valuestr, 0, 0);
    genTorqueQB = calllib('QBladeDLL','getCustomData_at_num',valuestr2, 0, 0);
    TSR = calllib('QBladeDLL','getCustomData_at_num',valuestr3, 0, 0);
    Azimuth1 = calllib('QBladeDLL','getCustomData_at_num', Azimu1, 0, 0);
    Azimuth2 = calllib('QBladeDLL','getCustomData_at_num', Azimu2, 0, 0);
    Azimuth3 = calllib('QBladeDLL','getCustomData_at_num', Azimu3, 0, 0);
    Pitch1 = calllib('QBladeDLL','getCustomData_at_num', Pit1, 0, 0);
    Pitch2 = calllib('QBladeDLL','getCustomData_at_num', Pit2, 0, 0);
    Pitch3 = calllib('QBladeDLL','getCustomData_at_num', Pit3, 0, 0);
    Power = calllib('QBladeDLL','getCustomData_at_num', PowerVar, 0, 0);
    Cp = calllib('QBladeDLL','getCustomData_at_num', CpVar, 0, 0);
    Moop1 = calllib('QBladeDLL','getCustomData_at_num', Moop1Var, 0, 0);
    Mip1 = calllib('QBladeDLL','getCustomData_at_num', Mip1Var, 0, 0);

    % Define transform matrix 
    invMBC = [1 cosd(Azimuth1+AzimuthOffset) sind(Azimuth1+AzimuthOffset);
              1 cosd(Azimuth2+AzimuthOffset) sind(Azimuth2+AzimuthOffset);
              1 cosd(Azimuth3+AzimuthOffset) sind(Azimuth3+AzimuthOffset)];
    invR_helix = [cos(omega_e*t(i)) sin(omega_e*t(i)); 
                  -sin(omega_e*t(i)) cos(omega_e*t(i))];

    % ==================== LiDAR data sampling (Circle) 
    windspeed = Circle_LiDAR_Parallel(LiDAR_x, LiDAR_y, LiDAR_z, D_NREL5MW, LiDAR_num_sample); 
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
    centerZ = wakeCenter(1) - meanZ;
    centerY = wakeCenter(2) - meanY;
    center_e = invR_helix * [centerZ; centerY];
    [HF_helixCenter_filtered(i, 1), filterState3] = filter(b_fir, 1, center_e(1), filterState3);
    [HF_helixCenter_filtered(i, 2), filterState4] = filter(b_fir, 1, center_e(2), filterState4);
    % Sign change because of opposite model
    HF_helixCenter_filtered(i, :) = HF_helixCenter_filtered(i, :) * [-1 0; 0 1]; 

    % ==================== Control
    % I. Torque control to maintain optimal TSR of 9 
    omega_g = omega*N;                      % rotor to generator
    genTorque = K.*(omega_g*(2*pi/60))^2;

    % II. Wake mixing
    if i < Trigger
        % Normal Helix Control
        u(i, :) = [sigTilt_e(i) sigYaw_e(i)];
    else
        % Activate CL Control
        % Update controller
%         x_Kbuf = A_K * xk(i, :)' + B_K * y(i-1, :)'; % y / yc
%         xk(i+1, :) = x_Kbuf';
%         y_Kbuf = C_K * xk(i, :)' + D_K * y(i-1, :)';
%         yk(i, :) = y_Kbuf';

        x_Kbuf = A_K * xk(i, :)' + B_K * yc(i-1, :)'; % y / yc
        xk(i+1, :) = x_Kbuf';
        y_Kbuf = C_K * xk(i, :)' + D_K * yc(i-1, :)';
        yk(i, :) = y_Kbuf';
        
        % Get error / input of the plant
        u(i, :) = r(i, :) - yk(i, :);
    end

    % 1. Get tilt and yaw signals
    beta_tilt_e = u(i, 1);
    beta_yaw_e = u(i, 2);
    % 2. Inverse MBC & Blade pitch signal
    betaTiltYaw = invR_helix * [beta_tilt_e; 
                                beta_yaw_e];    
    betaBlade_Helix = invMBC * [0; 
                                betaTiltYaw(1); 
                                betaTiltYaw(2)];
    % Send control signal to qblade
    calllib('QBladeDLL','setControlVars_at_num',[genTorque 0 ...
        betaBlade_Helix(1) betaBlade_Helix(2) betaBlade_Helix(3)],0)
    
    % Building Smith Predictor
    % Internal Model
    u_curr = [beta_tilt_e beta_yaw_e];
    xM_curr = xM(i, :); % 1*4
    yM_curr = (decoupled_sys.C * xM_curr' + decoupled_sys.D * u_curr')'; % 1*2
    yM_curr = yM_curr * [1 0; 0 1]; 
    xM_next = (decoupled_sys.A * xM_curr' + decoupled_sys.B * u_curr')'; % 1*4
    xM(i+1, :) = xM_next;
    y(i, :) = yM_curr;
    % Delayed Model 
    xMd_curr = xMd(i, :); % 1*n
    yMd_curr = (decoupled_delayed_sys.C * xMd_curr' + decoupled_delayed_sys.D * u_curr')'; % 1*2
    yMd_curr = yMd_curr * [1 0; 0 1]; 
    xMd_next = (decoupled_delayed_sys.A * xMd_curr' + decoupled_delayed_sys.B * u_curr')'; % 1*n
    xMd(i+1, :) = xMd_next;  
    ytilda(i, :) = yMd_curr;
    % Wind turbine activation 
    ym(i, :) = [HF_helixCenter_filtered(i, 1) HF_helixCenter_filtered(i, 2)];
%     ym(i, :) = ym(i, :) * [-1 0; 0 1];

    % Adaptive filter check
    bufy_error(i, :) = ym(i, :) - ytilda(i, :);
    [ybuf_fir(i, 1), filterState_adpFIR1] = filter(SP_adpFIR, 1, bufy_error(i, 1), filterState_adpFIR1);
    [ybuf_fir(i, 2), filterState_adpFIR2] = filter(SP_adpFIR, 1, bufy_error(i, 2), filterState_adpFIR2);
    % Combine output
    yc(i, :) = ybuf_fir(i, :) + y(i, :);

    % ==================== Store values 
%     omega_store(i,:) = omega;
%     genTorqueQB_store(i,:) = genTorqueQB;
%     genTorque_store(i,:) = genTorque;
    TSR_store(i) = TSR;
    Power_store(i) = Power;
    Cp_store(i) = Cp;
    Moop1_store(i) = Moop1;
    Mip1_store(i) = Mip1;
    Mflap1_store(i) = Moop1*cosd(Pitch1) + Mip1*sind(Pitch1);
    Medge1_store(i) = -Moop1*sind(Pitch1) + Mip1*cosd(Pitch1);
    FF_beta(i,:) = [betaTiltYaw(1) betaTiltYaw(2)];
    HF_beta(i,:) = [beta_tilt_e beta_yaw_e];
%     AzimuthAngles(i,:) = [Azimuth1 Azimuth2 Azimuth3];
    PitchAngles(i,:) = [Pitch1 Pitch2 Pitch3];
    FF_helixCenter(i, :) = [wakeCenter(1) wakeCenter(2)]; % Z(tilt), Y(yaw)
    HF_helixCenter(i, :) = [center_e(1) center_e(2)];   % Ze(tilt), Ye(yaw) 
    LiDAR_data(i) = windspeed;

    waitbar(i/simTime, f, sprintf('Simulation Running: %.1f%%', (i/simTime)*100));

end
close(f)
calllib('QBladeDLL','storeProject', [turbineName caseName QprName])
calllib('QBladeDLL','closeInstance')
% save([turbineName caseName fileName], 'LiDAR_data', ...
%                                       'FF_helixCenter', ...
%                                       'FF_helixCenter_filtered', ...
%                                       'HF_helixCenter', ...
%                                       'HF_helixCenter_filtered', ...
%                                       'FF_beta', ...
%                                       'HF_beta');
% save([turbineName caseName fileName], 'FF_helixCenter', ...
%                                       'FF_helixCenter_filtered', ...
%                                       'HF_helixCenter', ...
%                                       'HF_helixCenter_filtered', ...
%                                       'FF_beta', ...
%                                       'HF_beta', ...
%                                       'u', ...
%                                       'e', ...
%                                       'r', ...
%                                       'y', ...
%                                       'ym', ...
%                                       'ytilda', ...
%                                       'yc');
save([turbineName caseName fileName], 'Power_store', ...
                                      'Cp_store', ...
                                      'Moop1_store', ...
                                      'Mip1_store', ...
                                      'Mflap1_store', ...
                                      'Medge1_store', ...
                                      'PitchAngles');
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

% % Comparison between Delayed Model and Wind Turbine Real Output
% figure('Name', 'Output Comparison', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
% plot((1:length(ytilda)) * timeStep, ytilda(:, 1),'m','LineWidth', 1)
% hold on
% plot((1:length(ytilda)) * timeStep, ytilda(:, 2),'b', 'LineWidth', 1)
% % plot((1:length(ytilda)) * timeStep, ym(:, 1),'m','LineWidth', 1)
% % hold on
% % plot((1:length(ytilda)) * timeStep, ym(:, 2),'b', 'LineWidth', 1)
% plot((1:length(ym)) * timeStep, ym(:, 1),'m--', 'LineWidth', 1)
% plot((1:length(ym)) * timeStep, ym(:, 2),'b--', 'LineWidth', 1)
% plot((1:length(sigYaw_e)) * timeStep, sigYaw_e,'m:', 'LineWidth', 0.5)
% plot((1:length(sigTilt_e)) * timeStep, sigTilt_e,'b:', 'LineWidth', 0.5)
% xline(trigger_time, '--k', 'Activate CL Ctrl', 'LabelOrientation', 'horizontal', 'LineWidth', 1);
% yline(0, '--', 'LineWidth', 1)
% hold off;
% xlabel('Time [s]')
% ylabel('Magnitude')
% title('Model Percision Check')
% legend('z_{Delay}','y_{Delay}','z_{WTm}','y_{WTm}','\beta_{yaw}','\beta_{tilt}')

% % check different errors
% figure('Name', 'Error', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
% plot((1:length(e)) * timeStep, e(:, 1),'m','LineWidth', 1)
% hold on
% plot((1:length(e)) * timeStep, e(:, 2),'b','LineWidth', 1)
% plot((1:length(e)) * timeStep, r(:, 1)-y(:, 1),'m--','LineWidth', 1)
% plot((1:length(e)) * timeStep, r(:, 2)-y(:, 2),'b--','LineWidth', 1)
% plot((1:length(e)) * timeStep, r(:, 1)-yc(:, 1),'m:','LineWidth', 1)
% plot((1:length(e)) * timeStep, r(:, 2)-yc(:, 2),'b:','LineWidth', 1)
% yline(0, '--', 'LineWidth', 1)
% xline(trigger_time, '--k', 'Activate CL Ctrl', 'LabelOrientation', 'horizontal', 'LineWidth', 1);
% hold off
% xlabel('Time [s]')
% ylabel('Magnitude')
% title('Error check')
% legend('e_{r,1}', 'e_{r,2}','e_{y,1}', 'e_{y,2}', 'e_{yc,1}', 'e_{yc,2}')

% % Compare wind turbine real output to the reference
% figure('Name', 'Controller performance', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
% plot((1:length(ym)) * timeStep, ym(:, 1),'m','LineWidth', 1)
% hold on
% plot((1:length(ym)) * timeStep, ym(:, 2),'b','LineWidth', 1)
% plot((1:length(r)) * timeStep, delayseq(r(:, 1), DeadtimeDelay),'m--','LineWidth', 1)
% plot((1:length(r)) * timeStep, delayseq(r(:, 2), DeadtimeDelay),'k--','LineWidth', 1)
% yline(0, '--', 'LineWidth', 1)
% xline(trigger_time, '--k', 'Activate CL Ctrl', 'LabelOrientation', 'horizontal', 'LineWidth', 1);
% hold off
% xlabel('Time [s]')
% ylabel('Magnitude')
% title('Controller Performance Check')
% legend('y_{WTm1}','y_{WTm2}','r_z','r_y')

% % Adaptive filter check
% figure
% plot((1:length(bufy_error)) * timeStep, bufy_error(:, 1),'m','LineWidth', 1)
% hold on
% plot((1:length(bufy_error)) * timeStep, bufy_error(:, 2),'b','LineWidth', 1)
% plot((1:length(ybuf_fir)) * timeStep, ybuf_fir(:, 1),'m--','LineWidth', 1)
% plot((1:length(ybuf_fir)) * timeStep, ybuf_fir(:, 2),'b--','LineWidth', 1)
% yline(0, '--', 'LineWidth', 1)
% hold off
% xlabel('Time [s]')
% ylabel('Magnitude')
% title('Adaptive Filter Check')
% legend('preFir_1','preFir_1','aftFir_1','aftFir_2')

% % See which output is dominate
% figure('Name', 'Output Component Check', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
% plot((1:length(y)) * timeStep, yc(:, 1), 'm','LineWidth', 1)
% hold on
% plot((1:length(y)) * timeStep, yc(:, 2), 'b','LineWidth', 1)
% plot((1:length(y)) * timeStep, y(:, 1), 'm--','LineWidth', 1)
% plot((1:length(y)) * timeStep, y(:, 2), 'b--','LineWidth', 1)
% xline(trigger_time, '--k', 'Activate CL Ctrl', 'LabelOrientation', 'horizontal', 'LineWidth', 1);
% yline(0, '--', 'LineWidth', 1)
% hold off
% title('Ouput Component Check')
% legend('y_{c1}','y_{c2}','y_{1}','y_{2}')

%% Unload Library 
% unloadlibrary 'QBladeDLL'