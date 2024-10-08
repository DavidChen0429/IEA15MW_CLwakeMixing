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
simFile = [SourcePath 'NREL5MW_Torque_Helix.sim'];
addpath('.\Functions');

loadlibrary(DllPath,'QBladeLibInclude.h','alias','QBladeDLL') 
m = libfunctions('QBladeDLL');
if isempty(m)
    fprintf('Error')
end

%% Data file 
fileName = 'train_120min_1bw_noise5%_AzimuthOffset.mat';   % Fixed Frame 'train_30min_1bw.mat'
turbineName = '.\Data\NREL5MW\';
caseName = 'Str0.3_U10_1Dd_10Hz_CCW\sysIDE\';

%% Load project and Initialize simulation
%this is setup using relative path and depends on the location of this file
calllib('QBladeDLL','createInstance',2,64)  % 64 for ring
calllib('QBladeDLL','setLibraryPath',DllPath)   % set lib path
calllib('QBladeDLL','loadSimDefinition',simFile)
calllib('QBladeDLL','initializeSimulation')
simTime = 72000;   % in timestep, actual time is simTime*timestep(Q-blade define)
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

%% Set Turbulent Wind
U_inflow = 10;        % Inflow wind speed, same with the Q-blade setting
D_NREL5MW = 126;     % Rotor diameter
Hub_NREL5MW = 90;   % Hub height
Wind_Height = Hub_NREL5MW;
dimension = D_NREL5MW;     % span dim*dim meters
grid_point = 50;     % sqaure grid
Turb_time = 10;      % Simulation length of the windfield in seconds
Turb_dt = timeStep;  % Temporal resolution of the windfield
Turb_class = 'A';    % A, B, C
Turb_type = 'NTM';   % NTM, ETM, etc   
seed = 43;
vertInf = 0;         % Vertical inflow angle in degrees
horInf = 0;          % Horizontal inflow angle in degrees
% calllib('QBladeDLL', 'addTurbulentWind', ...
%     U_inflow,Hub_IEA15MW,Hub_IEA15MW,dimension,grid_point, ...
%     Turb_time,Turb_dt,Turb_class,Turb_type,seed,vertInf,horInf,1)

%% Defining Torque Control Setting
% This need to be changed when inflow windspeed is varied
K = 2.24;
N = 97;          % Gearbox ratio

%% Signals for system IDE
% ======== Train data
bw = 0.0175;            % estimated bandwidth
N_signal = simTime;     % signal length [s] simLen
AMPL_signal = 1;        % amplitude
% === Pseudoradom Binary
N_prbn = simTime;       % signal length [s] simLen
AMPL_prbn = 1;          % amplitude
Ts_prbn = timeStep;     % sampling time [s] timeStep
F_prbn = 1*bw;          % cutoff frequency [Hz] 2*bandwidth (0.0175)
Fstop_prbn = 1*bw;      % band-stop filtered around this frequency
T0_prbn = 0;            % starting time [s]
P_prbn = 2;             % number of channels
IDEsig = idprbs(N_prbn,AMPL_prbn,Ts_prbn,F_prbn,Fstop_prbn,T0_prbn,P_prbn);
ns_prbn = floor((length(IDEsig)-N_prbn)/2);
sigTilt_e = IDEsig(ns_prbn+1:N_prbn+ns_prbn,1);   % tailor length
sigYaw_e = IDEsig(ns_prbn+1:N_prbn+ns_prbn,2);    % tailor length

% === Chirp
% signal_length = simTime;      
% t = 0:timeStep:simLen;  
% f0 = 0.5*bw;                   
% f1 = 1.5*bw;                  
% sigTilt_e = chirp(t, f0, simLen, f1) * AMPL_signal;
% sigYaw_e = chirp(t, f0, simLen, f1) * AMPL_signal;

% === Add disturbances (Gaussian noise)
disturbance = randn(N_signal, 2);                   % noise
noise_level = AMPL_signal * 0.05; % Adjust the noise level as needed
noise_tilt = noise_level * randn(size(sigTilt_e));
noise_yaw = noise_level * randn(size(sigYaw_e));
sigTilt_e = sigTilt_e + noise_tilt;
sigYaw_e = sigYaw_e + noise_yaw;

% [u1s,Du1,u2s,Du2] = sigscale(sigTilt_e,sigYaw_e); % signal scaling

% % ======== Test data
% Helix_amplitude = 1;
% steps = [0*ones(1, simTime/5) Helix_amplitude*ones(1, simTime/5) 0*ones(1, simTime/5) Helix_amplitude*ones(1, simTime/5) 0*ones(1, simTime/5)];
% sigTilt_e = steps;                  % 0 * ones(simTime, 1)
% sigYaw_e = steps;                   % 0 * ones(simTime, 1)
% step_size = 1000;
% num_steps = simTime / step_size;
% Helix_amplitude = 1;
% steps = repmat([0 Helix_amplitude], 1, num_steps/2);
% original_time = linspace(0, 1, length(steps));
% new_time = linspace(0, 1, simTime);
% stretched_signal = interp1(original_time, steps, new_time, 'previous');
% steps = stretched_signal;
% sigTilt_e = steps;   % 0 * ones(simTime, 1)
% sigYaw_e = 0 * ones(simTime, 1);    % 0 * ones(simTime, 1)

% Power Spectrum Density
[M1,F1] = pwelch(sigTilt_e,[],[],[],1/timeStep);
[M2,F2] = pwelch(sigYaw_e,[],[],[],1/timeStep);
figure
semilogx(F1,mag2db(M1),'k','LineWidth',1)
hold on
semilogx(F2,mag2db(M2),'r','LineWidth',1)
yline(0, '--', 'LineWidth', 1)
hold off
xlabel('Frequency [Hz]');
ylabel('Amplitude [dB]');
legend('\beta^e_{tilt}', '\beta^e_{yaw}')
title('Input PSD')

%% Helix Setting
Str = 0.3;                          % Strouhal number
Helix_amplitude = 1;                % Helix amplitude                
Freq = Str*U_inflow/D_NREL5MW;      % From Str, in Hz
omega_e = Freq*2*pi;
t = linspace(1, simLen, simTime);
AzimuthOffset = 0; % optimal 8 but with wrong relation

% !!! If input signals has been generated from 'idprbs', then note below
% line 
% sigTilt_e = 0 * ones(simTime, 1);                 % basic
% sigYaw_e = -Helix_amplitude * ones(simTime, 1);   % basic

%% Defining LiDAR sampling 
% When you change this, don't forget to change the name of data.mat
LiDAR_x = 1*D_NREL5MW;   % Definition of x is pointing downwind
LiDAR_y = 0;
LiDAR_z = Hub_NREL5MW;   % Wind height
LiDAR_num_sample = 80;

%% Simulation
% pre-define array to speed up code
TSR_store = zeros(simTime, 1);
FF_beta = zeros(simTime, 2);
HF_beta = zeros(simTime, 2);
PitchAngles = zeros(simTime, 3);
FF_helixCenter_filtered = zeros(simTime, 2);
HF_helixCenter_filtered = zeros(simTime, 2);
FF_helixCenter = zeros(simTime, 2);
HF_helixCenter = zeros(simTime, 2);
templateStruct = struct('x', [], 'y', [], 'z', [], 'u_x', [], 'u_y', [], 'u_z', [], 'u_norm', [], 'u_los', []);
LiDAR_data(simTime, 1) = templateStruct;

% Sliding window
ws_filter = 100;
ws_centering = ceil(1/(Freq * timeStep));

% Deadtime Delay
timeDelay = LiDAR_x / U_inflow;

%% Low pass filter property
Fs = 1/timeStep;
Fc = 0.05;
Wn = Fc / (Fs / 2);

% Finite Impulse Response LPF (small phase lag in real-time)
n = 50; % Filter order
b_fir = fir1(n, Wn, 'low');
bufferSize = 50;
buffer = zeros(bufferSize, 2);
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
    omega = calllib('QBladeDLL','getCustomData_at_num',valuestr, 0, 0);
    genTorqueQB = calllib('QBladeDLL','getCustomData_at_num',valuestr2, 0, 0);
    TSR = calllib('QBladeDLL','getCustomData_at_num',valuestr3, 0, 0);
    Azimuth1 = calllib('QBladeDLL','getCustomData_at_num', Azimu1, 0, 0);
    Azimuth2 = calllib('QBladeDLL','getCustomData_at_num', Azimu2, 0, 0);
    Azimuth3 = calllib('QBladeDLL','getCustomData_at_num', Azimu3, 0, 0);
    Pitch1 = calllib('QBladeDLL','getCustomData_at_num', Pit1, 0, 0);
    Pitch2 = calllib('QBladeDLL','getCustomData_at_num', Pit2, 0, 0);
    Pitch3 = calllib('QBladeDLL','getCustomData_at_num', Pit3, 0, 0);
    
    % Control
    % Torque control to maintain optimal TSR of 9 
    omega_g = omega*N;                      % rotor to generator
    genTorque = K.*(omega_g*(2*pi/60))^2;

    % Helix control
    % 1. Get tilt and yaw signals
    beta_tilt_e = sigTilt_e(i);
    beta_yaw_e = sigYaw_e(i);
    % 2. Inverse MBC 
%     invMBC = [1 cosd(Azimuth1) sind(Azimuth1);
%               1 cosd(Azimuth2) sind(Azimuth2);
%               1 cosd(Azimuth3) sind(Azimuth3)];
    invMBC = [1 cosd(Azimuth1+AzimuthOffset) sind(Azimuth1+AzimuthOffset);
              1 cosd(Azimuth2+AzimuthOffset) sind(Azimuth2+AzimuthOffset);
              1 cosd(Azimuth3+AzimuthOffset) sind(Azimuth3+AzimuthOffset)];
    invR_helix = [cos(omega_e*t(i)) sin(omega_e*t(i)); 
                  -sin(omega_e*t(i)) cos(omega_e*t(i))];
    % 3. Blade pitch signal
    betaTiltYaw = invR_helix * [beta_tilt_e; 
                                beta_yaw_e];    
    betaBlade_Helix = invMBC * [0; 
                                betaTiltYaw(1); 
                                betaTiltYaw(2)];    

    % Send control signal to qblade
    calllib('QBladeDLL','setControlVars_at_num',[genTorque 0 ...
        betaBlade_Helix(1) betaBlade_Helix(2) betaBlade_Helix(3)],0)

    % LiDAR data sampling (Ring) 
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
    centerZ = wakeCenter(1) - meanZ;  % 91.9411
    centerY = wakeCenter(2) - meanY;  % -3.1245
    center_e = invR_helix * [centerZ; centerY];
    [HF_helixCenter_filtered(i, 1), filterState3] = filter(b_fir, 1, center_e(1), filterState3);
    [HF_helixCenter_filtered(i, 2), filterState4] = filter(b_fir, 1, center_e(2), filterState4);

    % Store values 
%     omega_store(i,:) = omega;
%     genTorqueQB_store(i,:) = genTorqueQB;
%     genTorque_store(i,:) = genTorque;
    TSR_store(i) = TSR;
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
%calllib('QBladeDLL','storeProject','15MW_Helix_Uni-U8_Str3.qpr') 
calllib('QBladeDLL','closeInstance')
% save([turbineName caseName fileName], 'LiDAR_data', ...
%                                       'FF_helixCenter', ...
%                                       'FF_helixCenter_filtered', ...
%                                       'HF_helixCenter', ...
%                                       'HF_helixCenter_filtered', ...
%                                       'FF_beta', ...
%                                       'HF_beta');
save([turbineName caseName fileName], 'FF_helixCenter', ...
                                      'FF_helixCenter_filtered', ...
                                      'HF_helixCenter', ...
                                      'HF_helixCenter_filtered', ...
                                      'FF_beta', ...
                                      'HF_beta');
toc 

%% Visualization
% figure;
% plot(TSR_store)
% xticks(0:100:length(TSR_store));
% xticklabels(0:100*timeStep:length(TSR_store)*timeStep);
% legend('TSR')
% xlabel("Time (s)");
% ylabel("TSR");

% figure;
% plot(genTorqueQB_store)
% hold on
% plot(genTorque_store)
% xticks(0:100:length(genTorqueQB_store));
% xticklabels(0:100*timeStep:length(genTorqueQB_store)*timeStep);
% xlabel("Time (s)");
% ylabel("Torque (Nm)")
% legend('QB HSS Torque','K omega^2')

% figure;
% plot(PitchAngles(:,1))
% hold on
% plot(PitchAngles(:,2))
% plot(PitchAngles(:,3))
% xticks(0:100:length(PitchAngles));
% xticklabels(0:100*timeStep:length(PitchAngles)*timeStep);
% xlabel("Time (s)");
% ylabel("Angle (deg)");
% legend('Blade 1','Blade 2','Blade 3')

% figure;
% plot(AzimuthAngles(:,1))
% hold on
% plot(AzimuthAngles(:,2))
% plot(AzimuthAngles(:,3))
% xticks(0:100:length(AzimuthAngles));
% xticklabels(0:100*timeStep:length(AzimuthAngles)*timeStep);
% xlabel("Time (s)");
% ylabel("Angle (deg)");
% legend('Blade 1','Blade 2','Blade 3')

figure;
subplot(2, 2, 1)
plot(FF_beta(:, 1));
hold on;
plot(FF_beta(:, 2));
hold off;
title('\beta FF')
legend('\beta_{tilt}', '\beta_{yaw}')
subplot(2, 2, 3);
plot(HF_beta(:, 1));
hold on;
plot(HF_beta(: ,2));
hold off;
title('\beta_e HF')
legend('\beta^e_{tilt}', '\beta^e_{yaw}')
subplot(2, 2, 2)
plot(FF_helixCenter(:, 1));
hold on;
plot(FF_helixCenter(:, 2));
plot(FF_helixCenter_filtered(:, 1));
plot(FF_helixCenter_filtered(:, 2));
hold off;
title('Center FF')
legend('z', 'y', 'z2', 'y2')
subplot(2, 2, 4)
plot(HF_helixCenter(:, 1));
hold on;
plot(HF_helixCenter(:, 2));
plot(HF_helixCenter_filtered(:, 1));
plot(HF_helixCenter_filtered(:, 2));
hold off;
title('Center HF')
legend('z_e', 'y_e', 'z_e2', 'y_e2')

% figure;
% subplot(2, 2, 1)
% plot(FF_beta(:, 1));
% hold on;
% plot(FF_beta(:, 2));
% hold off;
% title('\beta FF')
% legend('\beta_{tilt}', '\beta_{yaw}')
% subplot(2, 2, 3);
% plot(HF_beta(:, 1));
% hold on;
% plot(HF_beta(: ,2));
% hold off;
% title('\beta_e HF')
% legend('\beta^e_{tilt}', '\beta^e_{yaw}')
% subplot(2, 2, 2)
% plot(FF_helixCenter_filtered(:, 1));
% hold on;
% plot(FF_helixCenter_filtered(:, 2));
% hold off;
% title('Center FF')
% legend('z_f', 'y_f')
% subplot(2, 2, 4)
% plot(HF_helixCenter_filtered(:, 1));
% hold on;
% plot(HF_helixCenter_filtered(:, 2));
% hold off;
% title('Center HF')
% legend('z_{f,e}', 'y_{f,e}')

% % Power Spectrum Density
% [M1,F1] = pwelch(HF_helixCenter_filtered(:, 1),[],[],[],1/Ts_prbn);
% [M2,F2] = pwelch(HF_helixCenter_filtered(:, 2),[],[],[],1/Ts_prbn);
% figure 
% semilogx(F1,mag2db(M1),'k','LineWidth',1)
% hold on
% semilogx(F2,mag2db(M2),'r','LineWidth',1)
% hold off
% xlabel('Frequency [Hz]');
% ylabel('Amplitude [dB]');
% legend('z_{f,e}', 'y_{f,e}');
% title('Output PSD');

% % SPA average between input and output
% [Ga,ws] = spa_avf(u,y,1,25,[],[],'hamming');
% Ga = frd(Ga,ws);

% ringVisualization(LiDAR_data, D_NREL5MW)
%% Unload Library 
% unloadlibrary 'QBladeDLL'