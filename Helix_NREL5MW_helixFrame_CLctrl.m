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
fileName = 'try.mat';   % Fixed Frame
turbineName = '.\Data\NREL5MW\';
caseName = 'Str0.3_U10_1Dd_10Hz_CCW\CLctrl\';

%% Load project and Initialize simulation
%this is setup using relative path and depends on the location of this file
calllib('QBladeDLL','createInstance',2,64)  % 64 for ring
calllib('QBladeDLL','setLibraryPath',DllPath)   % set lib path
calllib('QBladeDLL','loadSimDefinition',simFile)
calllib('QBladeDLL','initializeSimulation')
simTime = 5000;     % in timestep, actual time is simTime*timestep(Q-blade define)
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

%% Load internal model
buf_sys = load('Model\ModelOrder4_decoupled.mat');
buf_sys2 = load('Model\ModelOrder4_decoupled_delayed.mat');
decoupled_sys = buf_sys.decouple_sys;
decoupled_delayed_sys = buf_sys2.delayed_sys;

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

%% Defining Helix Control Setting
Str = 0.3;                          % Strouhal number
Helix_amplitude = 1;                % Helix amplitude                
Freq = Str*U_inflow/D_NREL5MW;      % From Str, in Hz
omega_e = Freq*2*pi;

t = linspace(1, simLen, simTime);
sigTilt_e = Helix_amplitude * ones(simTime, 1);  % basic
sigYaw_e = 0 * ones(simTime, 1);                 % basic

% % Step input to test basic properties
% steps = [0*ones(1, simTime/5) Helix_amplitude*ones(1, simTime/5) 0*ones(1, simTime/5) Helix_amplitude*ones(1, simTime/5) 0*ones(1, simTime/5)];
% % steps = [0*ones(1, simTime/10) Helix_amplitude*ones(1, simTime/10) 0*ones(1, simTime/10) Helix_amplitude*ones(1, simTime/10) 2*ones(1, simTime/10) -2*ones(1, simTime/10) 0*ones(1, simTime/10) Helix_amplitude*ones(1, simTime/10) -2*ones(1, simTime/10) 0*ones(1, simTime/10)];
% sigTilt_e = 0 * ones(simTime, 1);                  % 0 * ones(simTime, 1)
% sigYaw_e = steps;    % 0 * ones(simTime, 1)

% figure;
% plot(t, sigTilt_e);
% hold on
% plot(t, sigYaw_e);
% hold off
% legend('\beta_{tilt,e}', '\beta_{yaw,e}')

%% Define CL Ctrl setting
Ctrlers = [1 0 ; 0 1];      % very simple SISO controller
Trigger =100;      % Time that CL ctrl is triggered
Tilt_r = 2 * ones(simTime, 1);
Yaw_r = 1 * ones(simTime, 1);
r = [Tilt_r Yaw_r];         % reference signal
e = zeros(simTime, 2);      % error
u = zeros(simTime, 2);      % control input
y = zeros(simTime, 2);      % internal model output
ym = zeros(simTime, 2);     % WT measurement
ytilda = zeros(simTime, 2); % delayed sys output
ybuf_fir = zeros(simTime, 2);
yc = zeros(simTime, 2);     % combined output

% state space variables (add 1 due to the loop simulation)
xM = zeros(simTime+1, size(decoupled_sys.A, 1));
xMd = zeros(simTime+1, size(decoupled_delayed_sys.A, 1));

%% Defining LiDAR sampling 
% When you change this, don't forget to change the name of data.mat
LiDAR_x = 1*D_NREL5MW;   % Definition of x is pointing downwind
LiDAR_y = 0;
LiDAR_z = Hub_NREL5MW;   % Wind height
LiDAR_num_sample = 80;   % 5(ring) to speed up sampling, only 4 valid points

%% Simulation
% pre-define array to speed up code
TSR_store = zeros(simTime, 1);
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

% Deadtime Delay
timeDelay = LiDAR_x / U_inflow;

%% Real-time LPF
Fs = 1/timeStep;
Fc = 0.05;
Wn = Fc / (Fs / 2);

% Finite Impulse Response LPF (small phase lag in real-time)
n = 50; % Filter order
b_fir = fir1(n, Wn, 'low');
bufferSize = 50;   % 50
buffer = zeros(bufferSize, 2);
filterState1 = zeros(n, 1);
filterState2 = zeros(n, 1);
filterState3 = zeros(n, 1);
filterState4 = zeros(n, 1);

%% Adaptive filter for Smith Predictor
filter_order_adpFIR = 10;
DeadtimeDelay = 110;
omega_adpFIR = pi / (8 * DeadtimeDelay);
Wn_adpFIR = omega_adpFIR / (Fs / 2);
SP_adpFIR = fir1(filter_order_adpFIR, Wn_adpFIR, 'low');
filterState_adpFIR = zeros(filter_order_adpFIR, 1);

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
    centerZ = wakeCenter(1) - meanZ;  % 91.9411
    centerY = wakeCenter(2) - meanY;  % -3.1245
    center_e = invR_helix * [centerZ; centerY];
    [HF_helixCenter_filtered(i, 1), filterState3] = filter(b_fir, 1, center_e(1), filterState3);
    [HF_helixCenter_filtered(i, 2), filterState4] = filter(b_fir, 1, center_e(2), filterState4);

    % ====================  Control
    % I. Torque control to maintain optimal TSR of 9 
    omega_g = omega*N;                      % rotor to generator
    genTorque = K.*(omega_g*(2*pi/60))^2;

    % II. Wake mixing
    if i < Trigger
    % Normal Helix Control
    % 1. Get tilt and yaw signals
    beta_tilt_e = sigTilt_e(i);
    beta_yaw_e = sigYaw_e(i);
    % 2. Inverse MBC 
    invMBC = [1 cosd(Azimuth1) sind(Azimuth1);
              1 cosd(Azimuth2) sind(Azimuth2);
              1 cosd(Azimuth3) sind(Azimuth3)];
    invR_helix = [cos(omega_e*t(i)) -sin(omega_e*t(i)); 
                  sin(omega_e*t(i)) cos(omega_e*t(i))];
    % 3. Blade pitch signal
    betaTiltYaw = invR_helix * [beta_tilt_e; 
                                beta_yaw_e];    
    betaBlade_Helix = invMBC * [0; 
                                betaTiltYaw(1); 
                                betaTiltYaw(2)];
    % Send control signal to qblade
    calllib('QBladeDLL','setControlVars_at_num',[genTorque 0 ...
        betaBlade_Helix(1) betaBlade_Helix(2) betaBlade_Helix(3)],0)

    else
    % Activate CL control
    % 1. Get the error
    r_curr = r(i, :);
    yc_curr = yc(i, :);
    e_curr = r_curr - yc_curr;

    % 2. Get control input u based on error signal
    u_curr = e_curr * Ctrlers;  % 1*2

    % 3. Feed input to wind turbine and internal model
    % Internal Model 
    xM_curr = xM(i, :); % 1*4
    yM_curr = (decoupled_sys.C * xM_curr' + decoupled_sys.D * u_curr')'; % 1*2
    xM_next = (decoupled_sys.A * xM_curr' + decoupled_sys.B * u_curr')'; % 1*4
    xM(i+1, :) = xM_next;
    % Delayed Model 
    xMd_curr = xMd(i, :); % 1*n
    yMd_curr = (decoupled_delayed_sys.C * xMd_curr' + decoupled_delayed_sys.D * u_curr')'; % 1*2
    xMd_next = (decoupled_delayed_sys.A * xMd_curr' + decoupled_delayed_sys.B * u_curr')'; % 1*n
    xMd(i+1, :) = xMd_next;
    % Actual Wind Turbine
    % 1). Get tilt and yaw signal
    beta_tilt_e = u(1);
    beta_yaw_e = u(2);
    % 2). Inverse MBC 
    invMBC = [1 cosd(Azimuth1) sind(Azimuth1);
              1 cosd(Azimuth2) sind(Azimuth2);
              1 cosd(Azimuth3) sind(Azimuth3)];
    invR_helix = [cos(omega_e*t(i)) -sin(omega_e*t(i)); 
                  sin(omega_e*t(i)) cos(omega_e*t(i))];
    % 3). Blade pitch signal
    betaTiltYaw = invR_helix * [beta_tilt_e; 
                                beta_yaw_e];    
    betaBlade_Helix = invMBC * [0; 
                                betaTiltYaw(1); 
                                betaTiltYaw(2)];
    % Send control signal to qblade
    calllib('QBladeDLL','setControlVars_at_num',[genTorque 0 ...
        betaBlade_Helix(1) betaBlade_Helix(2) betaBlade_Helix(3)],0)
    
    % 4. Get different output ym, ytilda, y
    y_curr = yM_curr; % 1*2
    ytilda_curr = yMd_curr; % 1*2
    ym_curr = [center_e(1) center_e(2)]; % 1*2  MEASUREMENT!!!!!!!!!

    % 5. Adaptive filter & Combine outputs
    buf_y_curr = ym_curr - ytilda_curr;
    [ybuf_fir(i, 1), filterState_adpFIR] = filter(SP_adpFIR, 1, buf_y_curr, filterState_adpFIR);
    yc_curr2 = ybuf_fir + y_curr;

    % 6. Update variables array
    xM(i+1, :) = xM_next;
    xMd(i+1, :) = xMd_next;
    yc(i+1, :) = yc_curr2;
    y(i+1,:) = y_curr;              % Should it be i or i+1??????
    ytilda(i+1,:) = ytilda_curr;

    end
     
     

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
% save([turbineName caseName fileName], 'FF_helixCenter', ...
%                                       'FF_helixCenter_filtered', ...
%                                       'HF_helixCenter', ...
%                                       'HF_helixCenter_filtered', ...
%                                       'FF_beta', ...
%                                       'HF_beta');
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
% ylim([-1.5 1.5])
% xlim([0 300])
% xlabel("Time (s)");
% ylabel("Angle (deg)");
% title('Blade Pitch Signal')
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

% figure();
% % plot(HF_helixCenter(:, 1));
% % hold on;
% % plot(HF_helixCenter(:, 2));
% plot(HF_helixCenter_filtered(:, 1));
% hold on;
% plot(HF_helixCenter_filtered(:, 2));
% hold off;
% title('Center HF')
% legend('z_e', 'y_e')

% figure()
% plot(PitchAngles(:,1))
% hold on
% plot(PitchAngles(:,2))
% plot(PitchAngles(:,3))
% hold off
% xticks(0:100:length(PitchAngles));
% xticklabels(0:100*timeStep:length(PitchAngles)*timeStep);
% ylim([-1.25 1.25])
% % xlim([0 300])
% xlabel("Time (s)");
% ylabel("Angle (deg)");
% title('Blade Pitch Signal')
% legend('\beta_1','\beta_2','\beta_3')
% 
% figure()
% plot(FF_beta(:, 1))
% hold on
% plot(FF_beta(:, 2))
% hold off
% xticks(0:100:length(PitchAngles));
% xticklabels(0:100*timeStep:length(PitchAngles)*timeStep);
% ylim([-1.25 1.25])
% xlabel("Time (s)");
% ylabel("Angle (deg)");
% title('Rotor Disc Signal')
% legend('\beta_{tilt}', '\beta_{yaw}')

% ringVisualization(LiDAR_data, D_NREL5MW)
%% Unload Library 
% unloadlibrary 'QBladeDLL'