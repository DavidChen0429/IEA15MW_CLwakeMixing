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
fileName = 'HF_Uni_sysIDE.mat';   % Fixed Frame
turbineName = '.\Data\NREL5MW\';
caseName = 'Str0.3_U10_1Dd_10Hz_CCW\sysIDE\';

%% Load project and Initialize simulation
%this is setup using relative path and depends on the location of this file
calllib('QBladeDLL','createInstance',2,64)  % 64 for ring
calllib('QBladeDLL','setLibraryPath',DllPath)   % set lib path
calllib('QBladeDLL','loadSimDefinition',simFile)
calllib('QBladeDLL','initializeSimulation')
simTime = 1000;     % in timestep, actual time is simTime*timestep(Q-blade define)
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
n = simTime; 
seedA = [1 0 0 1]; 
seedB = [1 1 0 1]; 
tap = [1 0 0 1]; % polynomial x^4 + x + 1

% PRBN Signal A
stateA = seedA;
prbnA = zeros(1, n);
for i = 1:n
    prbnA(i) = stateA(end); % Output the last bit
    feedbackA = mod(sum(stateA .* tap), 2); % Calculate feedback using tap positions
    stateA = [feedbackA stateA(1:end-1)]; % Shift and insert feedback
end

% PRBN Signal B
stateB = seedB;
prbnB = zeros(1, n);
for i = 1:n
    prbnB(i) = stateB(end); % Output the last bit
    feedbackB = mod(sum(stateB .* tap), 2); % Calculate feedback using tap positions
    stateB = [feedbackB stateB(1:end-1)]; % Shift and insert feedback
end

% % Plot the PRBN sequences
% figure;
% subplot(2, 1, 1);
% stairs(prbnA, 'LineWidth', 1.5);
% title('Pseudo-Random Binary Noise (PRBN) Sequence A');
% xlabel('Sample Index');
% ylabel('Binary Value');
% grid on;
% 
% subplot(2, 1, 2);
% stairs(prbnB, 'LineWidth', 1.5);
% title('Pseudo-Random Binary Noise (PRBN) Sequence B');
% xlabel('Sample Index');
% ylabel('Binary Value');
% grid on;

%% Defining Helix Control Setting
Str = 0.3;                          % Strouhal number
Helix_amplitude = 3;                % Helix amplitude                
Freq = Str*U_inflow/D_NREL5MW;      % From Str, in Hz
omega_e = Freq*2*pi;

t = linspace(1, simLen, simTime);
sigTilt_e = 0 * ones(simTime, 1);                 % basic
sigYaw_e = -Helix_amplitude * ones(simTime, 1);   % basic
% sigTilt_e = [linspace(0, -6, simTime*9/20) linspace(-6, 0, simTime*9/20) 0*ones(1, simTime/10)];
% sigTilt_e = [4*ones(1, simTime/6) 3*ones(1, simTime/6) 2*ones(1, simTime/6) 1*ones(1, simTime/6) 0*ones(1, simTime/3)];
% sigYaw_e = [-6*ones(1, simTime/10) linspace(-6, 6, simTime*4/5) 6*ones(1, simTime/10)];
% sigYaw_e = [linspace(-3, 3, simTime*9/20) linspace(3, -3, simTime*9/20) -3*ones(1, simTime/10)];
% sigYaw_e = [-6*ones(1, simTime/6) -5*ones(1, simTime/6) -4*ones(1, simTime/6) -3*ones(1, simTime/6) -2*ones(1, simTime/3)];
% sigYaw_e = [-2*ones(1, simTime/10) linspace(-2, 0, simTime*2/5) linspace(0, -2, simTime*2/5) -2*ones(1, simTime/10)];

% for sysIDE
sigTilt_e = prbnA * 3;
sigYaw_e = prbnB * 3;

% figure;
% plot(t, sigTilt_e);
% hold on
% plot(t, sigYaw_e);
% hold off
% legend('\beta_{tilt,e}', '\beta_{yaw,e}')

%% Defining LiDAR sampling 
% When you change this, don't forget to change the name of data.mat
LiDAR_x = 1*D_NREL5MW;   % Definition of x is pointing downwind
LiDAR_y = 0;
LiDAR_z = Hub_NREL5MW;   % Wind height
LiDAR_num_sample = 80;   % 5(ring) to speed up sampling, only 4 valid points

%% Simulation
% pre-define array to speed up code
TSR_store = zeros(simTime, 1);
FF_theta = zeros(simTime, 2);
HF_theta = zeros(simTime, 2);
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
    theta_tilt_e = sigTilt_e(i);
    theta_yaw_e = sigYaw_e(i);
    % 2. Inverse MBC 
    invMBC = [1 cosd(Azimuth1) sind(Azimuth1);
              1 cosd(Azimuth2) sind(Azimuth2);
              1 cosd(Azimuth3) sind(Azimuth3)];
    invR_helix = [cos(omega_e*t(i)) -sin(omega_e*t(i)); 
                  sin(omega_e*t(i)) cos(omega_e*t(i))];
    % 3. Blade pitch signal
    thetaTiltYaw = invR_helix * [theta_tilt_e; 
                                 theta_yaw_e];    
    thetaBlade_Helix = invMBC * [0; 
                                 thetaTiltYaw(1); 
                                 thetaTiltYaw(2)];    

    % Send control signal to qblade
    calllib('QBladeDLL','setControlVars_at_num',[genTorque 0 ...
        thetaBlade_Helix(1) thetaBlade_Helix(2) thetaBlade_Helix(3)],0)

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
    FF_theta(i,:) = [thetaTiltYaw(1) thetaTiltYaw(2)];
    HF_theta(i,:) = [theta_tilt_e theta_yaw_e];
%     AzimuthAngles(i,:) = [Azimuth1 Azimuth2 Azimuth3];
    PitchAngles(i,:) = [Pitch1 Pitch2 Pitch3];
    FF_helixCenter(i, :) = [wakeCenter(1) wakeCenter(2)]; % Z(tilt), Y(yaw)
    HF_helixCenter(i, :) = [center_e(1) center_e(2)];   % Ze(tilt), Ye(yaw) 
    LiDAR_data(i) = windspeed;

%     if mod(i, 1/timeStep) == 0
% %         fprintf('%d seconds.\n', i*timeStep);
%         LiDAR_data = [LiDAR_data; windspeed];   % 1Hz sampling
%     end
    waitbar(i/simTime, f, sprintf('Simulation Running: %.1f%%', (i/simTime)*100));

end
close(f)
%calllib('QBladeDLL','storeProject','15MW_Helix_Uni-U8_Str3.qpr') 
calllib('QBladeDLL','closeInstance')
% save([turbineName caseName fileName], 'LiDAR_data', ...
%                                       'FF_helixCenter', ...
%                                       'HF_helixCenter', ...
%                                       'FF_theta', ...
%                                       'HF_theta');
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
plot(FF_theta(:, 1));
hold on;
plot(FF_theta(:, 2));
hold off;
title('\beta FF')
legend('\theta_{tilt}', '\theta_{yaw}')
subplot(2, 2, 3);
plot(HF_theta(:, 1));
hold on;
plot(HF_theta(: ,2));
hold off;
title('\beta_e HF')
legend('\theta^e_{tilt}', '\theta^e_{yaw}')
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

% ringVisualization(LiDAR_data, D_NREL5MW)
%% Unload Library 
% unloadlibrary 'QBladeDLL'