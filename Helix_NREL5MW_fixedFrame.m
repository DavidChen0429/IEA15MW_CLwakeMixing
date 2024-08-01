%% Helix of IEA15MW in script
clear
close all 
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
fileName = 'try.mat';
dataPath = '.\Data\MAT\LiDAR_sampling\';
caseName = 'Uni\Str0.3_U8_1Dd_10Hz_CCW\';

%% Load project and Initialize simulation
%this is setup using relative path and depends on the location of this file
calllib('QBladeDLL','createInstance',2,64)  
calllib('QBladeDLL','setLibraryPath',DllPath)   % set lib path
calllib('QBladeDLL','loadSimDefinition',simFile)
calllib('QBladeDLL','initializeSimulation')
simTime = 3000;     % in timestep, actual time is simTime*timestep(Q-blade define)
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
%     U_inflow,Hub_NREL5MW,Hub_NREL5MW,dimension,grid_point, ...
%     Turb_time,Turb_dt,Turb_class,Turb_type,seed,vertInf,horInf,1)

%% Defining Torque Control Setting
K = 2.24;
N = 97;          % Gearbox ratio

%% Defining Helix Control Setting
Str = 0.3;                          % Strouhal number
Helix_amplitude = 3;                % Helix amplitude                
Freq = Str*U_inflow/D_NREL5MW;      % From Str, in Hz
omega_e = Freq*2*pi;

% Genereate the tilt and yaw signal 
% t = linspace(timeStep, simLen, simTime);
t = linspace(1, simLen, simTime);
% cutPoint = simLen/4;
% cutPoint2 = simLen*3/4;
% amplitudeTilt = Helix_amplitude*(t<=cutPoint) + 2*Helix_amplitude*(t>cutPoint) - Helix_amplitude*(t>=cutPoint2);

cutPoint = simLen/6;
sigTilt = Helix_amplitude * sin(2*pi*Freq*t);          
% sigTilt(t >= cutPoint) = sigTilt(t >= cutPoint) + 2;
sigYaw = Helix_amplitude * sin(2*pi*Freq*t - pi/2);  % CCW
% sigYaw(t >= cutPoint) = sigYaw(t >= cutPoint) - 5;

% sigTilt = amplitudeTilt .* sin(2*pi*Freq*t);  
% sigYaw = amplitudeTilt .* sin(2*pi*Freq*t - pi/2);  % CCW
% figure()
% plot(t, sigTilt)
% hold on
% plot(t, sigYaw)
% hold off
% xlabel('Time [s]')
% ylabel('Magnitude')
% legend('\beta_{tilt}','\beta_{yaw}')

%% Defining LiDAR sampling 
% When you change this, don't forget to change the name of data.mat
LiDAR_x = 1*D_NREL5MW;   % Definition of x is pointing downwind
LiDAR_y = 0;
LiDAR_z = Hub_NREL5MW;   % Wind height
LiDAR_num_sample = 80;   % 5(ring) to speed up sampling, only 4 valid points

%% Preparation for Simulation 
% pre-define array to speed up code
TSR_store = zeros(simTime, 1);
thetaTiltYaw_fixedFrame = zeros(simTime, 2);
thetaTiltYaw_helixFrame = zeros(simTime, 2);
PitchAngles = zeros(simTime, 3);
helixCenter = zeros(simTime, 2);
helixCenter_e = zeros(simTime, 2);
templateStruct = struct('x', [], 'y', [], 'z', [], 'u_x', [], 'u_y', [], 'u_z', [], 'u_norm', [], 'u_los', []);
LiDAR_data(simTime, 1) = templateStruct;

% Sliding window
result = zeros(simTime, 2);
result_e = zeros(simTime, 2);
ws_filter = 500;
ws_centering = 1/(Freq * timeStep);

% Deadtime Delay
timeDelay = LiDAR_x / U_inflow;

% Low pass filter property
Fs = 1/timeStep;
Fc = 0.03;

%% Simulation
% start simulation
tic
f = waitbar(0,'Initializing Simulation');
for i = 1:1:simTime
    calllib('QBladeDLL','advanceTurbineSimulation')
  
    % Get current value
    omega = calllib('QBladeDLL','getCustomData_at_num',valuestr, 0.5, 0);
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
    theta_tilt = sigTilt(i);
    theta_yaw = sigYaw(i);
    theta_col = 0;
    % 2. Inverse MBC 
    invMBC = [1 cosd(Azimuth1) sind(Azimuth1);
              1 cosd(Azimuth2) sind(Azimuth2);
              1 cosd(Azimuth3) sind(Azimuth3)];
    % 3. Blade pitch signal
    thetaBlade_Helix = invMBC * [theta_col; theta_tilt; theta_yaw];   

    % 4. Corresponding helix frame value
    R_helix = [cos(omega_e*t(i)) sin(omega_e*t(i)); 
               -sin(omega_e*t(i)) cos(omega_e*t(i))];
    invR_helix = [cos(omega_e*t(i)) -sin(omega_e*t(i)); 
                  sin(omega_e*t(i)) cos(omega_e*t(i))];
    thetaTiltYaw_helix = R_helix * [theta_tilt; theta_yaw]; 

    % Send control signal to qblade
    calllib('QBladeDLL','setControlVars_at_num',[genTorque 0 ...
        thetaBlade_Helix(1) thetaBlade_Helix(2) thetaBlade_Helix(3)],0)

    % LiDAR data sampling (Ring)
    windspeed = Circle_LiDAR_Parallel(LiDAR_x, LiDAR_y, LiDAR_z, D_NREL5MW, LiDAR_num_sample);
    wakeCenter = HelixCenter(windspeed, U_inflow, D_NREL5MW);
    helixCenter(i, :) = [wakeCenter(1) wakeCenter(2)]; % Z(tilt), Y(yaw)

    % Get the helix center from the helix frame
    % LPF the single element
    if i > ws_filter
        buf_f1 = lowpassFilter(helixCenter(i-ws_filter:i, 1), Fs, Fc);
        buf_f2 = lowpassFilter(helixCenter(i-ws_filter:i, 2), Fs, Fc);
        wakeCenterZ_f = buf_f1(end);
        wakeCenterY_f = buf_f2(end);
        result(i, :) = [wakeCenterZ_f, wakeCenterY_f];
    end

    % Get the mean
    meanZ = Hub_NREL5MW;
    meanY = 0;
    if i > ws_centering
        meanZ = mean(helixCenter(i-ws_filter:i, 1));
        meanY = mean(helixCenter(i-ws_filter:i, 2));
    end
    
    % Low pass filter
    % Centering
    centerZ = wakeCenter(1) - meanZ;  % 91.9411
    centerY = wakeCenter(2) - meanY;  % -3.1245
    center_e = invR_helix * [centerZ; centerY];
%     center_e2 = invR_helix * [result(i, 1)-91.9411; result(i, 2)+3.1245];

    % Store values 
%     omega_store(i,:) = omega;
%     genTorqueQB_store(i,:) = genTorqueQB;
%     genTorque_store(i,:) = genTorque;
    TSR_store(i) = TSR;
    thetaTiltYaw_fixedFrame(i,:) = [theta_tilt, theta_yaw];
    thetaTiltYaw_helixFrame(i,:) = [thetaTiltYaw_helix(1), thetaTiltYaw_helix(2)];
%     AzimuthAngles(i,:) = [Azimuth1 Azimuth2 Azimuth3];
    PitchAngles(i,:) = [Pitch1 Pitch2 Pitch3];
    helixCenter(i, :) = [wakeCenter(1) wakeCenter(2)]; % Z(tilt), Y(yaw)
    helixCenter_e(i, :) = [center_e(1) center_e(2)];   % Ze(tilt), Ye(yaw) 
%     result_e(i, :) = [center_e2(1) center_e2(2)];
    LiDAR_data(i) = windspeed;

%     if mod(i, 1/timeStep) == 0
% %         fprintf('%d seconds.\n', i*timeStep);
%         LiDAR_data(i) = windspeed;   % 1Hz sampling
%     end
    waitbar(i/simTime, f, sprintf('Simulation Running: %.1f%%', (i/simTime)*100));

end
close(f)
%calllib('QBladeDLL','storeProject','15MW_Helix_Uni-U8_Str3.qpr') 
calllib('QBladeDLL','closeInstance')
save([dataPath caseName fileName], 'LiDAR_data', ...
                                   'helixCenter', ...    
                                   'thetaTiltYaw_fixedFrame', ...
                                   'thetaTiltYaw_helixFrame');
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

figure;
subplot(2, 2, 1)
plot(thetaTiltYaw_fixedFrame(:, 1));
hold on;
plot(thetaTiltYaw_fixedFrame(:, 2));
hold off;
title('\beta FF')
legend('\theta_{tilt}', '\theta_{yaw}')
subplot(2, 2, 2);
plot(thetaTiltYaw_helixFrame(:, 1));
hold on;
plot(thetaTiltYaw_helixFrame(: ,2));
hold off;
title('\beta_e HF')
legend('\theta^e_{tilt}', '\theta^e_{yaw}')
subplot(2, 2, 3)
plot(helixCenter(:, 1));
hold on;
plot(helixCenter(:, 2));
hold off;
title('Center FF')
legend('z', 'y')
subplot(2, 2, 4)
plot(helixCenter_e(:, 1));
hold on;
plot(helixCenter_e(:, 2));
hold off;
title('Center HF')
legend('z_e', 'y_e')

% figure;
% subplot(2, 1, 1)
% plot(helixCenter(:, 1));
% hold on;
% plot(helixCenter(:, 2));
% plot(result(:, 1));
% plot(result(:, 2));
% legend('z', 'y', 'z_f', 'y_f')
% hold off;
% subplot(2, 1, 2)
% plot(helixCenter_e(:, 1));
% hold on;
% plot(helixCenter_e(:, 2));
% plot(result_e(:, 1));
% plot(result_e(:, 2));
% legend('z', 'y', 'z_f', 'y_f')

% ringVisualization(LiDAR_data, D_NREL5MW)

%% Unload Library 
% unloadlibrary 'QBladeDLL'