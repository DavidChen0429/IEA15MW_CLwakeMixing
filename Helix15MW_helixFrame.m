%% Helix of IEA15MW in script
clear
close all 
addpath('.\Functions');
%clc

%% Define paths
UserPath = 'C:\Users\DAVID CHEN\Desktop\TU_Delft\Thesis\IEA15MW_CLwakeMixing\'; 
QBladePath = 'C:\Users\DAVID CHEN\Desktop\TU_Delft\Thesis\QBladeEE_2.0.6.4\'; 
SourcePath = [UserPath 'Source\'];
DllPath = [QBladePath 'QBladeEE_2.0.6.dll'];
simFile = [SourcePath 'IEA15MW_torque_Helix.sim'];
addpath('.\Functions');

loadlibrary(DllPath,'QBladeLibInclude.h','alias','QBladeDLL') 
m = libfunctions('QBladeDLL');
if isempty(m)
    fprintf('Error')
end

%% Data file 
fileName = '600s_Center_HF_tilt_s,d.mat';
dataPath = '.\Data\MAT\LiDAR_sampling\';
caseName = 'Uni\Str0.3_U8_1Dd_10Hz_CCW\';

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

%% Set Turbulent Wind
U_inflow = 8;        % Inflow wind speed, same with the Q-blade setting
D_IEA15MW = 240;     % Rotor diameter
Hub_IEA15MW = 150;   % Hub height
Wind_Height = Hub_IEA15MW;
dimension = D_IEA15MW;     % span dim*dim meters
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
K_list = [0.33358007e6 
          0.344363428e6 
          0.344363428e6 
          0.344363428e6]; % k accord with u_in = [7 8 9 10]
K = K_list(2);   % gain for k-omega^2, T/RS^2
N = 1;          % gearbox ratio of IEA15MW direct drive

%% Defining Helix Control Setting
Str = 0.3;                          % Strouhal number
Helix_amplitude = 4;                % Helix amplitude                
Freq = Str*U_inflow/D_IEA15MW;      % From Str, in Hz
omega_e = Freq*2*pi;

% sigTilt_e = 0 * ones(simTime, 1);
% sigTilt_e = [linspace(0, 8, simTime*9/20) linspace(8, 0, simTime*9/20) 0*ones(1, simTime/10)];
sigTilt_e = [4*ones(1, simTime/6) 3*ones(1, simTime/6) 2*ones(1, simTime/6) 1*ones(1, simTime/6) 0*ones(1, simTime/3)];
% sigTilt_e = [0*ones(1, simTime/10) linspace(0, 2, simTime*2/5) linspace(2, 0, simTime*2/5) 0*ones(1, simTime/10)];
sigYaw_e = -2 * ones(simTime, 1);
% sigYaw_e = [linspace(-4, -12, simTime*9/20) linspace(-12, -4, simTime*9/20) -4*ones(1, simTime/10)];
% sigYaw_e = [-6*ones(1, simTime/6) -5*ones(1, simTime/6) -4*ones(1, simTime/6) -3*ones(1, simTime/6) -2*ones(1, simTime/3)];
% sigYaw_e = [-2*ones(1, simTime/10) linspace(-2, 0, simTime*2/5) linspace(0, -2, simTime*2/5) -2*ones(1, simTime/10)];
t = linspace(timeStep, simLen, simTime);

%% Defining LiDAR sampling 
% When you change this, don't forget to change the name of data.mat
LiDAR_x = 1*D_IEA15MW;   % Definition of x is pointing downwind
LiDAR_y = 0;
LiDAR_z = Wind_Height;   % Wind height
LiDAR_num_sample = 60;   % 5(ring) to speed up sampling, only 4 valid points

%% Simulation
% pre-define array to speed up code
TSR_store = zeros(simTime, 1);
thetaTiltYaw_fixedFrame = zeros(simTime, 2);
thetaTiltYaw_helixFrame = zeros(simTime, 2);
templateStruct = struct('x', [], 'y', [], 'z', [], 'u_x', [], 'u_y', [], 'u_z', [], 'u_los', [], 'centerY', [], 'centerZ', []);
LiDAR_data(simTime, 1) = templateStruct;

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
    omega_g = omega*N;  % rotor to generator
    genTorque = K.*(omega_g)^2;

    % Helix control
    % 1. Get tilt and yaw signals
    theta_tilt_e = sigTilt_e(i);
    theta_yaw_e = sigYaw_e(i);
    % 2. Inverse MBC 
    invMBC = [1 cosd(Azimuth1) sind(Azimuth1);
              1 cosd(Azimuth2) sind(Azimuth2);
              1 cosd(Azimuth3) sind(Azimuth3)];
    invR_helix = [1 0 0;
                  0 cos(omega_e*t(i)) -sin(omega_e*t(i)); 
                  0 sin(omega_e*t(i)) cos(omega_e*t(i))];
    % 3. Blade pitch signal
    thetaTiltYaw = invR_helix * [0; theta_tilt_e; theta_yaw_e];    
    thetaBlade_Helix = invMBC * [thetaTiltYaw(1); thetaTiltYaw(2); thetaTiltYaw(3)];    

    % Send control signal to qblade
    calllib('QBladeDLL','setControlVars_at_num',[genTorque 0 ...
        thetaBlade_Helix(1) thetaBlade_Helix(2) thetaBlade_Helix(3)],0)

    % LiDAR data sampling (Ring)
    %windspeed = ZXTM_lidar(LiDAR_x, LiDAR_y, LiDAR_z, LiDAR_num_sample);    
%     windspeed = Circle_LiDAR_Parallel(LiDAR_x, LiDAR_y, LiDAR_z, LiDAR_num_sample);    
    windspeed = Circle_LiDAR_Parallel_WakeCenter(LiDAR_x, LiDAR_y, LiDAR_z, LiDAR_num_sample, U_inflow);    % LiDAR with center calculation

    % Store values 
%     omega_store(i,:) = omega;
%     genTorqueQB_store(i,:) = genTorqueQB;
%     genTorque_store(i,:) = genTorque;
    TSR_store(i) = TSR;
%     AzimuthAngles(i,:) = [Azimuth1 Azimuth2 Azimuth3];
%     PitchAngles(i,:) = [Pitch1 Pitch2 Pitch3];
    thetaTiltYaw_helixFrame(i,:) = [theta_tilt_e theta_yaw_e];
    thetaTiltYaw_fixedFrame(i,:) = [thetaTiltYaw(2) thetaTiltYaw(3)];
    LiDAR_data(i) = windspeed;
%     if mod(i, 1/timeStep) == 0
% %         fprintf('%d seconds.\n', i*timeStep);
%         LiDAR_data = [LiDAR_data; windspeed];   % 1Hz sampling
%     end
    waitbar(i/simTime, f, sprintf('Simulation Running: %.2f%%', (i/simTime)*100));

end
close(f)
%calllib('QBladeDLL','storeProject','15MW_Helix_Uni-U8_Str3.qpr') 
calllib('QBladeDLL','closeInstance')
save([dataPath caseName fileName], 'LiDAR_data', ...
                                   'thetaTiltYaw_fixedFrame', ...
                                   'thetaTiltYaw_helixFrame');
toc 

%% Visualization
figure;
plot(TSR_store)
xticks(0:100:length(TSR_store));
xticklabels(0:100*timeStep:length(TSR_store)*timeStep);
legend('TSR')
xlabel("Time (s)");
ylabel("TSR");

% figure;
% plot(genTorqueQB_store)
% hold on
% plot(genTorque_store)
% xticks(0:100:length(genTorqueQB_store));
% xticklabels(0:100*timeStep:length(genTorqueQB_store)*timeStep);
% xlabel("Time (s)");
% ylabel("Torque (Nm)")
% legend('QB HSS Torque','K omega^2')
% 
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
subplot(2, 1, 1)
plot(thetaTiltYaw_fixedFrame(:, 1));
hold on;
plot(thetaTiltYaw_fixedFrame(:, 2));
hold off;
title('Fixed Frame')
legend('\theta_{tilt}', '\theta_{yaw}')
subplot(2, 1, 2);
plot(thetaTiltYaw_helixFrame(:, 1));
hold on;
plot(thetaTiltYaw_helixFrame(: ,2));
hold off;
title('Helix Frame')
legend('\theta^e_{tilt}', '\theta^e_{yaw}')

% ringVisualization(LiDAR_data)

%% Unload Library 
% unloadlibrary 'QBladeDLL'