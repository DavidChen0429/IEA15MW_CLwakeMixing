%% Helix of IEA15MW in script
clear
close all 
%clc

%% Define paths
UserPath = 'C:\Users\DAVID CHEN\Desktop\TU_Delft\Thesis\IEA15MW_CLwakeMixing\'; 
QBladePath = 'C:\Users\DAVID CHEN\Desktop\TU_Delft\Thesis\QBladeEE_2.0.6.4\'; 
SourcePath = [UserPath 'Source\'];
DllPath = [QBladePath 'QBladeEE_2.0.6.dll'];
simFile = [SourcePath 'IEA15MW_torque_Helix.sim'];
addpath('.\Functions');

loadlibrary(DllPath,'.\QBladeDLLFunctions.h','alias','QBladeDLL') 
m = libfunctions('QBladeDLL');
if isempty(m)
    fprintf('Error')
end

%% Load project and Initialize simulation
%this is setup using relative path and depends on the location of this file
calllib('QBladeDLL','createInstance',2,64)  % 64 for ring
calllib('QBladeDLL','setLibraryPath',DllPath)   % set lib path
calllib('QBladeDLL','loadSimDefinition',simFile)
calllib('QBladeDLL','initializeSimulation')
simTime = 12000;     % in timestep, actual time is simTime*timestep(Q-blade define)
% simTime = 500;     % test GPU speed 
timeStep = 0.05;    % same with the Q-blade setting
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
Turb_class = 'C';    % A, B, C
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

% Genereate the tilt and yaw signal 
t = linspace(timeStep, simLen, simTime);
sigTilt = Helix_amplitude * sin(2*pi*Freq*t);          
sigYaw = Helix_amplitude * sin(2*pi*Freq*t + pi/2);  % CCW

%% Defining LiDAR sampling 
% When you change this, don't forget to change the name of data.mat
LiDAR_x = 1*D_IEA15MW;   % Definition of x is pointing downwind
LiDAR_y = 0;
LiDAR_z = Wind_Height;   % Wind height
LiDAR_num_sample = 26;   % 5(ring) to speed up sampling, only 4 valid points
LiDAR_data = [];         % Array that store the windspeed struct 

%% Simulation
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
    theta_tilt = sigTilt(i);
    theta_yaw = sigYaw(i);
    theta_col = 0;
    % 2. Inverse MBC 
    invMBC = [1 cosd(Azimuth1) sind(Azimuth1);
              1 cosd(Azimuth2) sind(Azimuth2);
              1 cosd(Azimuth3) sind(Azimuth3)];
    % 3. Blade pitch signal
    thetaBlade_Helix = invMBC * [theta_col; theta_tilt; theta_yaw];    

    % Send control signal to qblade
    calllib('QBladeDLL','setControlVars_at_num',[genTorque 0 ...
        thetaBlade_Helix(1) thetaBlade_Helix(2) thetaBlade_Helix(3)],0)

    % LiDAR data sampling (Ring)
    %windspeed = ZXTM_lidar(LiDAR_x, LiDAR_y, LiDAR_z, LiDAR_num_sample);    
    windspeed = CircleLiDAR(LiDAR_x, LiDAR_y, LiDAR_z, LiDAR_num_sample);    
    

    % Store values 
    omega_store(i,:) = omega;
    genTorqueQB_store(i,:) = genTorqueQB;
    genTorque_store(i,:) = genTorque;
    TSR_store(i,:) = TSR;
    AzimuthAngles(i,:) = [Azimuth1 Azimuth2 Azimuth3];
    PitchAngles(i,:) = [Pitch1 Pitch2 Pitch3];
    thetaTilt_store(i,:) = theta_tilt;
    thetaYaw_store(i,:) = theta_yaw;

    % When change the store frequency, change the Fs in FFT 
    % and the name of the data.mat accordingly
    if mod(i, 1/timeStep) == 0
        fprintf('%d seconds.\n', i*timeStep);
        LiDAR_data = [LiDAR_data; windspeed];   % store every second
    end
    %LiDAR_data = [LiDAR_data; windspeed];   % 50Hz, same as ZXTM-LiDAR

    waitbar(i/simTime,f,'Simulation Running')

end
close(f)
%calllib('QBladeDLL','storeProject','15MW_Helix_Uni-U8_Str3.qpr') 
calllib('QBladeDLL','closeInstance')
save('.\Data\MAT\IEA15_Helix_CCW_Str0.3_U8_Uni_600s_1Dd_1Hz_Circle150_windspeedData.mat', 'LiDAR_data');
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

figure;
plot(PitchAngles(:,1))
hold on
plot(PitchAngles(:,2))
plot(PitchAngles(:,3))
xticks(0:100:length(PitchAngles));
xticklabels(0:100*timeStep:length(PitchAngles)*timeStep);
xlabel("Time (s)");
ylabel("Angle (deg)");
legend('Blade 1','Blade 2','Blade 3')

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

% figure;
% plot(thetaTilt_store, 'r');
% hold on;
% plot(thetaYaw_store, 'b');
% xticks(0:100:length(thetaTilt_store));
% xticklabels(0:100*timeStep:length(thetaTilt_store)*timeStep);
% legend('\theta_{tilt}', '\theta_{yaw}')
% xlabel("Time (s)");
% title("Helix Signal")
