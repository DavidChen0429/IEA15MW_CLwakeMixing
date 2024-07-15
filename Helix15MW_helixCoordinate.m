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

% loadlibrary(DllPath,'.\QBladeDLLFunctions.h','alias','QBladeDLL') 
loadlibrary(DllPath,'QBladeLibInclude.h','alias','QBladeDLL') 
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

%% debug
simTime = 18000;     % in timestep, actual time is simTime*timestep(Q-blade define)
% simTime = 500;     % test GPU speed 
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
clear 
close all
U_inflow = 8;        % Inflow wind speed, same with the Q-blade setting
D_IEA15MW = 240;     % Rotor diameter
simTime = 18000;     % in timestep, actual time is simTime*timestep(Q-blade define)
% simTime = 500;     % test GPU speed 
timeStep = 0.1;    % same with the Q-blade setting
simLen = simTime * timeStep; % seconds
Str = 0.3;                          % Strouhal number
Helix_amplitude = 4;                % Helix amplitude                
Freq = Str*U_inflow/D_IEA15MW;      % From Str, in Hz

% Tranform matrix to Helix coordinate frame
% R_helix = [1 0 0; 
%            0 cos(Freq) -sin(Freq); 
%            0 sin(Freq) cos(Freq)];
% inverseR_helix = [1 0 0; 
%                   0 cos(Freq) sin(Freq); 
%                   0 -sin(Freq) cos(Freq)];

% Genereate the tilt and yaw signal (fixed frame)
t = linspace(timeStep, simLen, simTime);
midPoint = simLen / 2;
amplitudeTilt = Helix_amplitude * (t <= midPoint) + 2 * Helix_amplitude * (t > midPoint);
% sigTilt = amplitudeTilt .* sin(2*pi*Freq*t);  
sigTilt = Helix_amplitude * sin(2*pi*Freq*t);          
sigYaw = Helix_amplitude * sin(2*pi*Freq*t + pi/2);  % CCW
% sigYaw = amplitudeTilt .* sin(2*pi*Freq*t + pi/2);  % CCW

% Transfer to helix frame
thetaTilt_helixFrame_store = zeros(simTime, 1);
thetaYaw_helixFrame_store = zeros(simTime, 1);
thetaTilt_fixFrame_store = zeros(simTime, 1);
thetaYaw_fixFrame_store = zeros(simTime, 1);
% f = waitbar(0,'Runing');

for i = 1:1:simTime
    R_helix = [1 0 0; 
           0 cos(2*pi*Freq*t(i)) -sin(2*pi*Freq*t(i)); 
           0 sin(2*pi*Freq*t(i)) cos(2*pi*Freq*t(i))];
    invR_helix = [1 0 0; 
                  0 cos(2*pi*Freq*t(i)) sin(2*pi*Freq*t(i)); 
                  0 -sin(2*pi*Freq*t(i)) cos(2*pi*Freq*t(i))];
    theta_tilt = sigTilt(i);
    theta_yaw = sigYaw(i);
    theta_col = 0;
    MtMy_helix = R_helix * [theta_col; theta_tilt; theta_yaw]; 
    thetaTilt_helixFrame_store(i) = MtMy_helix(2);
    thetaYaw_helixFrame_store(i) = MtMy_helix(3);

    MtMy_fix = invR_helix * MtMy_helix; 
    thetaTilt_fixFrame_store(i) = MtMy_fix(2);
    thetaYaw_fixFrame_store(i) = MtMy_fix(3);
%     waitbar(i/endTime, f, sprintf('Convert Running: %.2f%%', (i/endTime)*100));
end

figure()
subplot(2, 1, 1)
plot(t, thetaTilt_helixFrame_store)
hold on
plot(t, thetaYaw_helixFrame_store)
hold off
xlabel('Time [s]')
ylabel('Magnitude')
% ylim([-1 5])
legend('M^e_{tilt}','M^e_{yaw}')

subplot(2, 1, 2)
plot(t, sigTilt)
hold on
plot(t, sigYaw)
hold off
xlabel('Time [s]')
ylabel('Magnitude')
legend('M_{tilt}','M_{yaw}')

%% Helix frame tilt and yaw reference
clear
close all
U_inflow = 8;        % Inflow wind speed, same with the Q-blade setting
D_IEA15MW = 240;     % Rotor diameter
simTime = 18000;     % in timestep, actual time is simTime*timestep(Q-blade define)
timeStep = 0.1;    % same with the Q-blade setting
simLen = simTime * timeStep; % seconds
Str = 0.3;                          % Strouhal number
Helix_amplitude = 4;                % Helix amplitude                
Freq = Str*U_inflow/D_IEA15MW;      % From Str, in Hz

% Tranform matrix to Helix coordinate frame
% R_helix = [1 0 0; 
%            0 cos(Freq) -sin(Freq); 
%            0 sin(Freq) cos(Freq)];
% inverseR_helix = [1 0 0; 
%                   0 cos(Freq) sin(Freq); 
%                   0 -sin(Freq) cos(Freq)];

% Genereate the tilt and yaw signal (fixed frame)
t = linspace(timeStep, simLen, simTime);
midPoint = simLen / 3;
amplitudeTilt = Helix_amplitude * (t <= midPoint) + 2 * Helix_amplitude * (t > midPoint);
sigTilt = amplitudeTilt .* sin(2*pi*Freq*t);  
% sigTilt = Helix_amplitude * sin(2*pi*Freq*t);          
sigYaw = Helix_amplitude * sin(2*pi*Freq*t - pi/2);  % CCW

% sigTilt_e = 0 * ones(simTime, 1);
sigTilt_e = [0 * ones(simTime/3, 1); 3 * ones(simTime*2/3, 1)];
sigYaw_e = 4 * ones(simTime, 1);

% Transfer to helix frame
% thetaTilt_fixFrame_store = zeros(simTime, 1);
% thetaYaw_fixFrame_store = zeros(simTime, 1);
thetaTilt_helixFrame_store = zeros(simTime, 1);
thetaYaw_helixFrame_store = zeros(simTime, 1);

for i = 1:1:simTime
    R_helix = [1 0 0; 
           0 cosd(2*pi*Freq*t(i)) -sind(2*pi*Freq*t(i)); 
           0 sind(2*pi*Freq*t(i)) cosd(2*pi*Freq*t(i))];
%     invR_helix = [1 0 0; 
%                   0 cos(2*pi*Freq*t(i)) sin(2*pi*Freq*t(i)); 
%                   0 -sin(2*pi*Freq*t(i)) cos(2*pi*Freq*t(i))];
%     theta_tilt_e = sigTilt_e(i);
%     theta_yaw_e = sigYaw_e(i);
    theta_tilt = sigTilt(i);
    theta_yaw = sigYaw(i);
    theta_col = 0;
%     MtMy_fix = invR_helix * [theta_col; theta_tilt_e; theta_yaw_e]; 
    MtMy_helix = R_helix * [theta_col; theta_tilt; theta_yaw]; 

%     thetaTilt_fixFrame_store(i) = MtMy_fix(2);
%     thetaYaw_fixFrame_store(i) = MtMy_fix(3);
    thetaTilt_helixFrame_store(i) = MtMy_helix(2);
    thetaYaw_helixFrame_store(i) = MtMy_helix(3);
end

figure()
subplot(2, 1, 1)
% plot(t, sigTilt_e)
plot(t, sigTilt)
hold on
% plot(t, sigYaw_e)
plot(t, sigYaw)
hold off
% ylim([-1 10])
xlabel('Time [s]')
ylabel('Magnitude')
% legend('\beta_{tilt,e}','\beta_{yaw,e}')
legend('\beta_{tilt}','\beta_{yaw}')

subplot(2, 1, 2)
% plot(t, thetaTilt_fixFrame_store)
plot(t, thetaTilt_helixFrame_store)
hold on
% plot(t, thetaYaw_fixFrame_store)
plot(t, thetaYaw_helixFrame_store)
hold off
xlabel('Time [s]')
ylabel('Magnitude')
% legend('\beta_{tilt}','\beta_{yaw}')
legend('\beta_{tilt,e}','\beta_{yaw,e}')

%% Defining LiDAR sampling 
% When you change this, don't forget to change the name of data.mat
LiDAR_x = 1*D_IEA15MW;   % Definition of x is pointing downwind
LiDAR_y = 0;
LiDAR_z = Wind_Height;   % Wind height
LiDAR_num_sample = 60;   % 5(ring) to speed up sampling, only 4 valid points 

%% Simulation
TSR_store = zeros(1, simTime);
LiDAR_data = zeros(1, simTime);

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
    windspeed = Circle_LiDAR_Parallel(LiDAR_x, LiDAR_y, LiDAR_z, LiDAR_num_sample);    
%     windspeed = Circle_LiDAR_Parallel_WakeCenter(LiDAR_x, LiDAR_y, LiDAR_z, LiDAR_num_sample, U_inflow);    % LiDAR with center calculation

    % Store values 
%     omega_store(i,:) = omega;
%     genTorqueQB_store(i,:) = genTorqueQB;
%     genTorque_store(i,:) = genTorque;
    TSR_store(i,:) = TSR;
%     AzimuthAngles(i,:) = [Azimuth1 Azimuth2 Azimuth3];
%     PitchAngles(i,:) = [Pitch1 Pitch2 Pitch3];
%     thetaTilt_store(i,:) = theta_tilt;
%     thetaYaw_store(i,:) = theta_yaw;

    if mod(i, 1/timeStep) == 0
%         fprintf('%d seconds.\n', i*timeStep);
        LiDAR_data = [LiDAR_data; windspeed];   % 1Hz sampling
    end

    waitbar(i/simTime, f, sprintf('Simulation Running: %.2f%%', (i/simTime)*100));

end
close(f)
%calllib('QBladeDLL','storeProject','15MW_Helix_Uni-U8_Str3.qpr') 
calllib('QBladeDLL','closeInstance')
save('.\Data\MAT\LiDAR_sampling\Uni\Str0.3_U8_1Dd_1Hz\Point2724_Timestep0.1_1800s_Parallel_changeMTilt.mat', 'LiDAR_data');
% save('.\Data\MAT\LiDAR_sampling\Uni\Str0.3_U8_1Dd_1Hz\Point2724_Timestep0.1_600s_Parallel_Center.mat', 'LiDAR_data');
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

% figure;
% plot(thetaTilt_store, 'r');
% hold on;
% plot(thetaYaw_store, 'b');
% xticks(0:100:length(thetaTilt_store));
% xticklabels(0:100*timeStep:length(thetaTilt_store)*timeStep);
% legend('\theta_{tilt}', '\theta_{yaw}')
% xlabel("Time (s)");
% title("Helix Signal")

%% Unload Library 
% unloadlibrary 'QBladeDLL'